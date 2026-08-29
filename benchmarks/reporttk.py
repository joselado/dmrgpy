"""LaTeX + plot report generation for run_benchmarks.py.

Kept separate from run_benchmarks.py (which owns problem construction,
timing and CLI) and from calctk.py (which owns the catalogue of
calculations) so that a report can be regenerated from a previously saved
results.json (e.g. after tweaking a table/plot) without rerunning the
DMRG sweep -- see run_benchmarks.py's --from-json option.

All presentation metadata for a calculation -- its section title, its
sweep axis, and any caveat the reader needs -- lives here, so adding a
calculation means one calctk.CALCS entry plus one line in each dict
below.
"""
import os
import math
import subprocess

import matplotlib
matplotlib.use("Agg")  # no display needed to render report plots
import matplotlib.pyplot as plt
import numpy as np

# Section titles. These are LaTeX source (some carry math), so they are
# written into the document verbatim and never passed through
# tex_escape() -- unlike the backend labels, which are plain text.
CALC_LABELS = {
    "gs": "Ground state energy",
    "static": "Static correlator",
    "dynamic": "Dynamical correlator (KPM)",
    "dynamic_td": "Dynamical correlator (real time + Fourier)",
    "excited": "First excited state",
    "entropy": "Entanglement entropy",
    "evolution": "Real-time evolution after a quench",
    "sector": "Ground state in the $S_z=0$ sector",
    "idmrg_gs": "Infinite chain: iDMRG energy density",
    "idmrg_vev": "Infinite chain: iDMRG $\\langle S^z\\rangle$",
    "vumps_gs": "Infinite chain: VUMPS energy density",
    "vumps_vev": "Infinite chain: VUMPS $\\langle S^z\\rangle$",
}

# The order sections appear in, and the order the summary aggregates in.
CALC_ORDER = ["gs", "static", "dynamic", "dynamic_td", "excited", "entropy",
        "evolution", "sector", "idmrg_gs", "idmrg_vev", "vumps_gs", "vumps_vev"]

# Sweep axis of each calculation: chain length for the finite family,
# bond dimension for the infinite one (an infinite chain has no length).
CALC_AXIS = {c: ("chi" if c.startswith(("idmrg_", "vumps_")) else "N")
        for c in CALC_ORDER}
AXIS_MATH = {"N": "$N$", "chi": "$\\chi$"}

# Plain-text stand-ins for the math in the labels above, used for PDF
# bookmarks: hyperref cannot put math in one, and warns about every
# section title that carries some unless it is given an alternative.
PDF_LABEL_SUBS = [("$\\langle S^z\\rangle$", "<Sz>"), ("$S_z=0$", "Sz=0")]
AXIS_PLOT = {"N": "chain length $N$", "chi": "bond dimension $\\chi$"}

# What each accuracy column is measured against, per reference kind.
REF_HEADERS = {
    "gs": "ED energy", "static": None, "excited": "ED energy",
    "entropy": None, "evolution": None,
    "density": "exact density", "zero": "exact $\\langle S^z\\rangle$",
}

CALC_NOTES = {
    "dynamic":
        "No independent ED cross-check is run here (a KPM spectral function "
        "is more involved to compare pointwise against an exact spectrum); "
        "physics correctness for this path is covered by "
        "\\texttt{tests/test\\_dynamical\\_correlator.py}. This section is "
        "timing-only.",
    "dynamic_td":
        "The same spectral function as the previous section, computed "
        "instead by real-time evolution followed by a Fourier transform "
        "(\\texttt{submode=\"TD\"}). On \\texttt{itensor\\_version} 3 and "
        "\\texttt{\"python\"} that evolution is TDVP; on 2 it is the "
        "MPO-Taylor propagator, so that column times a different algorithm. "
        "Timing-only, for the same reason as above.",
    "evolution":
        "Real-time evolution of $\\langle S^z_0\\rangle(t)$ after a quench "
        "from the Neel ground state of a staggered field to the Heisenberg "
        "Hamiltonian. \\texttt{itensor\\_version} 3 and \\texttt{\"python\"} "
        "run TDVP here; \\texttt{itensor\\_version=2} has no TDVP and uses "
        "the 2nd-order MPO-Taylor propagator instead "
        "(\\texttt{timedependent.py}) -- a different algorithm with a "
        "different accuracy/cost tradeoff, so read that column as such and "
        "not as a like-for-like comparison.",
    "sector":
        "The same ground state as the first section, confined to the "
        "$S_z=0$ quantum-number sector "
        "(\\texttt{set\\_conserved\\_sector(Sz=0)}), which is where the "
        "Heisenberg ground state lives -- so the energy must match the "
        "dense run and only the time differs. Only "
        "\\texttt{itensor\\_version=3} and \\texttt{\"python\"} implement "
        "quantum numbers at all; the other columns are blank by "
        "construction rather than by failure.",
    "idmrg_vev":
        "$\\langle S^z\\rangle$ vanishes by symmetry on this model, so the "
        "accuracy column is simply the deviation from zero. It is worth "
        "measuring separately from the energy: a correct energy density "
        "computed from a wrongly gauged unit cell shows up here and nowhere "
        "else.",
    "vumps_vev":
        "As above: the exact answer is zero by symmetry.",
    "idmrg_gs":
        "Swept over bond dimension $\\chi$, not chain length -- this is the "
        "thermodynamic limit (\\texttt{infinitechain.py}), and $\\chi$ is "
        "the knob that controls both cost and accuracy. The reference is "
        "the Bethe-ansatz energy density $1/4-\\ln 2$.",
    "vumps_gs":
        "Same model and same reference as the iDMRG section, reached by a "
        "different algorithm (\\texttt{gs\\_method=\"vumps\"}).",
}


# ---------------------------------------------------------------------
# JSON (de)serialization helpers
# ---------------------------------------------------------------------

def to_jsonable(obj):
    """Recursively convert numpy scalars/arrays (and complex numbers)
    into plain JSON-serializable Python types."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, complex):
        return obj.real
    return obj


# ---------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------

def tex_escape(s):
    s = str(s)
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
            "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
            "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(repl.get(c, c) for c in s)


def pdf_label(label):
    """The label with its math replaced by ASCII, for PDF bookmarks."""
    for math, plain in PDF_LABEL_SUBS:
        label = label.replace(math, plain)
    return label


def fmt_time(t):
    if t is None:
        return "--"
    if t < 1.0:
        return f"{t*1000:.1f}\\,ms"
    return f"{t:.3f}\\,s"


def fmt_err(e):
    if e is None:
        return "--"
    return f"{e:.2e}"


def short_err(err, n=80):
    s = "" if err is None else str(err)
    return (s[:n] + "...") if len(s) > n else s


# ---------------------------------------------------------------------
# Result lookups
# ---------------------------------------------------------------------

def _boxed(tabular_tex):
    """Wrap a tabular in adjustbox so it shrinks to fit the text width
    when there are enough backend columns to overflow it (e.g. all four
    of v2/v3/python/julia_live with their accuracy-column headers), and
    is left at its natural size otherwise."""
    return "\\begin{adjustbox}{max width=\\linewidth}\n%s\n\\end{adjustbox}" % tabular_tex


def _find(results, calc, backend, n):
    for r in results:
        if r["calc"] == calc and r["backend"] == backend and r["n"] == n:
            return r
    return None


def _series(results, calc, backend):
    """Sorted (n, time) pairs for successful runs of one (calc, backend)."""
    rows = [r for r in results if r["calc"] == calc and r["backend"] == backend and r["ok"]]
    rows.sort(key=lambda r: r["n"])
    return [r["n"] for r in rows], [r["time"] for r in rows]


def _deviation(value, ref):
    """|value - ref| for a scalar result, max|value - ref| for a vector
    one. Returns None if the shapes don't line up (a backend that
    returned a different number of points is a failure to report, not a
    number to average)."""
    if value is None or ref is None:
        return None
    if np.ndim(ref) == 0:
        try:
            return abs(float(np.real(value)) - float(np.real(ref)))
        except (TypeError, ValueError):
            return None
    a, b = np.asarray(value, dtype=float), np.asarray(ref, dtype=float)
    if a.shape != b.shape:
        return None
    return float(np.max(np.abs(a - b)))


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

def make_plots(results, backends, outdir):
    """One time-vs-size plot per calculation type actually present in
    `results`. Returns {calc: basename} for the calcs that got a plot
    (a calc with zero successful runs anywhere is skipped)."""
    plot_files = {}
    for calc in CALC_ORDER:
        if not any(r["calc"] == calc for r in results):
            continue
        fig, ax = plt.subplots(figsize=(5.2, 3.6))
        any_line = False
        for key, label, _iv in backends:
            ns, ts = _series(results, calc, key)
            if not ns:
                continue
            ax.plot(ns, ts, marker="o", label=label)
            any_line = True
        if not any_line:
            plt.close(fig)
            continue
        ax.set_xlabel(AXIS_PLOT[CALC_AXIS[calc]])
        ax.set_ylabel("wall time [s]")
        ax.set_yscale("log")
        ax.set_title(CALC_LABELS[calc])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fname = f"plot_{calc}.pdf"
        fig.savefig(os.path.join(outdir, fname))
        plt.close(fig)
        plot_files[calc] = fname
    return plot_files


# ---------------------------------------------------------------------
# LaTeX tables
# ---------------------------------------------------------------------

def time_table(results, backends, sizes, calc):
    axis = AXIS_MATH[CALC_AXIS[calc]]
    sizes = [n for n in sizes if any(_find(results, calc, k, n) for k, _, _ in backends)]
    cols = "l" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            axis + " & " + " & ".join(tex_escape(l) for _, l, _ in backends) + r" \\",
            r"\midrule"]
    for n in sizes:
        cells = []
        for key, _label, _iv in backends:
            r = _find(results, calc, key, n)
            if r is None:
                cells.append("--")
            elif not r["ok"]:
                cells.append(r"\textit{fail}")
            else:
                cells.append(fmt_time(r["time"]))
        lines.append(f"{n} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines))


def relative_speed_table(results, backends, sizes, calc):
    """Each backend's time as a multiple of the fastest backend at that
    (calc, size) -- the table that most directly answers "what's worth
    optimizing": 1.00x is the current best, everything else is how much
    slower that mode is on the same problem."""
    axis = AXIS_MATH[CALC_AXIS[calc]]
    sizes = [n for n in sizes if any(_find(results, calc, k, n) and _find(results, calc, k, n)["ok"]
                                      for k, _, _ in backends)]
    if not sizes:
        return None, {}
    cols = "l" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            axis + " & " + " & ".join(tex_escape(l) for _, l, _ in backends) + r" \\",
            r"\midrule"]
    ratios_by_backend = {k: [] for k, _, _ in backends}
    for n in sizes:
        times = {}
        for key, _label, _iv in backends:
            r = _find(results, calc, key, n)
            if r is not None and r["ok"]:
                times[key] = r["time"]
        if not times:
            continue
        best = min(times.values())
        cells = []
        for key, _label, _iv in backends:
            if key in times:
                ratio = times[key] / best if best > 0 else float("inf")
                ratios_by_backend[key].append(ratio)
                cells.append(f"{ratio:.2f}x")
            else:
                cells.append("--")
        lines.append(f"{n} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines)), ratios_by_backend


def accuracy_table(results, backends, refs, sizes, calc, ref_kind):
    """Deviation of each backend's answer from the reference at every
    size that has one: |dE| for a scalar result, max|d| over the whole
    array for a vector one. A scalar reference is printed in its own
    column, so the reader can see what the deviation is relative to."""
    per_size = refs.get(calc, {})
    sizes = [n for n in sizes if n in per_size]
    if not sizes:
        return None
    header = REF_HEADERS.get(ref_kind)
    show_ref = header is not None and np.ndim(per_size[sizes[0]]) == 0
    axis = AXIS_MATH[CALC_AXIS[calc]]
    cols = "l" + ("r" if show_ref else "") + "r" * len(backends)
    head = [axis]
    if show_ref:
        head.append(header)
    head += [f"{tex_escape(l)} $|\\Delta|$" for _, l, _ in backends]
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            " & ".join(head) + r" \\", r"\midrule"]
    for n in sizes:
        ref = per_size[n]
        cells = []
        if show_ref:
            cells.append(f"{float(np.real(ref)):.6f}")
        for key, _label, _iv in backends:
            r = _find(results, calc, key, n)
            dev = _deviation(r["value"], ref) if (r is not None and r["ok"]) else None
            cells.append(fmt_err(dev))
        lines.append(f"{n} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines))


def sector_ratio_table(results, backends, sizes):
    """Sector-mode time divided by dense-mode time at the same N, per
    backend. This is the number the conserved-sector feature is actually
    judged on, and it is strongly size-dependent (block sparsity scales,
    its per-block bookkeeping does not), which is exactly why it belongs
    in a table rather than in a single quoted speedup."""
    have = [(k, l, iv) for k, l, iv in backends
            if any(r["calc"] == "sector" and r["backend"] == k and r["ok"]
                   for r in results)]
    if not have:
        return None
    rows = []
    for n in sizes:
        cells = []
        for key, _l, _iv in have:
            a = _find(results, "sector", key, n)
            b = _find(results, "gs", key, n)
            if a is not None and b is not None and a["ok"] and b["ok"] and b["time"] > 0:
                cells.append("%.2fx" % (a["time"] / b["time"]))
            else:
                cells.append("--")
        if any(c != "--" for c in cells):
            rows.append((n, cells))
    if not rows:
        return None
    cols = "l" + "r" * len(have)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "$N$ & " + " & ".join(tex_escape(l) for _, l, _ in have) + r" \\",
            r"\midrule"]
    for n, cells in rows:
        lines.append(f"{n} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines))


def failures_list(results, backends):
    label_of = {k: l for k, l, _ in backends}
    fails = [r for r in results if not r["ok"]]
    if not fails:
        return None
    lines = [r"\begin{itemize}"]
    for r in fails:
        label = label_of.get(r["backend"], r["backend"])
        lines.append(r"\item %s, %s, %s$=%d$: \texttt{%s}" % (
                CALC_LABELS.get(r["calc"], tex_escape(r["calc"])),
                tex_escape(label), AXIS_MATH[CALC_AXIS.get(r["calc"], "N")],
                r["n"], tex_escape(short_err(r["err"]))))
    lines.append(r"\end{itemize}")
    return "\n".join(lines)


def coverage_table(results, backends, calc_sizes):
    """Which calculation ran on which backend at all -- the answer to
    "what does this benchmark actually cover", including the entries that
    are blank because a backend does not implement that capability."""
    cols = "l" + "c" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "Calculation & " + " & ".join(tex_escape(l) for _, l, _ in backends) + r" \\",
            r"\midrule"]
    for calc in CALC_ORDER:
        if calc not in calc_sizes:
            continue
        cells = []
        for key, _l, _iv in backends:
            rows = [r for r in results if r["calc"] == calc and r["backend"] == key]
            if not rows:
                cells.append("--")
            elif all(r["ok"] for r in rows):
                cells.append(r"\checkmark")
            elif any(r["ok"] for r in rows):
                cells.append("partial")
            else:
                cells.append(r"\textit{fail}")
        lines.append(CALC_LABELS[calc] + " & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines))


# ---------------------------------------------------------------------
# Auto-generated recommendation paragraph
# ---------------------------------------------------------------------

def recommendation_text(results, backends):
    """Geometric-mean relative slowdown per backend, aggregated over
    every (calc, size) point it completed successfully, expressed
    relative to whichever backend was fastest at that same point.
    Highlights the backend most worth optimizing next."""
    ratios_all = {k: [] for k, _, _ in backends}
    for calc in CALC_ORDER:
        ns = sorted(set(r["n"] for r in results if r["calc"] == calc))
        for n in ns:
            times = {r["backend"]: r["time"] for r in results
                    if r["calc"] == calc and r["n"] == n and r["ok"]}
            if not times:
                continue
            best = min(times.values())
            if best <= 0:
                continue
            for key, t in times.items():
                ratios_all[key].append(t / best)

    label_of = {k: l for k, l, _ in backends}
    geo = {}
    for key, ratios in ratios_all.items():
        if ratios:
            geo[key] = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
    if not geo:
        return ("No backend completed enough runs to compare -- see the "
                "failure list below.")

    ordered = sorted(geo.items(), key=lambda kv: kv[1])
    fastest_key, _ = ordered[0]
    parts = [f"Averaged (geometric mean) across every problem size and "
             f"calculation type it completed, \\textbf{{%s}} was the "
             f"fastest backend, normalized to 1.00x at each point it "
             f"was fastest." % tex_escape(label_of[fastest_key])]
    if len(ordered) > 1:
        slowest_key, slowest_ratio = ordered[-1]
        parts.append(
            f"\\textbf{{%s}} was, on average, %.2fx slower than the "
            f"fastest backend on the same problem -- the strongest "
            f"candidate for further optimization work."
            % (tex_escape(label_of[slowest_key]), slowest_ratio))
    parts.append("Ranking (fastest to slowest, geometric-mean relative time): " +
            ", ".join(f"{tex_escape(label_of[k])} ({v:.2f}x)" for k, v in ordered) + ".")
    parts.append("Note that the aggregate mixes calculations that not "
            "every backend implements (the conserved sector and the "
            "infinite chain run on two backends only), so a backend is "
            "ranked on the subset of points it actually completed; the "
            "per-section tables above are the ones to read before acting "
            "on any single number.")
    return " ".join(parts)


# ---------------------------------------------------------------------
# Document assembly
# ---------------------------------------------------------------------

PREAMBLE = r"""\documentclass[11pt]{article}

\usepackage[margin=1in]{geometry}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{adjustbox}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{parskip}

\title{DMRGPY Backend Benchmark Report}
\author{generated by benchmarks/run\_benchmarks.py}
\date{%s}

\sloppy

\begin{document}
\maketitle
"""


def _params_block(params, backends, calc_sizes):
    es = np.asarray(params["es"])
    emin, emax, npts = float(np.min(es)), float(np.max(es)), len(es)
    lines = [
        r"\section*{Run parameters}",
        r"\begin{itemize}",
        r"\item Backends: " + ", ".join(tex_escape(l) for _, l, _ in backends),
        r"\item DMRG: maxm=%d, nsweeps=%d, cutoff=%.1e" % (
                params["maxm"], params["nsweeps"], params["cutoff"]),
        r"\item KPM: n=%d moments; spectral window [%.2f, %.2f] (%d points), "
        r"shared with the real-time correlator (broadening $\delta=%.2f$)" % (
                params["kpm_n"], emin, emax, npts, params["td_delta"]),
        r"\item Real-time evolution: %d steps of $dt=%.3f$" % (
                params["nt"], params["dt"]),
        r"\item Infinite chain: maxiter=%d, etol=%.1e" % (
                params["inf_maxiter"], params["inf_etol"]),
        r"\item Repeats per timing (minimum kept): %d" % params["repeats"],
        r"\end{itemize}",
        r"\subsection*{Sweep ranges}",
        r"\begin{itemize}",
    ]
    for calc in CALC_ORDER:
        if calc not in calc_sizes or not calc_sizes[calc]:
            continue
        axis = AXIS_MATH[CALC_AXIS[calc]]
        lines.append(r"\item %s: %s = %s" % (
                CALC_LABELS[calc], axis,
                ", ".join(str(x) for x in calc_sizes[calc])))
    lines.append(r"\end{itemize}")
    if any(k == "julia_live" for k, _, _ in backends):
        lines.append(
            r"\textit{Note: the Julia backend was warmed up with one "
            r"untimed run before the sweep below, since juliacall "
            r"JIT-compiles each Julia method the first time it is called "
            r"in a process -- an unavoidable, one-time cost of tens of "
            r"seconds that has nothing to do with the algorithm itself. "
            r"All Julia timings below are steady-state, post-JIT numbers.}")
    return "\n".join(lines)


def _calc_section(calc, results, backends, refs, plot_files, use_sizes, ref_kind):
    parts = [r"\section{\texorpdfstring{%s}{%s}}" % (
            CALC_LABELS[calc], pdf_label(CALC_LABELS[calc]))]
    if calc in CALC_NOTES:
        parts.append(r"\textit{%s}" % CALC_NOTES[calc])

    parts += [r"\subsection*{Wall time}",
              time_table(results, backends, use_sizes, calc)]

    if calc in plot_files:
        parts.append(r"\begin{center}\includegraphics[width=0.7\linewidth]{%s}\end{center}"
                % plot_files[calc])

    rel, _ = relative_speed_table(results, backends, use_sizes, calc)
    if rel is not None:
        axis = AXIS_MATH[CALC_AXIS[calc]]
        parts += [r"\subsection*{Relative speed (fastest backend = 1.00x at each %s)}"
                  % axis, rel]

    if calc == "sector":
        ratio = sector_ratio_table(results, backends, use_sizes)
        if ratio is not None:
            parts += [r"\subsection*{Sector vs. dense, same backend "
                      r"(below 1.00x means the sector is faster)}", ratio]

    if ref_kind is not None:
        acc = accuracy_table(results, backends, refs, use_sizes, calc, ref_kind)
        if acc is not None:
            parts += [r"\subsection*{Accuracy vs. reference}", acc]
    return "\n\n".join(parts)


def write_report(results, refs, params, backends, calc_sizes, ref_kinds, outdir,
        generated_at=None):
    """Build plots + report.tex under outdir and return the .tex path.

    `refs` is {calc: {size: reference value}}, `calc_sizes` is
    {calc: [sizes swept]}, and `ref_kinds` is {calc: reference kind}
    naming what that calculation is compared against (see REF_HEADERS);
    a calculation absent from `ref_kinds` gets a timing-only section."""
    os.makedirs(outdir, exist_ok=True)
    plot_files = make_plots(results, backends, outdir)

    sections = []
    for calc in CALC_ORDER:
        if not any(r["calc"] == calc for r in results):
            continue
        use_sizes = calc_sizes.get(calc) or sorted(
                set(r["n"] for r in results if r["calc"] == calc))
        sections.append(_calc_section(calc, results, backends, refs, plot_files,
                use_sizes, ref_kinds.get(calc)))

    summary = [r"\section{Summary}", recommendation_text(results, backends)]
    cov = coverage_table(results, backends, calc_sizes)
    if cov is not None:
        summary += [r"\subsection*{Coverage}", cov]
    fails = failures_list(results, backends)
    if fails is not None:
        summary += [r"\subsection*{Failed runs}", fails]

    doc = [PREAMBLE % (generated_at or ""), _params_block(params, backends, calc_sizes)]
    doc += sections
    doc += summary
    doc.append(r"\end{document}")

    tex_path = os.path.join(outdir, "report.tex")
    with open(tex_path, "w") as f:
        f.write("\n\n".join(doc) + "\n")
    return tex_path


def compile_pdf(tex_path):
    """Compile report.tex with pdflatex (twice, for stable cross-refs),
    run from the report's own directory so \\includegraphics{plot_*.pdf}
    resolves. Returns (ok, message)."""
    outdir = os.path.dirname(os.path.abspath(tex_path))
    base = os.path.basename(tex_path)
    cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", base]
    last = None
    for _ in range(2):
        last = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True)
        if last.returncode != 0:
            break
    pdf_path = os.path.join(outdir, os.path.splitext(base)[0] + ".pdf")
    if last is not None and last.returncode == 0 and os.path.exists(pdf_path):
        return True, f"Compiled {pdf_path}"
    tail = (last.stdout or "")[-3000:] if last is not None else ""
    return False, tail
