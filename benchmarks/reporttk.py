"""LaTeX + plot report generation for run_benchmarks.py.

Kept separate from run_benchmarks.py (which owns problem construction,
timing and CLI) so that a report can be regenerated from a previously
saved results.json (e.g. after tweaking a table/plot) without rerunning
the DMRG sweep -- see run_benchmarks.py's --from-json option.
"""
import os
import math
import subprocess

import matplotlib
matplotlib.use("Agg")  # no display needed to render report plots
import matplotlib.pyplot as plt
import numpy as np

CALC_LABELS = {"gs": "Ground state energy", "static": "Static correlator",
        "dynamic": "Dynamical correlator (KPM)"}
CALC_ORDER = ["gs", "static", "dynamic"]


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


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

def make_plots(results, backends, outdir):
    """One time-vs-N plot per calculation type actually present in
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
        ax.set_xlabel("chain length $N$")
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
    sizes = [n for n in sizes if any(_find(results, calc, k, n) for k, _, _ in backends)]
    cols = "l" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "$N$ & " + " & ".join(tex_escape(l) for _, l, _ in backends) + r" \\",
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
    (calc, N) -- the table that most directly answers "what's worth
    optimizing": 1.00x is the current best, everything else is how much
    slower that mode is at the same problem size."""
    sizes = [n for n in sizes if any(_find(results, calc, k, n) and _find(results, calc, k, n)["ok"]
                                      for k, _, _ in backends)]
    if not sizes:
        return None, {}
    cols = "l" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "$N$ & " + " & ".join(tex_escape(l) for _, l, _ in backends) + r" \\",
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


def gs_accuracy_table(results, backends, ed_refs, sizes):
    sizes = [n for n in sizes if ed_refs.get(n, {}).get("ok")]
    if not sizes:
        return None
    cols = "lr" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "$N$ & ED energy & " + " & ".join(f"{tex_escape(l)} $|\\Delta E|$" for _, l, _ in backends) + r" \\",
            r"\midrule"]
    for n in sizes:
        e_ed = ed_refs[n]["energy"]
        cells = [f"{e_ed:.6f}"]
        for key, _label, _iv in backends:
            r = _find(results, "gs", key, n)
            if r is not None and r["ok"]:
                cells.append(fmt_err(abs(r["energy"] - e_ed)))
            else:
                cells.append("--")
        lines.append(f"{n} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return _boxed("\n".join(lines))


def static_accuracy_table(results, backends, ed_refs, sizes):
    sizes = [n for n in sizes if ed_refs.get(n, {}).get("ok") and ed_refs[n].get("correlator") is not None]
    if not sizes:
        return None
    cols = "l" + "r" * len(backends)
    lines = [r"\begin{tabular}{%s}" % cols, r"\toprule",
            "$N$ & " + " & ".join(f"{tex_escape(l)} $\\max|\\Delta|$" for _, l, _ in backends) + r" \\",
            r"\midrule"]
    for n in sizes:
        corr_ed = ed_refs[n]["correlator"]
        cells = []
        for key, _label, _iv in backends:
            r = _find(results, "static", key, n)
            if r is not None and r["ok"]:
                diffs = [abs(a - b) for a, b in zip(r["value"], corr_ed)]
                cells.append(fmt_err(max(diffs)))
            else:
                cells.append("--")
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
        lines.append(r"\item %s, %s, $N=%d$: \texttt{%s}" % (
                tex_escape(CALC_LABELS.get(r["calc"], r["calc"])),
                tex_escape(label), r["n"], tex_escape(short_err(r["err"]))))
    lines.append(r"\end{itemize}")
    return "\n".join(lines)


# ---------------------------------------------------------------------
# Auto-generated recommendation paragraph
# ---------------------------------------------------------------------

def recommendation_text(results, backends):
    """Geometric-mean relative slowdown per backend, aggregated over
    every (calc, N) point it completed successfully, expressed relative
    to whichever backend was fastest at that same point. Highlights the
    backend most worth optimizing next."""
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
    parts = [f"Averaged (geometric mean) across every chain length and "
             f"calculation type it completed, \\textbf{{%s}} was the "
             f"fastest backend, normalized to 1.00x at each point it "
             f"was fastest." % tex_escape(label_of[fastest_key])]
    if len(ordered) > 1:
        slowest_key, slowest_ratio = ordered[-1]
        parts.append(
            f"\\textbf{{%s}} was, on average, %.2fx slower than the "
            f"fastest backend at the same problem size -- the strongest "
            f"candidate for further optimization work."
            % (tex_escape(label_of[slowest_key]), slowest_ratio))
    parts.append("Ranking (fastest to slowest, geometric-mean relative time): " +
            ", ".join(f"{tex_escape(label_of[k])} ({v:.2f}x)" for k, v in ordered) + ".")
    return " ".join(parts)


# ---------------------------------------------------------------------
# Document assembly
# ---------------------------------------------------------------------

PREAMBLE = r"""\documentclass[11pt]{article}

\usepackage[margin=1in]{geometry}
\usepackage{amsmath}
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


def _params_block(params, backends, sizes, dyn_sizes):
    es = params["kpm_es"]
    emin, emax, npts = float(np.min(es)), float(np.max(es)), len(es)
    lines = [
        r"\section*{Run parameters}",
        r"\begin{itemize}",
        r"\item Backends: " + ", ".join(tex_escape(l) for _, l, _ in backends),
        r"\item Chain lengths $N$: " + ", ".join(str(n) for n in sizes),
        r"\item Dynamical-correlator subset: " + ", ".join(str(n) for n in dyn_sizes),
        r"\item DMRG: maxm=%d, nsweeps=%d, cutoff=%.1e" % (
                params["maxm"], params["nsweeps"], params["cutoff"]),
        r"\item KPM: n=%d moments, energy window [%.2f, %.2f] (%d points)" % (
                params["kpm_n"], emin, emax, npts),
        r"\item Repeats per timing (minimum kept): %d" % params["repeats"],
        r"\end{itemize}",
    ]
    if any(k == "julia_live" for k, _, _ in backends):
        lines.append(
            r"\textit{Note: the Julia backend was warmed up with one "
            r"untimed run before the sweep below, since juliacall "
            r"JIT-compiles each Julia method the first time it is called "
            r"in a process -- an unavoidable, one-time cost of tens of "
            r"seconds that has nothing to do with the algorithm itself. "
            r"All Julia timings below are steady-state, post-JIT numbers.}")
    return "\n".join(lines)


def _calc_section(calc, results, backends, sizes, ed_refs, plot_files, use_sizes):
    title = CALC_LABELS[calc]
    parts = [r"\section{%s}" % title]

    tbl = time_table(results, backends, use_sizes, calc)
    parts += [r"\subsection*{Wall time}", tbl]

    if calc in plot_files:
        parts.append(r"\begin{center}\includegraphics[width=0.7\linewidth]{%s}\end{center}"
                % plot_files[calc])

    rel, _ = relative_speed_table(results, backends, use_sizes, calc)
    if rel is not None:
        parts += [r"\subsection*{Relative speed (fastest backend = 1.00x at each $N$)}", rel]

    if calc == "gs":
        acc = gs_accuracy_table(results, backends, ed_refs, use_sizes)
        if acc is not None:
            parts += [r"\subsection*{Accuracy vs. exact diagonalization}", acc]
    elif calc == "static":
        acc = static_accuracy_table(results, backends, ed_refs, use_sizes)
        if acc is not None:
            parts += [r"\subsection*{Accuracy vs. exact diagonalization}", acc]
    elif calc == "dynamic":
        parts.append(
            r"\textit{No independent ED cross-check is run here (a KPM "
            r"spectral function is more involved to compare pointwise "
            r"against an exact spectrum); physics correctness for this "
            r"path is covered by \texttt{tests/test\_dynamical\_correlator.py}. "
            r"This section is timing-only.}")
    return "\n\n".join(parts)


def write_report(results, ed_refs, params, backends, sizes, dyn_sizes, outdir, generated_at=None):
    """Build plots + report.tex under outdir and return the .tex path."""
    os.makedirs(outdir, exist_ok=True)
    plot_files = make_plots(results, backends, outdir)

    sections = []
    for calc in CALC_ORDER:
        use_sizes = dyn_sizes if calc == "dynamic" else sizes
        if not any(r["calc"] == calc for r in results):
            continue
        sections.append(_calc_section(calc, results, backends, sizes, ed_refs, plot_files, use_sizes))

    summary = [r"\section{Summary}", recommendation_text(results, backends)]
    fails = failures_list(results, backends)
    if fails is not None:
        summary += [r"\subsection*{Failed runs}", fails]

    doc = [PREAMBLE % (generated_at or ""), _params_block(params, backends, sizes, dyn_sizes)]
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
