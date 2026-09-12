"""Regression tests for the examples/ findings of the 2026-09 audit.

Findings #24, #33 and #34 of `docs/audit_2026_09_hole_hunt.md`. All three
are about `examples/*/*/main.py` scripts rather than about library code,
so most of what has to be pinned is a property of the script itself:

* #24 -- `boson_models/v2_VS_v3_boson` is an assert-carrying regression
  script whose assert fired on roughly 40% of clean-tree runs, because it
  left DMRG at the library defaults (nsweeps=15, maxm=30) and both
  compiled backends start from a random MPS. The fix is a pinned sweep
  schedule inside its own `get_energy()`. What has to stay pinned is that
  schedule, so the test reads it back out of the source; the *mechanism*
  (converged DMRG at those settings really does reach ED) is checked
  separately on the single hardest U point rather than by re-running the
  ~25s script.
* #33 -- `readme_examples/energy_VS_length` was an empty stub: no
  gs_energy() call, no length loop, no output at all.
* #34 -- `groundstate/energy_fluctuation` and
  `utilities/multioperator_density` each computed a full sequence,
  printed it, and never plotted it, against CLAUDE.md's "examples should
  plot" rule.

The "does it plot" tests mirror the audit's own repro, which was a grep
for a script that imports matplotlib and never calls a drawing method,
and additionally execute the two cheapest of the scripts to confirm the
drawing calls actually produce populated axes rather than merely being
present in the text.
"""

import os
import re
import runpy

import numpy as np
import pytest

from dmrgpy import bosonchain


EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "examples")

# the audit's own "is anything ever drawn" regex, verbatim
DRAW_CALL = re.compile(
    r"\.(plot|scatter|imshow|errorbar|semilogy|semilogx|loglog|bar|barh|"
    r"fill_between|pcolormesh|hist|step|contourf|axhline|axvline|stem|"
    r"plot_surface)\(")


def example_source(*parts):
    with open(os.path.join(EXAMPLES, *parts, "main.py")) as f:
        return f.read()


def run_example(monkeypatch, *parts):
    """Execute an example's main.py with its own directory as cwd (its
    sys.path preamble is relative to it) and with plt.show() disabled, and
    return the axes it drew on."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    plt.close("all")
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    monkeypatch.chdir(os.path.join(EXAMPLES, *parts))
    runpy.run_path("main.py", run_name="__main__")
    axes = plt.gcf().axes
    return axes


# --------------------------------------------------- #24, the flaky guard

def test_boson_v2_VS_v3_example_pins_its_sweep_schedule():
    """The script must not go back to the library defaults.

    At nsweeps=15/maxm=30 its U-sweep assert failed 3 of 6 clean-tree runs
    (the audit's reviewer measured 2 of 5, the hunter 5 of 7), with the
    outlier landing on v2 in some runs and on v3 in others -- shared DMRG
    under-convergence from a random start, not a backend divergence. At
    nsweeps=80/maxm=100 it passed 16 of 16, with a worst cross-backend
    disagreement of 1.7e-7 over those runs -- which is what the tightened
    tolerance rests on. maxm must stay at or above 4**3=64, the exact MPS
    bond dimension of a 6-site 4-level boson chain, so nothing is
    truncated at all.
    """
    src = example_source("boson_models", "v2_VS_v3_boson")
    nsweeps = int(re.search(r"bc\.nsweeps\s*=\s*(\d+)", src).group(1))
    maxm = int(re.search(r"bc\.maxm\s*=\s*(\d+)", src).group(1))
    tol = float(re.search(r"^tol\s*=\s*([0-9.e-]+)", src, re.M).group(1))
    assert nsweeps >= 80
    assert maxm >= 4 ** 3
    # a guard at 1e-2 could not distinguish a real v2/v3 divergence from
    # its own background failure rate, which is the whole finding
    assert tol <= 1e-5


def test_boson_dmrg_reaches_ed_at_the_examples_pinned_schedule():
    """The mechanism behind the fix, on the hardest point of the sweep.

    U=0.2 is where the under-convergence tail lived (every clean-tree
    failure the audit recorded was at U=0.1, 0.2, 0.4 or 0.5, and U=0.2
    was the worst in repeated measurement). Only itensor_version=3 is run
    here; the example itself covers v2/ED/python as well, and the finding
    was explicitly *not* backend-specific.

    The tolerance deliberately matches the example's own tol=1e-5 rather
    than being tighter: the worst deviation observed over 16 full runs of
    the example was 1.7e-7, so 1e-6 would leave only 6x of margin on a
    quantity whose whole problem is a random-start tail.
    """
    n, U = 6, 0.2

    def energy(mode, itensor_version=3):
        bc = bosonchain.Bosonic_Chain(n)
        bc.setup_cpp(itensor_version)
        bc.nsweeps, bc.maxm = 80, 100  # the example's pinned schedule
        np.random.seed(11)
        h = 0
        for i in range(n - 1):
            h = h + np.random.random() * (bc.Adag[i] * bc.A[i + 1]
                                          + bc.Adag[i + 1] * bc.A[i])
        for i in range(n):
            den = bc.Adag[i] * bc.A[i]
            h = h + U * den * den
        bc.set_hamiltonian(h)
        return bc.gs_energy(mode=mode)

    assert energy("DMRG") == pytest.approx(energy("ED"), abs=1e-5)


# ------------------------------------------------------ #33, the empty stub

def test_energy_VS_length_example_sweeps_length_and_plots(monkeypatch):
    """It used to build one fixed 30-site chain and stop: no gs_energy()
    call, no loop, no print, no plot, exit 0 with no output at all."""
    src = example_source("readme_examples", "energy_VS_length")
    assert "gs_energy(" in src
    assert DRAW_CALL.search(src) is not None

    axes = run_example(monkeypatch, "readme_examples", "energy_VS_length")
    # the swept curve, plus the Bethe-ansatz axhline reference
    lines = [l for l in axes[0].lines if len(l.get_xdata()) > 2]
    assert len(lines) >= 1
    ns, es = lines[0].get_xdata(), lines[0].get_ydata()
    assert len(ns) >= 5  # an actual sweep, not a single point
    # E/n is plotted, and for an open Heisenberg chain it falls towards
    # the Bethe-ansatz density 1/4-ln2 = -0.4431 from above with n
    assert np.all(np.diff(es) < 0)
    assert es[-1] < 0.25 - np.log(2.) + 0.05


# ------------------------------------------ #34, sequences that never plotted

@pytest.mark.parametrize("parts", [
    ("groundstate", "energy_fluctuation"),
    ("utilities", "multioperator_density"),
])
def test_sequence_examples_draw_what_they_compute(parts):
    """Both scripts imported pyplot at the top and never touched it,
    which is exactly the signature CLAUDE.md describes."""
    src = example_source(*parts)
    assert "matplotlib" in src
    assert DRAW_CALL.search(src) is not None


def test_multioperator_density_plots_both_backends(monkeypatch):
    """The site-resolved density profile is computed twice, once under
    mode='DMRG' and once under mode='ED'; plotting them overlaid is what
    makes the backend agreement visible instead of leaving the reader to
    diff two printed lists.

    The example seeds its random hopping matrix, so this comparison is
    deterministic: an unseeded draw can put a single-particle level at
    zero and leave the many-body ground state nearly degenerate, at which
    point the two solvers legitimately return different density profiles
    at the same energy (measured: 0.14 apart at seed 55, and more sweeps
    do not cure it). See the comment in the example itself."""
    axes = run_example(monkeypatch, "utilities", "multioperator_density")
    lines = axes[0].lines
    assert len(lines) == 2  # DMRG and ED
    dmrg, ed = lines[0].get_ydata(), lines[1].get_ydata()
    assert len(dmrg) == 6 and len(ed) == 6  # the chain has 6 sites
    assert dmrg == pytest.approx(ed, abs=1e-6)
