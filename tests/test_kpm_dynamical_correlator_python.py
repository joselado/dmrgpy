"""Regression guards for the KPM dynamical correlator on the pure-Python
backend (pyitensor/chain.py's kpm_dynamical_correlator/_kpm_moments*),
written so that a GPU port of that path can be validated against them
(docs/pyitensor_gpu_port_plan.md).

test_dynamical_correlator.py already checks that the KPM peak lands near
the exact excitation gap on all three DMRG backends, and its docstring
explains why a pointwise comparison against ED is the wrong assertion for
KPM: it is a controlled approximation with its own effective resolution,
so raising the moment count does not shrink the pointwise gap to ED. That
leaves untested exactly the properties a *port* has to preserve, which is
what this file covers:

* cross-backend agreement. The compiled v2/v3 sessions and pyitensor run
  the same Chebyshev recursion with the same truncation controls, so their
  spectra agree far more tightly than either agrees with ED -- measured at
  5.4e-12 max absolute deviation on the chain below. That makes the
  compiled backends a hard reference for any change to the Python one,
  moving its arrays onto a device included.
* the zeroth-moment sum rule. Integrating S_AA(omega) over a window wide
  enough to hold the whole spectrum must return mu_0 = <GS|A A|GS>, which
  for A = Sz or Sx on a spin-1/2 site is exactly 1/4 (A^2 = 1/4
  identically) whatever the chain, coupling or bond dimension. Unlike the
  peak position this is sensitive to every moment in the recursion, not
  just the dominant pole.
* run-to-run reproducibility. Two freshly built chains must agree to
  ~1e-14 even though each starts DMRG from its own random MPS -- without
  that, no GPU-vs-CPU diff at 1e-10 would mean anything.

Deliberately small (n=6, kpmmaxm=40): every property here is either exact
or backend-to-backend, so none of it needs a large chain.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain

from _helpers import setup_backend as _setup_backend

N_SITES = 6
DELTA = 0.1
ES = np.linspace(-1.0, 5.0, 400)
# Wide enough to contain this chain's entire Sz/Sx spectral weight, which is
# what makes the sum rule a statement about mu_0 rather than about where the
# window was cut.
ES_WIDE = np.linspace(-8.0, 8.0, 3000)

_COMPILED_VERSIONS = [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2),
        reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
]


def _heisenberg(itensor_version):
    """Uniform S=1/2 Heisenberg chain, the same model the rest of the
    dynamical-correlator tests use."""
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(N_SITES)])
    _setup_backend(sc, itensor_version)
    h = 0
    for i in range(N_SITES - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1]
        h = h + sc.Sy[i] * sc.Sy[i + 1]
        h = h + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm = 30
    sc.kpmmaxm = 40
    sc.nsweeps = 8
    return sc


def _kpm(sc, op, es=ES):
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                        name=(op, op), es=es, delta=DELTA)
    return np.array(x), np.array(y).real


@pytest.mark.parametrize("version", _COMPILED_VERSIONS)
def test_python_kpm_spectrum_matches_compiled_backend(version):
    """pyitensor's Chebyshev recursion against the compiled one, pointwise.

    Same algorithm and same kpmmaxm/cutoff, so this is a far sharper probe
    than the existing peak-position test: measured deviation is ~5e-12,
    asserted at 1e-8 to leave room for each backend's own iterative
    convergence to the ground state the recursion starts from.
    """
    sc_py = _heisenberg("python")
    sc_cpp = _heisenberg(version)
    xp, yp = _kpm(sc_py, sc_py.Sz[0])
    xc, yc = _kpm(sc_cpp, sc_cpp.Sz[0])
    assert np.allclose(xp, xc)
    assert np.max(np.abs(yp - yc)) < 1e-8


@pytest.mark.parametrize("opname", ["Sz", "Sx"])
def test_python_kpm_obeys_the_zeroth_moment_sum_rule(opname):
    """int S_AA(omega) domega == mu_0 == <GS|A A|GS> == 1/4 exactly, for
    A = Sz or Sx on a spin-1/2 site.

    Sensitive to every moment rather than just the dominant pole, so a
    port that misnormalizes, drops or reorders moments fails here even
    when the peak still lands in the right place.
    """
    sc = _heisenberg("python")
    x, y = _kpm(sc, getattr(sc, opname)[0], es=ES_WIDE)
    assert np.trapezoid(y, x) == pytest.approx(0.25, abs=1e-4)


def test_python_kpm_is_reproducible_across_fresh_chains():
    """Two independently built chains, each starting DMRG from its own
    random MPS, must give the same KPM spectrum. Measured at ~2e-14; the
    1e-10 bound is what makes a later GPU-vs-CPU comparison at that level
    meaningful rather than lost in run-to-run noise."""
    sc1, sc2 = _heisenberg("python"), _heisenberg("python")
    _, y1 = _kpm(sc1, sc1.Sz[0])
    _, y2 = _kpm(sc2, sc2.Sz[0])
    assert np.max(np.abs(y1 - y2)) < 1e-10


def test_python_kpm_spectrum_is_nonnegative_and_peaked_in_window():
    """S_AA(omega) >= 0 is a property of the spectral function itself, and
    KPM's damped reconstruction preserves it up to the kernel's own small
    negative ringing. Cheap guard against a port that breaks the
    moment-to-spectrum reconstruction -- losing the real part, or
    mismatching the energy rescaling -- in a way the sum rule's single
    number could average out."""
    sc = _heisenberg("python")
    x, y = _kpm(sc, sc.Sz[0])
    assert y.min() > -1e-3 * y.max()
    assert ES[0] < x[np.argmax(y)] < ES[-1]
