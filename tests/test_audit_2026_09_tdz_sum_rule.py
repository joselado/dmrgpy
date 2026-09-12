"""Regression for the 2026-09 audit's finding #7, the one fix in that audit
that lives in C++ (`mpscpp3/chain_session.h::Chain::tdvp_step`).

`tdvp_step` used to pass `{"DoNormalize",true}` unconditionally. exp(-i H dt)
preserves the norm only for REAL dt; along the complex-time contour that
`tdz.py`'s submode="TDZ" walks, the decay of ||psi|| is the physics, and
forcing unit norm after every step deletes it. The symptom was 36% too much
spectral weight on `itensor_version=3` (the audit recorded 0.340676 against an
exact sum rule of 0.25) while `itensor_version="python"` -- whose own stepper
never normalized -- was right.

What this pins is the SUM RULE, not a golden number, and that choice is the
point: `tests/test_dynamical_correlator.py::test_tdz_dynamical_correlator_peak_matches_exact_gap`
already existed and did not catch this, because it checks the peak POSITION
and a uniform rescaling does not move a peak. The audit's own Status note
recorded this missing assertion as outstanding; this file is it.

    integral dw C_AB(w) = <GS|A B|GS>,  which for A = B = Sz_0 on S=1/2 is 1/4.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, cppext


def _chain(backend, n=6):
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    if backend == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(version=backend)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + 0.3 * sc.Sz[i]  # breaks the Sz symmetry, as the audit's chain does
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 12
    return sc


@pytest.mark.parametrize("backend", [3, "python"])
def test_tdz_integrates_to_the_exact_sum_rule(backend):
    """<Sz0 Sz0> = 1/4 exactly for S=1/2, so the TDZ spectral weight must
    integrate to 0.25 regardless of backend. Pre-fix, itensor_version=3
    returned 0.3407 here (1.36x) and "python" returned 0.2495."""
    if backend == 3 and not cppext.available(3):
        pytest.skip("ITensor v3 extension not compiled")
    sc = _chain(backend)
    es = np.linspace(-20, 20, 4001)
    x, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), submode="TDZ",
                                       es=es, delta=0.1)
    assert np.trapezoid(y.real, x) == pytest.approx(0.25, abs=2e-3)


def test_tdz_agrees_between_the_compiled_and_pure_python_backends():
    """The two implement the same contour by different code (ITensor's TDVP
    vs pyitensor's). Before the fix they disagreed by 36%; agreement is the
    sharper statement, since nothing forces it but both being right."""
    if not cppext.available(3):
        pytest.skip("ITensor v3 extension not compiled")
    es = np.linspace(-20, 20, 4001)
    out = []
    for backend in (3, "python"):
        sc = _chain(backend)
        out.append(sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                               submode="TDZ", es=es, delta=0.1)[1])
    assert np.max(np.abs(out[0] - out[1])) < 1e-6


def test_real_time_evolution_is_unaffected_by_the_normalization_change():
    """The other half of the fix: DoNormalize stays on for a REAL dt, so the
    real-time callers (which restore the input norm themselves) are untouched.
    A quench trajectory must still track ED."""
    if not cppext.available(3):
        pytest.skip("ITensor v3 extension not compiled")
    from dmrgpy import timedependent
    ts_ref = cs_ref = None
    for backend in (3, "ED"):
        sc = _chain(backend if backend != "ED" else 3, n=5)
        if backend == "ED":
            sc.mode = "ED"
        ts, cs = timedependent.evolution_DC(sc, mode=sc.get_mode(),
                                            name=(sc.Sz[0], sc.Sz[0]), nt=6, dt=0.2)
        if ts_ref is None:
            ts_ref, cs_ref = ts, cs
        else:
            assert np.max(np.abs(np.array(cs) - np.array(cs_ref))) < 1e-4
