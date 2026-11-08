"""Regression tests for the `fluctuation` cluster of the third 2026-09-24
hole hunt (docs/audit_2026_09_24c_hole_hunt.md, findings 6 and 7).

   6. gs_energy(maxde=...) returned the energy of its first, unrefined
      solve (-3.279373 against a stored -3.374932 on the 8-site chain at
      maxm=3), and the next dynamical correlator re-solved at the original
      maxm, discarding the refinement, because the send cache stayed keyed
      on the doubled maxm. It returns the refined energy now, and the
      refined state survives.
   7. gs_energy_fluctuation() computed <H^2>-<H>^2 with H|psi> truncated to
      the chain's maxm, so the ordinary solve-then-measure workflow
      under-reported the variance 10 to 51 times and a state wider than
      maxm over-reported it at order one, and maxde= stopped early. It is
      ||(H-<H>)|psi>||^2 with an uncapped application now, on every mode.

The anchor below full bond dimension is vev(h*h), an exact H^2 MPO with no
MPS truncation, which is accurate wherever the variance is far above its
roundoff floor; at full bond dimension it is ED.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import cppext, spinchain


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


from _helpers import julia_live_param

BACKENDS = [_backend("python"), _backend(3), _backend(2)]
N = 8
E0_EXACT = -3.37493260


def _chain(version, maxm):
    sc = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=version)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(N-1))
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, 10
    return sc, h


def _exact_h2(sc, h):
    e = sc.vev(h).real
    return np.sqrt(abs(sc.vev(h*h).real - e*e))


# ------------------------------------------------------------ finding 7

@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("maxm", [3, 6])
def test_fluctuation_is_the_variance_of_the_state(version, maxm):
    """Solved and measured at the same maxm, the old route under-reported:
    6.725e-03 against 3.436e-01 on "python" at maxm=3."""
    sc, h = _chain(version, maxm)
    np.random.seed(0)
    sc.gs_energy()
    true = _exact_h2(sc, h)
    assert true > 5e-2
    assert sc.gs_energy_fluctuation() == pytest.approx(true, rel=1e-6)


@pytest.mark.parametrize("version", BACKENDS)
def test_fluctuation_of_a_state_wider_than_maxm(version):
    """An essentially exact state (solved at full bond dimension 16)
    measured at maxm=3 reported 0.98 ("python", v2) and 1.64 (v3)."""
    sc, h = _chain(version, 16)
    np.random.seed(0)
    sc.gs_energy()
    wf = sc.get_gs()
    sc.maxm = 3
    assert sc.gs_energy_fluctuation(wf=wf) < 1e-9


def test_fluctuation_on_ed_is_the_same_quantity():
    """<(H-<H>)^2> on ED as well, rather than <H^2>-<H>^2, whose roundoff
    floor (1.1e-7 here) sits far above a converged state's fluctuation."""
    sc, h = _chain("python", 16)
    assert sc.gs_energy_fluctuation(mode="ED") < 1e-10


def test_fluctuation_rejects_unknown_keywords():
    sc, h = _chain("python", 8)
    with pytest.raises(TypeError, match="maxmm"):
        sc.gs_energy_fluctuation(maxmm=4)


# ------------------------------------------------------------ finding 6

@pytest.mark.parametrize("version", BACKENDS)
def test_maxde_returns_the_refined_energy_and_keeps_it(version):
    """Returned value, stored energy and ED agree, and a KPM correlator
    afterwards measures the refined state rather than re-solving at
    maxm=3 (overlap 0.856 on "python", 0.765 on v2/v3, before)."""
    sc, h = _chain(version, 3)
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ret = sc.gs_energy(maxde=1e-4)
    assert sc.maxm == 3
    assert ret == pytest.approx(sc.gs_energy(), abs=1e-12)
    assert ret == pytest.approx(E0_EXACT, abs=1e-6)
    wf = sc.get_gs()
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    assert abs(wf.dot(sc.get_gs()))**2 == pytest.approx(1.0, abs=1e-10)
    assert sc.gs_energy() == pytest.approx(ret, abs=1e-12)


@pytest.mark.parametrize("version", BACKENDS)
def test_maxde_meets_its_tolerance(version):
    """maxde is a fluctuation per site; the old loop read the truncated
    route, stopped early and returned states 2.3 times over it."""
    sc, h = _chain(version, 3)
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.gs_energy(maxde=1e-3)
    assert _exact_h2(sc, h)/N <= 1e-3


def test_fluctuation_runs_on_julia_live():
    """julia_live's vev took no npow, so the old route raised TypeError
    there; it is ||(H-<H>)|psi>||^2 with an uncapped application now, as
    on the session backends."""
    from _helpers import julia_available
    ok, reason = julia_available()
    if not ok:
        pytest.skip(reason)
    from dmrgpy import multioperator
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="julia_live")
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(5))
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 8, 10
    sc.gs_energy()
    assert sc.gs_energy_fluctuation() < 1e-6 # full bond dimension: exact


@pytest.mark.parametrize("version", [julia_live_param()])
def test_julia_fluctuation_of_a_state_that_is_not_an_eigenstate(version):
    """On julia_live the variance of a state that is not an eigenstate is
    the ED one: the state is (1+0.5 Sx0)|GS> built on both sides from each
    one's own ground state, which differ only by a phase. And on a 10-site
    chain at maxm=4 the uncapped route equals vev(h*h)-<h>^2, both exact,
    for the same state."""
    from dmrgpy import multioperator
    n = 6
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    ed = spinchain.Spin_Chain(["S=1/2"]*n)
    for c in (sc, ed):
        c.set_hamiltonian(sum(c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1]
                              + c.Sz[i]*c.Sz[i+1] for i in range(n-1)))
        c.maxm, c.nsweeps = 30, 10
    np.random.seed(1)
    g = sc.get_gs()
    x = g + 0.5*(sc.Sx[0]*g)
    gE = ed.get_gs(mode="ED")
    xE = gE + 0.5*(ed.Sx[0]*gE)
    xE = xE*(1.0/np.sqrt(np.real(xE.dot(xE))))
    hE = ed.hamiltonian
    eE = np.real(xE.dot(hE*xE))
    var_ed = np.real(xE.dot(hE*(hE*xE))) - eE**2
    assert var_ed > 1e-2
    assert sc.gs_energy_fluctuation(wf=x)**2 == pytest.approx(var_ed, rel=1e-8)
    # 10 sites, maxm=4: a state far from converged, two exact routes
    m = 10
    sk = spinchain.Spin_Chain(["S=1/2"]*m, itensor_version=version)
    h = sum(sk.Sx[i]*sk.Sx[i+1] + sk.Sy[i]*sk.Sy[i+1] + sk.Sz[i]*sk.Sz[i+1]
            for i in range(m-1))
    sk.set_hamiltonian(h)
    sk.maxm, sk.nsweeps = 4, 10
    np.random.seed(2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sk.gs_energy()
    f = sk.gs_energy_fluctuation()
    assert f > 1e-2
    assert f == pytest.approx(_exact_h2(sk, h), rel=1e-6)
