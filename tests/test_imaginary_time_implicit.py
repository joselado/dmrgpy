"""Imaginary-time evolution by implicit steps on a logarithmic grid
(arXiv:2606.02930): imaginarytime.py, pyitensor/implicit.py,
mpscpp3's Chain::implicit_step and mpsjulialive/implicit.jl, and the
METTS route built on them."""
from math import factorial

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from dmrgpy import spinchain, imaginarytime
from dmrgpy import cppext
from _helpers import julia_live_param

BACKENDS = ["python"] + ([3] if cppext.available(3) else [])
# every julia_live test carries the julia_live id, so -k "not julia_live"
# skips them; they reuse the 8-site chain of _neel_chain
JL = pytest.mark.parametrize("backend", [julia_live_param()])


def test_pade_factors_are_the_diagonal_approximant():
    """The product of the factors of order 2s is the (s,s) Pade
    approximant of exp(-z): it matches exp(-z) to O(z^(2s+1)) and has
    modulus at most one on z >= 0."""
    assert imaginarytime.gauss_legendre_roots(1) == pytest.approx([-2.0])
    for order in (2, 4, 6):
        z = 1e-2
        R = np.prod([(1-b*z)/(1+b*z) for b in imaginarytime.step_factors(1.0, order)])
        assert abs(R-np.exp(-z)) < 10*z**(order+1)/factorial(order)
        for z in np.geomspace(1e-3, 1e4, 50):
            R = np.prod([(1-b*z)/(1+b*z) for b in imaginarytime.step_factors(1.0, order)])
            assert abs(R) <= 1+1e-12
            assert abs(R.imag) < 1e-12
    with pytest.raises(ValueError):
        imaginarytime.step_factors(1.0, 3)


def test_log_grid():
    steps = imaginarytime.log_grid(100., 0.01, 4)
    assert steps.sum() == pytest.approx(100.)
    assert steps[0] <= 0.01/4 + 1e-15
    assert len(steps) == 4*(int(np.ceil(np.log2(1e4)))+1)
    assert steps[-1] == pytest.approx(50./4)
    assert len(imaginarytime.log_grid(0.005, 0.01, 3)) == 3


def _neel_chain(backend, n=8, cutoff=0.0, maxm=64):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=backend,
                              maxm=maxm, nsweeps=10, cutoff=cutoff)
    sc.set_hamiltonian(sum((-1)**i*sc.Sz[i] for i in range(n)))
    wf, wfed = sc.get_gs(), sc.get_gs(mode="ED")
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    return sc, h, wf, wfed


def _rational(sc, h, v, taus, order, shift):
    """The exact rational propagator of the steps, by sparse solves."""
    Hm = sp.csc_matrix(sc.get_ED_obj().MO2matrix(h))
    Id = sp.identity(Hm.shape[0], format="csc")
    v = v/np.linalg.norm(v)
    out = [np.vdot(v, Hm@v).real]
    for dt in np.diff(taus):
        for b in imaginarytime.step_factors(dt, order):
            v = spsolve(Id+b*(Hm-shift*Id), (Id-b*(Hm-shift*Id))@v)
        v = v/np.linalg.norm(v)
        out.append(np.vdot(v, Hm@v).real)
    return np.array(out)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("order", [2, 4])
def test_sweeps_solve_the_rational_steps_exactly(backend, order):
    """With no truncation and tight local solves, the MPS steps are the
    exact rational propagator: what is left is the time step alone."""
    sc, h, wf, wfed = _neel_chain(backend)
    taus, e = sc.evolve_imaginary_time_and_measure(h, wf, 5., nsteps=2,
            order=order, nsweeps=3, tol=1e-13, krylov=30, maxiter=1000)
    er = _rational(sc, h, np.array(wfed.v, dtype=complex), taus, order,
                   sc.gs_energy())
    assert np.max(np.abs(e.real-er)) < 1e-10


@pytest.mark.parametrize("backend", BACKENDS)
def test_trapezoid_is_second_order_and_reaches_the_ground_state(backend):
    sc, h, wf, wfed = _neel_chain(backend, cutoff=1e-12)
    e0 = sc.gs_energy(mode="ED")
    errs = []
    for n in (2, 4):
        taus, e = sc.evolve_imaginary_time_and_measure(h, wf, 40., nsteps=n)
        _, eex = sc.evolve_imaginary_time_and_measure(h, wfed, 40., nsteps=n,
                                                      mode="ED")
        errs.append(np.max(np.abs(e-eex)))
        assert e[-1].real-e0 < 1e-8
    assert 3. < errs[0]/errs[1] < 5.


def test_v3_agrees_with_python():
    if 3 not in BACKENDS: pytest.skip("itensor_version=3 is not compiled")
    out = []
    for backend in ("python", 3):
        sc, h, wf, wfed = _neel_chain(backend, cutoff=1e-12)
        out.append(sc.evolve_imaginary_time_and_measure(
            [h, sc.Sz[0]*sc.Sz[1]], wf, 10., nsteps=2, order=4)[1])
    assert np.max(np.abs(out[0]-out[1])) < 1e-7


@JL
def test_julia_live_solves_the_rational_steps_exactly(backend):
    """test_sweeps_solve_the_rational_steps_exactly on julia_live, at
    order=4: complex betas and the mixed environments <x|H|psi>"""
    sc, h, wf, wfed = _neel_chain(backend)
    taus, e = sc.evolve_imaginary_time_and_measure(h, wf, 5., nsteps=2,
            order=4, nsweeps=3, tol=1e-13, krylov=30, maxiter=1000)
    er = _rational(sc, h, np.array(wfed.v, dtype=complex), taus, 4,
                   sc.gs_energy())
    assert np.max(np.abs(e.real-er)) < 1e-10


@JL
def test_julia_live_agrees_with_python(backend):
    """test_v3_agrees_with_python on julia_live. The shift is the
    ground-state energy, as by default, solved once on "python" and handed
    to both, since a julia_live solve varies run to run by solver noise"""
    out, shift = [], None
    for b in ("python", backend):
        sc, h, wf, wfed = _neel_chain(b, cutoff=1e-12)
        if shift is None: shift = float(np.real(sc.gs_energy()))
        out.append(sc.evolve_imaginary_time_and_measure(
            [h, sc.Sz[0]*sc.Sz[1]], wf, 10., nsteps=2, order=4,
            shift=shift)[1])
    assert np.max(np.abs(out[0]-out[1])) < 1e-7


@JL
def test_julia_live_keeps_tdvp_in_metts(backend):
    """METTS on julia_live stays on TDVP: the automatic plan is None and
    itime="implicit" is refused, with no solve on the chain"""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=backend)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(5))
    sc.set_hamiltonian(h)
    assert imaginarytime.metts_plan(sc, 0.05) is None
    with pytest.raises(NotImplementedError):
        imaginarytime.metts_plan(sc, 0.05, itime="implicit")
    with pytest.raises(NotImplementedError):
        sc.metts_vev(h, 0.05, itime="implicit", nsamples=2, nwarmup=1, seed=1)


@pytest.mark.parametrize("backend", BACKENDS)
def test_state_and_running_shift(backend):
    sc, h, wf, wfed = _neel_chain(backend, cutoff=1e-12)
    e0 = sc.gs_energy()
    out = sc.evolve_imaginary_time(wf, 30., shift="energy")
    assert out.dot(out).real == pytest.approx(1.0, abs=1e-10)
    assert out.aMb(h, out).real == pytest.approx(e0, abs=1e-6)
    out2 = sc.evolve_imaginary_time(wf, 30., shift=e0)
    assert out2.aMb(h, out2).real == pytest.approx(e0, abs=1e-6)
    # a shift well below E0 stops the projection: at long steps the
    # trapezoid rule sends every mode to -1, the ground state included
    out3 = sc.evolve_imaginary_time(wf, 30., shift=e0-0.5)
    assert out3.aMb(h, out3).real > e0+1e-2


@pytest.mark.parametrize("backend", BACKENDS)
def test_refusals(backend):
    sc, h, wf, wfed = _neel_chain(backend, n=6)
    with pytest.raises(ValueError):
        sc.evolve_imaginary_time(wf, 1., order=3)
    with pytest.raises(ValueError):
        sc.evolve_imaginary_time(wf, 1., shift="lowest")
    with pytest.raises(ValueError):
        sc.evolve_imaginary_time(wf, -1.)
    sc.set_hamiltonian(h+0.3j*sc.Sz[0])
    with pytest.raises(ValueError):
        sc.evolve_imaginary_time(wf, 1.)


@pytest.mark.parametrize("backend", BACKENDS)
def test_metts_implicit_steps_follow_the_tdvp_ones(backend):
    """The same Markov chain (same seed) evolved by TDVP and by implicit
    steps: the two evolutions differ only by their time steps, so the
    samples, and so the averages, agree to that accuracy."""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=backend, maxm=30,
                              cutoff=1e-10)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(5))
    sc.set_hamiltonian(h)
    kw = dict(nsamples=6, nwarmup=1, seed=4)
    m1, _ = sc.metts_vev(h, 0.1, itime="tdvp", **kw)
    m2, _ = sc.metts_vev(h, 0.1, itime="implicit", **kw)
    m3, _ = sc.metts_vev(h, 0.1, **kw)    # the default takes implicit here
    assert abs(m1-m2) < 1e-3
    assert m2 == m3
    assert imaginarytime.metts_plan(sc, 1.0) is None   # and TDVP here
    assert imaginarytime.metts_plan(sc, 0.1) is not None
    ts, c1, _ = sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 0.1,
            nt=3, dt=0.1, itime="tdvp", **kw)
    ts, c2, _ = sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 0.1,
            nt=3, dt=0.1, itime="implicit", **kw)
    assert np.max(np.abs(c1-c2)) < 1e-3
    with pytest.raises(TypeError):
        sc.metts_vev(h, 0.1, itime="implicit", itime_options=dict(nstep=2), **kw)
    with pytest.raises(ValueError):
        sc.metts_vev(h, 0.1, itime="tdvp", itime_options=dict(nsteps=2), **kw)


def test_ed_reference_does_not_overflow_on_long_steps():
    """e^{-dtau H} on a step of order 100 overflows for E0 < 0 unless H is
    shifted by its ground-state energy first."""
    sc, h, wf, wfed = _neel_chain("python", n=6)
    taus, e = sc.evolve_imaginary_time_and_measure(h, wfed, 2000., mode="ED")
    assert np.all(np.isfinite(e))
    assert e[-1].real == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-10)


@pytest.mark.parametrize("backend", BACKENDS)
def test_shift_is_the_lower_edge_after_a_set_state(backend):
    """After set_gs() of a state that is not the ground state the chain's
    energy is that state's own, 1.6 above E0 here; the shift must still be
    E0, or the long steps meet the pole of the trapezoid factor."""
    sc, h, wf, wfed = _neel_chain(backend, cutoff=1e-12)
    e0 = sc.gs_energy()
    sc.set_gs(wf)
    assert sc.gs_energy() > e0+1.
    assert imaginarytime.metts_plan(sc, 0.05)["shift"] == pytest.approx(e0, abs=1e-8)
    out = sc.evolve_imaginary_time(wf, 30.)
    assert out.aMb(h, out).real == pytest.approx(e0, abs=1e-6)


def test_metts_implicit_steps_in_parallel_chains():
    """The plan reaches the spawned workers of njobs>1 ("python")."""
    sc = spinchain.Spin_Chain(["S=1/2"]*5, itensor_version="python", maxm=20,
                              cutoff=1e-10)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(4))
    sc.set_hamiltonian(h)
    kw = dict(nsamples=4, nwarmup=1, seed=5, njobs=2)
    m1, _ = sc.metts_vev(h, 0.1, itime="tdvp", **kw)
    m2, _ = sc.metts_vev(h, 0.1, itime="implicit", **kw)
    assert abs(m1-m2) < 1e-3
