"""Functional coverage for dmrgpy.timeevolution.evolve_WF, distilled
from examples/time_evolution/time_evolution_wavefunction down to a
small, fast, deterministic system.

Unlike the other test files here, this doesn't do an ED/v2/v3
cross-backend comparison: the DMRG-side time evolution
(mpsalgebra.exponential_dmrg) is a Trotterized/fitted approximation,
while the ED side is an exact matrix exponential, so the two are not
expected to agree to DMRG_TOL-style precision without also tuning the
Trotter step count -- exercising exact ED evolution here is enough to
cover the feature.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy import spinchain
from dmrgpy import timedependent
from dmrgpy.timeevolution import evolve_WF


def test_time_evolution_conserves_norm_and_matches_reference_trajectory():
    """Two-site chain: prepare a Neel-like initial state as the ground
    state of a staggered field, then time-evolve under the Heisenberg
    Hamiltonian. The evolution must be norm-conserving (it's a unitary
    exponential), and <Sz_0>(t) must match a golden regression
    trajectory."""
    n = 2
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf0 = sc.get_gs(mode="ED")

    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    h = h + h.get_dagger()

    ts = np.linspace(0., 4., 5)
    wfs = evolve_WF(h, wf0, ts=ts)

    norms = np.array([w.dot(w) for w in wfs])
    assert norms.real == pytest.approx(np.ones(len(ts)), abs=1e-8)
    assert norms.imag == pytest.approx(np.zeros(len(ts)), abs=1e-8)

    sz0 = np.array([w.dot(sc.Sz[0] * w).real for w in wfs])
    expected = np.array([
        -0.5,
        0.20807341827357126,
        0.32682181043180597,
        -0.48008514332518304,
        0.07275001690430655,
    ])
    assert sz0 == pytest.approx(expected, abs=1e-6)


def test_tebd_matches_ed_spin_chain():
    """pyitensor/tebd.py's TEBDEvolver (self.tevol_method="TEBD",
    itensor_version="python" only) against exact diagonalization: prepare
    the GS of a Neel-favoring staggered field, quench to the (strictly
    nearest-neighbor) Heisenberg Hamiltonian, and require <Sz_0>(t) to
    match ED closely -- modeled on
    examples/time_evolution/tdvp_VS_ED_time_evolution, the same pattern
    used to validate TDVP."""
    n = 4
    spins = [2 for _ in range(n)]  # S=1/2
    sc = spinchain.Spin_Chain(spins)
    sc.setup_python()
    sc.tevol_method = "TEBD"

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]

    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    wfED = sc.get_gs(mode="ED")
    sc.set_hamiltonian(h1)

    nt, dt = 50, 0.05
    (_ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
    (_tsED, szED) = timedependent.evolve_and_measure(
            sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wfED, mode="ED")

    assert np.array(sz).real == pytest.approx(np.array(szED).real, abs=1e-4)


def test_tebd_matches_ed_fermion_chain():
    """Same cross-check as test_tebd_matches_ed_spin_chain, but for a
    spinless-fermion nearest-neighbor hopping Hamiltonian, to exercise
    TEBD's Jordan-Wigner handling (bond_hamiltonians() extracts each
    term's two adjacent-site matrices directly from HTerm.resolve(),
    relying on JW strings never needing to extend past an already-
    adjacent bond -- this is the case that would break if that reasoning
    were wrong)."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.tevol_method = "TEBD"

    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)

    nt, dt = 40, 0.05
    (_x, y) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="ED", A=fc.Cdag[0], B=fc.N[0])
    (_x1, y1) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="DMRG", A=fc.Cdag[0], B=fc.N[0])

    assert np.array(y1) == pytest.approx(np.array(y), abs=1e-4)


def test_tebd_rejects_non_nearest_neighbor_hamiltonian():
    """bond_hamiltonians() must fail loudly, not silently drop the
    long-range piece, when the Hamiltonian isn't strictly nearest-
    neighbor -- TEBD has no representation for a term spanning 3+ sites."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.tevol_method = "TEBD"

    h = fc.Cdag[0] * fc.C[2]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)

    with pytest.raises(NotImplementedError):
        timedependent.evolution_ABA(fc, nt=5, dt=0.05, mode="DMRG", A=fc.Cdag[0], B=fc.N[0])


def test_tebd_dynamical_correlator_matches_ed_spin_chain():
    """timedependent.evolution_DC (the raw time-domain quench correlator
    C(t)=<GS|A(t)B(0)|GS> behind get_dynamical_correlator(submode="TD"),
    dispatching to Chain.quench_tebd() when tevol_method="TEBD") must
    match exact diagonalization on a small nearest-neighbor system. This
    exercises quench_tebd() specifically -- unlike
    test_tebd_matches_ed_spin_chain above, which only covers
    evolve_and_measure_tebd() (a single operator measured along one
    evolved state, no ground-state-energy-shift/two-state-overlap logic).
    Compared directly in the time domain, not after the submode="TD"
    Fourier transform, since ED's own get_dynamical_correlator() uses an
    unrelated (KPM-moment-based) construction with a different
    normalization convention -- comparing post-FFT spectra would conflate
    that pre-existing convention mismatch (present identically for
    tevol_method="TDVP") with a genuine TEBD correctness check."""
    n = 4
    spins = [2 for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    sc.setup_python()
    sc.tevol_method = "TEBD"

    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)

    name = (sc.Sz[0], sc.Sz[0])
    nt, dt = 60, 0.05
    sc.get_gs()
    (_ts_ed, cs_ed) = timedependent.evolution_DC(sc, mode="ED", name=name, nt=nt, dt=dt)
    (_ts_tebd, cs_tebd) = timedependent.evolution_DC(sc, mode="DMRG", name=name, nt=nt, dt=dt)

    assert np.array(cs_tebd) == pytest.approx(np.array(cs_ed), abs=1e-4)


def test_tebd_dynamical_correlator_matches_ed_fermion_chain():
    """Same cross-check as test_tebd_dynamical_correlator_matches_ed_spin_chain,
    but for spinless-fermion nearest-neighbor hopping, to exercise
    quench_tebd()'s Jordan-Wigner handling on the two-operator (A,B)
    quench correlator rather than the single-operator evolve_and_measure
    path already covered by test_tebd_matches_ed_fermion_chain."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.tevol_method = "TEBD"

    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)

    name = (fc.Cdag[0], fc.C[0])
    nt, dt = 40, 0.05
    fc.get_gs()
    (_ts_ed, cs_ed) = timedependent.evolution_DC(fc, mode="ED", name=name, nt=nt, dt=dt)
    (_ts_tebd, cs_tebd) = timedependent.evolution_DC(fc, mode="DMRG", name=name, nt=nt, dt=dt)

    assert np.array(cs_tebd) == pytest.approx(np.array(cs_ed), abs=1e-4)


def test_tdvp_lanczos_falls_back_when_stemr_fails_to_converge():
    """pyitensor/tdvp.py's _eigh_tridiagonal_robust() must fall back to
    lapack_driver="stebz" and still return a correct eigendecomposition
    when the default "stemr" driver raises LinAlgError -- a real failure
    mode hit while benchmarking TDVP's dynamical correlator on a
    20-orbital native-Hubbard chain (Chain.quench_tdvp's per-step Lanczos
    expm-multiply, niter=50): scipy's eigh_tridiagonal("stemr") does not
    converge for some tridiagonal matrices with (near-)degenerate
    eigenvalues, aborting the whole quench with an uncaught
    numpy.linalg.LinAlgError. Simulated here by monkeypatching
    eigh_tridiagonal to always raise on the default driver, so the test
    doesn't depend on reproducing the exact degenerate matrix that
    triggered this originally."""
    from scipy.linalg import eigh_tridiagonal as real_eigh_tridiagonal
    from dmrgpy.pyitensor import tdvp as tdvp_mod

    rng = np.random.default_rng(0)
    n = 6
    alphas = rng.normal(size=n)
    betas = np.abs(rng.normal(size=n - 1)) + 0.1

    def flaky_eigh_tridiagonal(a, b, lapack_driver="stemr"):
        if lapack_driver == "stemr":
            raise np.linalg.LinAlgError("stemr (eigh_tridiagonal) did not converge")
        return real_eigh_tridiagonal(a, b, lapack_driver=lapack_driver)

    original = tdvp_mod.eigh_tridiagonal
    tdvp_mod.eigh_tridiagonal = flaky_eigh_tridiagonal
    try:
        evals, evecs = tdvp_mod._eigh_tridiagonal_robust(alphas, betas)
    finally:
        tdvp_mod.eigh_tridiagonal = original

    expected_evals, expected_evecs = real_eigh_tridiagonal(alphas, betas, lapack_driver="stebz")
    assert evals == pytest.approx(expected_evals, abs=1e-10)
    dense = np.diag(alphas) + np.diag(betas, 1) + np.diag(betas, -1)
    assert dense @ evecs == pytest.approx(evecs @ np.diag(evals), abs=1e-8)


def test_mps_copy_is_independent_under_python_backend():
    """mps.MPS.copy() must return a wavefunction whose subsequent
    evolution doesn't affect the original. itensor_version="python"'s
    cpp_handle is a mutable pyitensor.mpscontainer.MPS that TDVP/TEBD
    mutate in place via set_A() (unlike the C++ backends' handles, which
    are immutable by convention -- every Chain method there takes wf by
    const reference and returns a brand-new MPS) -- so a naive shallow
    copy of the outer MPS wrapper aliases the same underlying tensor
    list, and evolving the "copy" silently corrupts the original too.
    This is a regression test for that exact bug, found while
    benchmarking TEBD against TDVP: comparing the two methods starting
    from the same wf, or using evolve_and_measure_dmrg's own documented
    forward+backward round-trip pattern, both silently produced wrong
    results before mps.py's copy() was fixed to special-case this."""
    from dmrgpy.pyitensor.autompo import AutoMPO
    from dmrgpy.pyitensor.mpobuilder import to_mpo
    from dmrgpy.pyitensor.mpsalgebra import inner

    n = 4
    sc = spinchain.Spin_Chain([2 for _ in range(n)])
    sc.setup_python()

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()

    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h1)

    op = sc.Sz[0]
    Aop = to_mpo(AutoMPO.from_terms(sc._session.sites, op.to_terms()))
    val_before = inner(wf.cpp_handle, Aop, wf.cpp_handle)

    sc.tevol_method = "TEBD"
    timedependent.evolve_and_measure(sc, operator=op, nt=10, dt=0.05, wf=wf.copy())

    val_after = inner(wf.cpp_handle, Aop, wf.cpp_handle)
    assert val_after == pytest.approx(val_before, abs=1e-10)
