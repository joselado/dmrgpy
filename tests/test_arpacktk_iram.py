"""Regression coverage for algebra/arpacktk.py's Implicitly Restarted
Arnoldi Method (IRAM, ported from ARPACK), the default MPS Arnoldi
solver used by get_excited_states (non-Hermitian, n>=2) and by
degeneracy.py's eigenvalue_degeneracy. Before this file existed, none of
that code path had any pytest coverage at all -- a real crash/silent-
wrong-answer bug in mpsiram_shift_invert would have shipped silently
through run_tests.py staying green (found and fixed via manual stress
testing, see the PR history).
"""
import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain
from dmrgpy.algebra import arnolditk
from dmrgpy.algebra.arpacktk import mpsiram, mpsiram_shift_invert
from dmrgpy.algebra.arpacktk import mpsiram_generalized, generalized_excited_states


@pytest.fixture(autouse=True)
def _isolate_shared_global_state():
    """Two pieces of shared, module/process-global state this file must
    not leak or depend on:

    arnolditk.arnoldimode -- the module-level flag arpacktk.py's
    random_state()/toMPO() read to pick ED vs DMRG. Tests here set it to
    exercise both modes; without saving/restoring it, whatever the last
    test in this file left it set to would leak into every other test
    file that runs after it (this file has no control over execution
    order).

    numpy's global RNG -- edtk/edchain.py's State.random_state() (the ED
    IRAM start-vector generator) draws from np.random.random() directly,
    not a seeded local generator, so how many prior np.random calls
    happened anywhere earlier in the process (other tests, or even
    earlier iterations of the same IRAM run) affects the exact numerical
    trajectory and hence convergence precision -- confirmed directly: a
    recursive multi-state test passed in isolation but failed when run
    after other tests in this file consumed different amounts of the
    global stream. Seeding it fresh before every test makes each test's
    outcome depend only on its own random draws, not on execution order."""
    original = arnolditk.arnoldimode
    np.random.seed(12345)
    yield
    arnolditk.arnoldimode = original


def _heisenberg_chain(n=6, seed=1):
    rng = np.random.RandomState(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * sc.Sx[i] * sc.Sx[i + 1]
        h = h + rng.random() * sc.Sy[i] * sc.Sy[i + 1]
        h = h + rng.random() * sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h


def _staggered_fermion_chain(n=6, seed=4):
    """Interacting fermion chain with a staggered imaginary potential --
    non-Hermitian, complex spectrum, exact complex-conjugate degenerate
    ground state (used elsewhere in this repo, e.g.
    examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi)."""
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    mh = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        mh[i, i + 1] = 1.0
        mh[i + 1, i] = 1.0
    for i in range(n):
        mh[i, i] = 1j * (-1) ** i
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + mh[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


def test_mpsiram_hermitian_ground_state_ed():
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=1)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=1, tol=1e-8)
    assert abs(es[0].real - es_ed[0]) < 1e-6
    wf = wfs[0]
    assert abs(wf.dot(wf) - 1.0) < 1e-8 # normalized


def test_mpsiram_simultaneous_multi_eigenvalue_ed():
    """A genuine nev=3 simultaneous block search (not one-at-a-time)."""
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=3)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=3, tol=1e-6)
    assert len(es) == 3
    for e_iram, e_ed in zip(sorted(es.real), es_ed):
        assert abs(e_iram - e_ed) < 1e-4


def test_mpsiram_non_hermitian_conjugate_pair_ed():
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    es, wfs = mpsiram(fc, h, which="SR", nev=1, tol=1e-6)
    assert abs(es[0].real - es_ed[0].real) < 1e-4
    assert min(abs(es[0] - x) for x in es_ed) < 1e-4


def test_mpsiram_shift_invert_finds_nearest_eigenvalue_ed():
    sc, h = _heisenberg_chain(n=6, seed=1)
    arnolditk.arnoldimode = "ED"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=6)).real)
    target = es_ed[2] + 0.01
    es, wfs = mpsiram_shift_invert(sc, h, e=target, nev=1, delta=1e-6, maxn=40)
    closest = es_ed[np.argmin(np.abs(es_ed - target))]
    assert abs(es[0].real - closest) < 1e-3


def test_mpsiram_shift_invert_deflated_second_state_is_a_genuine_eigenvalue_ed():
    """A recursive (deflated) second shift-invert search must still land
    on *some* genuine eigenvalue of H, not a spurious combination --
    even though it is not guaranteed to land on the exact conjugate
    partner of the first (see module note on deflation and non-Hermitian
    operators: right eigenvectors of a non-Hermitian H need not be
    mutually orthogonal under the plain Hermitian inner product -- e.g.
    confirmed directly for this model's ground-state conjugate pair,
    |<v1|v2>| ~= 0.40 -- so .dot()-based deflation, used identically by
    arnolditk's own wfskip mechanism, is only a heuristic here, not a
    rigorous one)."""
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    e0 = es_ed[0].real
    es, wfs = mpsiram_shift_invert(fc, h, e=1.2 * e0, nev=2, delta=1e-3, maxn=40)
    for e in es:
        assert min(abs(e - x) for x in es_ed) < 1e-2


def test_excited_states_non_hermitian_default_matches_ed():
    """get_excited_states' non-Hermitian route (excited.py's
    excited_states_non_hermitian, IRAM by default) must agree with ED."""
    fc, h = _staggered_fermion_chain(n=6, seed=4)
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=3), key=lambda x: x.real))
    es, wfs = fc.get_excited_states(n=3, purify=False)
    es_sorted = np.array(sorted(es, key=lambda x: x.real))
    for e_iram, e_ed in zip(es_sorted, es_ed):
        assert abs(e_iram.real - e_ed.real) < 1e-3


def test_excited_states_non_hermitian_accepts_legacy_nkry_kwargs():
    """excited_states_non_hermitian's old arnolditk-based signature had
    nkry_min/nkry_max as explicit named parameters; the arpacktk swap
    replaced them with ncv, but any caller still passing nkry_min/
    nkry_max explicitly (e.g. via get_excited_states(..., nkry_min=...))
    would otherwise hit a TypeError several calls deep once those stray
    kwargs reached mpsiram (no **kwargs catch-all there) -- confirmed via
    direct reproduction. Both must still be accepted (and produce a
    correct result), even though IRAM has no direct equivalent of
    arnolditk's adaptive Krylov-size range."""
    fc, h = _staggered_fermion_chain(n=4, seed=3)
    fc.setup_python()
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=3), key=lambda x: x.real))
    es, wfs = fc.get_excited_states(n=2, nkry_min=2)
    es_sorted = np.array(sorted(es, key=lambda x: x.real))
    for e_iram, e_ed in zip(es_sorted, es_ed[:2]):
        assert abs(e_iram.real - e_ed.real) < 1e-2
    # nkry_max together with nkry_min must also still work
    es2, wfs2 = fc.get_excited_states(n=2, nkry_min=3, nkry_max=8)
    es2_sorted = np.array(sorted(es2, key=lambda x: x.real))
    for e_iram, e_ed in zip(es2_sorted, es_ed[:2]):
        assert abs(e_iram.real - e_ed.real) < 1e-2


def test_eigenvalue_degeneracy_non_degenerate_target_ed():
    """degeneracy.py's eigenvalue_degeneracy (shift-invert IRAM, called
    recursively with deflation) targeted at a genuinely non-degenerate
    eigenvalue must report a degeneracy close to 1.

    Note: this deliberately avoids this model's own exact
    complex-conjugate ground-state pair as a target -- right
    eigenvectors of a non-Hermitian H need not be mutually orthogonal
    under the plain Hermitian inner product (confirmed directly here,
    |<v1|v2>| ~= 0.40 for that pair), so the .dot()-based deflation this
    routine relies on (identical to arnolditk's own wfskip mechanism) is
    only a heuristic, not a rigorous one, for resolving non-orthogonal
    degenerate clusters -- a pre-existing characteristic of the
    deflation approach itself, not specific to this solver."""
    from dmrgpy.degeneracy import eigenvalue_degeneracy
    fc, h = _staggered_fermion_chain(n=4, seed=3)
    arnolditk.arnoldimode = "ED"
    es_ed = np.array(sorted(fc.get_excited(mode="ED", n=6), key=lambda x: x.real))
    # the "-0.25" level in this model is an isolated, non-degenerate
    # eigenvalue (confirmed via the ED spectrum)
    target = es_ed[np.argmin(np.abs(es_ed.real - (-0.25)))].real
    deg = eigenvalue_degeneracy(fc, h, target, n=1, delta=5e-2)
    assert 0.8 < deg < 1.5


@pytest.mark.parametrize("itensor_version", ["python"])
def test_mpsiram_dmrg_mode_matches_ed(itensor_version):
    """Real MPS/MPO (pyitensor backend, no compiled extension needed)
    ground state via IRAM must agree with ED."""
    sc, h = _heisenberg_chain(n=4, seed=1)
    sc.setup_python()
    arnolditk.arnoldimode = "DMRG"
    es_ed = np.sort(np.array(sc.get_excited(mode="ED", n=1)).real)
    es, wfs = mpsiram(sc, h, which="SR", nev=1, tol=1e-6)
    assert abs(es[0].real - es_ed[0]) < 1e-4


def _generalized_fermion_problem(n=5, seed=2):
    """A*x = lambda*M*x: A a Hermitian interacting fermion Hamiltonian,
    M = 1 + 0.5*sum(N_i) a Hermitian positive-definite weight operator
    (diagonal, eigenvalues in [1, 1+0.5n], never zero)."""
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    a = 0
    for i in range(n - 1):
        a = a + rng.random() * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    for i in range(n - 1):
        a = a + rng.random() * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    m = 1
    for i in range(n):
        m = m + 0.5 * fc.N[i]
    fc.set_hamiltonian(a)
    return fc, a, m


def _generalized_ground_truth(fc, a, m):
    import scipy.linalg as sla
    Amat = fc.get_ED_obj().MO2matrix(a)
    Mmat = fc.get_ED_obj().MO2matrix(m)
    Amat = Amat.toarray() if hasattr(Amat, "toarray") else np.array(Amat)
    Mmat = Mmat.toarray() if hasattr(Mmat, "toarray") else np.array(Mmat)
    return np.sort(sla.eigh(Amat, Mmat, eigvals_only=True))


def test_mpsiram_generalized_ground_state_ed():
    """mode 2 (A*x=lambda*M*x, OP=inv(M)*A, B=M) ground state must match
    scipy.linalg.eigh's generalized Hermitian-definite eigensolver."""
    fc, a, m = _generalized_fermion_problem(n=5, seed=2)
    arnolditk.arnoldimode = "ED"
    w = _generalized_ground_truth(fc, a, m)
    es, wfs = mpsiram_generalized(fc, a, m, which="SR", nev=1, tol=1e-8)
    assert abs(es[0].real - w[0]) < 1e-5
    wf = wfs[0]
    assert abs(wf.dot(wf) - 1.0) < 1e-8 # plain-normalized, not M-normalized


def test_mpsiram_generalized_simultaneous_multi_eigenvalue_ed():
    fc, a, m = _generalized_fermion_problem(n=5, seed=2)
    arnolditk.arnoldimode = "ED"
    w = _generalized_ground_truth(fc, a, m)
    es, wfs = mpsiram_generalized(fc, a, m, which="SR", nev=3, tol=1e-6)
    for e_iram, e_true in zip(sorted(es.real), w[:3]):
        assert abs(e_iram - e_true) < 1e-4


def test_generalized_excited_states_recursive_ed():
    fc, a, m = _generalized_fermion_problem(n=5, seed=2)
    arnolditk.arnoldimode = "ED"
    w = _generalized_ground_truth(fc, a, m)
    es, wfs = generalized_excited_states(fc, a, m, nwf=3, which="SR",
            recursive=True, tol=1e-6)
    for e_iram, e_true in zip(sorted(es.real), w[:3]):
        assert abs(e_iram - e_true) < 1e-4


@pytest.mark.parametrize("itensor_version", ["python"])
def test_mpsiram_generalized_dmrg_mode_matches_ed(itensor_version):
    """Real MPS/MPO (pyitensor backend) mode-2 ground state must agree
    with the generalized-eigenvalue ED ground truth."""
    fc, a, m = _generalized_fermion_problem(n=4, seed=2)
    w = _generalized_ground_truth(fc, a, m)
    fc.setup_python()
    arnolditk.arnoldimode = "DMRG"
    es, wfs = mpsiram_generalized(fc, a, m, which="SR", nev=1, tol=1e-6)
    assert abs(es[0].real - w[0]) < 1e-4
