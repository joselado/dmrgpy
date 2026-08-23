"""The ED path must return every member of a degenerate level.

Above `algebra.maxsize` the ED path goes through ARPACK, which stops once
the Ritz pairs it holds have converged -- a condition ONE copy of a
degenerate eigenvalue satisfies. The remaining copies were then silently
never produced, so `mode="ED"` appeared to return distinct levels while
`mode="DMRG"` returned multiplet members, and DMRG got blamed for
duplicating states (see docs/ed_vs_dmrg_degenerate_multiplets.md).

The model is a FERROMAGNETIC spin-1/2 Heisenberg chain, chosen because its
answer is analytic: the ground state is the fully polarized maximal-spin
multiplet, so with H = -sum_i S_i.S_{i+1} the ground level is exactly
-(N-1)/4 and exactly (N+1)-fold degenerate. That makes the assertions
independent of any reference diagonalization at the sizes where a dense one
would be expensive.
"""
import numpy as np
import pytest
import scipy.linalg as dlg

from dmrgpy import spinchain
from dmrgpy.algebra import algebra


def _ferro(n, maxm=40):
    """H = -sum S_i.S_{i+1}: ground level -(n-1)/4, degeneracy n+1."""
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    sc.maxm, sc.nsweeps = maxm, 20
    h = 0
    for i in range(n - 1):
        h = h - (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                 + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(h)
    return sc


def _random_field(n, seed=1):
    """A chain with no degeneracy at all: the other side of the contract."""
    rng = np.random.default_rng(seed)
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + 0.7 * sc.Sy[i] * sc.Sy[i + 1] \
            + 1.3 * sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + float(rng.normal()) * sc.Sz[i] \
            + 0.3 * float(rng.normal()) * sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc


@pytest.fixture
def force_arpack(monkeypatch):
    """Route even a tiny Hilbert space through ARPACK, so the path under
    test can be checked against a dense reference that is cheap to
    compute."""
    monkeypatch.setattr(algebra, "maxsize", 8)


@pytest.mark.parametrize("n_ask", [4, 8, 13])
def test_degenerate_ground_multiplet_is_returned_in_full(n_ask):
    """N=12 (dim 4096, above maxsize, so genuinely the ARPACK path): the
    ground level is 13-fold, so every one of the first n_ask <= 13 energies
    must equal it. Before the fix these came back as 1, 3 and 4 copies."""
    n = 12
    es = np.array(_ferro(n).get_excited(n=n_ask, mode="ED")).real
    assert len(es) == n_ask
    assert es == pytest.approx(-(n - 1) / 4.0, abs=1e-8)


def test_asking_past_the_multiplet_finds_the_next_level():
    """The degeneracy must not be padded out forever either: level 14 of the
    N=12 chain is a genuinely higher one."""
    n = 12
    es = np.array(_ferro(n).get_excited(n=15, mode="ED")).real
    e0 = -(n - 1) / 4.0
    assert es[:13] == pytest.approx(e0, abs=1e-8)
    assert es[13] > e0 + 1e-3
    assert np.all(np.diff(es) > -1e-8)  # sorted


def test_matches_dense_elementwise_on_the_arpack_path(force_arpack):
    """Not just the ground level: the whole returned list must equal a dense
    diagonalization entry by entry, multiplicities included."""
    sc = _ferro(8)
    h = sc.get_ED_obj().get_hamiltonian()
    ref = np.sort(dlg.eigvalsh(np.array(h.todense())))
    es, _ = algebra.lowest_states(h, n=12)
    assert np.array(es).real == pytest.approx(ref[:12], abs=1e-8)


def test_a_nondegenerate_spectrum_is_unaffected(force_arpack):
    """Deflation must not perturb the ordinary case it is not there for."""
    sc = _random_field(8)
    h = sc.get_ED_obj().get_hamiltonian()
    ref = np.sort(dlg.eigvalsh(np.array(h.todense())))
    es, _ = algebra.lowest_states(h, n=12)
    assert np.array(es).real == pytest.approx(ref[:12], abs=1e-8)


def test_returned_vectors_are_orthonormal_eigenvectors(force_arpack):
    """Copies recovered by deflation must be real eigenvectors, not Ritz
    vectors taken on faith, and mutually orthogonal."""
    sc = _ferro(8)
    h = sc.get_ED_obj().get_hamiltonian()
    es, vs = algebra.lowest_states(h, n=10)
    v = np.array([np.asarray(x).reshape(-1) for x in vs])
    assert v @ np.conjugate(v.T) == pytest.approx(np.eye(len(vs)), abs=1e-8)
    for e, w in zip(es, v):
        r = np.asarray(h @ w).reshape(-1) - e * w
        assert np.linalg.norm(r) < 1e-7


def test_ed_and_dmrg_agree_index_by_index():
    """The point of the whole thing: on a degenerate spectrum the two modes
    now return the same list, so validating DMRG against ED positionally --
    the natural thing to do -- no longer reports a spurious defect."""
    sc = _ferro(12)
    ed = np.array(sc.get_excited(n=6, mode="ED")).real
    dmrg = np.array(sc.get_excited(n=6, mode="DMRG")).real
    assert dmrg == pytest.approx(ed, abs=1e-4)
