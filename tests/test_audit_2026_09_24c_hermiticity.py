"""Regression tests for the `hermiticity` cluster of the third 2026-09-24
hole hunt (docs/audit_2026_09_24c_hole_hunt.md, findings 12 and 18).

  12. The DMRG Hermiticity probe compared ||(A-A^dagger)w||^2, on a unit
      witness but for an unnormalized A, against an absolute 1e-4, so any
      anti-Hermitian part below about 1e-2 in the units the operator was
      written in was called Hermitian: gs_energy() dropped a weak decay rate
      whole on v2 and v3 and returned a real part 3 to 18 per cent off on
      "python", and disentangle_manifold diagonalized the Hermitian part of
      a small non-Hermitian operator. The probe now decides on A/max|coef|
      at a roundoff tolerance, and the two ED checks share one relative
      test.
  18. disentangle._is_hermitian fell back to the bare proof only when a
      state's MBO was None, but an ED State carries its EDchain, which had
      no is_hermitian, so every ED manifold raised AttributeError. EDchain
      now decides exactly on the operator's matrix.

The anchors are dense diagonalizations of the ED matrix and exact
eigenvectors of the representation matrix, so no convergence enters beyond
the pinned sweep schedules.
"""

import warnings

import numpy as np
import pytest
import scipy.linalg as dlg
import scipy.sparse as sparse

from dmrgpy import cppext, fermionchain, mpsalgebra, spinchain
from dmrgpy.algebra import algebra
from dmrgpy.mpsalgebratk.disentangle import get_representation


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3), _backend(2)]


def _heisenberg(sc, n):
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h


def _exact_lowest(sc):
    """Eigenvalue of lowest real part of the ED matrix, and its conjugate
    (a real non-Hermitian matrix has its complex levels in conjugate pairs,
    and NH-DMRG may land on either member)."""
    ev = np.linalg.eigvals(sc.get_ED_obj().get_hamiltonian().toarray())
    return ev[np.argsort(ev.real)][0]


# ------------------------------------------------------------ finding 12

@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("s", [1.0, 1e-2, 1e-3, 1e-6])
def test_probe_verdict_does_not_depend_on_units(version, s):
    """H = Heis4 + 0.3*Sz_tot + 1j*Sz0 is non-Hermitian at every scale; the
    old probe called it Hermitian from s=1e-2 down. An unprovable Hermitian
    operator (1j*Sx0*Sy0 is exactly -Sz0/2) stays Hermitian at every scale."""
    n = 4
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    base = _heisenberg(sc, n) + 0.3*sum(sc.Sz[i] for i in range(n))
    np.random.seed(1)
    assert not sc.is_hermitian(s*(base + 1j*sc.Sz[0]))
    np.random.seed(1)
    assert sc.is_hermitian(s*(base + 1j*sc.Sx[0]*sc.Sy[0]))


@pytest.mark.parametrize("version", BACKENDS)
def test_weak_loss_keeps_its_decay_rate(version):
    """Heis4 + 1.5*Sz_tot + 1j*5e-3*Sz0 in natural units: the ground state is
    polarized, so Im E0 is first order in the loss, -0.002134. v2 and v3
    returned the real part alone and "python" a real number 3 to 18 per cent
    off; the non-Hermitian solver is exact."""
    n = 4
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    h = _heisenberg(sc, n) + 1.5*sum(sc.Sz[i] for i in range(n)) + 1j*5e-3*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    np.random.seed(2)
    assert not sc.is_hermitian(sc.hamiltonian)
    e0 = np.complex128(sc.gs_energy())
    ex = _exact_lowest(sc)
    assert abs(ex.imag) == pytest.approx(0.002134, abs=1e-6)
    assert e0.real == pytest.approx(ex.real, abs=1e-6)
    assert abs(e0.imag) == pytest.approx(abs(ex.imag), abs=1e-6)


@pytest.mark.parametrize("version", BACKENDS)
def test_rescaled_non_hermitian_hamiltonian_is_solved_as_one(version):
    """The same non-Hermitian model at s=1e-2 (the old probe's knife edge,
    where the verdict flipped between runs) solves to s*E0 with E0 =
    -1.242697 +- 0.291419i, a conjugate pair."""
    n, s = 4, 1e-2
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    h = s*(_heisenberg(sc, n) + 0.3*sum(sc.Sz[i] for i in range(n)) + 1j*sc.Sz[0])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    np.random.seed(3)
    e0 = np.complex128(sc.gs_energy())/s
    ex = _exact_lowest(sc)/s
    assert ex.real == pytest.approx(-1.242697, abs=1e-6)
    assert e0.real == pytest.approx(ex.real, abs=1e-5)
    assert abs(e0.imag) == pytest.approx(abs(ex.imag), abs=1e-5)


@pytest.fixture(scope="module")
def spin_manifold_3():
    """The 8 eigenstates of a generic 3-spin Hamiltonian at full bond
    dimension (the misc hunter's script 01)."""
    sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
    h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
        + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 8, 20
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        es, wfs = sc.get_excited_states(n=8)
    return sc, wfs


def _eigen_residual(wfs, out, A):
    """max over output states of ||(ma - lambda) v|| / ||v|| in the input
    manifold's coordinates, relative to max|ma|"""
    ma = get_representation(wfs, A)
    C = np.array([[a.dot(b) for b in out] for a in wfs])
    worst = 0.
    for j in range(C.shape[1]):
        v = C[:, j]
        lam = (v.conj() @ ma @ v)/(v.conj() @ v)
        worst = max(worst, np.linalg.norm(ma @ v - lam*v)/np.linalg.norm(v))
    return worst/np.max(np.abs(ma))


@pytest.mark.parametrize("eps", [1e-2, 1e-3, 1e-4])
def test_small_non_hermitian_operator_is_disentangled_on_its_eigenvectors(spin_manifold_3, eps):
    """A = eps*(Sz0 + S+0) + 0.3*eps*Sz1 is non-Hermitian at every eps; from
    eps=1e-2 down the old probe sent it to eigh of its Hermitian part, whose
    vectors are 0.702 of max|ma| off being eigenvectors of A."""
    sc, wfs = spin_manifold_3
    A = eps*(sc.Sz[0] + sc.Sx[0] + 1j*sc.Sy[0]) + 0.3*eps*sc.Sz[1]
    np.random.seed(7)
    assert not sc.is_hermitian(A)
    out = mpsalgebra.disentangle_manifold(wfs, A)
    assert _eigen_residual(wfs, out, A) < 1e-10


def test_ed_hermiticity_checks_agree_and_are_relative():
    """algebra.is_hermitian (read by the ED correlator dispatch) and
    algebra.ishermitian (read by lowest_states) used to be an absolute 1e-8
    on ||h-h^dag||_F^2 and an absolute 1e-6 on max|h-h^dag|, and disagreed
    at c=3e-6. They are one relative test now."""
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    ed = sc.get_ED_obj()
    h = ed.MO2matrix(_heisenberg(sc, 4))
    d = ed.MO2matrix(1j*sc.Sz[0])
    for c in (1.0, 3e-6, 1e-8):
        m = h + c*d
        assert algebra.is_hermitian(m) is False
        assert algebra.ishermitian(m) is False
        assert algebra.ishermitian(m.toarray()) is False
    for s in (1.0, 1e-6, 1e-12):
        assert algebra.is_hermitian(s*h) and algebra.ishermitian(s*h)
    assert algebra.is_hermitian(sparse.csc_matrix((4, 4)))


# ------------------------------------------------------------ finding 18

@pytest.fixture(scope="module")
def ed_spin_manifold():
    sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
    h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
        + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
    sc.set_hamiltonian(h)
    es, wfs = sc.get_excited_states(n=8, mode="ED")
    return sc, wfs


@pytest.fixture(scope="module")
def ed_fermion_manifold():
    fc = fermionchain.Fermionic_Chain(3, itensor_version="python")
    hf = sum(fc.Cdag[i]*fc.C[i+1] for i in range(2))
    hf = hf + hf.get_dagger() + 0.3*fc.N[0] + 0.7*fc.N[0]*fc.N[1]
    fc.set_hamiltonian(hf)
    es, wfs = fc.get_excited_states(n=8, mode="ED")
    return fc, wfs


def _gram_and_offdiag(out, A):
    G = np.array([[a.dot(b) for b in out] for a in out])
    ma = get_representation(out, A)
    return (np.max(np.abs(G - np.eye(len(out)))),
            np.max(np.abs(ma - np.diag(np.diag(ma)))))


@pytest.mark.parametrize("which", ["Sz0", "1j*Sx0*Sy0"])
def test_ed_spin_manifold_hermitian_operator_is_disentangled(ed_spin_manifold, which):
    """Both raised AttributeError on an ED manifold after 867e2b4; the
    unprovable one gave a basis 0.340 off orthonormal before it."""
    sc, wfs = ed_spin_manifold
    A = sc.Sz[0] if which == "Sz0" else 1j*sc.Sx[0]*sc.Sy[0]
    assert wfs[0].MBO.is_hermitian(A)
    gram, off = _gram_and_offdiag(mpsalgebra.disentangle_manifold(wfs, A), A)
    assert gram < 1e-12 and off < 1e-12


def test_ed_spin_manifold_non_hermitian_operator_takes_eig(ed_spin_manifold):
    sc, wfs = ed_spin_manifold
    A = sc.Sx[0] + 1j*sc.Sy[0]
    assert not wfs[0].MBO.is_hermitian(A)
    out = mpsalgebra.disentangle_manifold(wfs, A)
    assert _eigen_residual(wfs, out, A) < 1e-10


@pytest.mark.parametrize("which", ["N0", "C0+Cdag0"])
def test_ed_fermion_manifold_is_disentangled(ed_fermion_manifold, which):
    fc, wfs = ed_fermion_manifold
    A = fc.N[0] if which == "N0" else fc.C[0] + fc.Cdag[0]
    gram, off = _gram_and_offdiag(mpsalgebra.disentangle_manifold(wfs, A), A)
    assert gram < 1e-12 and off < 1e-12
