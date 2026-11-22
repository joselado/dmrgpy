"""The Gaussian basis of an impurity model (gaussianbasis.py, after
arXiv:2212.09798): bath modes rotated by the Fishman-White circuit of the
ground state without the interaction, the impurity left as it is.

The rotation is a unitary of the single-particle modes, so the rotated
Hamiltonian must have the spectrum of the original one exactly, and an
observable rotated the same way the same value; ED checks both to
1e-12. A complex hopping is used, since a real one passes with a
transposed or conjugated V."""

import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy.gaussianmpstk import basis as gbasis

from _helpers import julia_live_param


def _cpp_param(version):
    from dmrgpy import cppext
    return pytest.param(version, id="v%d" % version, marks=pytest.mark.skipif(
        not cppext.available(version),
        reason="the v%d extension is not compiled" % version))


BACKENDS = [_cpp_param(3), pytest.param("python", id="python"),
            julia_live_param()]


def _impurity_hopping(n, V=0.8, phase=0.4, imp=0):
    """An impurity orbital coupled by V to the end of a uniform bath with
    a Peierls phase, the impurity on orbital imp (0 or n-1)."""
    h = np.zeros((n, n), dtype=complex)
    order = list(range(n)) if imp == 0 else list(range(n-1, -1, -1))
    for k in range(n-1):
        i, j = order[k], order[k+1]
        h[i, j] = -(V if k == 0 else np.exp(1j*phase))
        h[j, i] = np.conj(h[i, j])
    h[order[2], order[2]] = 0.3  # breaks particle-hole symmetry
    return h


def _spinful_hamiltonian(fc, h, U=2.0, B=0.3):
    n = h.shape[0]
    H = sum(h[i, j]*(fc.Cdagup[i]*fc.Cup[j]+fc.Cdagdn[i]*fc.Cdn[j])
            for i in range(n) for j in range(n) if h[i, j] != 0)
    H = H + U*(fc.Nup[0]-0.5)*(fc.Ndn[0]-0.5)
    return H + B*fc.Sz[n-1]  # a field on a bath orbital, a rotated N


def _free_entropies(L):
    """Entanglement entropy of a Slater determinant at every cut."""
    out = []
    for k in range(1, L.shape[0]):
        nu = np.clip(np.linalg.eigvalsh(L[:k, :k]), 1e-15, 1-1e-15)
        out.append(float(-np.sum(nu*np.log(nu)+(1-nu)*np.log(1-nu))))
    return np.array(out)


def _uniform_correlation(n, imp):
    h = _impurity_hopping(n, V=1.0, imp=imp)
    e, u = np.linalg.eigh(h)
    phi = u[:, :n//2]
    return phi.conj() @ phi.T


def test_rotation_is_unitary_and_leaves_the_impurity():
    n = 12
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    b = fc.get_gaussian_basis(h=_impurity_hopping(n), nf=n, fixed=[0])
    V = b.V
    assert np.max(np.abs(V @ V.conj().T-np.eye(2*n))) < 1e-13
    for m in (0, 1):  # both spins of the impurity orbital
        e = np.eye(2*n)[m]
        assert np.max(np.abs(V[m]-e)) == 0.0
        assert np.max(np.abs(V[:, m]-e)) == 0.0
    assert np.max(np.abs(V[0::2, 1::2])) == 0.0  # the spins do not mix
    assert np.max(np.abs(V[1::2, 0::2])) == 0.0


@pytest.mark.parametrize("imp", [0, "end"])
def test_reference_is_a_product_away_from_the_impurity(imp):
    """In the new modes the reference state is a product state except at
    the impurity and its partner orbital, whichever end the impurity is
    on: the circuit runs toward it."""
    n = 40
    i0 = 0 if imp == 0 else n-1
    L = _uniform_correlation(n, i0)
    V, info = gbasis.rotation(L, fixed=[i0])
    S0 = _free_entropies(L)
    S1 = _free_entropies(V @ L @ V.conj().T)
    cut = 0 if i0 == 0 else n-2  # the cut that isolates the impurity
    p = np.real(L[i0, i0])
    Simp = -(p*np.log(p)+(1-p)*np.log(1-p))  # entropy of the impurity
    assert S1[cut] == pytest.approx(Simp, abs=1e-6)
    near = [0, 1] if i0 == 0 else [n-3, n-2]  # the impurity and its partner
    assert np.max(np.delete(S1, near)) < 1e-6
    assert np.max(S0) > 0.8  # the site basis is entangled everywhere


def test_two_baths_around_one_impurity():
    n = 21
    h = np.zeros((n, n), dtype=complex)
    for i in range(n-1):
        h[i, i+1] = -np.exp(0.3j)
        h[i+1, i] = np.conj(h[i, i+1])
    e, u = np.linalg.eigh(h)
    L = u[:, :10].conj() @ u[:, :10].T
    V, info = gbasis.rotation(L, fixed=[10])
    assert [d["modes"] for d in info] == [(0, 10), (11, 21)]
    assert np.max(np.abs(V @ V.conj().T-np.eye(n))) < 1e-13
    assert np.max(np.abs(V[10]-np.eye(n)[10])) == 0.0
    # the two baths are entangled with each other through the impurity,
    # which a rotation of each bath on its own cannot remove: it stays on
    # the few orbitals next to the impurity
    S1 = _free_entropies(V @ L @ V.conj().T)
    far = [k for k in range(n-1) if abs(k+1-10) > 6]  # cut after k+1 modes
    assert np.max(S1[far]) < 1e-6
    assert np.min(_free_entropies(L)[far]) > 0.3


def test_ed_energy_and_observables_are_the_same_in_both_bases():
    n = 4
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = _impurity_hopping(n)
    H = _spinful_hamiltonian(fc, h)
    obs = [fc.Nup[0]*fc.Ndn[0], fc.Sz[0]*fc.Sz[2], fc.Cdagup[1]*fc.Cup[3],
           fc.Ntot[n-1], fc.Cdagdn[0]*fc.Cdn[2]]
    fc.set_hamiltonian(H)
    e0 = fc.gs_energy(mode="ED")
    wf0 = fc.get_gs(mode="ED")
    v0 = np.array([fc.vev(o, wf=wf0, mode="ED") for o in obs])
    b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
    fc.set_hamiltonian(b.rotate(H))
    e1 = fc.gs_energy(mode="ED")
    wf1 = fc.get_gs(mode="ED")
    v1 = np.array([fc.vev(b.rotate(o), wf=wf1, mode="ED") for o in obs])
    assert e1 == pytest.approx(e0, abs=1e-12)
    assert np.max(np.abs(v1-v0)) < 1e-10
    # the impurity operators are not touched by the rotation
    assert np.max(np.abs(v1[0]-fc.vev(obs[0], wf=wf1, mode="ED"))) < 1e-10


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_dmrg_in_the_gaussian_basis_against_ed(itensor_version):
    n = 6
    fc = fermionchain.Spinful_Fermionic_Chain(n,
                                              itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 30
    h = _impurity_hopping(n)
    H = _spinful_hamiltonian(fc, h)
    b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
    fc.set_hamiltonian(H)
    e0 = fc.gs_energy(mode="ED")
    d0 = fc.vev(fc.Nup[0]*fc.Ndn[0], mode="ED").real
    fc.set_hamiltonian(b.rotate(H))
    assert fc.gs_energy() == pytest.approx(e0, abs=1e-8)
    assert fc.vev(fc.Nup[0]*fc.Ndn[0]).real == pytest.approx(d0, abs=1e-6)


def test_spinless_chain_with_the_impurity_at_the_end():
    n = 6
    fc = fermionchain.Fermionic_Chain(n)
    h = _impurity_hopping(n, imp=n-1)
    H = sum(h[i, j]*fc.Cdag[i]*fc.C[j] for i in range(n) for j in range(n)
            if h[i, j] != 0)
    H = H + 1.5*fc.N[n-1]*fc.N[n-2]  # an interaction on the two fixed modes
    fc.set_hamiltonian(H)
    e0 = fc.gs_energy(mode="ED")
    b = fc.get_gaussian_basis(h=h, nf=3, fixed=[n-2, n-1])
    assert b.info[0]["modes"] == (0, n-2)
    fc.set_hamiltonian(b.rotate(H))
    assert fc.gs_energy(mode="ED") == pytest.approx(e0, abs=1e-12)


def test_rotated_hamiltonian_keeps_the_spin():
    """Each spin is rotated on its own, so no term of the rotated
    Hamiltonian moves an electron from one spin to the other."""
    n = 8
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = _impurity_hopping(n)
    b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
    for term in b.rotate(_spinful_hamiltonian(fc, h)).op:
        ladders = [(name, site) for (name, site) in term[1:]
                   if name in ("C", "Cdag")]
        if len(ladders) == 2:
            assert ladders[0][1] % 2 == ladders[1][1] % 2


def test_mixed_reference_is_accepted():
    """A correlation matrix of a mixed state (a thermal one) gives a basis
    too, unitary and with the impurity left as it is."""
    n = 10
    h = _impurity_hopping(n)
    e, u = np.linalg.eigh(h)
    f = 1.0/(np.exp(e/0.2)+1.0)
    cm = (u.conj()*f) @ u.T
    fc = fermionchain.Fermionic_Chain(n)
    b = fc.get_gaussian_basis(cm=cm, fixed=[0])
    assert np.max(np.abs(b.V @ b.V.conj().T-np.eye(n))) < 1e-13
    assert np.max(np.abs(b.V[0]-np.eye(n)[0])) == 0.0


def test_bad_input_raises():
    n = 6
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = _impurity_hopping(n)
    b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
    with pytest.raises(ValueError, match="interaction must act on fixed"):
        b.rotate(fc.Nup[2]*fc.Ndn[2])
    with pytest.raises(ValueError, match="cannot be rotated"):
        b.rotate(fc.A[4])
    with pytest.raises(ValueError, match="h \\(with nf\\) or cm"):
        fc.get_gaussian_basis(fixed=[0])
    with pytest.raises(ValueError, match="out of range"):
        fc.get_gaussian_basis(h=h, nf=n, fixed=[n])
    with pytest.raises(ValueError, match="outside \\[0, 1\\]"):
        fc.get_gaussian_basis(cm=2.0*np.eye(n), fixed=[0])
    sc = fermionchain.Fermionic_Chain(9)
    hs = _impurity_hopping(9)
    with pytest.raises(NotImplementedError, match="between two fixed"):
        sc.get_gaussian_basis(h=hs, nf=4, fixed=[2, 6])


def test_ed_vev_takes_the_state_of_get_gs():
    """get_gs(mode="ED") of a fermion chain returns a Fermionic_State, a
    subclass of the ED State, and vev(wf=) took it for a foreign object
    and failed inside numpy."""
    fc = fermionchain.Fermionic_Chain(4)
    h = _impurity_hopping(4)
    fc.set_hamiltonian(sum(h[i, j]*fc.Cdag[i]*fc.C[j] for i in range(4)
                           for j in range(4) if h[i, j] != 0))
    wf = fc.get_gs(mode="ED")
    assert type(wf).__name__ == "Fermionic_State"
    assert fc.vev(fc.N[1], wf=wf, mode="ED") == pytest.approx(
        fc.vev(fc.N[1], mode="ED"), abs=1e-12)
