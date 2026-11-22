"""Fermionic Gaussian states as MPS (gaussianmps.py, the Fishman-White
circuit of arXiv:1504.07701).

A state whose correlation matrix <Cdag_i C_j> has eigenvalues 0 and 1 is
exactly one Slater determinant, so comparing the whole complex matrix
against the requested one is a complete check of the state. A complex
hopping is used throughout, since a real one passes with a transposed or
conjugated matrix."""

import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy.gaussianmpstk import circuit
from dmrgpy.gaussianmpstk.gates import SiteArrays

from _helpers import julia_live_param, setup_backend


def _cpp_param(version):
    from dmrgpy import cppext
    return pytest.param(version, id="v%d" % version, marks=pytest.mark.skipif(
        not cppext.available(version),
        reason="the v%d extension is not compiled" % version))


BACKENDS = [_cpp_param(3), pytest.param("python", id="python"),
            julia_live_param()]


def _hopping(n, seed=0, phase=0.4):
    """A first-neighbor chain with a Peierls phase, random onsite energies
    and a weak complex third-neighbor hopping."""
    rng = np.random.default_rng(seed)
    h = np.diag(rng.uniform(-0.5, 0.5, n)).astype(complex)
    for i in range(n-1):
        h[i, i+1] = -np.exp(1j*phase)
    for i in range(n-3):
        h[i, i+3] = 0.2j
    return h+np.triu(h, 1).conj().T


def _correlations(fc, wf, mode=None):
    n = fc.ns
    kw = {} if mode is None else {"mode": mode}
    return np.array([[fc.vev(fc.Cdag[i]*fc.C[j], wf=wf, **kw)
                      for j in range(n)] for i in range(n)])


def _dense_slater(phi):
    """The Slater determinant as a dense vector, with site 0 leftmost and
    the local basis (empty, occupied)."""
    n = phi.shape[0]
    a, Z, I2 = np.array([[0, 1], [0, 0]]), np.diag([1, -1]), np.eye(2)
    cs = []
    for i in range(n):
        m = np.array([[1.0]])
        for k in range(n):
            m = np.kron(m, Z if k < i else (a if k == i else I2))
        cs.append(m)
    psi = np.zeros(2**n, dtype=complex)
    psi[0] = 1.0
    for k in range(phi.shape[1]):
        psi = sum(phi[i, k]*cs[i].T for i in range(n)) @ psi
    return psi/np.linalg.norm(psi)


@pytest.mark.parametrize("n,nf,maxblocksize", [(8, 3, None), (9, 5, 4),
                                               (6, 0, None), (6, 6, None)])
def test_circuit_is_the_slater_determinant(n, nf, maxblocksize):
    """The rotations diagonalize Lambda, and the gates on the product state
    give the dense Slater determinant, with nothing truncated."""
    rng = np.random.default_rng(n+nf)
    h = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    e, u = np.linalg.eigh(h+h.conj().T)
    phi = u[:, :nf]
    L = phi.conj() @ phi.T
    ns, rots, err = circuit.correlation_matrix_to_circuit(
            L, eigval_cutoff=1e-12, maxblocksize=maxblocksize)
    V = np.eye(n, dtype=complex)
    for (k, c, s) in rots:
        g = np.eye(n, dtype=complex)
        g[k:k+2, k:k+2] = [[c, s], [-np.conj(s), c]]
        V = g @ V
    D = V @ L @ V.conj().T
    assert np.max(np.abs(D-np.diag(np.diag(D)))) < max(1e-10, 10*err)
    st = SiteArrays.product([[1, 0] if round(x) == 0 else [0, 1] for x in ns])
    for (k, c, s) in reversed(rots):
        st.apply(circuit.gate_matrix(c, s), k, 10**6, 0.0)
    v = st.A[0]
    for A in st.A[1:]:
        v = np.tensordot(v, A, axes=(v.ndim-1, 0))
    ov = abs(np.vdot(_dense_slater(phi), v.reshape(-1)))
    assert ov == pytest.approx(1.0, abs=max(1e-10, 10*err))


def test_ed_route_has_the_correlation_matrix():
    n, nf = 7, 3
    fc = fermionchain.Fermionic_Chain(n)
    fc.set_hamiltonian(sum(fc.N))
    h = _hopping(n)
    e, u = np.linalg.eigh(h)
    cm = u[:, :nf].conj() @ u[:, :nf].T
    wf = fc.get_gaussian_state(h=h, nf=nf, mode="ED")
    assert np.max(np.abs(_correlations(fc, wf, mode="ED")-cm)) < 1e-12
    wf2 = fc.get_gaussian_state(cm=cm, mode="ED")
    assert abs(wf.dot(wf2)) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_mps_has_the_correlation_matrix(itensor_version):
    n, nf = 10, 4
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, itensor_version)
    fc.maxm = 64
    h = _hopping(n)
    e, u = np.linalg.eigh(h)
    cm = u[:, :nf].conj() @ u[:, :nf].T
    H = sum(h[i, j]*fc.Cdag[i]*fc.C[j] for i in range(n) for j in range(n)
            if h[i, j] != 0)
    fc.set_hamiltonian(H)
    wf = fc.get_gaussian_state(h=h, nf=nf, eigval_cutoff=1e-12)
    assert abs(wf.dot(wf)) == pytest.approx(1.0, abs=1e-10)
    assert np.max(np.abs(_correlations(fc, wf)-cm)) < 1e-8
    assert fc.vev(H, wf=wf).real == pytest.approx(np.sum(e[:nf]), abs=1e-8)


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_spinful_gutzwiller_against_ed(itensor_version):
    """The Gutzwiller-projected Fermi sea of a spinful chain: the same
    correlations on DMRG as on ED, and no double occupancy."""
    n, nf = 4, 4
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    setup_backend(fc, itensor_version)
    fc.maxm = 64
    fc.set_hamiltonian(sum(fc.N))
    h = _hopping(n, seed=3)
    wed = fc.get_gaussian_state(h=h, nf=nf, gutzwiller=True, mode="ED")
    wf = fc.get_gaussian_state(h=h, nf=nf, gutzwiller=True)
    ref = _correlations(fc, wed, mode="ED")
    assert np.max(np.abs(_correlations(fc, wf)-ref)) < 1e-8
    docc = sum(fc.Nup[i]*fc.Ndn[i] for i in range(n))
    assert abs(fc.vev(docc, wf=wf)) < 1e-10
    assert abs(fc.vev(sum(fc.N), wf=wf)-nf) < 1e-10


@pytest.mark.parametrize("itensor_version", [_cpp_param(3),
                                             pytest.param("python",
                                                          id="python")])
def test_dmrg_starts_from_the_gaussian_state(itensor_version):
    """A Hubbard chain solved from the U=0 Fermi sea gives the ED energy."""
    n, nf = 4, 4
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    setup_backend(fc, itensor_version)
    fc.maxm = 40
    fc.nsweeps = 6
    h = _hopping(n, seed=1)
    H = sum(h[i, j]*(fc.Cdagup[i]*fc.Cup[j]+fc.Cdagdn[i]*fc.Cdn[j])
            for i in range(n) for j in range(n) if h[i, j] != 0)
    H = H + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(n))
    H = H + 5.0*(sum(fc.N)-nf)*(sum(fc.N)-nf)
    fc.set_hamiltonian(H)
    e_ed = fc.gs_energy(mode="ED")
    wf = fc.get_gaussian_state(h=h, nf=nf)
    fc.set_initial_wf(wf)  # taken as it is: its own energy
    assert fc.gs_energy() == pytest.approx(fc.vev(H, wf=wf).real, abs=1e-10)
    fc.set_initial_wf_guess(wf)  # swept from it
    assert fc.gs_energy() == pytest.approx(e_ed, abs=1e-6)


def test_truncation_warns():
    n, nf = 12, 6
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.set_hamiltonian(sum(fc.N))
    with pytest.warns(UserWarning, match="discarded"):
        fc.get_gaussian_state(h=_hopping(n), nf=nf, maxm=2)


def test_bad_input_raises():
    fc = fermionchain.Fermionic_Chain(6)
    fc.setup_python()
    fc.set_hamiltonian(sum(fc.N))
    h = _hopping(6)
    with pytest.raises(ValueError, match="nf"):
        fc.get_gaussian_state(h=h, nf=7)
    with pytest.raises(ValueError, match="either h"):
        fc.get_gaussian_state(h=h, nf=2, cm=np.eye(6))
    with pytest.raises(ValueError, match="pure Gaussian"):
        fc.get_gaussian_state(cm=0.5*np.eye(6))
    with pytest.raises(ValueError, match="Hermitian"):
        fc.get_gaussian_state(h=h+np.triu(np.ones((6, 6)), 1), nf=2)
    with pytest.raises(ValueError, match="Spinful"):
        fc.get_gaussian_state(h=h, nf=2, gutzwiller=True)
    with pytest.raises(TypeError):
        fc.get_gaussian_state(h=h, nf=2, eigvalcutoff=1e-9)
    fc.set_conserved_sector(Nf=2)
    with pytest.raises(ValueError, match="conserved-sector"):
        fc.get_gaussian_state(h=h, nf=2)
    fc.set_conserved_sector()
    fc.hamiltonian = None
    with pytest.raises(ValueError, match="set_hamiltonian"):
        fc.get_gaussian_state(h=h, nf=2, mode="ED")


def test_degenerate_filling_warns():
    n = 6
    h = np.diag(np.arange(n, dtype=float))
    h[2, 2] = h[3, 3] = 1.0  # levels 0, 1, 1, 1, 4, 5
    fc = fermionchain.Fermionic_Chain(n)
    fc.set_hamiltonian(sum(fc.N))
    with pytest.warns(UserWarning, match="degenerate"):
        fc.get_gaussian_state(h=h, nf=2, mode="ED")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fc.get_gaussian_state(h=h, nf=4, mode="ED")


def test_ed_vev_reads_an_explicit_state():
    """vev(op, wf=x, mode="ED") is <x|op|x>, as on DMRG; it used to drop
    wf= and answer for the ground state."""
    fc = fermionchain.Fermionic_Chain(4)
    fc.set_hamiltonian(sum(fc.N))  # ground state: the vacuum
    h = _hopping(4)
    wf = fc.get_gaussian_state(h=h, nf=2, mode="ED")
    assert fc.vev(sum(fc.N), wf=wf, mode="ED").real == pytest.approx(2.0,
                                                                    abs=1e-12)
    assert fc.vev(sum(fc.N), mode="ED").real == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(TypeError, match="ED state"):
        fc.vev(sum(fc.N), wf=object(), mode="ED")


@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_gutzwiller_of_a_full_band_raises(mode):
    """Two electrons in every orbital leave nothing after the projection."""
    fc = fermionchain.Spinful_Fermionic_Chain(3)
    fc.setup_python()
    fc.set_hamiltonian(sum(fc.N))
    with pytest.raises(ValueError, match="annihilates"):
        fc.get_gaussian_state(h=_hopping(3), nf=6, gutzwiller=True, mode=mode)


def _kitaev(n, mu=0.4, delta=0.7, seed=5):
    """A Kitaev chain with a complex hopping, a complex pairing and a weak
    pairing to the second neighbor, in the trivial phase |mu| > 0 or the
    topological one: the bulk gap is what matters here, with open ends
    the zero modes of the topological phase are split by n."""
    rng = np.random.default_rng(seed)
    h = np.diag(rng.uniform(-0.3, 0.3, n)-mu).astype(complex)
    D = np.zeros((n, n), dtype=complex)
    for i in range(n-1):
        h[i, i+1] = -np.exp(0.3j)
        h[i+1, i] = np.conj(h[i, i+1])
        D[i, i+1] = delta*np.exp(0.5j*i)
        D[i+1, i] = -D[i, i+1]
    for i in range(n-2):
        D[i, i+2] = 0.15
        D[i+2, i] = -0.15
    return h, D


def _generalized(fc, wf, mode=None):
    """<gamma_a^dag gamma_b> on gamma = (C_0, Cdag_0, C_1, ...)."""
    kw = {} if mode is None else {"mode": mode}
    gam, gamd = [], []
    for i in range(fc.ns):
        gam += [fc.C[i], fc.Cdag[i]]
        gamd += [fc.Cdag[i], fc.C[i]]
    return np.array([[fc.vev(gamd[a]*gam[b], wf=wf, **kw)
                      for b in range(2*fc.ns)] for a in range(2*fc.ns)])


def test_bdg_circuit_is_the_dense_ground_state():
    """The BdG gates on the product state give the ground state of the
    dense Jordan-Wigner BdG Hamiltonian, with no correlation matrix in the
    reference."""
    n = 6
    h, D = _kitaev(n)
    fc = fermionchain.Fermionic_Chain(n)
    from dmrgpy import gaussianmps
    G = gaussianmps.bdg_correlation_matrix(fc, h, D)
    ns, gates, err = circuit.bdg_correlation_matrix_to_circuit(
            G, eigval_cutoff=1e-12)
    assert any(g[0] == "pair" for g in gates)
    st = SiteArrays.product([[1, 0] if round(x) == 0 else [0, 1] for x in ns])
    for (kind, k, c, s) in reversed(gates):
        M = (circuit.gate_matrix(c, s) if kind == "hop"
             else circuit.pair_gate_matrix(c, s))
        st.apply(M, k, 10**6, 0.0)
    v = st.A[0]
    for A in st.A[1:]:
        v = np.tensordot(v, A, axes=(v.ndim-1, 0))
    a, Z, I2 = np.array([[0, 1], [0, 0]]), np.diag([1, -1]), np.eye(2)
    cs = []
    for i in range(n):
        m = np.array([[1.0]])
        for k in range(n):
            m = np.kron(m, Z if k < i else (a if k == i else I2))
        cs.append(m)
    H = sum(h[i, j]*cs[i].T @ cs[j] for i in range(n) for j in range(n))
    H = H + 0.5*sum(D[i, j]*cs[i].T @ cs[j].T
                    + np.conj(D[i, j])*cs[j] @ cs[i]
                    for i in range(n) for j in range(n))
    psi = np.linalg.eigh(H)[1][:, 0]
    assert abs(np.vdot(psi, v.reshape(-1))) == pytest.approx(1.0, abs=1e-10)


def test_bdg_ed_route_has_the_generalized_correlation_matrix():
    n = 5
    h, D = _kitaev(n)
    fc = fermionchain.Fermionic_Chain(n)
    fc.set_hamiltonian(sum(fc.N))
    from dmrgpy import gaussianmps
    G = gaussianmps.bdg_correlation_matrix(fc, h, D)
    wf = fc.get_gaussian_state(h=h, delta=D, mode="ED")
    assert np.max(np.abs(_generalized(fc, wf, mode="ED")-G)) < 1e-10
    wf2 = fc.get_gaussian_state(cm=G, mode="ED")
    assert abs(wf.dot(wf2)) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_bdg_mps_has_the_generalized_correlation_matrix(itensor_version):
    n = 8
    h, D = _kitaev(n)
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, itensor_version)
    fc.maxm = 64
    fc.set_hamiltonian(sum(fc.N))
    from dmrgpy import gaussianmps
    G = gaussianmps.bdg_correlation_matrix(fc, h, D)
    wf = fc.get_gaussian_state(h=h, delta=D, eigval_cutoff=1e-12)
    assert abs(wf.dot(wf)) == pytest.approx(1.0, abs=1e-10)
    assert np.max(np.abs(_generalized(fc, wf)-G)) < 1e-8


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_spinful_swave_gutzwiller_against_ed(itensor_version):
    """An s-wave BdG state of a spinful chain, projected: the same
    generalized correlations on DMRG as on ED, and no double occupancy."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    setup_backend(fc, itensor_version)
    fc.maxm = 64
    fc.set_hamiltonian(sum(fc.N))
    h = _hopping(n, seed=4)
    delta = 0.6*np.eye(n)+0.2*(np.eye(n, k=1)+np.eye(n, k=-1))
    for g in (False, True):
        wed = fc.get_gaussian_state(h=h, delta=delta, gutzwiller=g,
                                    mode="ED")
        wf = fc.get_gaussian_state(h=h, delta=delta, gutzwiller=g)
        ref = _generalized(fc, wed, mode="ED")
        assert np.max(np.abs(ref[0::2, 1::2])) > 1e-2  # pairing is there
        assert np.max(np.abs(_generalized(fc, wf)-ref)) < 1e-8
    docc = sum(fc.Nup[i]*fc.Ndn[i] for i in range(n))
    assert abs(fc.vev(docc, wf=wf)) < 1e-10


def test_bdg_bad_input_raises():
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.set_hamiltonian(sum(fc.N))
    h, D = _kitaev(n)
    with pytest.raises(ValueError, match="antisymmetric"):
        fc.get_gaussian_state(h=h, delta=np.ones((n, n)))
    with pytest.raises(ValueError, match="grand canonical"):
        fc.get_gaussian_state(h=h, delta=D, nf=2)
    G = np.eye(2*n)*0.5
    with pytest.raises(ValueError, match="pure Gaussian"):
        fc.get_gaussian_state(cm=G)
    sf = fermionchain.Spinful_Fermionic_Chain(2)
    with pytest.raises(ValueError, match="symmetric"):
        sf.get_gaussian_state(h=np.eye(2), delta=np.array([[0, 1], [0, 0]]))


@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_bdg_exact_zero_modes_give_a_ground_state(mode):
    """At the ideal Kitaev point, mu=0 and t=Delta with open ends, the two
    end Majoranas give an exact zero mode: the state is one of the two
    ground states, with the ground-state energy, and a warning says so."""
    n = 6
    h = -np.eye(n, k=1)-np.eye(n, k=-1)
    D = np.eye(n, k=1)-np.eye(n, k=-1)
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    H = sum(h[i, j]*fc.Cdag[i]*fc.C[j] for i in range(n) for j in range(n)
            if h[i, j] != 0)
    H = H + 0.5*sum(D[i, j]*fc.Cdag[i]*fc.Cdag[j]
                    + np.conj(D[i, j])*fc.C[j]*fc.C[i]
                    for i in range(n) for j in range(n) if D[i, j] != 0)
    fc.set_hamiltonian(H)
    e0 = fc.gs_energy(mode="ED")
    with pytest.warns(UserWarning, match="zero"):
        wf = fc.get_gaussian_state(h=h, delta=D, mode=mode)
    kw = {"mode": "ED"} if mode == "ED" else {}
    assert fc.vev(H, wf=wf, **kw).real == pytest.approx(e0, abs=1e-8)
