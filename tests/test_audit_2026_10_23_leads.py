"""Regressions for the reviewed leads of docs/audit_2026_10_23_hole_hunt.md
("The leads, reviewed"), one test group per finding."""
import numpy as np
import pytest

from dmrgpy import spinchain, entropy

BACKENDS = [2, 3, "python"]


def xxz(n, v, jz=0.7):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + jz*sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    sc.maxm = 40
    sc.nsweeps = 10
    return sc


def exact_entropy(v, n, keep):
    """Von Neumann entropy of the sites `keep` of the ED vector v"""
    t = np.asarray(v).ravel().reshape([2]*n)
    rest = [k for k in range(n) if k not in keep]
    t = np.transpose(t, list(keep)+rest).reshape(2**len(keep), -1)
    p = np.linalg.svd(t, compute_uv=False)**2
    p = p/p.sum()
    p = p[p > 1e-300]
    return -np.sum(p*np.log(p))


# Finding 6: the pair entropy and the mutual information did not normalize

@pytest.mark.parametrize("v", BACKENDS)
@pytest.mark.parametrize("mode", ["DMRG", "ED"])
def test_pair_entropy_and_mutual_information_are_scale_invariant(v, mode):
    sc = xxz(6, v)
    wf = sc.get_gs(mode=mode)
    for c in (2.0, 1.0002):  # 1.0002 used to pass the trace check, 1.1e-4 off
        w = c*wf
        assert sc.get_pair_entropy(w, 1, 3) == pytest.approx(
            sc.get_pair_entropy(wf, 1, 3), abs=1e-10)
        assert sc.get_mutual_information(w, 1, 3) == pytest.approx(
            sc.get_mutual_information(wf, 1, 3), abs=1e-10)
        assert sc.get_site_entropy(w, 1) == pytest.approx(
            sc.get_site_entropy(wf, 1), abs=1e-10)


def test_entropy_of_a_zero_state_raises():
    sc = xxz(4, "python")
    wf = sc.get_gs()
    with pytest.raises(ValueError, match="zero state"):
        sc.get_pair_entropy(0.*wf, 1, 2)
    with pytest.raises(ValueError, match="zero state"):
        sc.get_site_entropy(0.*wf, 1)


# Finding 7: the absolute 1e-6 eigenvalue floor of entropy_dm

@pytest.mark.parametrize("p", [5.6e-7, 1e-6, 1e-9])
def test_entropy_dm_keeps_small_weights(p):
    dm = np.diag([1-p, p]).astype(complex)
    exact = -(1-p)*np.log(1-p) - p*np.log(p)
    assert entropy.entropy_dm(dm) == pytest.approx(exact, rel=1e-10)


def test_entropy_dm_ignores_roundoff_eigenvalues():
    dm = np.diag([0.5, 0.5, 1e-18, -1e-17]).astype(complex)
    assert entropy.entropy_dm(dm) == pytest.approx(np.log(2), abs=1e-14)


def test_site_entropy_of_a_nearly_product_state_against_ed():
    # staggered field and a weak XY coupling: site 0 has a Schmidt weight
    # of 5.6e-7, and the site entropy was 5.6e-7 for an exact 8.66e-6
    n, eps = 4, 0.003
    sc = spinchain.Spin_Chain([2]*n, itensor_version=3)
    h = sum(eps*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]) for i in range(n-1))
    h = h + sum((-1)**i*sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h)
    wf = sc.get_gs(mode="ED")
    v = np.asarray(wf.v).ravel()
    assert sc.get_site_entropy(wf, 0) == pytest.approx(
        exact_entropy(v, n, [0]), rel=1e-6)
    assert sc.get_pair_entropy(wf, 1, 2) == pytest.approx(
        exact_entropy(v, n, [1, 2]), rel=1e-6)


def test_compute_entropy_with_a_geometric_mean_is_gone():
    assert not hasattr(entropy, "compute_entropy")


# Finding 8: an unconverged NH-DMRG state was reused silently

def pt_chain(n, g, v):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] for i in range(n-1))
    h = h + sum(1j*g*(-1)**i*sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("v", BACKENDS)
def test_unconverged_nhdmrg_state_warns_on_every_reuse(v, capsys):
    sc = pt_chain(4, 0.3, v)
    sc.gs_energy(tol=1e-30, ntries=1)  # a tolerance no solve can reach
    assert "did not reach the residual tolerance" in capsys.readouterr().out
    for read in (sc.gs_energy, sc.get_gs, lambda: sc.vev(sc.Sz[0])):
        read()
        assert "did not converge" in capsys.readouterr().out
    sc.set_gs(sc.wf0.copy())  # a state the caller sets is not the solve's
    sc.vev(sc.Sz[0])
    assert "did not converge" not in capsys.readouterr().out


@pytest.mark.parametrize("v", BACKENDS)
def test_converged_nhdmrg_state_is_reused_quietly(v, capsys):
    sc = pt_chain(4, 0.3, v)
    sc.gs_energy()
    capsys.readouterr()
    sc.gs_energy()
    sc.vev(sc.Sz[0])
    assert "converge" not in capsys.readouterr().out


# Finding 9: gs_energy_generalized with a metric that is not Hermitian
# positive definite

def heisenberg(n, v):
    return xxz(n, v, jz=1.0)


@pytest.mark.parametrize("v", [3, "python"])
def test_generalized_refuses_a_non_hermitian_metric(v):
    sc = heisenberg(6, v)
    with pytest.raises(ValueError, match="not Hermitian"):
        sc.gs_energy_generalized(1 + 1j*sc.Sz[0])


@pytest.mark.parametrize("v", [3, "python"])
@pytest.mark.parametrize("which", ["negative", "singular"])
def test_generalized_refuses_a_metric_that_is_not_positive(v, which):
    # A=-1 returned the largest eigenvalue, 2.4936, and on v3 the singular
    # A=Sz_0+1/2 aborted the interpreter
    sc = heisenberg(6, v)
    A = -1 + 0*sc.Sz[0] if which == "negative" else sc.Sz[0] + 0.5
    with pytest.raises(RuntimeError, match="not positive"):
        sc.gs_energy_generalized(A)


@pytest.mark.parametrize("v", [3, "python"])
def test_generalized_with_a_positive_metric_matches_ed(v):
    import scipy.linalg
    sc = heisenberg(6, v)
    sc.nsweeps = 20
    A = 1 + 0.5*sc.Sz[0]
    ed = sc.get_ED_obj()
    Hm = np.asarray(ed.MO2matrix(sc.hamiltonian).todense())
    Am = np.asarray(ed.MO2matrix(A).todense())
    exact = np.min(scipy.linalg.eigh(Hm, Am, eigvals_only=True))
    assert sc.gs_energy_generalized(A) == pytest.approx(exact, abs=1e-6)


# Finding 10: METTS with a diagonal basis whose own sum H does not conserve

def hubbard_native(n, v):
    from dmrgpy import fermionchain
    sc = fermionchain.Spinful_Fermionic_Chain_Native(n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h - (sc.Cdagup[i]*sc.Cup[i+1] + sc.Cdagdn[i]*sc.Cdn[i+1])
        h = h - (sc.Cdagup[i+1]*sc.Cup[i] + sc.Cdagdn[i+1]*sc.Cdn[i])
    h = h + sum(4.0*sc.Nup[i]*sc.Ndn[i] for i in range(n))
    sc.set_hamiltonian(h)
    return sc


def bose_hubbard(n, v):
    from dmrgpy import bosonchain
    sc = bosonchain.Bosonic_Chain(n, maxnb=[4]*n, itensor_version=v)
    h = sum(sc.Adag[i]*sc.A[i+1] + sc.Adag[i+1]*sc.A[i] for i in range(n-1))
    h = h + sum(2.0*sc.N[i]*sc.N[i] for i in range(n))
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("v", [3, "python"])
@pytest.mark.parametrize("chain,basis", [(hubbard_native, ("Nupdn",)),
                                         (hubbard_native, ("F",)),
                                         (hubbard_native, ("Nupdn", "Sz")),
                                         (bose_hubbard, ("N2",)),
                                         (bose_hubbard, ("N0", "N1"))])
def test_metts_refuses_a_diagonal_basis_on_a_charge_conserving_chain(v, chain, basis):
    # <Ntot> came out 6, 4, 2 by seed with stderr 0 for an exact 2.74
    sc = chain(4, v)
    ntot = sum(sc.Ntot) if hasattr(sc, "Ntot") else sum(sc.N)
    with pytest.raises(ValueError, match="conserves"):
        sc.metts_vev(ntot, T=1.0, nsamples=5, nwarmup=2, seed=1,
                     basis_ops=basis)


# Finding 11: "python" DMRG kept the sector of its first local solve

def two_spins_and_bosons(v, hop=1.0, nb=4):
    from dmrgpy import bosonchain
    b = bosonchain.SpinBoson_Chain(["S=1/2"]*2 + ["B"]*nb, itensor_version=v)
    h = b.Sx[0]*b.Sx[1] + b.Sy[0]*b.Sy[1] + b.Sz[0]*b.Sz[1]
    h = h + (b.Sx[0] + b.Sx[1])*(b.A[2] + b.Adag[2])
    h = h + sum((b.N[i] - 1.0)*(b.N[i] - 1.0) for i in range(2, 2+nb))
    h = h + sum(hop*b.Adag[i]*b.A[i+1] for i in range(2, 1+nb))
    b.set_hamiltonian(h + h.get_dagger())
    b.maxm = 20
    b.nsweeps = 10
    return b


def dense_two_spins_and_bosons(hop=1.0, nb=4):
    """The same Hamiltonian as a dense matrix, built from scratch"""
    sx = np.array([[0, 1], [1, 0]])/2
    sy = np.array([[0, -1j], [1j, 0]])/2
    sz = np.diag([0.5, -0.5])
    a = np.diag(np.sqrt(np.arange(1, 4)), 1)
    n = np.diag(np.arange(4.))
    dims = [2, 2] + [4]*nb
    D = int(np.prod(dims))

    def emb(o, k):
        m = np.eye(1)
        for j, dj in enumerate(dims):
            m = np.kron(m, o if j == k else np.eye(dj))
        return m
    H = sum(emb(o, 0) @ emb(o, 1) for o in (sx, sy, sz))
    H = H + (emb(sx, 0) + emb(sx, 1)) @ (emb(a, 2) + emb(a.T, 2))
    H = H + sum((emb(n, i) - np.eye(D)) @ (emb(n, i) - np.eye(D)) for i in range(2, 2+nb))
    H = H + sum(hop*emb(a.T, i) @ emb(a, i+1) for i in range(2, 1+nb))
    return np.linalg.eigvalsh(H + H.conj().T)[0]


@pytest.mark.parametrize("hop", [1.0, 0.0])
def test_python_dmrg_leaves_the_sector_of_its_first_solve(hop):
    # two equivalent spins: -4.470573 for an exact -5.917370 in 11 of 12 runs
    exact = dense_two_spins_and_bosons(hop=hop, nb=3)
    for _ in range(4):
        sc = two_spins_and_bosons("python", hop=hop, nb=3)
        assert sc.gs_energy() == pytest.approx(exact, abs=1e-5)
