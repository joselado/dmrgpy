"""julia_live MPOs: parity-odd fermion operators (a single C or Cdag),
and the algebra of already-built MPOs (product, sum, scalar, dagger,
trace, <a|M|b>).

ITensorMPS's OpSum refuses a parity-odd fermionic term, so
mpsjulialive/mpo.py sends such an operator in the Jordan-Wigner-dressed
form MultiOperator.to_terms() gives (A/Adag/F names, bosonic for OpSum,
see mpsjulialive/fermion_ops.jl). Everything is checked against
mode="ED", and the MPO algebra also against itensor_version="python".
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, timedependent
from _helpers import julia_live_param, setup_backend

BACKENDS = [julia_live_param(), pytest.param("python", id="python")]


def _fermion_chain(itensor_version, n=4):
    """Spinless fermions with first-neighbor hopping and a weak staggered
    onsite energy, so the ground state is not degenerate"""
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    for i in range(n):
        h = h + 0.2*(-1)**i*fc.N[i] + 0.1*i*fc.N[i]
    fc.set_hamiltonian(h)
    setup_backend(fc, itensor_version)
    fc.maxm = 30
    fc.nsweeps = 10
    return fc


def _sp(sc, i):
    """S+ on site i"""
    return sc.Sx[i] + 1j*sc.Sy[i]


def _spin_chain(itensor_version, n=4):
    """Heisenberg chain with a field on the first two sites only, which
    breaks the SU(2) symmetry, so that <S+> is not zero"""
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.6*sc.Sx[0] + 0.3*sc.Sz[1]
    sc.set_hamiltonian(h)
    setup_backend(sc, itensor_version)
    sc.maxm = 30
    sc.nsweeps = 10
    return sc


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_single_fermion_operators_reproduce_ed(itensor_version):
    """<Cdag_i C_j> three ways: the even MultiOperator, the product of the
    two odd MPOs, and the overlap of two odd operators applied to the
    ground state. The pair (2,0) checks the sign of the string."""
    fc = _fermion_chain(itensor_version)
    gs = fc.get_gs()
    for (i, j) in [(0, 2), (2, 0), (1, 2), (0, 3)]:
        ref = fc.vev(fc.Cdag[i]*fc.C[j], mode="ED")
        assert abs(ref) > 1e-2 # a pair that is not zero by symmetry
        assert fc.vev(fc.Cdag[i]*fc.C[j]) == pytest.approx(ref, abs=1e-6)
        ov = (fc.C[i]*gs).dot(fc.C[j]*gs)
        assert ov == pytest.approx(ref, abs=1e-6)
        assert fc.toMPO(fc.Cdag[i]*fc.C[j]).aMb(gs, gs) == pytest.approx(
                ref, abs=1e-6)


def _cpp_param(version):
    from dmrgpy import cppext
    return pytest.param(version, id="v%d" % version, marks=pytest.mark.skipif(
        not cppext.available(version),
        reason="the v%d extension is not compiled" % version))


# toMPO(A)*toMPO(B) is the operator A*B on every backend, as for
# MultiOperators and mode="ED"
ORDER_BACKENDS = [julia_live_param(), pytest.param("python", id="python"),
                  _cpp_param(3), _cpp_param(2)]


@pytest.mark.parametrize("itensor_version", ORDER_BACKENDS)
def test_product_of_odd_mpos_is_ordered(itensor_version):
    """toMPO(Cdag_i)*toMPO(C_j) is Cdag_i C_j. In the reversed order the
    two odd operators anticommute and the sign flips"""
    fc = _fermion_chain(itensor_version)
    gs = fc.get_gs()
    for (i, j) in [(0, 2), (2, 0), (1, 2), (0, 3)]:
        ref = fc.vev(fc.Cdag[i]*fc.C[j], mode="ED")
        M = fc.toMPO(fc.Cdag[i])*fc.toMPO(fc.C[j])
        assert M.aMb(gs, gs) == pytest.approx(ref, abs=1e-6)
    ref = fc.vev(fc.C[0]*fc.C[1], mode="ED")
    M = fc.toMPO(fc.C[0])*fc.toMPO(fc.C[1])
    assert M.aMb(gs, gs) == pytest.approx(ref, abs=1e-6)


@pytest.mark.parametrize("submode", ["KPM", "CVM"])
def test_fermion_dynamical_correlator_julia_matches_ed(submode):
    """The single-particle Green's function on julia_live, which needs the
    odd MPOs of Cdag_0 and C_2, against mode="ED". CVM agrees with ED to
    1e-6. KPM on every DMRG backend differs from the ED KPM by about 2e-4
    here, from its own estimate of the spectral window, so julia_live KPM
    is pinned to "python" KPM to 1e-6 and to ED to 1e-3 (a wrong string
    sign would move it by the height of the peak, 0.5)."""
    pytest.importorskip("juliacall")
    fc = _fermion_chain("julia_live")
    fp = _fermion_chain("python")
    es = np.linspace(-3., 3., 13)
    for (A, B) in [(0, 2), (1, 1)]:
        for (Aname, Bname) in [("Cdag", "C"), ("C", "Cdag")]:
            kw = dict(es=es, delta=0.2, submode=submode)
            name = (getattr(fc, Aname)[A], getattr(fc, Bname)[B])
            x, y = fc.get_dynamical_correlator(mode="DMRG", name=name, **kw)
            xe, ye = fc.get_dynamical_correlator(mode="ED", name=name, **kw)
            assert np.allclose(x, xe)
            assert np.max(np.abs(ye)) > 1e-2
            tol = 1e-6 if submode == "CVM" else 1e-3
            assert np.max(np.abs(np.array(y) - np.array(ye))) < tol
            if submode == "KPM":
                namep = (getattr(fp, Aname)[A], getattr(fp, Bname)[B])
                _xp, yp = fp.get_dynamical_correlator(mode="DMRG",
                                                      name=namep, **kw)
                assert np.max(np.abs(np.array(y) - np.array(yp))) < 1e-6


def test_fermion_evolution_aba_julia_matches_ed():
    """<GS|C_0 e^{iHt} N_2 e^{-iHt} Cdag_0|GS> on julia_live against ED"""
    pytest.importorskip("juliacall")
    fc = _fermion_chain("julia_live")
    kw = dict(nt=40, dt=0.05, A=fc.Cdag[0], B=fc.N[2])
    _t, y = timedependent.evolution_ABA(fc, mode="DMRG", **kw)
    _te, ye = timedependent.evolution_ABA(fc, mode="ED", **kw)
    assert np.max(np.abs(np.array(y) - np.array(ye))) < 1e-5


@pytest.mark.parametrize("itensor_version", ORDER_BACKENDS)
def test_mpo_product_is_operator_product(itensor_version):
    """toMPO(A)*toMPO(B) is the operator A*B, B acting first: Sp and Sz
    on the same site do not commute, so the order shows"""
    sc = _spin_chain(itensor_version)
    gs = sc.get_gs()
    for (A, B) in [(_sp(sc, 0), sc.Sz[0]), (sc.Sz[0], _sp(sc, 0)),
                   (sc.Sx[1], sc.Sz[2]*sc.Sy[3])]:
        ref = sc.vev(A*B, mode="ED")
        M = sc.toMPO(A)*sc.toMPO(B)
        assert M.aMb(gs, gs) == pytest.approx(ref, abs=1e-6)
        assert (sc.toMPO(A)*B).aMb(gs, gs) == pytest.approx(ref, abs=1e-6)
        assert (A*sc.toMPO(B)).aMb(gs, gs) == pytest.approx(ref, abs=1e-6)
    assert abs(sc.vev(_sp(sc, 0)*sc.Sz[0], mode="ED")
               - sc.vev(sc.Sz[0]*_sp(sc, 0), mode="ED")) > 1e-2


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_mpo_dagger_sum_and_scalars(itensor_version):
    sc = _spin_chain(itensor_version)
    gs = sc.get_gs()
    w = (sc.Sx[0] + 0.5*sc.Sz[1])*gs # a second state, not the ground state
    A = _sp(sc, 0)*sc.Sz[1]
    B = sc.Sx[2] + 0.3j*sc.Sy[3]
    M, N = sc.toMPO(A), sc.toMPO(B)
    assert abs(M.aMb(gs, w)) > 1e-2
    # dagger: <a|M^dagger|b> = conj(<b|M|a>)
    assert M.get_dagger().aMb(gs, w) == pytest.approx(
            np.conj(M.aMb(w, gs)), abs=1e-6)
    assert M.get_dagger().aMb(gs, gs) == pytest.approx(
            np.conj(M.aMb(gs, gs)), abs=1e-6)
    assert M.get_dagger().aMb(gs, w) == pytest.approx(
            sc.toMPO(A.get_dagger()).aMb(gs, w), abs=1e-6)
    # sums, differences and scalars against the MultiOperator algebra
    for (got, op) in [(M + N, A + B), (M - N, A - B), (2*M, 2*A),
                      (M*(1j), 1j*A), (M/2, A/2), (-M, -A)]:
        assert got.aMb(gs, w) == pytest.approx(
                sc.toMPO(op).aMb(gs, w), abs=1e-6)
    # aMb is the overlap with the applied state
    assert M.aMb(gs, w) == pytest.approx(gs.dot(M*w), abs=1e-6)


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_mpo_trace(itensor_version):
    """Full traces: Tr(Sz_0 Sz_0) = 2^L/4 on S=1/2, and Tr(H) against the
    ED matrix, on a spin chain and on a fermion chain"""
    L = 4
    sc = _spin_chain(itensor_version, n=L)
    assert sc.trace(sc.Sz[0]*sc.Sz[0]) == pytest.approx(2**(L-2), abs=1e-10)
    for chain in (sc, _fermion_chain(itensor_version, n=L)):
        H = chain.hamiltonian + 0.7*chain.hamiltonian*chain.hamiltonian
        ref = np.trace(chain.toMPO(H, mode="ED").SO.toarray()
                       if hasattr(chain.toMPO(H, mode="ED").SO, "toarray")
                       else chain.toMPO(H, mode="ED").SO)
        assert chain.trace(H) == pytest.approx(ref, abs=1e-8)
        assert chain.toMPO(H).get_dagger().trace() == pytest.approx(
                np.conj(ref), abs=1e-8)


@pytest.mark.parametrize("itensor_version", BACKENDS)
def test_mpo_product_is_cut_at_mpomaxm_not_maxm(itensor_version):
    """The MPO product is truncated at the operator bond dimension
    mpomaxm, not at the wavefunction maxm. H^2 on a 6-site Heisenberg
    chain has bond dimension up to 25, so a cut at maxm=2 would move
    <H^2> far away from ED."""
    sc = spinchain.Spin_Chain(["1/2"]*6)
    h = 0
    for i in range(5):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    setup_backend(sc, itensor_version)
    sc.maxm = 30
    sc.nsweeps = 10
    gs = sc.get_gs()  # exact at maxm=30 on 6 sites
    ref = sc.vev(h*h, mode="ED")
    sc.maxm = 2  # mpomaxm stays at its default
    M = sc.toMPO(h)
    assert (M*M).aMb(gs, gs) == pytest.approx(ref, abs=1e-8)
