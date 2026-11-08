"""Chain-level MPS algebra on itensor_version="julia_live": overlap, aMb,
applyoperator, summps, scale_mps, conjugate, random_state(orthogonal=),
the delegating MPS methods, exponential() and applyinverse(), and the
three correlators built on applyinverse (CVM_explicit, the analytic
continuation behind CVMimag, and the non-Hermitian CVM).

Every check is against mode="ED" on the same chain, and against
itensor_version="python" where the two MPS backends run the same
algorithm (the Taylor steps of exponential, the BiCGSTAB of
applyinverse). julia_live DMRG starts from a random MPS, so only
quantities that do not depend on the global phase of the ground state
are compared across backends."""

import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain
from dmrgpy import cvm
from dmrgpy.nonhermitian import dynamics as nhdynamics

from _helpers import julia_live_param, julia_available

JULIA = [julia_live_param()]


def _heisenberg(n, version, maxm=40, extra=0, field=0.):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    if extra: h = h + extra*sc.Sz[0]
    if field: h = h + field*sum(sc.Sz)
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = 12
    # a cutoff far below the tolerances checked here, so that the sums and
    # MPO applications inside the solvers do not set the floor
    sc.cutoff = 1e-16
    return sc, h


def _fermions(n, version, maxm=40):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
        h = h + 0.5*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    h = h + 0.2*fc.N[0]
    fc.set_hamiltonian(h)
    fc.maxm = maxm
    fc.nsweeps = 12
    fc.cutoff = 1e-16
    return fc, h


def _norm(c, w):
    return np.sqrt(abs(c.overlap(w, w).real))


def _algebra_numbers(c, h, ops, mode=None):
    """Phase-independent numbers of the chain-level primitives on the ground
    state of c. mode="ED" computes the same numbers with the operations of
    an ED state, as the reference (the chain-level ED calls rewrap their
    argument as a spin State, which a Fermionic_State does not survive)."""
    A, B, C = ops
    if mode is None:
        gs = c.get_gs()
        ov, amb, app = c.overlap, c.aMb, c.applyoperator
        add, scale = c.summps, c.scale_mps
        amb_wf = lambda a, M, b: a.aMb(M, b)
    else:
        gs = c.get_gs(mode=mode)
        ov = lambda a, b: a.dot(b)
        amb = lambda a, M, b: a.dot(M*b)
        amb_wf = amb
        app = lambda M, w: M*w
        add = lambda a, b: a + b
        scale = lambda x, w: x*w
    nrm = lambda w: np.sqrt(abs(ov(w, w).real))
    out = {}
    out["e0"] = amb(gs, h, gs).real
    out["AB"] = amb(gs, A*B, gs)
    out["wfAB"] = amb_wf(gs, A*B, gs)
    out["norm1"] = ov(gs, gs).real
    out["normsum"] = nrm(add(app(A, gs), app(B, gs)))
    out["scale"] = nrm(scale(2.+1j, gs))
    w = add(gs, scale(1j, app(C, gs)))
    wc = w.get_conjugate()
    out["ww"] = ov(w, w).real
    out["conjw"] = abs(ov(wc, w)) # bilinear w.w, phase^2 dropped by abs
    out["conjconj"] = ov(wc, wc).real
    return out


@pytest.mark.parametrize("version", JULIA)
@pytest.mark.parametrize("model", ["spin", "fermion"])
def test_chain_algebra_against_ed_and_python(version, model):
    n = 6
    build = _heisenberg if model == "spin" else _fermions
    cj, hj = build(n, version)
    cp, hp = build(n, "python")
    if model == "spin":
        opsj = (cj.Sz[0], cj.Sz[1], cj.Sx[0])
        opsp = (cp.Sz[0], cp.Sz[1], cp.Sx[0])
    else:
        opsj = (cj.Cdag[0]*cj.C[1]+cj.Cdag[1]*cj.C[0], cj.N[2], cj.N[0])
        opsp = (cp.Cdag[0]*cp.C[1]+cp.Cdag[1]*cp.C[0], cp.N[2], cp.N[0])
    rj = _algebra_numbers(cj, hj, opsj)
    rp = _algebra_numbers(cp, hp, opsp)
    re = _algebra_numbers(cp, hp, opsp, mode="ED")
    for k in rj:
        assert rj[k] == pytest.approx(re[k], abs=1e-6), (k, rj[k], re[k])
        assert rj[k] == pytest.approx(rp[k], abs=1e-6), (k, rj[k], rp[k])
    # the conjugate is not the state itself: |w.w| < <w|w> here
    assert rj["conjw"] < rj["ww"] - 0.1


@pytest.mark.parametrize("version", JULIA)
def test_random_state_orthogonal(version):
    c, h = _heisenberg(6, version, maxm=64)
    gs = c.get_gs()
    x1 = c.random_state(orthogonal=[gs])
    x2 = c.random_state(orthogonal=[gs, x1])
    for x, ws in ((x1, [gs]), (x2, [gs, x1])):
        assert abs(c.overlap(x, x) - 1.) < 1e-10
        for w in ws:
            assert abs(c.overlap(w, x)) < 1e-10


@pytest.mark.parametrize("version", JULIA)
def test_wrapper_methods_against_python(version):
    out = {}
    for v in (version, "python"):
        c, h = _heisenberg(8, v)
        gs = c.get_gs()
        out[v] = (gs.get_entropy(), gs.get_mutual_information(0, 3),
                  gs.get_pair_entropy(0, 1))
    for a, b in zip(out[version], out["python"]):
        assert a == pytest.approx(b, abs=1e-6)
    # central charge of the 20-site XX chain, c=1 in the scaling limit
    cc = {}
    for v in (version, "python"):
        sc = spinchain.Spin_Chain([2]*20, itensor_version=v)
        h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] for i in range(19))
        sc.set_hamiltonian(h)
        sc.maxm = 60
        sc.nsweeps = 14
        cc[v] = sc.get_gs().get_CFT_central_charge()
    assert cc[version] == pytest.approx(cc["python"], abs=1e-4)


def _evolved_sz(c, h, neel, tau, mode=None, nt0=None):
    """<Sz0> and <H> in exp(tau*h)|neel>/norm"""
    kw = {}
    if nt0 is not None: kw["nt0"] = nt0
    if mode is not None: kw["mode"] = mode
    w = c.exponential(tau*h, neel, **kw)
    nrm = c.overlap(w, w).real
    return (c.aMb(w, c.Sz[0], w)/nrm).real, (c.aMb(w, h, w)/nrm).real


def _taylor_dense(c, h, v, tau, nt0):
    """The same nt0 Taylor steps 1 + zH + z^2/2 H^2, z = tau/nt0, on the
    dense ED vector: what the MPS routes compute up to truncation"""
    ed = c.get_ED_obj()
    H = ed.MO2matrix(h)
    Sz = ed.MO2matrix(c.Sz[0])
    z = tau/nt0
    w = v.copy()
    for i in range(nt0):
        hw = H@w
        w = w + z*hw + 0.5*z*z*(H@hw)
    nrm = np.vdot(w, w).real
    return np.vdot(w, Sz@w).real/nrm, np.vdot(w, H@w).real/nrm


@pytest.mark.parametrize("version", JULIA)
def test_exponential_against_ed_and_python(version):
    n = 8
    cases = ((-0.5j, 200), (-1j, 200), (0.3, 200)) # exp(-iHt) and exp(+0.3H)
    res = {}
    for v in (version, "python"):
        c, h = _heisenberg(n, v)
        c.set_hamiltonian(sum((-1)**i*c.Sz[i] for i in range(n)))
        neel = c.get_gs() # the Neel state
        for tau, nt0 in cases:
            res[v, tau] = _evolved_sz(c, h, neel, tau, nt0=nt0)
        if v == version:
            neel_ed = c.get_gs(mode="ED")
            for tau, nt0 in cases:
                res["taylor", tau] = _taylor_dense(c, h, neel_ed.v, tau, nt0)
                res["ed", tau] = _evolved_sz(c, h, neel_ed, tau, mode="ED")
                res["big", tau] = _evolved_sz(c, h, neel, tau, nt0=20*nt0)
    for tau, nt0 in cases:
        for k in (0, 1):
            # the Taylor steps themselves, exactly
            assert res[version, tau][k] == pytest.approx(res["taylor", tau][k], abs=1e-6)
            # python truncates its step operator at the chain's cutoff, so
            # it sits about 1e-6 away from the exact Taylor steps
            assert res[version, tau][k] == pytest.approx(res["python", tau][k], abs=1e-5)
            # and the exact exponential once the steps are small
            assert res["big", tau][k] == pytest.approx(res["ed", tau][k], abs=1e-5)
    # exp(+0.3 H) raises <H> from the Neel value -(n-1)/4: the sign is +
    assert res[version, 0.3][1] > -0.25*(n-1) + 0.5


def _inverse_numbers(c, h, eps, mode=None):
    kw = {} if mode is None else {"mode": mode}
    gs = c.get_gs(**kw)
    e0 = c.gs_energy(**kw)
    b = c.applyoperator(eps*c.Sz[0], gs)
    A = h - (e0 + 1.0) + 0.1j
    if mode is None:
        x = c.applyinverse(A, b, delta=1e-8, maxn=400)
    else:
        x = c.applyinverse(A, b)
    r = c.summps(c.applyoperator(A, x), c.scale_mps(-1., b))
    return c.overlap(b, x), _norm(c, r), _norm(c, b)


@pytest.mark.parametrize("version", JULIA)
@pytest.mark.parametrize("eps", [1.0, 1e-2])
def test_applyinverse_against_ed(version, eps):
    c, h = _heisenberg(6, version, maxm=64)
    bx, res, nb = _inverse_numbers(c, h, eps)
    ce, he = _heisenberg(6, "python", maxm=64)
    bxe, _, _ = _inverse_numbers(ce, he, eps, mode="ED")
    assert res <= 1e-6*nb # BiCGSTAB stops on the bilinear residual
    assert abs(bx - bxe) <= 1e-6*abs(bxe)


ES = np.linspace(0., 4., 21)


def _pair(c, which):
    if which == "zz": return (c.Sz[0], c.Sz[0])
    return (c.Sx[0], c.Sy[1])


def _continuation_inputs(c, A, B, delta, es):
    """The resolvent values CVMimag continues from, the grid of
    cvm.dynamical_correlator_analytic_continuation"""
    wf = c.get_gs()
    wfa = A.get_dagger()*wf
    wfb = B*wf
    e0 = c.gs_energy()
    xz = np.linspace(min(es), max(es), 20) + delta*40*1j
    return np.array([wfb.dot(c.applyinverse(-c.hamiltonian+(e0+z), wfa))
                     for z in xz])


@pytest.mark.parametrize("version", JULIA)
@pytest.mark.parametrize("which", ["zz", "xy"])
def test_cvm_explicit_and_cvmimag(version, which):
    out = {}
    for v in (version, "python", "ED"):
        # the field splits the triplet, without it <Sx0 ... Sy1> is zero
        c, h = _heisenberg(6, "python" if v == "ED" else v, maxm=64,
                           field=0.3 if which == "xy" else 0.)
        if v == "ED": c.mode = "ED"
        c.cvm_tol = 1e-9
        A, B = _pair(c, which)
        es, d = nhdynamics.dynamical_correlator_cvm_explicit(c,
                name=(A, B), delta=0.1, es=ES)
        es, a = cvm.dynamical_correlator_analytic_continuation(c,
                name=(A, B), delta=0.1, es=ES)
        f = _continuation_inputs(c, A, B, 0.1, ES)
        out[v] = (np.array(d), np.array(a), f)
    peak = np.max(np.abs(out["ED"][0]))
    for v in (version, "python"):
        assert np.max(np.abs(out[v][0]-out["ED"][0])) < 1e-5*peak
    # CVMimag: the resolvent on the complex grid agrees with ED, and the
    # Pade continuation to the real axis runs. The continuation itself
    # amplifies a 1e-10 difference of its inputs to several percent of the
    # peak on every MPS backend ("python" included, 0.03 of a 0.23 peak
    # for (Sz0,Sz0)), so the continued curve is not compared.
    fe = out["ED"][2]
    assert np.max(np.abs(out[version][2]-fe)) < 1e-6*np.max(np.abs(fe))
    assert np.all(np.isfinite(out[version][1]))


@pytest.mark.parametrize("version", JULIA)
def test_non_hermitian_cvm(version):
    out = {}
    for v in (version, "ED"):
        c, h = _heisenberg(4, "python" if v == "ED" else v, maxm=32,
                           extra=0.3j)
        if v == "ED": c.mode = "ED"
        c.cvm_tol = 1e-9
        es, d = nhdynamics.dynamical_correlator_non_hermitian(c,
                name=(c.Sz[1], c.Sz[1]), delta=0.1, es=ES)
        out[v] = np.array(d)
    peak = np.max(np.abs(out["ED"]))
    assert np.max(np.abs(out[version]-out["ED"])) < 1e-5*peak


@pytest.mark.parametrize("version", JULIA)
def test_cvm_explicit_of_a_vanishing_pair(version):
    """(Sx0,Sy1) on the Heisenberg singlet is exactly zero. At the default
    cutoff, sums by the density-matrix algorithm inside BiCGSTAB broke the
    recursion down and this came back as about 5e9 at omega=1.8"""
    c, h = _heisenberg(6, version, maxm=64)
    c.cutoff = 1e-12 # the default
    c.cvm_tol = 1e-9
    es, d = nhdynamics.dynamical_correlator_cvm_explicit(c,
            name=(c.Sx[0], c.Sy[1]), delta=0.1, es=ES)
    assert np.max(np.abs(d)) < 1e-6
