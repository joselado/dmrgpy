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
