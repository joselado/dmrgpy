"""The exact arithmetic of the 2026-11-21 hole hunt, findings 16, 17, 20 and
21 (docs/audit_2026_11_21_hole_hunt.md).

mpsalgebra.check_bicgstab formed A x and A x - b by the default routes of
the sessions, the density-matrix apply and sum, whose floor is set by
||A|| ||x||: next to a level (||x||/||b|| = 9e7 on 8 sites) it read 1.0 on
"python" and 0.54 on v3 for a true residual of 1.4e-7 and 2.4e-7, and the
public applyinverse warned on a converged v3 solve (2.58e-3 for 1.19e-4 at
delta = 1e-4). It now takes the session's residual_norm (julia_live:
mps_residual_norm): A x and the difference with nothing truncated, the norm
read off the orthogonality centre.

sum_mps_weighted took the density-matrix sum on "python" and a tree of
truncated pairwise sums on v3, so a sum whose coefficients cancel lost its
small part: [b, a] with [1/eps, -1/eps], b = a + eps c, was 4.7e-3 off c at
eps = 1e-6 on "python" and 6.2e-4 on v3, and [A, B, A] with [M, 1, -M] was
1.4e-2 off B at M = 1e5 on v3. method="exact" is the direct sum truncated
once, on both, and nhlevels._combine passes it. The default is unchanged,
since the ESR Krylov route passes hundreds of states.

The anchors are dense vectors and matrices read off the site arrays (8 and
10 sites), and a residual planted by construction.
"""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, multioperator, nhlevels, spinchain
from dmrgpy import mps as dmps
from dmrgpy import mpsalgebra as dma
from dmrgpy.pyitensor import mpsalgebra as ma
from _helpers import julia_live_param

HAVE3 = cppext.available(3)
NEEDS3 = pytest.mark.skipif(not HAVE3,
        reason="needs the compiled itensor_version=3 extension")
BACKENDS = ["python", pytest.param(3, marks=NEEDS3)]
N = 8
DELTA = 1e-5 # cvm_tol, the default delta of applyinverse


def heis(c, n=N):
    h = 0
    for i in range(n-1):
        h = h + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
    return h


def python_reference():
    """8 spins held exactly on "python" (bond dimension 16 at maxm=32), no
    truncation anywhere in the solves"""
    np.random.seed(2)
    sp = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version="python")
    hp = heis(sp)
    sp.set_hamiltonian(hp)
    sp.maxm, sp.nsweeps, sp.cutoff = 32, 12, 0.0
    sp.gs_energy()
    gs = sp.get_gs()
    gs = gs*(1./np.sqrt(abs(gs.dot(gs))))
    sp._session.set_sweep_params(sp.maxm, sp.nsweeps, sp.cutoff, sp.noise)
    E0 = sp.vev(hp, wf=gs).real
    return sp, hp, gs, E0


@pytest.fixture(scope="module")
def ref():
    return python_reference()


@pytest.fixture(scope="module")
def chain3():
    if not HAVE3:
        pytest.skip("needs the compiled itensor_version=3 extension")
    s3c = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=3)
    assert s3c.itensor_version == 3 and \
        getattr(s3c, "_session", None) is not None
    h3 = heis(s3c)
    s3c.set_hamiltonian(h3)
    s3c.maxm, s3c.nsweeps, s3c.cutoff = 32, 12, 0.0
    s3c._session.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff,
                                  s3c.noise)
    return s3c, h3


def exact_pair(sp, hp, gs, E0, shift):
    """A = H - e with e = E0 - shift, x = gs/shift + Sx_0 gs and b = A x,
    x and b formed by the exact routes, so ||A x - b|| is roundoff"""
    s = sp._session
    A = hp - (E0 - shift)*multioperator.identity()
    w1 = sp.Sx[0]*gs
    x_h = ma.sum_many([gs.cpp_handle*complex(1./shift), w1.cpp_handle*1.0],
                      cutoff=0.0, maxdim=10**6, method="exact")
    b_h = ma.applyMPO(s._mpo(A.to_terms()), x_h, cutoff=0.0, maxdim=None,
                      method="exact")
    b_h.noPrime("Site")
    return A, x_h, b_h


def to_v3(sp, s3c, handle_py):
    """A "python" MPS moved into the v3 session through its site arrays"""
    arrs = [np.ascontiguousarray(np.asarray(a), dtype=complex)
            for a in sp._session.mps_site_arrays(handle_py)]
    return dmps.MPS(MBO=s3c, cpp_handle=s3c._session.mps_from_site_arrays(arrs))


def vec_py(sp, h_):
    sites = sp._session.sites
    T = h_.A(1)
    for i in range(2, h_.length()+1):
        T = T*h_.A(i)
    arr = T.transpose_to([sites.si(i) for i in range(1, h_.length()+1)])
    return np.asarray(arr).reshape(-1)


def mpo_dense_py(sp, mpo):
    sites = sp._session.sites
    T = mpo.A(1)
    for i in range(2, mpo.length()+1):
        T = T*mpo.A(i)
    order = ([sites.si(i).prime(1) for i in range(1, N+1)]
             + [sites.si(i) for i in range(1, N+1)])
    return np.asarray(T.transpose_to(order)).reshape(2**N, 2**N)


def vec_arrays(arrs):
    arrs = [np.asarray(a) for a in arrs]
    v = arrs[0].reshape(arrs[0].shape[1], arrs[0].shape[2])
    for a in arrs[1:]:
        v = np.tensordot(v, a, axes=([1], [0])).reshape(-1, a.shape[2])
    return v.reshape(-1)


def mpo_dense_arrays(Ws):
    Ws = [np.asarray(W) for W in Ws]
    M = Ws[0][0] # (s', s, r)
    for W in Ws[1:]:
        T = np.tensordot(M, W, axes=([2], [0])) # (R, C, s', s, r)
        R, C = M.shape[0], M.shape[1]
        M = T.transpose(0, 2, 1, 3, 4).reshape(R*W.shape[1], C*W.shape[2],
                                               W.shape[3])
    return M[:, :, 0]


def dense_residual3(s3c, A3, x3, b3):
    s3 = s3c._session
    Ad = mpo_dense_arrays(s3.mpo_site_arrays(s3.build_operator(A3.to_terms())))
    vx = vec_arrays(s3.mps_site_arrays(x3.cpp_handle))
    vb = vec_arrays(s3.mps_site_arrays(b3.cpp_handle))
    return np.linalg.norm(Ad @ vx - vb)/np.linalg.norm(vb)


def dense_state(s, handle):
    """The amplitudes of an MPS of either session, from its site arrays"""
    arrs = [np.asarray(x) for x in s.mps_site_arrays(handle)]
    T = arrs[0][0] # (d, chi)
    for a in arrs[1:]: # (chi_l, d, chi_r)
        T = np.einsum("xa,asb->xsb", T, a).reshape(-1, a.shape[2])
    return T[:, 0]


def read_check(chain, A, x, b, delta=DELTA):
    """check_bicgstab's value and whether it warned, with BICGSTAB_RAISE
    lifted, so that a floor above 1 reads as a number"""
    nb = np.sqrt(abs(b.dot(b)))
    old = dma.BICGSTAB_RAISE
    dma.BICGSTAB_RAISE = np.inf
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = dma.check_bicgstab(chain, A, x, b, nb, delta)
    finally:
        dma.BICGSTAB_RAISE = old
    return res, any("BiCGSTAB" in str(m.message) for m in w)


@pytest.mark.parametrize("iv", BACKENDS)
@pytest.mark.parametrize("shift", [1e-2, 1e-4, 1e-6, 1e-8])
def test_check_is_not_above_the_dense_residual(ref, iv, shift, request):
    """An exactly formed pair (x, b = A x) next to a level: the check reads
    at most 10 times the dense residual, and does not warn. The default
    routes read 2.1e-9, 4.4e-5 and 1.00 on "python" for dense residuals of
    1e-11, 1.1e-9 and 1.4e-7 at shifts 1e-4, 1e-6 and 1e-8, and 3.6e-7,
    1.2e-2 and 0.54 on v3 for 6.9e-11, 2.5e-8 and 2.4e-7 (v3's own MPO
    against a b formed with the "python" one, so its dense residual is the
    mismatch of the two MPOs times ||x||)"""
    sp, hp, gs, E0 = ref
    A, x_h, b_h = exact_pair(sp, hp, gs, E0, shift)
    if iv == "python":
        chain = sp
        x = dmps.MPS(MBO=sp, cpp_handle=x_h)
        b = dmps.MPS(MBO=sp, cpp_handle=b_h)
        Ad = mpo_dense_py(sp, sp._session._mpo(A.to_terms()))
        vb = vec_py(sp, b_h)
        res_dense = np.linalg.norm(Ad @ vec_py(sp, x_h) - vb)/np.linalg.norm(vb)
    else:
        chain, h3 = request.getfixturevalue("chain3")
        A = h3 - (E0 - shift)*multioperator.identity()
        x, b = to_v3(sp, chain, x_h), to_v3(sp, chain, b_h)
        res_dense = dense_residual3(chain, A, x, b)
    res, warned = read_check(chain, A, x, b)
    assert res <= 10.*max(res_dense, 1e-14), (res, res_dense)
    assert not warned


@pytest.mark.parametrize("iv", BACKENDS + [julia_live_param()])
def test_check_reads_a_planted_residual(iv):
    """b2 = A x + eta c on top of the exactly formed pair at shift 1e-8
    (||x||/||b|| = 9e7): the check reads eta ||c||/||b2||, about 1e-3 on
    "python" and v3, to 1e-2 of itself, where the default routes read
    about 1.0 for the pair alone. On v3 the
    pair is moved from "python", and its MPO mismatch (2.4e-7) is below a
    percent of the planted part. On julia_live x comes from the public
    arithmetic and b = A x from the exact contraction of the same MPO."""
    shift, eta = 1e-8, 2e-3
    if iv == "julia_live":
        from dmrgpy.mpsjulialive.juliasession import Main
        np.random.seed(2)
        sc = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version="julia_live")
        h = heis(sc)
        sc.set_hamiltonian(h)
        sc.maxm, sc.nsweeps, sc.cutoff = 32, 12, 0.0
        sc.gs_energy()
        gs = sc.get_gs()
        gs = gs*(1./np.sqrt(abs(gs.dot(gs))))
        E0 = sc.vev(h, wf=gs).real
        A = h - (E0 - shift)*multioperator.identity()
        x = gs*(1./shift) + sc.Sx[0]*gs
        c = sc.Sz[N-1]*gs
        bj = Main.apply_exact(dma._julia_mpo(sc, A), x.jlmps)
        b2 = dma._julia_mps(sc, Main.sum_direct(bj, (c*eta).jlmps))
        chain = sc
    else:
        sp, hp, gs, E0 = python_reference()
        A, x_h, b_h = exact_pair(sp, hp, gs, E0, shift)
        c = sp.Sz[N-1]*gs
        b2_h = ma.sum_many([b_h, c.cpp_handle*complex(eta)], cutoff=0.0,
                           maxdim=10**6, method="exact")
        if iv == "python":
            chain = sp
            x = dmps.MPS(MBO=sp, cpp_handle=x_h)
            b2 = dmps.MPS(MBO=sp, cpp_handle=b2_h)
        else:
            chain = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=3)
            h3 = heis(chain)
            chain.set_hamiltonian(h3)
            chain.maxm, chain.nsweeps, chain.cutoff = 32, 12, 0.0
            A = h3 - (E0 - shift)*multioperator.identity()
            x, b2 = to_v3(sp, chain, x_h), to_v3(sp, chain, b2_h)
    expected = eta*np.sqrt(abs(c.dot(c)))/np.sqrt(abs(b2.dot(b2)))
    res, _ = read_check(chain, A, x, b2)
    assert abs(res - expected) <= 1e-2*expected, (res, expected)


@NEEDS3
def test_public_applyinverse_does_not_warn_on_a_converged_v3_solve():
    """v3's own BiCGSTAB of (H - e) x = gs + Sx_0 gs at shift 1e-6 and
    delta = 1e-4 converges (dense residual 1.19e-4, below 10 delta), and the
    public applyinverse is quiet; the density-matrix apply read 2.58e-3
    and warned"""
    np.random.seed(2)
    s3c = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=3)
    h3 = heis(s3c)
    s3c.set_hamiltonian(h3)
    s3c.maxm, s3c.nsweeps, s3c.cutoff = 32, 12, 0.0
    s3 = s3c._session
    s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
    s3c.gs_energy()
    gs = s3c.get_gs()
    gs = gs*(1./np.sqrt(abs(gs.dot(gs))))
    E0 = s3c.vev(h3, wf=gs).real
    b3 = gs + s3c.Sx[0]*gs
    nb3 = np.sqrt(abs(b3.dot(b3)))
    shift, delta = 1e-6, 1e-4
    A3 = h3 - (E0 - shift)*multioperator.identity()
    s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
    s3.set_mpomaxm(max(s3c.maxm, s3c.mpomaxm))
    x3 = dmps.MPS(MBO=s3c, cpp_handle=s3.apply_inverse(
        A3.to_terms(), b3.cpp_handle, delta*nb3, 400))
    rd3 = dense_residual3(s3c, A3, x3, b3)
    assert rd3 <= 10.*delta, "the solve itself did not converge: %g" % rd3
    res, warned = read_check(s3c, A3, x3, b3, delta)
    assert res <= 10.*rd3 and not warned, (res, rd3)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        s3c.applyinverse(A3, b3, delta=delta, maxn=400)
    assert not any("BiCGSTAB" in str(m.message) for m in w)


@pytest.mark.parametrize("iv", BACKENDS)
def test_exact_weighted_sum_keeps_what_survives_a_cancellation(iv):
    """[A, B, A] with [M, 1, -M] on 10 sites: method="exact" returns B, to
    the truncation at cutoff 1e-14, at every M (v3's default tree was
    8.0e-5, 1.1e-3 and 1.4e-2 off at M = 1e2, 1e4, 1e5); an unknown method
    raises"""
    n = 10
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=iv)
    sc.maxm, sc.nsweeps, sc.cutoff = 80, 12, 1e-14
    hA = 0
    for i in range(n-1):
        hA = hA + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        hA = hA + (0.3*(1.0 if i % 3 == 0 else -0.5))*sc.Sz[i]
    sc.set_hamiltonian(hA)
    A = sc.get_gs().copy()
    sc.set_hamiltonian(heis(sc, n))
    B = sc.get_gs().copy()
    A = A*(1./np.sqrt(A.dot(A).real))
    B = B*(1./np.sqrt(B.dot(B).real))
    s = sc._session
    s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    vB = dense_state(s, B.cpp_handle)
    assert abs(np.vdot(dense_state(s, A.cpp_handle), vB)) > 0.5 # a real cancellation
    handles = [A.cpp_handle, B.cpp_handle, A.cpp_handle]
    for M in (1e2, 1e4, 1e5):
        h = s.sum_mps_weighted(handles, [complex(M), 1.0+0j, complex(-M)],
                               method="exact")
        assert np.linalg.norm(dense_state(s, h) - vB) < 1e-6, M
    with pytest.raises(ValueError):
        s.sum_mps_weighted(handles, [1.0+0j]*3, method="exactly")


@pytest.mark.parametrize("iv", BACKENDS)
def test_nh_levels_combination_resolves_a_cancellation(iv):
    """nhlevels._combine, the Ritz vectors of get_nh_levels: [b, a] with
    [1/eps, -1/eps], b = a + eps c formed exactly, returns c (eps = 1e-6,
    8 sites); the default sums were 4.7e-3 ("python") and 6.2e-4 (v3) off,
    the exact one 4.4e-9"""
    eps = 1e-6
    np.random.seed(1)
    sp = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version="python")
    sp.set_hamiltonian(heis(sp))
    sp.maxm, sp.nsweeps, sp.cutoff = 32, 12, 1e-14
    sp.gs_energy()
    a = sp.get_gs()
    a = a*(1./np.sqrt(abs(a.dot(a))))
    sp._session.set_sweep_params(sp.maxm, sp.nsweeps, sp.cutoff, sp.noise)
    c = sp.Sz[0]*a
    b_h = ma.sum_many([a.cpp_handle, c.cpp_handle*complex(eps)], cutoff=0.0,
                      maxdim=10**6, method="exact")
    vc = dense_state(sp._session, c.cpp_handle)
    if iv == "python":
        chain = sp
        snaps = [dmps.MPS(MBO=sp, cpp_handle=b_h), a]
    else:
        chain = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=3)
        chain.set_hamiltonian(heis(chain))
        chain.maxm, chain.nsweeps, chain.cutoff = 32, 12, 1e-14
        chain._session.set_sweep_params(chain.maxm, chain.nsweeps,
                                        chain.cutoff, chain.noise)
        snaps = [to_v3(sp, chain, b_h), to_v3(sp, chain, a.cpp_handle)]
    y = nhlevels._combine(chain, snaps, [1./eps, -1./eps])
    vy = dense_state(chain._session, y.cpp_handle)
    assert np.linalg.norm(vy - vc)/np.linalg.norm(vc) < 1e-5
