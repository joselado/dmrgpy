"""Paths of the julia_live ports of 2026-11-15 that the other tests do not
run: the running shift and the two-site edge of mpsjulialive/implicit.jl
(tests/test_imaginary_time_implicit.py has the rest), CLIK-MPS
(mpsjulialive/clik.py and clik.jl, driven by cliksession.py) at a bond
dimension that truncates and with the exact residuals of a non-Hermitian
generator (tests/test_krylov_spectrum.py has the rest), and the direct
Hadamard product (mpsalgebratk/hadamard.py, mpsjulialive/arrays.jl) at a
bond dimension that truncates (tests/test_rsi_hadamard.py has the rest).

Every test carries the julia_live id, so -k "not julia_live" skips them.
Each new Julia signature costs about a minute of compilation, so the
chains are small and each test makes few Julia calls."""
import numpy as np
import pytest

from dmrgpy import imaginarytime, spinchain
from dmrgpy.krylovspectrum import _res, krylov_spectrum_dmrg
from dmrgpy.pyitensor import rsi
from _helpers import julia_live_param
from test_imaginary_time_implicit import _neel_chain
from test_nh_dmrg import nh_asymmetric_hopping_chain

JL = pytest.mark.parametrize("backend", [julia_live_param()])


def _dense_matrix(M):
    """A dense array of an ED operator, sparse or not"""
    return M.toarray() if hasattr(M, "toarray") else np.asarray(M)


def _site_arrays(sc, wf):
    """The site tensors of an MPS of either backend, as (left link, site,
    right link) arrays with length-1 links at the open ends"""
    if hasattr(wf, "jlmps"):
        from dmrgpy.mpsjulialive.juliasession import Main
        return [np.asarray(a) for a in Main.mps_site_arrays(wf.jlmps)]
    return [np.asarray(a) for a in sc._session.mps_site_arrays(wf.cpp_handle)]


def _bond_dimension(arrays):
    """The largest bond dimension of an MPS given by its site arrays"""
    return max(a.shape[2] for a in arrays)


def _tails(v, n, d, maxm):
    """eps_k^2 for the bonds k = 1..n-1 of a dense vector v of n sites of
    dimension d, site 1 the slowest index: the weight of the singular
    values beyond the first maxm when v is cut between sites k and k+1,
    meaning that eps_k is the error of the best rank-maxm approximation
    at that bond"""
    out = []
    for k in range(1, n):
        s = np.linalg.svd(np.reshape(v, (d**k, d**(n-k))), compute_uv=False)
        out.append(float(np.sum(s[maxm:]**2)))
    return np.array(out)


def _svd_cut(v, n, d, maxm):
    """The dense vector of the MPS made from v by one SVD sweep from the
    left that keeps at most maxm states per bond, site 1 the slowest
    index"""
    cores, rest, r = [], np.reshape(v, (1, -1)), 1
    for _ in range(n-1):
        u, s, vh = np.linalg.svd(np.reshape(rest, (r*d, -1)),
                                 full_matrices=False)
        m = min(maxm, len(s))
        cores.append(np.reshape(u[:, :m], (r, d, m)))
        rest, r = s[:m, None]*vh[:m], m
    cores.append(np.reshape(rest, (r, d, 1)))
    return rsi.dense_vector(cores)


# implicit_step on julia_live (mpsjulialive/implicit.jl)

@JL
def test_julia_live_running_shift_agrees_with_python(backend):
    """test_julia_live_agrees_with_python with shift="energy": each step
    takes its shift from <psi|H|psi>/<psi|psi> of its input state, the
    shift_energy branch of implicit_step. order=4 as there, since its
    complex factors keep the pole of a step off the real axis, which a
    running shift above E0 could otherwise meet. The flag goes to Julia
    with a shift of 0, so a julia_live that dropped the flag would follow
    the fixed shift 0 instead, and the python run at that shift says how
    far apart the two are."""
    kw = dict(nsteps=2, order=4)
    out = []
    for b in ("python", backend):
        sc, h, wf, wfed = _neel_chain(b, cutoff=1e-12)
        ops = [h, sc.Sz[0]*sc.Sz[1]]
        out.append(sc.evolve_imaginary_time_and_measure(ops, wf, 10.,
                   shift="energy", **kw)[1])
        if b == "python":
            fixed = sc.evolve_imaginary_time_and_measure(ops, wf, 10.,
                    shift=0., **kw)[1]
    branch = np.max(np.abs(out[0]-fixed))
    assert branch > 1e-5, ("on python the running shift and the fixed "
            "shift 0 give the same numbers (%.2e apart), so this test "
            "cannot see the shift_energy branch" % branch)
    diff = np.max(np.abs(out[0]-out[1]))
    assert diff < 1e-7, ("with shift='energy' julia_live differs from "
            "python by %.2e; julia_live is %.2e from python's fixed shift "
            "0, which it would match if it dropped the flag"
            % (diff, np.max(np.abs(out[1]-fixed))))


@JL
def test_julia_live_two_site_chain(backend):
    """evolve_imaginary_time on two sites, where the sweep has one bond:
    no environment is extended, and both half sweeps update the same
    bond. From the Neel state |dn,up> the dimer S_0.S_1 has the singlet
    (E = -3/4) and the Sz=0 triplet (E = 1/4) with equal weights. With the
    shift at E0 = -3/4 and order=2, a step of size dt keeps the singlet
    and multiplies the triplet by (1-dt/2)/(1+dt/2). After the steps of
    the grid the triplet over singlet ratio is the product rho of these
    factors, meaning that
        <H> = (-3/4 + rho^2/4)/(1 + rho^2),  <Sz_0> = -rho/(1 + rho^2).
    With one bond the local problem is the whole problem, so the default
    solver settings give these numbers to the GMRES tolerance on both
    backends. tau = 1 is about the inverse gap, so the state is still far
    from the singlet and the numbers depend on every step."""
    from dmrgpy.mpsalgebra import wavefunction_mode
    tau, nsteps, tau0 = 1.0, 2, 0.01
    rho = np.prod([(1-dt/2)/(1+dt/2)
                   for dt in imaginarytime.log_grid(tau, tau0, nsteps)])
    ref = np.array([1.0, (-0.75+0.25*rho**2)/(1+rho**2), -rho/(1+rho**2)])
    out = {}
    for b in ("python", backend):
        sc, h, wf, _ = _neel_chain(b, n=2, cutoff=1e-12)
        assert wavefunction_mode(wf) == "DMRG", ("the two-site start on %r "
                "is an ED state, so the sweep would never run" % (b,))
        psi = sc.evolve_imaginary_time(wf, tau, nsteps=nsteps, tau0=tau0,
                                       order=2, shift=-0.75)
        out[b] = np.array([psi.dot(psi).real, psi.aMb(h, psi).real,
                           psi.aMb(sc.Sz[0], psi).real])
    for b, vals in out.items():
        assert np.max(np.abs(vals-ref)) < 1e-8, ("%r on two sites: (norm, "
                "<H>, <Sz_0>) = %s, the exact steps give %s" % (b, vals, ref))
    diff = np.max(np.abs(out["python"]-out[backend]))
    assert diff < 1e-8, ("on two sites julia_live differs from python by "
                         "%.2e" % diff)


# CLIK-MPS on julia_live (mpsjulialive/clik.py, clik.jl)

@JL
def test_julia_live_krylov_spectrum_at_a_binding_maxm(backend):
    """get_krylov_spectrum on an 8-site Heisenberg chain at maxm=4, where
    the exact ground state needs 16 states at the middle bond: the random
    starts are drawn at maxm, and the global subspace expansion (cut back
    to maxm by clik_expand) and every TDVP step truncate.
    krylov_spectrum_dmrg is the DMRG half of get_krylov_spectrum, called
    with the same arguments, and it keeps the snapshots, whose bond
    dimension is checked here."""
    from dmrgpy.mpsjulialive.juliasession import Main
    n, maxm = 8, 4
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=backend,
                              maxm=maxm)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    ks = krylov_spectrum_dmrg(sc, -1*h)
    ev, res = ks.eigenvalues, _res(ks)
    assert np.all(np.isfinite(ev)), "non-finite levels at maxm=%d: %s" % (
        maxm, ev[~np.isfinite(ev)])
    assert np.all(np.isfinite(res)), ("non-finite eigenvalue changes at "
                                      "maxm=%d: %s" % (maxm, res))
    chi = max(_bond_dimension([np.asarray(a) for a in
                               Main.mps_site_arrays(s)]) for s in ks.states)
    assert chi <= maxm, ("a snapshot has bond dimension %d above maxm=%d: "
                         "the expansion or the TDVP step did not cut it"
                         % (chi, maxm))
    Hm = _dense_matrix(sc.get_ED_obj().MO2matrix(h))
    w, V = np.linalg.eigh(Hm)
    e0 = w[0]
    v = _svd_cut(V[:, 0], n, 2, maxm)
    e_cut = np.vdot(v, Hm@v).real/np.vdot(v, v).real
    assert e_cut-e0 > 1e-4, ("maxm=%d does not bind here (the cut ground "
            "state is %.2e above E0), so the test checks nothing about the "
            "truncation" % (maxm, e_cut-e0))
    # A guard, not a precision check. Every snapshot is an MPS of bond
    # dimension maxm, and the lowest Ritz value lies below the energy of
    # each of them, so the level is at most about the energy of the best
    # bond-maxm state the evolution reaches. e_cut, the energy of the exact
    # ground state cut to maxm by SVD, is the energy of one bond-maxm
    # state, and the variational optimum at maxm lies below it; twice its
    # distance to E0 leaves room for an evolution that stops short of that
    # optimum. From below the level cannot go under E0: the overlaps and
    # the matrix elements of G are exact contractions, so the projection
    # is a true Rayleigh-Ritz, and 1e-6 is room for its rounding.
    E = -ev[0].real
    assert E > e0-1e-6, ("the lowest level %.10f is below the exact E0 "
                         "%.10f at maxm=%d" % (E, e0, maxm))
    assert E-e0 < 2*(e_cut-e0), ("the lowest level is %.3e above E0 at "
            "maxm=%d, more than twice the %.3e of the SVD-cut ground state"
            % (E-e0, maxm, e_cut-e0))


@JL
def test_julia_live_exact_residuals_of_a_non_hermitian_generator(backend):
    """exact_residuals=True on the asymmetric fermion chain of
    test_krylov_spectrum.py, which is not normal: the residuals are read
    off the MPO of G^dag G, which the driver asks for as
    ops.product(G, ops.dagger(G)), G acting first. Taken in the other
    order, the quadratic form is ||(G^dag - conj(lambda)) r||^2, which is
    not small for a right eigenvector of a non-normal G. The order is
    also checked directly on two random states, against <G a|G b>."""
    from dmrgpy.mpsjulialive import clik as jclik
    sc, h = nh_asymmetric_hopping_chain(6, backend)
    assert not sc.is_hermitian(sc.hamiltonian)
    ev, res, _ = sc.get_krylov_spectrum(exact_residuals=True)
    res = np.asarray(res)
    assert np.all(np.isfinite(res)), "non-finite exact residuals: %s" % res
    assert np.all(res >= 0), "negative exact residuals: %s" % res
    k = int(np.argmin(res))
    assert res[k] < 1e-3, ("the best-converged level %s has a residual of "
                           "%.2e" % (ev[k], res[k]))
    assert np.max(res) > 10*res[k], ("the residuals do not tell converged "
            "from unconverged levels: %.2e to %.2e" % (res[k], np.max(res)))
    # Bauer-Fike: G = V diag(mu) V^-1 has an eigenvalue within
    # cond(V) ||G r - lambda r|| of lambda for a normalized r, so a
    # residual that is too small shows up here. 1e-7 of the operator
    # norm is room for the floor of the residuals, near 1e-8 of it.
    Gd = _dense_matrix(sc.get_ED_obj().get_operator(-1*h))
    mu, V = np.linalg.eig(Gd)
    V = V/np.linalg.norm(V, axis=0)
    kappa = np.linalg.cond(V)
    dist = np.min(np.abs(mu-ev[k]))
    floor = 1e-7*np.linalg.norm(Gd, 2)
    assert dist <= kappa*(res[k]+floor), ("the best-converged level %s is "
            "%.2e from the nearest exact eigenvalue, more than cond(V) = "
            "%.3g times its residual %.2e allows" % (ev[k], dist, kappa,
                                                    res[k]))
    # the product order, on the operations the driver uses
    ops = jclik.Session(sc, hermitian=False)
    G = ops.mpo(-1*h)
    Gdag = ops.dagger(G)
    W = ops.product(G, Gdag)   # as in cliksession.krylov_spectrum
    a, b = ops.random_starts(4, 0)
    ga, gb = ops.apply(G, a), ops.apply(G, b)
    da, db = ops.apply(Gdag, a), ops.apply(Gdag, b)
    gap = 0.
    for x, y, gx, gy, dx, dy in ((a, a, ga, ga, da, da),
                                 (a, b, ga, gb, da, db)):
        scale = ops.norm(gx)*ops.norm(gy)
        right = ops.inner(gx, gy)   # <x|G^dag G|y>
        wrong = ops.inner(dx, dy)   # <x|G G^dag|y>
        wxy = ops.sandwich(x, W, y)
        assert abs(wxy-right) <= 1e-6*scale, ("<x|W|y> = %s but "
                "<Gx|Gy> = %s and <G^dag x|G^dag y> = %s: the MPO of the "
                "residuals is not G^dag G" % (wxy, right, wrong))
        gap = max(gap, abs(right-wrong)/scale)
    assert gap > 1e-4, ("G G^dag and G^dag G agree to %.2e on these states, "
                        "so the check above cannot see the order" % gap)


# the direct Hadamard product on julia_live (mpsjulialive/arrays.jl)

def _field_heisenberg(n, backend):
    """The chain of test_rsi_hadamard.py's _heisenberg, with the cutoff
    set here since the bound below reads it"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=backend,
                              maxm=40, nsweeps=12, cutoff=1e-12)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h+0.1*sc.Sz[0])
    return sc


@JL
def test_julia_live_direct_hadamard_at_a_binding_maxm(backend):
    """|psi(s)|^2 of a 10-site ground state by method="direct" at maxm=6,
    where the exact product needs up to 32 states at the middle bond: the
    Kronecker product goes through mps_from_site_arrays and the cut
    through compress_mps. The state is multiplied by a global phase
    first, which |psi|^2 does not see and psi^2 would, so a conjugation
    lost on the way shows up.

    One SVD sweep from a canonical form is an orthogonal projection of the
    exact product p, with an error err between max_k eps_k (no state with
    bonds of at most maxm does better) and sqrt(sum_k eps_k^2) (the
    TT-SVD bound), eps_k the tail of p beyond maxm at bond k. The cutoff
    can drop a relative weight of at most cutoff more per bond, which is
    the pad in the bound."""
    n, maxm = 10, 6
    phase = complex(np.exp(0.7j))
    out = {}
    for b in ("python", backend):
        sc = _field_heisenberg(n, b)
        wf = sc.get_gs()*phase
        pex = np.abs(rsi.dense_vector(_site_arrays(sc, wf)))**2
        p = sc.hadamard_product(wf, wf, conj=(True, False),
                                method="direct", maxm=maxm)
        arrs = _site_arrays(sc, p)
        pv = rsi.dense_vector(arrs)
        chi = _bond_dimension(arrs)
        nex = np.linalg.norm(pex)
        eps2 = _tails(pex, n, 2, maxm)
        hi = np.sqrt(np.sum(eps2)+(n-1)*sc.cutoff*nex**2)
        err = np.linalg.norm(pv-pex)
        assert np.sqrt(np.max(eps2)) > 1e-6*nex, ("maxm=%d does not bind "
                "for |psi|^2, so the test checks nothing about the cut"
                % maxm)
        assert chi <= maxm, ("%r: the product has bond dimension %d above "
                             "maxm=%d" % (b, chi, maxm))
        assert err <= hi*(1+1e-8)+1e-14, ("%r: the cut product is %.3e "
                "from |psi|^2, above the SVD bound %.3e (a lost conjugation "
                "would put it near %.3e)" % (b, err, hi, abs(phase**2-1)*nex))
        assert abs(np.linalg.norm(pv)**2+err**2-nex**2) < 1e-10*nex**2, (
                "%r: the cut product is not an orthogonal projection of "
                "|psi|^2" % (b,))
        out[b] = (pex, pv, hi)
    (pp, vp, hp), (pj, vj, hj) = out["python"], out[backend]
    dpex = np.linalg.norm(pj-pp)
    assert dpex < 1e-3*np.linalg.norm(pp), ("the two ground states differ "
            "(|psi|^2 apart by %.2e), so the products cannot be compared"
            % dpex)
    # python cuts in one sweep from the left after a right-canonical sweep
    # (sum_mps), ITensorMPS's truncate in one sweep from the right after a
    # left-canonical one, so the two results are different vectors. Each
    # is within its bound of its own exact product, and the exact products
    # are dpex apart, so the two agree within the sum.
    d = np.linalg.norm(vj-vp)
    assert d <= (hp+hj+dpex)*(1+1e-8), ("julia_live and python products "
            "differ by %.3e, above the %.3e the two SVD bounds allow"
            % (d, hp+hj+dpex))
