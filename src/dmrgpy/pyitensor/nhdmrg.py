"""Non-Hermitian two-site DMRG (NH-DMRG) for the pure-Python backend.

Port of mpscpp3/chain_session.h's Chain::nhdmrg (the annotated original --
see its long comment for the algorithm, its ITensorNHDMRG.jl reference,
and the two deliberate deviations) to this package's engine:

- the same "onesided" local solver (independent Krylov-Schur solves of
  the projected two-site eigenproblems for H and H^T, smallest real
  part, with the left solve anchored to the right solve's eigenvalue and
  Re-degenerate Ritz values tie-broken toward the previous bond's), and
- the same "fidelity" truncation (one shared isometry from the hermitian
  average rho = (rho_l + rho_r)/2 of the left/right two-site reduced
  density matrices), which keeps psil and psir on identical site *and*
  link Index objects throughout -- so the environments hit the exact same
  bra/ket link-identity collision dmrg.py's self-overlap environments do,
  and reuse its _extend_left/_extend_right/_all_right_environments
  builders directly (their bra= parameter draws the bra tensors from the
  *other* MPS).

Inside a sweep the left state is the reference's own: the TRANSPOSED
left eigenvector wl, H^T|wl> = E|wl> (ProjNHMPO's second projector is
dag(swapprime(conj(H))) = H^T), not dmrgpy's adjoint one psil = conj(wl),
H^dagger|psil> = conj(E)|psil>. _nhdmrg_one_sweep conjugates psil on the
way in and back on the way out, so its callers keep the adjoint
convention. The shared isometry has to hold the right state and the
sweep's left state, and the two conventions ask it for different spaces:
{r, wl} or {r, conj(wl)}. For a complex-symmetric H (every imaginary
potential, H^T = H) wl is r itself, so the transposed sweep keeps r's
own Schmidt basis while the adjoint one had to hold r and conj(r) at
once. Measured on the 10-site PT-symmetric XX chain with i 0.3 (-1)^i
Sz_i, whose four lowest eigenvalues share one real part
(docs/audit_2026_11_15_hole_hunt.md, finding 3): the adjoint sweep had a
fixed point that is not an eigenpair, relative residual 0.29 to 0.44 with
the old local Arnoldi and still 1.3e-2 with a converged one at 10, 30 and
60 sweeps, Re E up to 0.03 below the spectrum; the transposed sweep with
the Krylov-Schur solver gives 3.4e-6 on every seed, the residual of the
1e-10 cutoff, as ITensorNHDMRG.jl does on the same chain. Either change
alone does not: the transposed sweep with the single-vector restart stops
at 8e-4 to 4e-2.

Everything bond-local is done on flat numpy arrays (via dmrg.py's
two_site_heff matvec factory): the biorthogonal normalization and the
fidelity truncation (a dense eigh of the small
rho matrix, truncated with svd.py's own eigh_truncate() -- shared with
gse.py's Krylov-subspace basis enrichment, which reduces to the same
step). Unlike dmrg.py
-- which deliberately skips noise-term support -- the reference's
noise/perturbation term IS implemented here (hermitized, exactly as in
the C++ ports), because for NH-DMRG it demonstrably matters: measured on
the C++ backend's PT-symmetric test chain, noise=0.1 converged 5/5 runs
where noise=0 converged 2/5 (see tests/test_nh_dmrg.py's models).

Like the C++ ports, every run starts from a fresh random MPS: the
non-Hermitian "energy" is not a variational bound, so a stalled run can
only be detected by the caller's eigen-residual certificate and cured by
re-running (see dmrgpy's nhdmrg.py retry loop, shared by all three
backends).

nhdmrg_generalized() below solves the *generalized* eigenproblem
H|psi_R>=lambda*A|psi_R> (A Hermitian positive definite, H possibly
non-Hermitian) by reusing the same per-sweep machinery
(_nhdmrg_one_sweep, factored out of nhdmrg()'s own loop body) against a
complex-shifted pair (H-lambda*A, HA-conj(lambda)*A), self-consistently
updating lambda between sweeps via the biorthogonal generalized Rayleigh
quotient -- see its own docstring for the derivation. The non-Hermitian
counterpart of dmrg.py's dmrg_generalized(), and (like
mpsiram_generalized in algebra/arpacktk.py, ARPACK mode 2) needs no
Hermiticity assumption on H itself, only on the metric A.
"""

import numpy as np

from . import backend
from .dmrg import (_all_right_environments, _extend_left, _extend_right,
                   lam0_is_unset, two_site_heff)
from .index import Index
from .mpsalgebra import inner
from .mpsalgebra import sum as mps_sum
from .svd import eigh_truncate
from .tensor import ITensor, dag


def _select_ritz(evals, sel, target):
    """Index of the Ritz value selected by `sel` ("SR" = smallest real
    part, "closest" = closest to `target`, "SRTieBreak" = smallest real
    part with Re-degenerate candidates tie-broken toward `target`). This
    global-min-then-candidates formulation is the reference semantics all
    three backends implement (the C++ ports mirror it exactly), so the
    backends cannot silently target different members of a degenerate
    Ritz cluster.

    The tie window is measured against the Ritz spectrum's own spread,
    1e-6*(max Re - min Re), plus a roundoff floor of 100 eps times the
    largest |Ritz value| (the two members of a genuinely Re-degenerate
    pair are split by at least that much), and compared with <=, so a
    lone Ritz value or an all-equal set is its own candidate. It used to
    be 1e-6*(1+|remin|): absolute at small |E|, so for s*H it was
    1e-6/s wide in units of H, and growing with |E|, so a constant
    offset c widened it to 1e-6*c. Once that exceeded the real-part gap
    the sweep followed the previous bond onto an excited eigenpair, a
    converged one that no residual certificate can reject: 0.482 off ED
    at s=1e-6 and at c=5e5 on a 6-site Heisenberg + 0.3j*Sz0 chain,
    1.06 at s=1e-7 (2026-09-25b audit, finding 25). The spread is
    shift-invariant and scale-covariant, so neither trigger remains.
    What it does not buy is headroom at large L: the spread grows like
    the local bandwidth, ~L*J, so a gapless chain whose gap is ~J/L
    still meets the window once L^2 is of order 1e6, the same regime
    where the old form failed through |E0| ~ L*J."""
    if sel == "closest":
        return int(np.argmin(np.abs(evals - target)))
    if sel == "SRTieBreak":
        re = evals.real
        remin = re.min()
        degtol = (1e-6 * (re.max() - remin)
                  + 100.0 * np.finfo(float).eps * np.max(np.abs(evals)))
        cand = np.flatnonzero(re <= remin + degtol)
        return int(cand[np.argmin(np.abs(evals[cand] - target))])
    return int(np.argmin(evals.real))


def _ritz_order(evals, sel, target):
    """Every Ritz index, most wanted first under `sel`; the first entry is
    _select_ritz(evals, sel, target) exactly (stable sorts, so a tie
    resolves to the lowest index, as np.argmin does). "SR" sorts by real
    part, "closest" by the distance to `target`, and "SRTieBreak" puts
    _select_ritz's candidate set first, sorted by the distance to
    `target`, and the rest after it by real part. The thick restart of
    _krylov_schur_smallest_real keeps a leading slice of this order;
    chain_session.h's nh_ritz_order is the same rule."""
    evals = np.asarray(evals)
    idx = np.arange(evals.size)
    if sel == "closest":
        return idx[np.argsort(np.abs(evals - target), kind="stable")]
    re = evals.real
    if sel == "SRTieBreak":
        remin = re.min()
        degtol = (1e-6 * (re.max() - remin)
                  + 100.0 * np.finfo(float).eps * np.max(np.abs(evals)))
        incand = re <= remin + degtol
        cand = idx[incand]
        cand = cand[np.argsort(np.abs(evals[cand] - target), kind="stable")]
        rest = idx[~incand]
        rest = rest[np.argsort(re[rest], kind="stable")]
        return np.concatenate([cand, rest])
    return idx[np.argsort(re, kind="stable")]


def _orthonormal_columns(Y, droptol=1e-8):
    """An orthonormal basis of the span of Y's columns, by modified
    Gram-Schmidt with one reorthogonalization pass, in column order: a
    column left with less than `droptol` of its norm after the projection
    is dropped, so a nearly defective cluster of Ritz vectors gives fewer
    columns rather than an ill-conditioned basis. The first column is
    always kept."""
    cols = []
    for c in range(Y.shape[1]):
        v = np.array(Y[:, c], dtype=complex)
        nv0 = np.linalg.norm(v)
        for _ in range(2):
            for u in cols:
                v = v - np.vdot(u, v) * u
        nv = np.linalg.norm(v)
        if nv > droptol * nv0 or not cols:
            cols.append(v / nv)
    return np.column_stack(cols)


def _krylov_schur_smallest_real(matvec, x0, krylovdim, restarts,
                                sel="SR", target=0.0):
    """Krylov-Schur (thick-restart) Arnoldi on a matrix-free operator:
    `restarts` thick restarts, so restarts+1 Krylov builds of at most
    `krylovdim` vectors (KrylovKit's maxiter=3 at the default 2), returning the
    Ritz pair selected by `sel` (_select_ritz). Between builds it keeps
    the 3/5 most wanted Ritz vectors of the build (_ritz_order), not only
    the selected one: Y those Ritz vectors of the small matrix h, Z an
    orthonormal basis of their span (_orthonormal_columns), the new basis
    Q Z, the new leading block Z^dag h Z and the coupling row
    h[m, m-1] Z[m-1, :] onto the old residual vector, which stays the next
    Krylov vector. span(Y) is invariant under h, so A Q Z = Q Z (Z^dag h
    Z) + q_m (h[m, m-1] Z[m-1, :]) exactly and the Arnoldi relation
    survives the restart: the residual bound |h[m, m-1]| |y[m-1]| of the
    selected pair stays valid, and it stops the builds early below 1e-10
    (1 + |lambda|). The kept count is KrylovKit's (3 krylovdim + 2
    converged)/5 with nothing converged; KrylovKit's eigsolve is what
    ITensorNHDMRG.jl calls.

    It replaced a restart on the single selected Ritz vector, which loses
    the rest of a cluster of eigenvalues with one real part at every
    restart: on the 10-site PT chain of finding 3 the center bonds' local
    residuals stayed at 1e-2 with it. The Ritz vectors of a non-normal h
    are where Schur reordering usually comes in; an eigendecomposition of
    h plus Gram-Schmidt on m-vectors needs nothing but eig, which is what
    chain_session.h's nh_krylov_schur ports.

    GPU: the Krylov basis lives in a backend.zeros store written with
    backend.setrow and orthogonalized with backend.project_rows; only
    the overlaps, the small matrix h and the m-vectors of the restart go
    to the host. The returned vector is a host array, as the sweep's
    bond-local algebra is."""
    x0 = backend.asarray(x0)
    d = int(x0.size)
    nx = float(np.linalg.norm(backend.to_host(x0)))
    if nx < 1e-14:
        return 0.0 + 0.0j, backend.to_host(x0)
    m = min(int(krylovdim), d)
    Q = backend.setrow(backend.zeros((m + 1, d)), 0, x0 / nx)
    h = np.zeros((m + 1, m), dtype=complex)
    k = 0  # leading vectors carried over by the last thick restart
    lam, x = 0.0 + 0.0j, backend.to_host(x0)
    nbuilds = max(0, int(restarts)) + 1
    for build in range(nbuilds):
        built = m
        breakdown = False
        for j in range(k, m):
            w = matvec(Q[j])
            for _ in range(2):  # classical Gram-Schmidt, then once more
                c, w = backend.project_rows(Q, w, j + 1)
                h[:j + 1, j] += backend.to_host(c)[:j + 1]
            nw = float(np.linalg.norm(backend.to_host(w)))
            h[j + 1, j] = nw
            if nw < 1e-13:  # happy breakdown: an invariant subspace
                built = j + 1
                breakdown = True
                break
            Q = backend.setrow(Q, j + 1, w / nw)
        evals, Y = np.linalg.eig(h[:built, :built])
        order = _ritz_order(evals, sel, target)
        i = int(order[0])
        lam = complex(evals[i])
        x = backend.to_host(backend.asarray(Y[:, i]) @ Q[:built])
        resid = abs(h[built, built - 1]) * abs(Y[built - 1, i])
        if (breakdown or built == d or build == nbuilds - 1
                or resid < 1e-10 * (1.0 + abs(lam))):
            break
        Z = _orthonormal_columns(Y[:, order[:max(1, (3 * built) // 5)]])
        k = Z.shape[1]
        T = Z.conj().T @ h[:built, :built] @ Z
        b = h[built, built - 1] * Z[built - 1, :]
        Qk = backend.asarray(Z.T) @ Q[:built]
        qnext = Q[built]
        Qn = backend.zeros((m + 1, d))
        for c in range(k):
            Qn = backend.setrow(Qn, c, Qk[c])
        Q = backend.setrow(Qn, k, qnext)
        h = np.zeros((m + 1, m), dtype=complex)
        h[:k, :k] = T
        h[k, :k] = b
    return lam, x


def _nhdmrg_one_sweep(psil, psir, H, HA, maxdim, cutoff, noise,
                      krylovdim, restarts, energy_hint, have_energy_hint):
    """One full forward-then-backward NH-DMRG sweep against a fixed
    (H, HA=H^dagger) pair, mutating psil/psir in place. Builds its own
    left/right environments from scratch (self-contained, like dmrg.py's
    _dmrg_one_sweep) rather than assuming any carried over from a
    previous sweep -- needed by nhdmrg_generalized() below, whose H/HA
    change every outer iteration, so stale environments from a previous
    (H,HA) pair would be wrong. Returns the final bond's local eigenvalue
    (nhdmrg()'s own within-sweep tie-break anchor; see _select_ritz),
    starting the SRTieBreak anchor from (energy_hint, have_energy_hint)
    instead of always (0, False) so a caller can continue anchoring
    across sweeps that share the same (H, HA) -- nhdmrg() itself does,
    nhdmrg_generalized() below deliberately does not (see its own
    docstring).

    psil comes in and goes out as the adjoint left state (H^dagger|psil>
    = conj(E)|psil>, the callers' convention); in between it holds
    conj(psil), the transposed left state wl the sweep works with (module
    docstring), solved against H^T = conj(HA) with its eigenvalue
    anchored to the right solve's E itself."""
    n = psir.length()
    energy = energy_hint
    have_energy = have_energy_hint
    _conj_in_place(psil)  # adjoint -> transposed left state
    HT = _conj_chain(HA)  # H^T, from the adjoint's terms
    # (env tensor, dangling bra link) per boundary, one family per
    # projected problem: H sandwiched between <wl| and |psir>, and
    # H^T sandwiched between <psir| and |wl>. Built once per sweep and
    # then maintained incrementally as the sweep itself progresses (each
    # backward half-sweep refreshes every right entry, each forward
    # half-sweep every left entry), like the C++ ports -- a per-bond
    # rebuild would be redundant work.
    right_h = _all_right_environments(H, psir, bra=psil)
    right_a = _all_right_environments(HT, psil, bra=psir)
    left_h = {0: (None, None)}
    left_a = {0: (None, None)}

    for ha in (1, 2):
        bonds = range(1, n) if ha == 1 else range(n - 1, 0, -1)
        for b in bonds:
            Lh, Lbra_h = left_h[b - 1]
            Rh, Rbra_h = right_h[b + 2]
            La, Lbra_a = left_a[b - 1]
            Ra, Rbra_a = right_a[b + 2]
            mv_r, order_in, shape, x0r = two_site_heff(
                Lh, Lbra_h, H, psir, b, Rh, Rbra_h)
            mv_l, _, _, x0l = two_site_heff(
                La, Lbra_a, HT, psil, b, Ra, Rbra_a)
            er, thr = _krylov_schur_smallest_real(
                mv_r, x0r, krylovdim, restarts,
                sel="SRTieBreak" if have_energy else "SR", target=energy)
            # only the eigenvector is needed from the transposed solve
            # (its eigenvalue is er itself, H^T having H's spectrum)
            _, thl = _krylov_schur_smallest_real(
                mv_l, x0l, krylovdim, restarts,
                sel="closest", target=er)
            energy = er
            have_energy = True

            thl = thl / np.linalg.norm(thl)
            thr = thr / np.linalg.norm(thr)
            # rescale so <thl|thr> = 1, split between the two states,
            # with inner()'s conjugated overlap as the reference
            # normalizes its (wl, r) pair; the scale and phase only set
            # the size of the center tensors and the phase of the noise
            # term, rho_l and rho_r keep equal weight either way, and
            # nhdmrg() makes <psil|psir> = 1 once at the end (separate
            # real branch: see chain_session.h's note on sqrt's branch
            # cut for real negative overlaps)
            ov = np.vdot(thl, thr)
            if abs(ov) > 1e-12:
                if abs(ov.imag) < 1e-14 * abs(ov):
                    sq = np.sqrt(abs(ov.real))
                    thl = thl / sq
                    thr = thr / (-sq if ov.real < 0 else sq)
                else:
                    thl = thl / np.sqrt(np.conj(ov))
                    thr = thr / np.sqrt(ov)

            # fidelity truncation: keep-indices are a contiguous slice
            # of order_in ((left_link?, s_b) for ha==1, (s_{b+1},
            # right_link?) for ha==2 -- the boundary link is present
            # exactly when the boundary element of order_in is not a
            # Site index), so a transpose+reshape turns each theta
            # into a (kept x rest) matrix
            if ha == 1:
                nk = 1 if order_in[0].hastags("Site") else 2
                keep_inds = list(order_in[:nk])
                perm = list(range(len(order_in)))
            else:
                nk = 1 if order_in[-1].hastags("Site") else 2
                keep_inds = list(order_in[-nk:])
                perm = (list(range(len(order_in) - nk, len(order_in))) +
                        list(range(len(order_in) - nk)))
            rest_inds = [ind for ind in order_in if ind not in keep_inds]
            kd = int(np.prod([ind.dim for ind in keep_inds], dtype=int))
            Ml = thl.reshape(shape).transpose(perm).reshape(kd, -1)
            Mr = thr.reshape(shape).transpose(perm).reshape(kd, -1)
            rho = 0.5 * (Ml @ Ml.conj().T + Mr @ Mr.conj().T)
            if noise > 0:
                # reference's noiseterm() (cross term between the left
                # and right solutions through the half-contracted MPO),
                # hermitized before the hermitian eigensolver -- the
                # same construction as the C++ ports, done on matrices
                if ha == 1:
                    X = H.A(b) if Lh is None else Lh * H.A(b)
                else:
                    X = H.A(b + 1) if Rh is None else H.A(b + 1) * Rh
                thl_t = ITensor(tuple(order_in), thl.reshape(shape))
                thr_t = ITensor(tuple(order_in), thr.reshape(shape))
                # kept-side output legs of X, in the same order as
                # keep_inds' basis: (Lbra?, s_b') / (s_{b+1}', Rbra?)
                if ha == 1:
                    keep_out = ([Lbra_h] if Lbra_h is not None else []) \
                        + [keep_inds[-1].prime(1)]
                else:
                    keep_out = [keep_inds[0].prime(1)] \
                        + ([Rbra_h] if Rbra_h is not None else [])
                Xl = X * thl_t
                Xr = X * thr_t
                rest_out = [ind for ind in Xl.inds
                            if ind not in keep_out]
                Al = Xl.transpose_to(keep_out + rest_out).reshape(kd, -1)
                Ar = Xr.transpose_to(keep_out + rest_out).reshape(kd, -1)
                nt = Al @ Ar.conj().T
                rho = rho + (noise / 2.0) * (nt + nt.conj().T)
            Uk, keep = eigh_truncate(rho, cutoff, maxdim, mindim=1)

            link = Index(keep, tags="Link")
            keep_shape = tuple(ind.dim for ind in keep_inds)
            site_T = ITensor(tuple(keep_inds) + (link,),
                             Uk.reshape(keep_shape + (keep,)))
            rest_shape = tuple(ind.dim for ind in rest_inds)
            Cl = (Uk.conj().T @ Ml).reshape((keep,) + rest_shape)
            Cr = (Uk.conj().T @ Mr).reshape((keep,) + rest_shape)
            center_l = ITensor((link,) + tuple(rest_inds), Cl)
            center_r = ITensor((link,) + tuple(rest_inds), Cr)
            if ha == 1:
                psil.set_A(b, site_T)
                psir.set_A(b, site_T)
                psil.set_A(b + 1, center_l)
                psir.set_A(b + 1, center_r)
                psil.center = psir.center = b + 1
                left_h[b] = _extend_left(Lh, Lbra_h, H, psir, b, bra=psil)
                left_a[b] = _extend_left(La, Lbra_a, HT, psil, b, bra=psir)
            else:
                psil.set_A(b + 1, site_T)
                psir.set_A(b + 1, site_T)
                psil.set_A(b, center_l)
                psir.set_A(b, center_r)
                psil.center = psir.center = b
                right_h[b + 1] = _extend_right(Rh, Rbra_h, H, psir, b + 1, bra=psil)
                right_a[b + 1] = _extend_right(Ra, Rbra_a, HT, psil, b + 1, bra=psir)

    _conj_in_place(psil)  # transposed -> adjoint left state
    return energy


def _conj_chain(chain):
    """A copy of an MPS or MPO with every tensor complex-conjugated
    (tensor.dag, which has no Index arrows to reverse here): the adjoint
    MPO's conjugate is the transpose H^T, and an MPS's is the conjugate
    state on the same Index objects."""
    out = chain.copy()
    _conj_in_place(out)
    return out


def _conj_in_place(chain):
    """Conjugate every tensor of `chain` in place (its center stays)."""
    for i in range(1, chain.length() + 1):
        chain.set_A(i, dag(chain.A(i)))


def _biorthonormalize(psil, psir):
    """Scale the pair so that <psil|psir> = 1, splitting the factor
    between the two the way the per-bond normalization does (the real
    branch for a real negative overlap). Both end a sweep with their center
    at site 1, so the factor goes on that tensor. Inside the sweep the
    pair is normalized as the reference normalizes its (wl, r) pair, which
    leaves <psil|psir> = wl^T r rather than 1 (_nhdmrg_one_sweep). A
    vanishing overlap is the left and right states of two different
    eigenvalues, which no scale repairs; it raises RuntimeError, which
    nhdmrg.py's retry loop treats as a bad draw."""
    ov = complex(inner(psil, psir))
    if not (abs(ov) > 1e-14):
        raise RuntimeError(
            "nhdmrg: <psil|psir> collapsed to ~0 (or NaN): the left and "
            "right states belong to different eigenvalues")
    if abs(ov.imag) < 1e-14 * abs(ov):
        sq = np.sqrt(abs(ov.real))
        fl, fr = 1.0 / sq, 1.0 / (-sq if ov.real < 0 else sq)
    else:
        sq = np.sqrt(ov)
        fl, fr = 1.0 / np.conj(sq), 1.0 / sq
    psil.set_A(1, psil.A(1) * fl)
    psir.set_A(1, psir.A(1) * fr)


def nhdmrg(H, HA, psi0, sweeps, krylovdim=20, restarts=2, quiet=True):
    """Run NH-DMRG for the non-Hermitian MPO H (HA must be its adjoint,
    built from MultiOperator.get_dagger()'s terms). psi0 seeds both the
    right state and the sweep's transposed left state wl, as the
    reference starts from psil0 = psir0 = psi (so the adjoint psil starts
    as conj(psi0)). Returns (energy, psil, psir) with <psil|psir> = 1."""
    psir = psi0
    psir.position(1)
    psir.normalize()
    psil = _conj_chain(psir)

    energy = 0.0 + 0.0j
    have_energy = False
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, _niter = sweeps.at(sweep_i)
        energy = _nhdmrg_one_sweep(psil, psir, H, HA, maxdim, cutoff, noise,
                krylovdim, restarts, energy, have_energy)
        have_energy = True
        if not quiet:
            print("NH-DMRG sweep {}: energy = {}".format(sweep_i, energy))

    _biorthonormalize(psil, psir)
    # definitive energy: the biorthogonal Rayleigh quotient of the final pair
    energy = inner(psil, H, psir) / inner(psil, psir)
    return energy, psil, psir


def nhdmrg_generalized(H, HA, A, psi0, sweeps, krylovdim=20, restarts=2,
                       quiet=True, lam0=None):
    """Non-Hermitian generalized-eigenvalue NH-DMRG: finds the smallest-
    real-part lambda solving H|psi_R>=lambda*A|psi_R> (equivalently
    H^dagger|psi_L>=conj(lambda)*A|psi_L>, since A is Hermitian) for a
    possibly non-Hermitian H and a Hermitian positive-definite metric MPO
    A. Same self-consistent Lagrange-multiplier trick as dmrg.py's
    dmrg_generalized(), generalized to a complex lambda and biorthogonal
    (rather than plain) expectation values: minimizing (in the NH-DMRG
    sense -- smallest real part, not a variational bound) subject to the
    metric biorthogonal normalization <psi_L|A|psi_R>=1 gives stationarity
    condition (H-lambda*A)|psi_R>=mu|psi_R>,
    (H-lambda*A)^dagger|psi_L>=conj(mu)|psi_L> -- the *ordinary*
    (<psi_L|psi_R>=1-normalized) NH-DMRG eigenproblem of the shifted pair
    (H-lambda*A, HA-conj(lambda)*A), exactly what one ordinary NH-DMRG
    sweep already finds (note the adjoint shift uses conj(lambda), not
    lambda: (lambda*A)^dagger=conj(lambda)*A^dagger=conj(lambda)*A since A
    is Hermitian). At mu=0 this is precisely H|psi_R>=lambda*A|psi_R>, so
    each outer iteration here (i) builds Heff=H-lambda*A and
    HAeff=HA-conj(lambda)*A from the current lambda estimate, (ii) runs
    one ordinary NH-DMRG sweep (_nhdmrg_one_sweep, unmodified) against
    them, then (iii) updates lambda to the freshly-swept pair's
    generalized biorthogonal Rayleigh quotient
    <psi_L|H|psi_R>/<psi_L|A|psi_R>. A=identity reduces this exactly to
    plain nhdmrg().

    Each outer sweep's own within-sweep SRTieBreak anchor
    (_nhdmrg_one_sweep's energy_hint/have_energy_hint) restarts fresh
    (target 0, no anchor) rather than carrying over from the previous
    outer iteration: (H,HA) changes every outer iteration (unlike
    nhdmrg()'s own sweep loop, which reuses the same fixed H/HA
    throughout and so *does* carry the anchor across sweeps), so an
    anchor built for a different lambda's shifted spectrum has no
    particular reason to still be meaningful for the next one.

    lam0 (optional): starting lambda estimate; None or a NaN real part
    both mean "unset" (see this function's own note on why NaN is
    accepted too), defaulting to the seed state's own (plain, since
    psi_L=psi_R=psi0 initially) generalized Rayleigh quotient
    <psi0|H|psi0>/<psi0|A|psi0>.
    """
    psir = psi0
    psir.position(1)
    psir.normalize()
    psil = _conj_chain(psir)  # wl = psi0 in the sweep, as in nhdmrg()

    # None or NaN (real part) both mean "unset" -- NaN accepted too since
    # Chain::nhdmrg_generalized's own C++/pybind11 binding has no None
    # equivalent and uses a NaN real part as its own "unset" sentinel.
    # lam0_is_unset() (dmrg.py) coerces via complex(lam0) rather than
    # gating on isinstance(lam0, complex): a caller passing a bare float
    # NaN here (a reasonable thing to do given this docstring's own "None
    # or a NaN real part" wording) would otherwise silently skip this
    # fallback and feed a live NaN into the first outer iteration instead
    # -- confirmed directly, found via code review.
    lam = lam0
    if lam0_is_unset(lam):
        # the seed's own plain Rayleigh quotient, <psi0| on both sides
        a0 = inner(psir, A, psir)
        h0 = inner(psir, H, psir)
        lam = h0 / a0 if abs(a0) > 1e-14 else 0.0 + 0.0j
    lam = complex(lam)

    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, _niter = sweeps.at(sweep_i)
        Heff = mps_sum(H, A * (-lam))
        HAeff = mps_sum(HA, A * (-np.conj(lam)))
        _nhdmrg_one_sweep(psil, psir, Heff, HAeff, maxdim, cutoff, noise,
                krylovdim, restarts, 0.0 + 0.0j, False)
        a_psi = inner(psil, A, psir)
        if not (abs(a_psi) > 1e-14):
            # "not (> tol)", not "< tol": see dmrg.py's identical guard
            # for why the negated form is needed to also catch NaN.
            raise RuntimeError(
                "nhdmrg_generalized: <psi_L|A|psi_R> collapsed to ~0 (or "
                "NaN) mid-iteration (A may not be positive definite, or "
                "the sweep drove the biorthogonal pair toward A's near-"
                "null-space) -- cannot form a meaningful generalized "
                "Rayleigh quotient")
        h_psi = inner(psil, H, psir)
        lam = h_psi / a_psi
        if not quiet:
            print("sweep {}: lambda = {}".format(sweep_i, lam))

    _biorthonormalize(psil, psir)
    return lam, psil, psir
