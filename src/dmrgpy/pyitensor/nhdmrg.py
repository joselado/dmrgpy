"""Non-Hermitian two-site DMRG (NH-DMRG) for the pure-Python backend.

Port of mpscpp3/chain_session.h's Chain::nhdmrg (the annotated original --
see its long comment for the algorithm, its ITensorNHDMRG.jl reference,
and the two deliberate deviations) to this package's engine:

- the same "onesided" local solver (independent Arnoldi solves of the
  projected two-site eigenproblems for H and Hdag, smallest real part,
  with the adjoint solve anchored to the right solve's eigenvalue and
  Re-degenerate Ritz values tie-broken toward the previous bond's), and
- the same "fidelity" truncation (one shared isometry from the hermitian
  average rho = (rho_l + rho_r)/2 of the left/right two-site reduced
  density matrices), which keeps psil and psir on identical site *and*
  link Index objects throughout -- so the environments hit the exact same
  bra/ket link-identity collision dmrg.py's self-overlap environments do,
  and reuse its _relabel_bra_local machinery, just with the bra tensors
  drawn from the *other* MPS.

Everything bond-local is done on flat numpy arrays (via dmrg.py's
two_site_heff matvec factory): the Arnoldi iteration, the biorthogonal
normalization, and the fidelity truncation (a dense eigh of the small
rho matrix, truncated with svd.py's own _truncate rule). Unlike dmrg.py
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
"""

import numpy as np

from .dmrg import (_relabel_bra_local, two_site_heff)
from .index import Index
from .mpsalgebra import inner
from .svd import _truncate
from .tensor import ITensor, contract_many


def _extend_left_nh(L, left_bra, W, ket, bra, i):
    """One more site of the two-sided left environment <bra|W|ket>: same
    as dmrg.py's _extend_left except the bra tensor comes from a different
    MPS than the ket (they share link Index objects, so the relabeling
    dance is identical)."""
    bra_piece, _, right_bra = _relabel_bra_local(bra.A(i), bra, i, left_bra, None)
    pieces = [p for p in (L, bra_piece, W.A(i), ket.A(i)) if p is not None]
    return contract_many(pieces), right_bra


def _extend_right_nh(R, right_bra, W, ket, bra, i):
    bra_piece, left_bra, _ = _relabel_bra_local(bra.A(i), bra, i, None, right_bra)
    pieces = [p for p in (R, bra_piece, W.A(i), ket.A(i)) if p is not None]
    return contract_many(pieces), left_bra


def _all_right_nh(W, ket, bra):
    n = ket.length()
    env = {n + 1: (None, None)}
    for i in range(n, 1, -1):
        R_next, bra_next = env[i + 1]
        env[i] = _extend_right_nh(R_next, bra_next, W, ket, bra, i)
    return env


def _arnoldi_smallest_real(matvec, x0, krylovdim, restarts,
                           sel="SR", target=0.0):
    """Restarted Arnoldi over flat numpy vectors; keeps the Ritz pair
    selected by `sel` ("SR" = smallest real part, "closest" = closest to
    `target`, "SRTieBreak" = smallest real part with Re-degenerate values
    tie-broken toward `target`). Mirrors chain_session.h's
    arnoldi_smallest_real (see the Sel rationale there)."""
    lam = 0.0 + 0.0j
    x0 = np.asarray(x0, dtype=complex)
    for _ in range(restarts):
        nx = np.linalg.norm(x0)
        if nx < 1e-14:
            break
        q = x0 / nx
        Q = [q]
        m = min(krylovdim, x0.size)
        h = np.zeros((m + 1, m), dtype=complex)
        built = 0
        for j in range(m):
            w = matvec(Q[j])
            for i in range(j + 1):
                c = np.vdot(Q[i], w)
                h[i, j] = c
                w = w - c * Q[i]
            for i in range(j + 1):  # one re-orthogonalization pass
                c = np.vdot(Q[i], w)
                h[i, j] += c
                w = w - c * Q[i]
            built = j + 1
            nw = np.linalg.norm(w)
            h[j + 1, j] = nw
            if nw < 1e-13:
                break  # happy breakdown: invariant subspace
            if j + 1 < m:
                Q.append(w / nw)
        evals, evecs = np.linalg.eig(h[:built, :built])
        if sel == "closest":
            k = int(np.argmin(np.abs(evals - target)))
        elif sel == "SRTieBreak":
            remin = evals.real.min()
            degtol = 1e-6 * (1.0 + abs(remin))
            cand = np.flatnonzero(evals.real < remin + degtol)
            k = int(cand[np.argmin(np.abs(evals[cand] - target))])
        else:
            k = int(np.argmin(evals.real))
        lam = complex(evals[k])
        x0 = np.column_stack(Q[:built]) @ evecs[:, k]
    return lam, x0


def nhdmrg(H, HA, psi0, sweeps, krylovdim=20, restarts=2, quiet=True):
    """Run NH-DMRG for the non-Hermitian MPO H (HA must be its adjoint,
    built from MultiOperator.get_dagger()'s terms). psi0 seeds both the
    left and right MPS (trivially biorthonormal, like the C++ ports).
    Returns (energy, psil, psir) with <psil|psir> = 1."""
    psir = psi0
    psir.position(1)
    psir.normalize()
    psil = psir.copy()
    n = psir.length()

    energy = 0.0 + 0.0j
    have_energy = False
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, _niter = sweeps.at(sweep_i)

        # (env tensor, dangling bra link) per boundary, one family per
        # projected problem: H sandwiched between <psil| and |psir>, and
        # Hdag sandwiched between <psir| and |psil>
        right_h = _all_right_nh(H, psir, psil)
        right_a = _all_right_nh(HA, psil, psir)
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
                mv_l, _order_l, _shape_l, x0l = two_site_heff(
                    La, Lbra_a, HA, psil, b, Ra, Rbra_a)
                er, thr = _arnoldi_smallest_real(
                    mv_r, x0r, krylovdim, restarts,
                    sel="SRTieBreak" if have_energy else "SR", target=energy)
                _el, thl = _arnoldi_smallest_real(
                    mv_l, x0l, krylovdim, restarts,
                    sel="closest", target=np.conj(er))
                energy = er
                have_energy = True

                thl = thl / np.linalg.norm(thl)
                thr = thr / np.linalg.norm(thr)
                # rescale so <thl|thr> = 1, split between the two states
                # (separate real branch: see chain_session.h's note on
                # sqrt's branch cut for real negative overlaps)
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
                w, U = np.linalg.eigh(rho)
                idx = np.argsort(w)[::-1]
                w = w[idx]
                U = U[:, idx]
                svals = np.sqrt(np.clip(w, 0.0, None))
                keep, _disc = _truncate(svals, cutoff, maxdim, 1)
                Uk = U[:, :keep]

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
                    left_h[b] = _extend_left_nh(Lh, Lbra_h, H, psir, psil, b)
                    left_a[b] = _extend_left_nh(La, Lbra_a, HA, psil, psir, b)
                else:
                    psil.set_A(b + 1, site_T)
                    psir.set_A(b + 1, site_T)
                    psil.set_A(b, center_l)
                    psir.set_A(b, center_r)
                    psil.center = psir.center = b
                    right_h[b + 1] = _extend_right_nh(Rh, Rbra_h, H, psir, psil, b + 1)
                    right_a[b + 1] = _extend_right_nh(Ra, Rbra_a, HA, psil, psir, b + 1)

        if not quiet:
            print("NH-DMRG sweep {}: energy = {}".format(sweep_i, energy))

    # definitive energy: the biorthogonal Rayleigh quotient of the final pair
    energy = inner(psil, H, psir) / inner(psil, psir)
    return energy, psil, psir
