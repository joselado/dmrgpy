"""Krylov global subspace expansion (GSE) for one-site TDVP -- the
Yang-White scheme (arXiv:2005.06104/Phys. Rev. B 102, 094315) that lets
one-site TDVP (which conserves bond dimension exactly on its own) grow
bond dimension the way two-site TDVP does via SVD. Mirrors
mpscpp3/TDVP/basisextension.h's addBasis()/addBasisWorker()/
denmatSumDecomp(): enrich phi's local bases, one bond at a time sweeping
right-to-left, with new directions spanned by a Krylov subspace
{H*phi, H^2*phi, ...} built via repeated MPO application -- but implemented
directly against pyitensor's dense-array ITensor representation (a
"combiner"/"plusser" in ITensor's own sparse-block formalism is just a
plain reshape/concatenate on the dense array here) rather than transliterating
ITensor's own sparse-tensor bookkeeping.

The state-preserving property (phi's own state comes out *exactly*
unchanged -- only its local bond bases get bigger, giving one-site TDVP
more room to grow into) rests on one linear-algebra fact, applied at every
bond: split phi's local frontier tensor via SVD as M = U1 @ S1 @ V1^dagger
(V1's rows orthonormal, spanning the "combined" physical(+already-enlarged
link) space phi's own state currently occupies). Any *new* basis
directions U2 appended orthogonal to V1's own row space automatically get
exactly-zero coefficient in phi's own representation, since
M @ [V1^dagger|U2] 's second block is U1@S1@V1^dagger@U2 = 0 (U2 built
orthogonal to V1 by construction). U2 is chosen as the dominant
eigenvectors (above `cutoff`) of the companions' local density matrix,
projected onto the orthogonal complement of V1's span -- exactly
denmatSumDecomp()'s rho2/proj2 construction, just spelled out with a
dense eigh instead of ITensor's diag_hermitian/combiner/plusser.

Where the arrays live
---------------------
Every array this sweep touches is a frontier tensor of phi or of one of
its Krylov companions, so it arrives already on whatever device
backend.py has selected and it must stay there. Three calls used to take
it off: `np.concatenate` of V1's rows with U2's new ones, which returns a
host array from device inputs with no error at all, so the enlarged
tensor res.A(b) was rebuilt from the host once per bond, and the two
`np.linalg.norm` calls that decide whether any new direction is worth
adding, each of which pulled the whole (combined, combined) matrix home
to produce one number. Measured on a 6-site spin chain at bond dimension
8, with the array type recorded at every ITensor construction: 5 of the
229 tensors a single expansion built came from a host array, exactly one
per bond, and the count is 0 now. A `backend.to_host` counter cannot see
any of this, since none of the three went through it; the array type is
what shows it.

What stays on the host is the same split the rest of the engine uses:
`combined` is index bookkeeping on Python ints, and the two norms cross
as scalars rather than as matrices. They cross one at a time rather than
in one stacked transfer, so the second norm is still left uncomputed when
the first is zero, which is what the host path has always done. Two
synchronizations per bond is affordable here in a way it would not be
inside a Krylov recursion: an expansion runs only for the leading
`tdvp_gse_sweeps` steps of an evolution (3 by default), not once per time
step.
"""

import numpy as np

from . import backend as bk
from .index import Index
from .mpsalgebra import applyMPO, _link_at
from .svd import svd, eigh_truncate
from .tensor import ITensor


def _gse_bond_step(B_phi, left_link, B_companions, right_inds, cutoff, bond_maxdim=None):
    """One bond of the GSE sweep (right-to-left). B_phi: phi's current
    frontier tensor, indices = left_link (the link to site b-1, about to
    be finalized) + right_inds (physical(b), plus -- for every bond but
    the first processed -- the already-enlarged link built by the
    previous iteration). B_companions: the same-shape frontier tensors
    for each Krylov companion (sharing the identical right_inds Index
    objects, but each with its own, generally different, left-link
    Index). Returns (new_res_b, new_B_phi, new_B_companions): new_res_b is
    the enlarged tensor to write as res.A(b); new_B_phi/new_B_companions
    are each frontier tensor projected into the new basis, ready for the
    caller to contract against site (b-1)'s own original tensor to build
    next bond's input."""
    _xp = bk.xp()
    # Index bookkeeping, on Python ints rather than on tensor data, so it
    # stays on the host on every backend.
    combined = int(np.prod([ind.dim for ind in right_inds])) if right_inds else 1

    # phi's own currently-occupied subspace of the combined (physical [+
    # already-enlarged-link]) space: V1's rows, orthonormal by
    # construction (cutoff=0/maxdim=None -> lossless split, only exactly-
    # zero singular values ever get dropped).
    _, S1, V1, spec1 = svd(B_phi, [left_link], cutoff=0.0, maxdim=None)
    bond_v = next(ind for ind in V1.inds if ind not in right_inds)
    m = bond_v.dim
    V1_mat = V1.transpose_to([bond_v] + right_inds).reshape(m, combined)
    # m is phi's RANK here, not the bond dimension svd() hands back, which
    # under backend.set_pad_bonds is the pad width K: the padded rows of V1
    # are zero, and counting them left `room = bond_maxdim - m` at 0 at the
    # recommended K=maxm, so no direction was ever added (2026-09-24 audit,
    # finding 16). The live rows come first, and the spectrum stops at the
    # true rank (svd.py), so this is a static slice on a host-side count.
    # Unpadded the two counts are equal and nothing here runs. The route
    # through Chain.global_subspace_expand suspends padding altogether;
    # this keeps a direct caller of this module right as well.
    m_true = len(spec1.eigs())
    if m_true < m:
        V1_mat = V1_mat[:m_true]
        m = m_true

    rho2 = bk.zeros((combined, combined))
    Bk_mats, Bk_lefts = [], []
    for Bk in B_companions:
        left_k = next(ind for ind in Bk.inds if ind not in right_inds)
        Bk_mat = Bk.transpose_to([left_k] + right_inds).reshape(left_k.dim, combined)
        Bk_mats.append(Bk_mat)
        Bk_lefts.append(left_k)
        rho2 += Bk_mat.conj().T @ Bk_mat

    # Project the companions' density matrix onto the orthogonal
    # complement of phi's own subspace, so any kept direction is
    # guaranteed new (and thus zero-weight in phi's own tensor below).
    proj = bk.eye(combined) - V1_mat.conj().T @ V1_mat
    rho2_proj = proj @ rho2 @ proj
    rho2_proj = 0.5 * (rho2_proj + rho2_proj.conj().T)  # enforce exact Hermiticity

    U2 = bk.zeros((combined, 0))
    # Two numbers per bond come home, the same split the rest of the engine
    # uses: the norms are computed where rho2 lives and only the scalars
    # cross. They are read one at a time rather than stacked into a single
    # transfer, so the second norm is still not computed when the first is
    # zero, which is the order the host path has always evaluated them in.
    norm_rho2 = float(bk.to_host(_xp.linalg.norm(rho2)))
    if norm_rho2 > 0 and float(bk.to_host(_xp.linalg.norm(rho2_proj))) / norm_rho2 >= 1e-12:
        # cap the *total* new_dim=m+keep2 at bond_maxdim, not just keep2
        # alone -- without this, only `cutoff` bounds how many new
        # directions get added at any one bond, and for a large-enough
        # system with genuinely entangled Krylov vectors that count can
        # keep exceeding the caller's own bond-dimension budget at every
        # bond of every GSE call with no ceiling (confirmed directly:
        # n=20, default settings, bond dimension ran to 370 on the v3
        # backend before being capped there the same way).
        room = None if bond_maxdim is None else max(0, bond_maxdim - m)
        U2, _keep2 = eigh_truncate(rho2_proj, cutoff, room, mindim=0)

    # new_res_mat: (new_dim, combined), ROWS orthonormal -- V1's own m rows
    # (phi's existing subspace, unchanged) followed by keep2 brand-new rows
    # (U2's columns, conjugate-transposed to sit alongside V1_mat's rows).
    # This mirrors svd.py's own "V"/"Vh" convention (a right-orthogonal
    # factor whose rows -- not columns -- are orthonormal), so every
    # *ordinary* (non-conjugating) `*` contraction elsewhere in this
    # package that chains an MPS tensor against its neighbor keeps working
    # unchanged: the "coefficient" factor propagated to the next bond is
    # `M @ new_res_mat^dagger` (the standard "U@S = M @ V^dagger" identity
    # for M = U@S@V), not `M @ new_res_mat` -- getting this backwards
    # (contracting via `basis` instead of `basis^dagger`, an earlier bug
    # here) silently corrupts the state by a phase/rotation that only
    # shows up once compounded with a *different* code path's own tensors
    # (e.g. TDVP's), since a single GSE call's own self-consistency check
    # can't detect a systematic conjugation mismatch.
    new_res_mat = _xp.concatenate([V1_mat, U2.conj().T], axis=0)  # (new_dim, combined)
    new_dim = new_res_mat.shape[0]
    new_link = Index(new_dim, tags="Link")

    new_res_arr = new_res_mat.reshape((new_dim,) + tuple(ind.dim for ind in right_inds))
    new_res_b = ITensor((new_link,) + tuple(right_inds), new_res_arr)

    new_res_dag = new_res_mat.conj().T  # (combined, new_dim)
    Bphi_mat = B_phi.transpose_to([left_link] + right_inds).reshape(left_link.dim, combined)
    new_B_phi = ITensor((left_link, new_link), Bphi_mat @ new_res_dag)

    new_B_companions = [ITensor((Bk_lefts[k], new_link), Bk_mats[k] @ new_res_dag)
                         for k in range(len(B_companions))]

    return new_res_b, new_B_phi, new_B_companions


def global_subspace_expand(H, phi, krylov_order, cutoff, maxdim=None, bond_maxdim=None):
    """Grow phi's bond dimension using a krylov_order-dimensional Krylov
    subspace {phi, H*phi, ..., H^(krylov_order-1)*phi} built from the MPO
    H, discarding directions below `cutoff` (maxdim caps each individual
    Krylov-vector applyMPO application, matching
    Chain::global_subspace_expand()'s own maxdim semantics; bond_maxdim
    hard-caps the *enlarged* bond dimension itself at every bond,
    matching what Chain::global_subspace_expand() passes as "MaxDim" to
    addBasis()'s own density-matrix truncation on the v3/mpscpp3 side --
    without it, only `cutoff` bounds bond growth, which for a
    large-enough system can keep exceeding any practical budget at every
    GSE call). krylov_order<=1 is a no-op. Returns a new MPS representing
    the *same* state as phi (see module docstring), with (generally)
    larger bond dimensions -- pair with one-site TDVP
    (tdvp_step(...,num_center=1)), which alone cannot grow bond
    dimension."""
    if krylov_order <= 1:
        return phi.copy()

    companions = []
    cur = phi
    for _ in range(krylov_order - 1):
        cur = applyMPO(H, cur, cutoff=cutoff, maxdim=maxdim)
        cur.noPrime("Site")  # applyMPO leaves the physical leg primed
        companions.append(cur)

    n = phi.length()
    res = phi.copy()
    res.position(n)
    for c in companions:
        c.position(n)

    B_phi = res.A(n)
    B_companions = [c.A(n) for c in companions]

    for b in range(n, 1, -1):
        left_link = _link_at(res, b - 1, b)
        right_inds = [ind for ind in B_phi.inds if ind != left_link]
        new_res_b, new_B_phi, new_B_companions = _gse_bond_step(
                B_phi, left_link, B_companions, right_inds, cutoff, bond_maxdim)
        res.set_A(b, new_res_b)
        B_phi = new_B_phi * res.A(b - 1)
        B_companions = [Ck * comp.A(b - 1)
                         for Ck, comp in zip(new_B_companions, companions)]

    res.set_A(1, B_phi)
    res.center = 1
    return res
