"""Tangent-space Krylov (TaSK) spectral functions on top of a ground-state MPS
(Kovalska, von Delft and Gleis, arXiv:2510.07279).

The idea
--------
A Lanczos recursion on MPS (Hallberg's continued fraction, or the global
Krylov of rootndmrg.py) has to compress after every H|v>, and the
compression destroys the orthogonality of the Krylov basis. TaSK runs the
recursion inside the one-site tangent space of the ground state instead,
which is an honest vector space: a vector is

    |Psi> = sum_l A_1 ... A_{l-1} T_l B_{l+1} ... B_L

with A_l (B_l) the left (right) isometries of the ground state and T_l in
the discarded space of A_l, sum_{a,s} conj(A_l[a,s,m]) T_l[a,s,r] = 0. The
gauge makes the L terms mutually orthogonal and orthogonal to |GS>, so
overlaps, sums and scalings act on the T_l alone, exactly, at O(L D^2 d).
The spectral function

    S(w) = <GS|A delta(w - H + E0) B|GS>

is approximated by projecting H - E0 and B|GS> onto that space,
H^{1perp} = P (H - E0) P and |b> = P B|GS>, and running an ordinary Lanczos
with full reorthogonalization. The Ritz pairs (w_alpha, |alpha>) give the
discrete poles and weights; no broadening enters the solve. The one
approximation is the projection P, which is exact only if the ground-state
bond dimension saturates the Hilbert space; otherwise an excitation that
needs more than a one-site change of the ground state (a multi-spinon
state, for instance) is represented through the tangent vectors closest
to it. Its size is estimated per pole by the two-site variance
||P^{2perp} H |alpha>||^2 (Hubig et al., Gleis et al.), the part of
H|alpha> that falls into the two-site discarded space.

The elastic term n = GS is not in the tangent space (which is orthogonal
to |GS>) and is returned separately as <GS|A|GS><GS|B|GS> at w = 0.

Padded bonds
------------
The size of the tangent space, sum_l (a_l d_l - b_l) b_l, is set by the
bond dimensions of the ground-state MPS, so a state of small Schmidt rank
has a small one, whatever its accuracy. Two cases where this is the whole
error: a product state (bond dimension 1, the J = 0 Kondo chain), where
the composite STM vertex, which flips the impurity spin and adds an
electron to the bath at once, is not a one-site change of |GS> and TASK
keeps a third of the weight; and a state solved in one quantum-number
sector and promoted to dense indices, whose bonds are those of one member
of a degenerate pair, smaller than those of the mixture a plain solve
returns. The paper notes the same limit for a state with an exact finite
bond dimension (AKLT) and leaves a tangent-space expansion as future work.

pad_bonds is the simplest such expansion. Every interior bond is padded
with zero slices on both neighbouring tensors, from chi_l to
max(chi_l, min(D, b_l)), b_l = min(d_1...d_l, d_{l+1}...d_L) being its
structural bound, or to b_l itself for pad_bonds="full". The state is
unchanged; the QR and LQ sweeps of mixed_gauge complete the rank-deficient
bonds with orthonormal directions, so A_l and B_l stay isometries, both
gauge relations hold, and every kept space still contains the support of
the state on its side of the bond (it contains the unpadded kept space
whenever that bond had full rank). The tangent space then contains
the unpadded one and the Lanczos is variational in a larger space. The
two-site variance is another matter: with "full" it is zero, since the
two-site discarded space is empty, but with an int below the bound the
Ritz vectors mix in the arbitrary padded directions and the variance is
inflated by them, so it is no longer a sharp error bar of the poles (on
d22, the S = 1 pair below at pad_bonds=16: a correlator exact to 1e-15
before and after padding goes from a largest pole variance of 2.6e-7 to
6.6e-4, and the one whose error falls from 8.7e-4 to 2.5e-10 keeps a
weight-averaged variance of 4.1e-6). With "full" every kept space is the
whole space of its side of the bond, the tangent space is the whole
Hilbert space orthogonal to |GS> (4^6 - 1 = 4095 on a six-site chain of
d = 4) and TASK is the Lanczos of exact diagonalization, which only fits
on short chains. With an int D below the bound the added directions are
the arbitrary completion of the QR (deterministic, but not chosen from
H|GS> or B|GS>), so their value on a long chain is not controlled; a
controlled expansion in the spirit of CBE would pick them from the
two-site discarded space.

Checked on d22 against ED CVM (tests/test_task_bond_padding.py): the
sector state of the degenerate pair S = 1, Dz = -0.05 on 4 bath sites
goes from a tangent space of 562 and an error of 8.7e-4 to 4095 and
4.9e-8 ("python"), 2.6e-8 (v3), and S = 3/2, Dz = -0.04 on 3 sites from
205 and 4.8e-2 to 4095 and 1.6e-9, 8.9e-12; padding to 16 already gives
4.9e-8 and 1.6e-9 there at a tangent space of 1791. At J = 0 TASK goes
from 1/3 of the weight (U = 0) and 0.405 (U = 0.3) to ED within 1.4e-14
("python") and 3.2e-12 (v3). A chain whose spectrum is already converged
moves by at most 4e-12, and one whose bonds are at their bound does not
move at all.

On a 21-site Wilson chain of a BCS band (v3, maxm 128, Delta = 1e-3,
Lambda = 4) padding to 128 and 256 does not help: at nkry = 100 the
subgap pole is 1.1e-5 above the sector energy difference unpadded and
4.2e-5 and 4.3e-5 above it padded, the four spin and particle-hole
copies of it split by 1e-6, at 2.4 and 15.6 times the cost; at nkry =
300 unpadded and 128 both agree with it to the rounding of the
reference (5e-9), padded at 2.6 times the cost, with the lowest
continuum pole at 1.154 Delta against 1.055 Delta unpadded. The
arbitrary directions slow the Krylov convergence and break the spin
symmetry of the tangent space, so an int pad_bonds is for short chains.

Expanded bonds
--------------
expand_bonds=k is the controlled expansion. The Lanczos from b = P B|GS>
loses, at its start, the part of B|GS> outside the tangent space, and at
its first step the part of H b outside it. Both have a two-site part,
P^{2perp} B|GS> and P^{2perp} H b, which lives in the two-site discarded
space of a bond: outside the kept columns of A_l on the left and the kept
rows of B_{l+1} on the right. On every bond the two-site tensors of these
vectors, each normalized over the chain, are stacked side by side and on
top of each other and reshaped to (a d) x (d' c); their leading left
singular vectors become at most k new columns of A_l, their leading right
singular vectors at most k new rows of B_{l+1}, and B_l and A_{l+1} get
zero columns and zero rows, C_l and C_{l+1} zeros, so that A_l and B_l
keep the same shapes. The state, E0 and both gauge relations do not
change, and the tangent space contains the unexpanded one and every
product of a new direction with anything on the neighbouring site, so
the two-site directions that carry the dropped weight are in it. This is
the two-site expansion of CBE (Gleis, Li and von Delft, arXiv:2207.14712)
fed by B|GS> and H P B|GS> instead of H|GS>. MPSKit.jl's OptimalExpand
is the closest reference code: in a left-to-right sweep it projects the
two-site update on both null spaces, takes a truncated SVD and appends the
right singular vectors to AR, completing AL by a zero slice and a gauge
step; here both sides get their own singular vectors, so that A_l and B_l
keep the same shapes; the shrewd selection of CBE, a cost shortcut, is
not needed at these bond dimensions, the full two-site tensor costs one
two-site variance per feed. The bonds are expanded in one left-to-right
sweep, and the two-site tensor of bond l is taken with the bras already
expanded on bonds 0..l-1 (the kets are the unexpanded state and b), so a
direction added on one bond is seen by the two-site tensor of the next.
Singular values below EXPANSION_FLOOR = 1e-8 of the largest weight of a
bond (from a first pass in the unexpanded gauge) are not used, a feed
whose two-site norm is below EXPANSION_FLOOR of the norm of its vector is
dropped (on v3 the rounding of a saturated bond is about 1e-10 of it,
and a looser cut added arbitrary directions), and the
largest possible dimension, from the shapes alone, is checked against
MAX_PADDED_DIMENSION before any SVD. On a real gauge a complex feed is
split into its real and imaginary parts, so the recursion stays real. The
SVDs and the assembly are xp calls on the active backend
(concatenations, no in-place write); only the singular values come back
to the host. "hamiltonian" (P^{2perp} H|GS>, the CBE-DMRG direction) is a
third feed of TangentSpace.expand, off by default. No charge blocks are
used: the SVD can mix degenerate directions of different Sz.

Checked on d22 against ED CVM (tests/test_task_expansion.py, run
short_all): the sector state of the degenerate pair S = 1, Dz = -0.05 on
4 bath sites goes from an error of 1.6e-3 to 1.3e-11 ("python") and
1.7e-11 (v3) at every k from 1 to 16, two of its four correlators growing
from a tangent space of 562 to 658 ("python") and 621 (v3) and the other
two having nothing to add; pad_bonds=16 gives 2.4e-8 at 1791. At J = 0
the whole weight comes back from k = 1 (a third of it unexpanded), but
the poles are not those of ED (0.39 off against 0.27 unexpanded): the
flipped impurity is added on the first bond, and the bath electron that
the Lanczos then moves along the chain needs it on the next bonds too,
which b and H b of the unexpanded state do not reach. expand_rounds
repeats the expansion from the expanded state, where b already holds the
flipped impurity: with two rounds at k = 8 the J = 0 spectrum is ED to
3.1e-14 ("python", run diag3) at a tangent space of 515, against 1023
for pad_bonds="full", and a third round adds nothing. On v3 two rounds
at k >= 4 give 2.8e-12, and k = 2 needs three rounds on both backends
(run rounds). Each round costs one first pass and one sweep, and the cap
is checked per round. On a converged chain both feeds are below the
rounding cut and nothing moves. The two-site variance stays an error
bar: on the S = 1
pair the largest variance of a pole of weight above 1e-3 of the largest
falls with the error, from 4.2e-2 to 8.7e-5 ("python") and from 4.2e-2
to 9.7e-5 (v3) at k = 4 (run diag3).

On a long chain the expansion does not help. On kondo_chain(20, J = 0.3
A_Lambda, Delta = 1e-3, Lambda = 4) on v3 at maxm 128 (run long16 on d22,
the four runs interleaved in one process, E_ysr = 2.4864516231e-4 from the
sector solves), expand_bonds=16 adds 192 directions, the tangent space
grows from 298751 to about 391000, and the subgap pole is 1.81e-5 above
E_ysr at nkry = 100 against 1.08e-5 unexpanded, at 1.9 times the time
(244 s against 129 s). At nkry = 300 both reach E_ysr, 9.3e-11 expanded
and 3.7e-11 unexpanded, at 760 s against 518 s, and the lowest continuum
pole (weight above 1e-6 of the largest) moves from 1.055 Delta to 1.103
Delta. There the unexpanded tangent space already holds the subgap state,
and what nkry = 100 lacks is Krylov convergence, which a larger space
cannot give. The four spin and particle-hole copies of the subgap pole
stay within 1.5e-8 of each other (1.7e-14 at nkry = 300), so the SVD does
not break the spin symmetry here. So expand_bonds is for a state whose
tangent space is too small, not for a slow Krylov.

Off-diagonal pairs
------------------
The paper treats A = B^dagger. For a general pair the weight of a pole is
<GS|A|alpha><alpha|B|GS> = conj(<alpha|a>) <alpha|b> with |a> = P A^dagger|GS>,
and <alpha|a> is evaluated by expanding |a> in the Krylov basis built from
|b>. That is one approximation more than the paper makes: the part of |a>
outside the Krylov space of |b> is dropped, which converges away with
nkry like any Krylov resolvent but is not controlled by the variance.

Everything here is dense arrays of the active pyitensor backend
(backend.xp(): NumPy, or JAX on a device) on (chi_left, d, chi_right) site
arrays and (w_left, d_out, d_in, w_right) MPO arrays; the ITensor objects
are only read at the entry point. The site and MPO arrays arrive on the
host (mpo_arrays, chain._mps_arrays_lpr, the v3 export) and real_if_exact
is the one place they are moved to the backend, so on JAX the environments,
the Krylov basis and every matvec stay on the device and only scalars and
the small Krylov-space matrices come back.

Only one array of the size of the Krylov basis, Q (one row per tangent
vector), is ever alive: H Q, the Ritz vectors U^T Q and O Q are never
formed, their overlaps with Q are accumulated column by column instead.
"""

import numpy as np

from . import backend


# ---------------------------------------------------------------- arrays

def mpo_arrays(H):
    """Every site tensor of the MPO `H` as a (w_left, d_out, d_in, w_right)
    array, d_out the primed (bra) leg and d_in the unprimed (ket) one, with
    length-1 axes at the open ends."""
    from .tensor import commonIndex
    n = H.length()
    out = []
    for i in range(1, n + 1):
        W = H.A(i)
        sout = next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 1)
        sin = next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 0)
        left = commonIndex(H.A(i - 1), W) if i > 1 else None
        right = commonIndex(W, H.A(i + 1)) if i < n else None
        order = [ind for ind in (left, sout, sin, right) if ind is not None]
        arr = np.asarray(W.transpose_to(order))
        if left is None:
            arr = arr.reshape((1,) + arr.shape)
        if right is None:
            arr = arr.reshape(arr.shape + (1,))
        out.append(np.ascontiguousarray(arr, dtype=complex))
    return out


def real_if_exact(arrays):
    """The arrays as float64 when every imaginary part is exactly zero,
    unchanged otherwise. A real Hamiltonian, a real ground state and a real
    operator then run the whole recursion in real arithmetic, about four
    times fewer flops per contraction; a complex start vector on real
    environments is still correct, numpy promotes it."""
    xp = backend.xp()
    arrays = [backend.to_host(a) for a in arrays]
    if all(not np.iscomplexobj(a) or not np.any(a.imag) for a in arrays):
        return [xp.asarray(np.ascontiguousarray(np.real(a), dtype=float))
                for a in arrays]
    return [xp.asarray(np.asarray(a, dtype=complex)) for a in arrays]


def real_up_to_phase(Ws, tol=1e-13):
    """(phase, real arrays) when the MPO arrays Ws are a global phase times
    a real MPO site by site, every site tensor a phase times a real tensor
    to tol of its largest entry; (1.0, Ws) unchanged otherwise.

    compress_mpo writes S_y = i times a real operator as complex site
    tensors with a phase on each (from its SVDs), so on a real basis every
    application of it promoted the real rows and ran in complex arithmetic.
    The product of the site phases is snapped to 1, -1, i or -i when it is
    one of them to rounding, so that i S_y comes out exactly real. The
    tensors are small (bond dimension 1 for a one-site operator) and are
    split on the host."""
    xp = backend.xp()
    host = [backend.to_host(W) for W in Ws]
    if all(not np.iscomplexobj(W) for W in host):
        return 1.0, list(Ws)
    phase, out = 1.0 + 0j, []
    for W in host:
        k = int(np.argmax(np.abs(W)))
        top = complex(W.flat[k])
        if top == 0.0:
            return 1.0, list(Ws)
        p = top / abs(top)
        R = W / p
        if float(np.max(np.abs(R.imag))) > tol * abs(top):
            return 1.0, list(Ws)
        phase *= p
        out.append(np.ascontiguousarray(R.real))
    for exact in (1.0, -1.0, 1j, -1j):
        if abs(phase - exact) < 1e-12:
            phase = exact
            break
    return phase, [xp.asarray(R) for R in out]


def real_mpo(Ws, tol=1e-12):
    """The MPO arrays of the same operator as Ws in a real gauge when that
    operator is real, None otherwise.

    The "python" AutoMPO writes S_y S_y with imaginary site tensors, and its
    compression mixes the channels, so the MPO of a real H is complex entry
    by entry. Bond by bond from the left, the left block M (written in the
    real basis of the bonds before) is replaced by a real orthonormal basis
    Q of the span of [Re M, Im M], and Q^T M is moved into the next tensor,
    which keeps the operator to tol. The span of a real operator's minimal
    left blocks is closed under conjugation, so Q has the bond dimension of
    M; a larger span (a complex operator, or a bond that is not minimal)
    returns None, and so does an imaginary part left on the last tensor,
    which the real bases above make the imaginary part of the operator.
    Done on the host: the tensors are those of an MPO."""
    xp = backend.xp()
    host = [np.asarray(backend.to_host(W)) for W in Ws]
    if all(not np.iscomplexobj(W) or not np.any(W.imag) for W in host):
        return [xp.asarray(np.ascontiguousarray(np.real(W))) for W in host]
    G = np.ones((1, 1))
    out = []
    for l, W in enumerate(host):
        X = np.tensordot(G, W, axes=(1, 0))
        a, d1, d2, b = X.shape
        if l == len(host) - 1:
            if float(np.max(np.abs(X.imag))) > tol * float(np.max(np.abs(X))):
                return None
            out.append(np.ascontiguousarray(X.real))
            break
        M = X.reshape(a * d1 * d2, b)
        U, S, _ = np.linalg.svd(np.concatenate([M.real, M.imag], axis=1),
                                full_matrices=False)
        r = max(1, int(np.sum(S > tol * S[0]))) if S[0] > 0 else 1
        if r > b:
            return None
        Q = U[:, :r]
        G = Q.T @ M
        out.append(np.ascontiguousarray(Q.reshape(a, d1, d2, r)))
    return [xp.asarray(W) for W in out]


def _exact_phase_form(M, phase, tol=1e-12):
    """M as phase times a real matrix with exact zeros, when phase is 1,
    -1, i or -i and M / phase is real to tol of its largest entry; M
    otherwise. On a real basis the matrix of i times a real operator is
    imaginary up to the rounding of its column, and conductance.py reads
    that form off exact zeros (_spin_triplet) to take its real route."""
    if phase not in (1.0, -1.0, 1j, -1j):
        return M
    R = M * np.conj(phase)  # exact for these four
    if float(np.max(np.abs(R.imag), initial=0.)) > tol * max(
            float(np.max(np.abs(R), initial=0.)), 1e-300):
        return M
    out = np.zeros(M.shape, dtype=complex)
    if phase in (1.0, -1.0):
        out.real = phase.real * R.real
    else:
        out.imag = phase.imag * R.real
    return out


# The largest tangent space pad_bonds may ask for. Measured on d22 (v3,
# 21-site BCS Wilson chain, nkry = 100): pad_bonds=256 gives 2621439 and a
# peak of about 1.9 GB in 2114 s; "full" there would be 4^21.
MAX_PADDED_DIMENSION = 4_000_000


def compress_mpo(Ws, tol=1e-13):
    """The same operator as the MPO arrays Ws with every bond at its exact
    rank: one SVD sweep from the right and one from the left, dropping the
    singular values below tol times the largest of their bond, which are
    the channels that carry nothing (rounding apart, the operator is
    unchanged).

    AutoMPO writes a one-site operator with bond dimension 2 (the channel
    "not placed yet" next to "placed"), and it comes out of here with bond
    dimension 1, so applying it to a tangent vector costs the contractions
    of a product operator and not half of those of H."""
    xp = backend.xp()
    Ws = list(Ws)
    L = len(Ws)

    def kept(S):
        top = float(S[0]) if S.shape[0] else 0.0
        return max(1, int(xp.sum(S > tol * top))) if top > 0.0 else 1

    for l in range(L - 1, 0, -1):
        wl, t, u, wr = Ws[l].shape
        U, S, Vh = xp.linalg.svd(Ws[l].reshape(wl, t * u * wr), full_matrices=False)
        k = kept(S)
        Ws[l] = Vh[:k].reshape(k, t, u, wr)
        Ws[l - 1] = xp.tensordot(Ws[l - 1], U[:, :k] * S[:k], axes=(3, 0))
    for l in range(L - 1):
        wl, t, u, wr = Ws[l].shape
        U, S, Vh = xp.linalg.svd(Ws[l].reshape(wl * t * u, wr), full_matrices=False)
        k = kept(S)
        Ws[l] = U[:, :k].reshape(wl, t, u, k)
        Ws[l + 1] = xp.tensordot(S[:k, None] * Vh[:k], Ws[l + 1], axes=(1, 0))
    return real_if_exact(Ws)


def check_pad_bonds(pad_bonds):
    """pad_bonds as given if it is None, "full" or a positive int; a
    ValueError otherwise."""
    if pad_bonds is None or (isinstance(pad_bonds, str) and pad_bonds == "full"):
        return pad_bonds
    if (isinstance(pad_bonds, bool) or not isinstance(pad_bonds, (int, np.integer))
            or pad_bonds < 1):
        raise ValueError("pad_bonds must be None, \"full\" or a positive "
                         "integer, got %r" % (pad_bonds,))
    return int(pad_bonds)


def structural_bounds(dims):
    """The largest useful dimension of every interior bond of an MPS with
    the physical dimensions `dims`: bond l (between sites l and l+1,
    counted from 0) has min(d_0...d_l, d_{l+1}...d_{L-1}). Python ints, so
    a long chain does not overflow."""
    out = []
    for l in range(len(dims) - 1):
        left, right = 1, 1
        for d in dims[:l + 1]:
            left *= int(d)
        for d in dims[l + 1:]:
            right *= int(d)
        out.append(min(left, right))
    return out


def pad_mps_bonds(Ms, pad_bonds):
    """The site arrays `Ms` (chi_left, d, chi_right) with every interior
    bond padded by zero slices on both of its tensors to
    max(chi_l, min(pad_bonds, b_l)), b_l the structural bound of the bond,
    or to max(chi_l, b_l) for pad_bonds="full". The state is unchanged. A
    bond already above b_l is left as it is (mixed_gauge's first two
    sweeps bring it down to b_l), and no bond is ever shrunk. pad_bonds
    None returns Ms itself.

    The padded tangent dimension, sum_l (chi'_{l-1} d_l - chi'_l) chi'_l
    with every bond clipped to its bound, is computed from the shapes
    first, and a ValueError is raised above MAX_PADDED_DIMENSION, since the
    Krylov basis holds nkry vectors of that length ("full" on a 21-site
    chain of d = 4 asks for 4^21).

    The padding is xp.pad on the arrays of the active backend, a new array
    per site and no in-place write; only the shapes are read on the host.
    This is unrelated to backend.set_pad_bonds, which freezes the bond
    shapes of a whole DMRG run for JAX."""
    pad_bonds = check_pad_bonds(pad_bonds)
    if pad_bonds is None:
        return Ms
    xp = backend.xp()
    L = len(Ms)
    dims = [int(M.shape[1]) for M in Ms]
    bounds = structural_bounds(dims)
    new = [1] + [max(int(Ms[l].shape[2]), b if pad_bonds == "full"
                     else min(pad_bonds, b))
                 for l, b in enumerate(bounds)] + [1]
    eff = [1] + [min(c, b) for c, b in zip(new[1:-1], bounds)] + [1]
    dim = sum((eff[l] * dims[l] - eff[l + 1]) * eff[l + 1] for l in range(L))
    if dim > MAX_PADDED_DIMENSION:
        raise ValueError(
            "pad_bonds=%r gives a tangent space of dimension %d, above the "
            "cap of %d (the Krylov basis holds nkry vectors of that length). "
            "Use a smaller integer pad_bonds, or a shorter chain."
            % (pad_bonds, dim, MAX_PADDED_DIMENSION))
    out = []
    for l, M in enumerate(Ms):
        pl, pr = new[l] - M.shape[0], new[l + 1] - M.shape[2]
        if l == 0:
            pl = 0
        if l == L - 1:
            pr = 0
        if pl > 0 or pr > 0:
            M = xp.pad(M, ((0, max(pl, 0)), (0, 0), (0, max(pr, 0))))
        out.append(M)
    return out


# The feeds of expand_bonds: P^{2perp} O|GS> and P^{2perp} H P O|GS>, what
# the Lanczos from P O|GS> loses at its start and at its first step.
EXPANSION_FEEDS = ("operator", "krylov")
# Bytes apply_many may hold for one stack of tangent vectors (the stacked
# environments and the largest intermediate, counted as complex).
BATCH_BYTES = 1 << 28

# Singular values of the stacked feeds below this fraction of the largest
# one over all bonds are not used: a bond whose discarded weight is
# rounding gets no direction.
EXPANSION_FLOOR = 1e-8


def check_expand_bonds(expand_bonds):
    """expand_bonds as given if it is None or a positive int; a ValueError
    otherwise."""
    if expand_bonds is None:
        return None
    if (isinstance(expand_bonds, bool)
            or not isinstance(expand_bonds, (int, np.integer))
            or expand_bonds < 1):
        raise ValueError("expand_bonds must be None or a positive integer, "
                         "got %r" % (expand_bonds,))
    return int(expand_bonds)


def check_expansion_keywords(pad_bonds, expand_bonds, expand_rounds=1):
    """The keywords validated; a ValueError when pad_bonds and expand_bonds
    are both given, since padding and the expansion together are not
    tested, and for expand_rounds other than a positive int, or above 1
    without expand_bonds."""
    pad_bonds = check_pad_bonds(pad_bonds)
    expand_bonds = check_expand_bonds(expand_bonds)
    if (isinstance(expand_rounds, bool)
            or not isinstance(expand_rounds, (int, np.integer))
            or expand_rounds < 1):
        raise ValueError("expand_rounds must be a positive integer, got %r"
                         % (expand_rounds,))
    if expand_rounds > 1 and expand_bonds is None:
        raise ValueError("expand_rounds=%d needs expand_bonds" % expand_rounds)
    if pad_bonds is not None and expand_bonds is not None:
        raise ValueError("pad_bonds and expand_bonds cannot be used "
                         "together; pass one of them")
    return pad_bonds, expand_bonds


def expansion_bounds(shapes, k):
    """The most directions expand_bonds=k can add to every interior bond of
    a mixed-gauge MPS whose isometries have the (a, d, b) shapes `shapes`,
    in the left-to-right sweep of TangentSpace.expand: min(k, (a_l +
    k_{l-1}) d_l - b_l, d_{l+1} c_{l+1} - b_l), the dimensions of the two
    discarded spaces of the bond once bond l-1 has grown by k_{l-1}."""
    out, prev = [], 0
    for l in range(len(shapes) - 1):
        a, d, b = shapes[l]
        _, d2, c = shapes[l + 1]
        prev = max(0, min(int(k), (a + prev) * d - b, d2 * c - b))
        out.append(prev)
    return out


def expanded_dimension(shapes, ks):
    """The tangent dimension sum_l (a'_l d_l - b'_l) b'_l once bond l has
    grown by ks[l], from the shapes alone (Python ints)."""
    L = len(shapes)
    add = [0] + list(ks) + [0]
    tot = 0
    for l, (a, d, b) in enumerate(shapes):
        a2, b2 = int(a) + add[l], int(b) + add[l + 1]
        tot += (a2 * int(d) - b2) * b2
    return tot


def mixed_gauge(Ms, pad_bonds=None):
    """Left isometries A_l, right isometries B_l and site-centre tensors
    C_l of the normalized MPS given by the arrays `Ms`, with A_l and B_l of
    the same shape, so that

        |GS> = A_1 ... A_{l-1} C_l B_{l+1} ... B_L   for every l.

    A left QR sweep then a right LQ sweep bring every bond to at most d
    times its neighbours, after which a second left QR sweep from the
    right-canonical form gives A_l of exactly the B_l bond dimensions (a
    rank-deficient bond keeps its dimension; the A_l there is completed
    arbitrarily by the QR, which leaves both gauge relations intact).

    pad_bonds (None, an int or "full", see pad_mps_bonds and the module
    docstring) pads the bonds with zero slices first; the completion above
    then fills them, and the tangent space grows without |GS> changing."""
    xp = backend.xp()
    Ms = list(real_if_exact(Ms))
    if pad_bonds is not None:
        Ms = pad_mps_bonds(Ms, pad_bonds)
    L = len(Ms)
    # left QR sweep (only to shrink structurally oversized bonds)
    for l in range(L - 1):
        a, d, b = Ms[l].shape
        Q, R = xp.linalg.qr(Ms[l].reshape(a * d, b))
        Ms[l] = Q.reshape(a, d, Q.shape[1])
        Ms[l + 1] = xp.tensordot(R, Ms[l + 1], axes=(1, 0))
    # right LQ sweep -> right isometries B
    Bs = [None] * L
    X = Ms[L - 1]
    for l in range(L - 1, 0, -1):
        a, d, b = X.shape
        Q, R = xp.linalg.qr(X.reshape(a, d * b).T)  # X = R^T Q^T
        Bs[l] = Q.T.reshape(Q.shape[1], d, b)
        X = xp.tensordot(Ms[l - 1], R.T, axes=(2, 0))
    nrm = float(xp.linalg.norm(X))
    if nrm == 0.0:
        raise ValueError("TaSK: the ground state has zero norm")
    Bs[0] = X / nrm  # the first site carries the norm; a (1,d,D1) unit vector
    # left QR sweep from the right-canonical form -> A, C
    As, Cs = [None] * L, [None] * L
    X = Bs[0]
    for l in range(L):
        Cs[l] = X
        if l == L - 1:
            As[l] = X  # a (D,d,1) unit vector, left-normalized
            break
        a, d, b = X.shape
        Q, R = xp.linalg.qr(X.reshape(a * d, b))
        As[l] = Q.reshape(a, d, b)
        X = xp.tensordot(R, Bs[l + 1], axes=(1, 0))
    Bs[0] = Bs[0] / xp.linalg.norm(Bs[0])
    return As, Bs, Cs


# ---------------------------------------------------------------- contractions
#
# Written as tensordot chains, not np.einsum: einsum runs complex
# contractions like these through its own loops rather than BLAS, measured
# 0.04 s per pairwise step at D=64, which made one Lanczos step on a
# 12-site chain take 18 s.

def _ext_left(E, bra, W, ket):
    """Left environment E[a,w,b] (bra, MPO, ket bond) grown by one site."""
    xp = backend.xp()
    X = xp.tensordot(E, ket, axes=(2, 0))  # a w s z
    X = xp.tensordot(X, W, axes=([1, 2], [0, 2]))  # a z t y
    X = xp.tensordot(bra.conj(), X, axes=([0, 1], [0, 2]))  # x z y
    return X.transpose(0, 2, 1)


def _ext_right(E, bra, W, ket):
    """Right environment E[x,y,z] (bra, MPO, ket bond) grown by one site."""
    xp = backend.xp()
    X = xp.tensordot(ket, E, axes=(2, 2))  # b s x y
    X = xp.tensordot(X, W, axes=([1, 3], [2, 3]))  # b x w t
    X = xp.tensordot(bra.conj(), X, axes=([1, 2], [3, 1]))  # a b w
    return X.transpose(0, 2, 1)


def _open1(Le, W, K, Re):
    """The site-l bra tensor of <left env| W |ket K> |right env>."""
    xp = backend.xp()
    X = xp.tensordot(Le, K, axes=(2, 0))  # a w s z
    X = xp.tensordot(X, W, axes=([1, 2], [0, 2]))  # a z t y
    return xp.tensordot(X, Re, axes=([1, 3], [2, 1]))  # a t x


def _open2(Le, W1, W2, K2, Re):
    """Two-site version of _open1 on the bond (l, l+1), K2[b,s,v,c]."""
    xp = backend.xp()
    X = xp.tensordot(Le, K2, axes=(2, 0))  # a w s v c
    X = xp.tensordot(X, W1, axes=([1, 2], [0, 2]))  # a v c t y
    X = xp.tensordot(X, W2, axes=([1, 4], [2, 0]))  # a c t u q
    return xp.tensordot(X, Re, axes=([1, 4], [2, 1]))  # a t u x


# The same contractions for a stack of tangent vectors, with one leading
# batch axis on the tensor that carries the vector (the ket K, or an
# environment that already holds it) and on the result. One tensordot then
# does the work of one per vector: at bond dimensions of a few tens the
# time of the unbatched calls was the call overhead, not the algebra
# (700 000 tensordot calls for the 1212 applications of a 24-site Kondo
# spectrum at D=40).

def _ext_left_k(E, bra, W, K):
    """_ext_left with a batched ket K[n,b,s,z]; returns [n,x,y,z]."""
    xp = backend.xp()
    X = xp.tensordot(E, K, axes=(2, 1))  # a w n s z
    X = xp.tensordot(X, W, axes=([1, 3], [0, 2]))  # a n z t y
    X = xp.tensordot(bra.conj(), X, axes=([0, 1], [0, 3]))  # x n z y
    return X.transpose(1, 0, 3, 2)


def _ext_left_e(E, bra, W, ket):
    """_ext_left with a batched environment E[n,a,w,b]."""
    xp = backend.xp()
    X = xp.tensordot(E, ket, axes=(3, 0))  # n a w s z
    X = xp.tensordot(X, W, axes=([2, 3], [0, 2]))  # n a z t y
    X = xp.tensordot(bra.conj(), X, axes=([0, 1], [1, 3]))  # x n z y
    return X.transpose(1, 0, 3, 2)


def _ext_right_k(E, bra, W, K):
    """_ext_right with a batched ket K[n,b,s,z]; returns [n,a,w,b]."""
    xp = backend.xp()
    X = xp.tensordot(K, E, axes=(3, 2))  # n b s x y
    X = xp.tensordot(X, W, axes=([2, 4], [2, 3]))  # n b x w t
    X = xp.tensordot(bra.conj(), X, axes=([1, 2], [4, 2]))  # a n b w
    return X.transpose(1, 0, 3, 2)


def _ext_right_e(E, bra, W, ket):
    """_ext_right with a batched environment E[n,x,y,z]."""
    xp = backend.xp()
    X = xp.tensordot(ket, E, axes=(2, 3))  # b s n x y
    X = xp.tensordot(X, W, axes=([1, 4], [2, 3]))  # b n x w t
    X = xp.tensordot(bra.conj(), X, axes=([1, 2], [4, 2]))  # a b n w
    return X.transpose(2, 0, 3, 1)


def _open1_k(Le, W, K, Re):
    """_open1 with a batched ket K[n,b,s,z]; returns [n,a,t,x]."""
    xp = backend.xp()
    X = xp.tensordot(Le, K, axes=(2, 1))  # a w n s z
    X = xp.tensordot(X, W, axes=([1, 3], [0, 2]))  # a n z t y
    X = xp.tensordot(X, Re, axes=([2, 4], [2, 1]))  # a n t x
    return X.transpose(1, 0, 2, 3)


def _open1_l(Le, W, K, Re):
    """_open1 with a batched left environment Le[n,a,w,b]."""
    xp = backend.xp()
    X = xp.tensordot(Le, K, axes=(3, 0))  # n a w s z
    X = xp.tensordot(X, W, axes=([2, 3], [0, 2]))  # n a z t y
    return xp.tensordot(X, Re, axes=([2, 4], [2, 1]))  # n a t x


def _open1_r(Le, W, K, Re):
    """_open1 with a batched right environment Re[n,x,y,z]."""
    xp = backend.xp()
    X = xp.tensordot(Le, K, axes=(2, 0))  # a w s z
    X = xp.tensordot(X, W, axes=([1, 2], [0, 2]))  # a z t y
    X = xp.tensordot(X, Re, axes=([1, 3], [3, 2]))  # a t n x
    return X.transpose(2, 0, 1, 3)


def _proj_left_many(X, A):
    """_proj_left on a stack X[n,a,d,b]."""
    xp = backend.xp()
    c = xp.tensordot(X, A.conj(), axes=([1, 2], [0, 1]))  # n b' b
    return X - xp.tensordot(c, A, axes=(2, 2)).transpose(0, 2, 3, 1)


def _proj_left(X, A):
    """Remove from X its component in the kept space of the left isometry
    A, on X's first two legs."""
    xp = backend.xp()
    c = xp.tensordot(A.conj(), X, axes=([0, 1], [0, 1]))
    return X - xp.tensordot(A, c, axes=(2, 0))


def _proj_right(X, B):
    """Remove from X its component in the kept space of the right isometry
    B, on X's last two legs."""
    xp = backend.xp()
    c = xp.tensordot(X, B.conj(), axes=([-2, -1], [1, 2]))
    return X - xp.tensordot(c, B, axes=(-1, 0))


def _pair(T1, T2):
    return backend.xp().tensordot(T1, T2, axes=(2, 0))


def _left_envs(As, Ws):
    L = len(As)
    out = [backend.xp().ones((1, 1, 1))]
    for l in range(L):
        out.append(_ext_left(out[-1], As[l], Ws[l], As[l]))
    return out  # out[l] holds sites 0..l-1


def _right_envs(Bs, Ws):
    L = len(Bs)
    out = [None] * (L + 1)
    out[L] = backend.xp().ones((1, 1, 1))
    for l in range(L - 1, -1, -1):
        out[l] = _ext_right(out[l + 1], Bs[l], Ws[l], Bs[l])
    return out  # out[l] holds sites l..L-1


# ---------------------------------------------------------------- tangent space

class TangentSpace:
    """The one-site tangent space of a ground-state MPS, orthogonal to it,
    with the Hamiltonian H projected onto it. pad_bonds enlarges it by
    padding the bonds of the MPS first (mixed_gauge, pad_mps_bonds);
    expand_bonds by at most that many directions per bond taken from the
    two-site discarded space of the operators expand_ops (MPO arrays), see
    expand(). The two cannot be combined."""

    def __init__(self, Ms, Ws, pad_bonds=None, expand_bonds=None,
                 expand_ops=(), expand_rounds=1):
        pad_bonds, expand_bonds = check_expansion_keywords(
            pad_bonds, expand_bonds, expand_rounds)
        self.As, self.Bs, self.Cs = mixed_gauge(Ms, pad_bonds=pad_bonds)
        self.Ws = Ws = real_if_exact(Ws)
        self.L = len(Ms)
        self.LH = _left_envs(self.As, Ws)
        self.RH = _right_envs(self.Bs, Ws)
        self.e0 = float(backend.xp().vdot(
            self.Cs[0], _open1(self.LH[0], Ws[0], self.Cs[0], self.RH[1])).real)
        self.added = [0] * (self.L - 1)
        if expand_bonds is not None:
            # every round sweeps the chain again from the expanded gauge,
            # so b and H b see the directions the previous round added; it
            # stops early once a round adds nothing
            for _ in range(int(expand_rounds)):
                add = self.expand(expand_bonds, list(expand_ops))
                self.added = [x + y for x, y in zip(self.added, add)]
                if sum(add) == 0:
                    break

    # vector-space operations on lists of T_l
    def dot(self, X, Y):
        xp = backend.xp()
        return sum(xp.vdot(x, y) for x, y in zip(X, Y))

    def norm(self, X):
        return float(np.sqrt(abs(self.dot(X, X))))

    def flat_size(self):
        return int(sum(a.size for a in self.As))

    def flat(self, X):
        """A tangent vector as one array; overlaps are then vdot."""
        return backend.xp().concatenate([x.ravel() for x in X])

    def unflat(self, v):
        out, k = [], 0
        for a in self.As:
            out.append(v[k:k + a.size].reshape(a.shape))
            k += a.size
        return out

    def gauge(self, X):
        return [_proj_left(x, a) for x, a in zip(X, self.As)]

    def dimension(self):
        return int(sum((a.shape[0] * a.shape[1] - a.shape[2]) * a.shape[2]
                       for a in self.As))

    def expectation(self, op):
        """<GS|O|GS> for an operator from operator()."""
        Wo, LO, RO = op
        y = _open1(LO[0], Wo[0], self.Cs[0], RO[1])
        return complex(backend.xp().vdot(self.Cs[0], y))

    def project_operator(self, Wo):
        """P O|GS> for the MPO arrays Wo, as a tangent vector. Exact for a
        sum of one-site operators, the projection otherwise."""
        Wo = real_if_exact(Wo)
        LO = _left_envs(self.As, Wo)
        RO = _right_envs(self.Bs, Wo)
        return self.gauge([_open1(LO[l], Wo[l], self.Cs[l], RO[l + 1])
                           for l in range(self.L)])

    def operator(self, Wo):
        """An MPO with its A-A and B-B environments, ready for apply()."""
        Wo = real_if_exact(Wo)
        return (Wo, _left_envs(self.As, Wo), _right_envs(self.Bs, Wo))

    def _mixed_envs(self, T, op=None):
        """The environments of the operator `op` (default H) whose ket
        carries the tangent tensor of T to their left (resp. right), A_l
        before it and B_l after it."""
        Ws, LE, RE = op if op is not None else (self.Ws, self.LH, self.RH)
        L, As, Bs = self.L, self.As, self.Bs
        zero = backend.xp().zeros((1, 1, 1))
        Lc = [zero]
        for l in range(L):
            Lc.append(_ext_left(Lc[-1], As[l], Ws[l], Bs[l])
                      + _ext_left(LE[l], As[l], Ws[l], T[l])
                      if l > 0 else
                      _ext_left(LE[l], As[l], Ws[l], T[l]))
        Rc = [None] * (L + 1)
        Rc[L] = zero
        for l in range(L - 1, -1, -1):
            Rc[l] = (_ext_right(RE[l + 1], Bs[l], Ws[l], T[l])
                     + (_ext_right(Rc[l + 1], Bs[l], Ws[l], As[l])
                        if l < L - 1 else 0.0))
        return Lc, Rc

    def apply(self, op, T):
        """P O P acting on the tangent vector T, for an operator from
        operator(); the part of O|T> along |GS> is dropped."""
        Ws, LE, RE = op
        L, As, Bs = self.L, self.As, self.Bs
        Lc, Rc = self._mixed_envs(T, op)
        out = []
        for l in range(L):
            y = _open1(LE[l], Ws[l], T[l], RE[l + 1])
            if l > 0:
                y = y + _open1(Lc[l], Ws[l], Bs[l], RE[l + 1])
            if l < L - 1:
                y = y + _open1(LE[l], Ws[l], As[l], Rc[l + 1])
            out.append(y)
        return self.gauge(out)

    def flat_many(self, X):
        """A stack of tangent vectors as an (n, N) array, one row each."""
        n = X[0].shape[0]
        return backend.xp().concatenate([x.reshape(n, -1) for x in X], axis=1)

    def unflat_many(self, V):
        """The rows of an (n, N) array as a stack of tangent vectors, each
        T_l with a leading batch axis."""
        out, k = [], 0
        for a in self.As:
            out.append(V[:, k:k + a.size].reshape((V.shape[0],) + a.shape))
            k += a.size
        return out

    def batch_size(self, op=None):
        """How many tangent vectors apply_many takes at once for the
        operator `op` (default H), so that the stacked right environments
        and the largest intermediate stay within BATCH_BYTES."""
        Ws = (op if op is not None else (self.Ws,))[0]
        envs = sum(a.shape[2] ** 2 * W.shape[3] for a, W in zip(self.As, Ws))
        work = max(a.size * max(W.shape[0], W.shape[3]) for a, W in zip(self.As, Ws))
        return int(max(1, BATCH_BYTES // (16 * (envs + 3 * work + 2 * self.flat_size()))))

    def apply_many(self, op, T):
        """apply() on a stack of tangent vectors, every T_l with a leading
        batch axis: P O P on each, as one contraction per site. The right
        environments are kept for the whole stack and the left ones are
        carried along the sweep, so the caller bounds the stack with
        batch_size()."""
        Ws, LE, RE = op
        L, As, Bs = self.L, self.As, self.Bs
        Rc = [None] * (L + 1)
        for l in range(L - 1, 0, -1):
            Rc[l] = _ext_right_k(RE[l + 1], Bs[l], Ws[l], T[l])
            if l < L - 1:
                Rc[l] = Rc[l] + _ext_right_e(Rc[l + 1], Bs[l], Ws[l], As[l])
        out, Lc = [], None
        for l in range(L):
            y = _open1_k(LE[l], Ws[l], T[l], RE[l + 1])
            if l > 0:
                y = y + _open1_l(Lc, Ws[l], Bs[l], RE[l + 1])
            if l < L - 1:
                y = y + _open1_r(LE[l], Ws[l], As[l], Rc[l + 1])
                new = _ext_left_k(LE[l], As[l], Ws[l], T[l])
                if l > 0:
                    new = new + _ext_left_e(Lc, As[l], Ws[l], Bs[l])
                Lc = new
            Rc[l + 1] = None
            out.append(_proj_left_many(y, As[l]))
        return out

    def apply_h_many(self, T):
        """P (H - E0) P on a stack of tangent vectors."""
        y = self.apply_many((self.Ws, self.LH, self.RH), T)
        return [a - self.e0 * t for a, t in zip(y, T)]

    def apply_h(self, T):
        """P (H - E0) P acting on the tangent vector T."""
        y = self.apply((self.Ws, self.LH, self.RH), T)
        return [a - self.e0 * t for a, t in zip(y, T)]

    def _two_site_discarded(self, T, op=None):
        """Bond by bond (l, l+1), the two-site tensor (a, d, d', c) of
        P^{2perp} O|T>, O the operator `op` from operator() (default H):
        the part of O|T> in the two-site discarded space of the bond, outside
        the kept spaces of A_l on the left and of B_{l+1} on the right."""
        L, As, Bs = self.L, self.As, self.Bs
        Ws, LE, RE = op if op is not None else (self.Ws, self.LH, self.RH)
        Lc, Rc = self._mixed_envs(T, op)
        for l in range(L - 1):
            Le, Re = LE[l], RE[l + 2]
            K = _pair(T[l], Bs[l + 1]) + _pair(As[l], T[l + 1])
            y = _open2(Le, Ws[l], Ws[l + 1], K, Re)
            if l > 0:
                KB = _pair(Bs[l], Bs[l + 1])
                y = y + _open2(Lc[l], Ws[l], Ws[l + 1], KB, Re)
            if l < L - 2:
                KA = _pair(As[l], As[l + 1])
                y = y + _open2(Le, Ws[l], Ws[l + 1], KA, Rc[l + 2])
            yield _proj_right(_proj_left(y, As[l]), Bs[l + 1])

    def _state_discarded(self, op):
        """Bond by bond, the two-site tensor of P^{2perp} O|GS> for an
        operator from operator(), from the two-site centre C_l B_{l+1}."""
        Wo, LO, RO = op
        for l in range(self.L - 1):
            K = _pair(self.Cs[l], self.Bs[l + 1])
            y = _open2(LO[l], Wo[l], Wo[l + 1], K, RO[l + 2])
            yield _proj_right(_proj_left(y, self.As[l]), self.Bs[l + 1])

    def two_site_variance(self, T):
        """||P^{2perp} H |T>||^2, the part of H|T> in the two-site discarded
        space, which no one-site tangent vector can represent."""
        tot = 0.0
        for y in self._two_site_discarded(T):
            tot += float(backend.xp().vdot(y, y).real)
        return tot

    def expand(self, expand_bonds, Wops, feeds=EXPANSION_FEEDS):
        """Enlarge the tangent space with directions taken from the two-site
        discarded space (the module docstring, "Expanded bonds"). Wops are
        the MPO arrays of the operators whose O|GS> feeds the expansion.
        feeds is a tuple of "operator" (P^{2perp} O|GS>), "krylov"
        (P^{2perp} H b, b = P O|GS> normalized) and "hamiltonian"
        (P^{2perp} H|GS>, the direction of CBE-DMRG). Returns the number of
        directions added to every bond. The state, E0 and both gauge
        relations do not change; LH and RH are rebuilt."""
        expand_bonds = check_expand_bonds(expand_bonds)
        if expand_bonds is None:
            return [0] * (self.L - 1)
        for f in feeds:
            if f not in EXPANSION_FEEDS + ("hamiltonian",):
                raise ValueError("expand: unknown feed %r" % (f,))
        xp = backend.xp()
        As, Bs, Cs, L = self.As, self.Bs, self.Cs, self.L
        kmax = expansion_bounds([a.shape for a in As], expand_bonds)
        dim = expanded_dimension([a.shape for a in As], kmax)
        if dim > MAX_PADDED_DIMENSION:
            raise ValueError(
                "expand_bonds=%r gives a tangent space of dimension up to %d, "
                "above the cap of %d (the Krylov basis holds nkry vectors of "
                "that length). Use a smaller expand_bonds."
                % (expand_bonds, dim, MAX_PADDED_DIMENSION))
        # A feed is O applied to a ket: |GS> ("operator", "hamiltonian") or
        # the tangent vector b ("krylov"). First pass, in the unexpanded
        # gauge: the weight of every feed in the two-site discarded space,
        # bond by bond, which normalizes the feeds and sets the floor.
        Ao, Bo, Co = self.As, self.Bs, self.Cs
        specs = []
        for Wo in Wops:
            Wv, LO, RO = self.operator(Wo)
            b = self.gauge([_open1(LO[l], Wv[l], Co[l], RO[l + 1])
                            for l in range(L)])
            nb = self.norm(b)
            if "operator" in feeds:
                specs.append(dict(
                    W=Wv, RE=RO, T=None,
                    ys=self._state_discarded((Wv, LO, RO)),
                    ref2=nb**2 + abs(self.expectation((Wv, LO, RO)))**2))
            if "krylov" in feeds and nb > 0.0:
                b = [x / nb for x in b]
                specs.append(dict(
                    W=self.Ws, RE=self.RH, T=b, Rc=self._mixed_envs(b)[1],
                    ys=self._two_site_discarded(b),
                    ref2=self.norm(self.apply_h(b))**2))
        if "hamiltonian" in feeds:
            specs.append(dict(
                W=self.Ws, RE=self.RH, T=None,
                ys=self._state_discarded((self.Ws, self.LH, self.RH)),
                ref2=abs(self.e0)**2))
        # a feed whose two-site norm is below EXPANSION_FLOOR of the norm of
        # its vector is dropped (on v3 the rounding of a saturated bond is
        # 1e-11 of it, and kept it would add arbitrary directions), the
        # others are scaled to unit weight
        kept, per_bond = [], np.zeros(L - 1)
        for sp in specs:
            n2l = np.array([float(np.asarray(backend.to_host(
                xp.vdot(y, y))).real) for y in sp.pop("ys")])
            n2 = float(np.sum(n2l))
            if n2 > 0.0 and n2 > EXPANSION_FLOOR**2 * (sp["ref2"] + n2):
                sp["scale"] = 1.0 / np.sqrt(n2)
                kept.append(sp)
                per_bond = per_bond + n2l / n2
        ks = [0] * (L - 1)
        if not kept:
            return ks
        floor = EXPANSION_FLOOR * float(np.sqrt(np.max(per_bond)))
        # Second pass, a left-to-right sweep: the two-site tensor of bond l
        # is taken with the bras already expanded on bonds 0..l-1 (kets
        # unexpanded), so a direction added on one bond is carried on to
        # the next one.
        from .cbe import _orthonormal_complement_columns
        A, B, C = list(Ao), list(Bo), list(Co)
        real = not np.iscomplexobj(np.zeros(0, dtype=Ao[0].dtype))
        for sp in kept:
            sp["LE"], sp["Lc"] = xp.ones((1, 1, 1)), None
        for l in range(L - 1):
            a2, d, b = A[l].shape
            _, d2, c = Bo[l + 1].shape
            bound = max(0, min(expand_bonds, a2 * d - b, d2 * c - b))
            k = 0
            if bound > 0:
                mats = []
                for sp in kept:
                    W, RE, T, Le = sp["W"], sp["RE"], sp["T"], sp["LE"]
                    if T is None:
                        y = _open2(Le, W[l], W[l + 1], _pair(Co[l], Bo[l + 1]),
                                   RE[l + 2])
                    else:
                        K = _pair(T[l], Bo[l + 1]) + _pair(Ao[l], T[l + 1])
                        y = _open2(Le, W[l], W[l + 1], K, RE[l + 2])
                        if l > 0:
                            y = y + _open2(sp["Lc"], W[l], W[l + 1],
                                           _pair(Bo[l], Bo[l + 1]), RE[l + 2])
                        if l < L - 2:
                            y = y + _open2(Le, W[l], W[l + 1],
                                           _pair(Ao[l], Ao[l + 1]),
                                           sp["Rc"][l + 2])
                    y = _proj_right(_proj_left(y, A[l]), Bo[l + 1])
                    mats.append(sp["scale"] * y.reshape(a2 * d, d2 * c))
                if real and any(np.iscomplexobj(np.zeros(0, dtype=m.dtype))
                                for m in mats):
                    mats = [p for m in mats for p in (m.real, m.imag)]
                U, sl, _ = xp.linalg.svd(xp.concatenate(mats, axis=1),
                                         full_matrices=False)
                _, sr, Vh = xp.linalg.svd(xp.concatenate(mats, axis=0),
                                          full_matrices=False)
                sl = np.asarray(backend.to_host(sl)).real
                sr = np.asarray(backend.to_host(sr)).real
                k = min(bound, int(np.count_nonzero(sl > floor)),
                        int(np.count_nonzero(sr > floor)))
            if k > 0:
                dt = Ao[l].dtype
                Ul = _orthonormal_complement_columns(
                    U[:, :k].astype(dt), A[l].reshape(a2 * d, b))
                Vr = _orthonormal_complement_columns(
                    Vh[:k].astype(dt).conj().T,
                    Bo[l + 1].reshape(b, d2 * c).conj().T)
                k = min(int(Ul.shape[1]), int(Vr.shape[1]))
            if k > 0:
                A[l] = xp.concatenate([A[l], Ul[:, :k].reshape(a2, d, k)],
                                      axis=2)
                B[l] = xp.concatenate([B[l], xp.zeros((a2, d, k), dtype=dt)],
                                      axis=2)
                C[l] = xp.pad(C[l], ((0, 0), (0, 0), (0, k)))
                A[l + 1] = xp.concatenate(
                    [A[l + 1], xp.zeros((k, d2, c), dtype=dt)], axis=0)
                B[l + 1] = xp.concatenate(
                    [B[l + 1], Vr[:, :k].conj().T.reshape(k, d2, c)], axis=0)
                C[l + 1] = xp.pad(C[l + 1], ((0, k), (0, 0), (0, 0)))
                ks[l] = k
            # environments of sites 0..l, bra expanded, ket unexpanded
            for sp in kept:
                W, T, Le = sp["W"], sp["T"], sp["LE"]
                if T is not None:
                    Lc = _ext_left(Le, A[l], W[l], T[l])
                    if l > 0:
                        Lc = Lc + _ext_left(sp["Lc"], A[l], W[l], Bo[l])
                    sp["Lc"] = Lc
                sp["LE"] = _ext_left(Le, A[l], W[l], Ao[l])
        self.As, self.Bs, self.Cs = A, B, C
        self.LH = _left_envs(self.As, self.Ws)
        self.RH = _right_envs(self.Bs, self.Ws)
        self.e0 = float(xp.vdot(
            self.Cs[0], _open1(self.LH[0], self.Ws[0], self.Cs[0],
                               self.RH[1])).real)
        return ks


# ---------------------------------------------------------------- Lanczos
#
# A Krylov basis is a preallocated store, one flattened tangent vector per
# row, so that every overlap with it is a matrix product over its first m
# rows: backend.row_overlaps and backend.project_rows, a view on NumPy and a
# loop over row chunks on a device, where a slice or a whole-store product
# would copy the basis.


def lanczos(ts, b, nkry, breakdown=1e-10):
    """Lanczos on P(H-E0)P from the tangent vector b, with full
    reorthogonalization done twice per step (the paper's step vii).

    Returns (alpha, beta, store, m): the tridiagonal elements and the
    orthonormal Krylov vectors, the first m rows of the flattened store.
    Stops early when the new vector's norm falls below `breakdown` times
    the norm of H|v_i> it came from, the Krylov space then being
    invariant."""
    xp = backend.xp()
    nb = ts.norm(b)
    if nb == 0.0:
        return np.zeros(0), np.zeros(0), None, 0
    # The basis is kept flattened, one row per vector, so that the
    # reorthogonalization is two matrix products rather than a Python sum
    # over sites per pair of vectors (a third of the time at L=32, D=64).
    v0 = ts.flat(b) / nb
    store = xp.zeros((nkry + 1, v0.size), dtype=v0.dtype)
    store = backend.setrow(store, 0, v0)
    alpha, beta = [], []
    m = 1
    for i in range(nkry + 1):
        w = ts.flat(ts.apply_h(ts.unflat(store[i])))
        if w.dtype != store.dtype:  # a complex P H P on a real start
            store = store.astype(w.dtype)
        scale = float(xp.linalg.norm(w))
        for rep in range(2):
            c, w = backend.project_rows(store, w, m)
            if rep == 0:
                alpha.append(float(c[i].real))
            w = ts.flat(ts.gauge(ts.unflat(w)))
        if i == nkry:
            break
        bnext = float(xp.linalg.norm(w))
        if bnext <= breakdown * max(scale, abs(alpha[-1]), 1e-300):
            break
        beta.append(bnext)
        store = backend.setrow(store, m, w / bnext)
        m += 1
    return np.array(alpha), np.array(beta), store, m


def poles(ts, b, nkry, a=None):
    """Discrete poles (w_alpha) and weights of <GS|A delta(w-H+E0) B|GS>
    projected to the tangent space, b = P B|GS> and a = P A^dagger|GS>
    (a = None means A = B^dagger). Also returns the Ritz vectors' Krylov
    coefficients and the basis (store, m), for the variance."""
    al, be, store, m = lanczos(ts, b, nkry)
    if m == 0:
        return np.zeros(0), np.zeros(0, dtype=complex), None, (store, m)
    Tm = np.diag(al[:m]) + np.diag(be[:m - 1], 1) + np.diag(be[:m - 1], -1)
    ws, U = np.linalg.eigh(Tm)
    bproj = U[0, :] * ts.norm(b)  # <alpha|b>
    if a is None:
        weights = np.abs(bproj) ** 2 + 0j
    else:
        ak = backend.to_host(backend.row_overlaps(store, ts.flat(a), m))[:m]  # <v_i|a>
        aproj = U.conj().T @ ak  # <alpha|a>
        weights = aproj.conj() * bproj
    return ws, weights, U, (store, m)


def real_starts(ts, starts, tol=1e-13):
    """The starting vectors with a global phase removed from each, as real
    arrays, when the state and H are real and every start is a real vector
    times a phase (S_y|GS> = i times a real vector, on a real ground
    state); the starts unchanged otherwise.

    The complex span of the starts, and so the Krylov space and every Ritz
    value, is the same. What changes is that the basis is then real: half
    the memory of the store, and real arithmetic in every application of H
    and every overlap, where one imaginary start out of three made all of
    them complex."""
    xp = backend.xp()
    if any(xp.iscomplexobj(x) for x in list(ts.As) + list(ts.Ws)):
        return starts
    out = []
    for w in starts:
        if not any(xp.iscomplexobj(x) for x in w):
            out.append(w)
            continue
        f = ts.flat(w)
        k = int(xp.argmax(xp.abs(f)))
        top = backend.scalar(f[k])
        if abs(top) == 0.0:
            out.append([x.real for x in w])
            continue
        phase = top / abs(top)
        f = f * phase.conjugate()
        if float(xp.max(xp.abs(f.imag))) > tol * abs(top):
            return starts
        out.append(ts.unflat(f.real))
    return out


def _overlap_columns(store, ys, n):
    """<q_i|y_k> for the first n rows q_i of the store and every row y_k of
    ys, as a host (n, len(ys)) array: one matrix product on the host, the
    chunked row_overlaps vector by vector on a device, where a product
    with the whole store as an operand is what ran out of memory."""
    if backend.backend_name() == "numpy":
        # written so that the store is never copied: its conjugate, or its
        # complex copy when it is real and ys is not, is as large as it is
        Qt = store[:n].T
        if np.iscomplexobj(store):
            return (ys.conj() @ Qt).conj().T
        if np.iscomplexobj(ys):
            return (ys.real @ Qt + 1j * (ys.imag @ Qt)).T
        return (ys @ Qt).T
    return np.stack([backend.to_host(backend.row_overlaps(store, y, n))[:n]
                     for y in ys], axis=1)


def block_krylov(ts, starts, nkry, breakdown=1e-10, residual=False):
    """Orthonormal basis of the block Krylov space of P(H-E0)P from the
    tangent vectors `starts`, nkry blocks beyond the starting one, and the
    projected Hamiltonian in it, Hq[i,j] = <q_i|P(H-E0)P|q_j>.

    The basis is kept flattened, one row per vector (ts.flat), so that
    every overlap is a matrix product. Every vector is orthogonalized twice
    against the whole basis, and a candidate is dropped when what is left
    of it is below `breakdown` times its norm before the
    orthogonalization, so the basis stops growing once the space reached
    from `starts` is invariant. The threshold is 1e-10 and not tighter:
    after two passes the remainder of a vector already in the space is
    rounding noise of relative size about 1e-13, and accepting it
    normalizes that noise, gauge violation included, up to a unit vector,
    which made P S P measurably non-Hermitian (0.4) on the basis.

    H Q is never stored: when P(H-E0)P|q_j> is computed, its overlaps with
    every vector already in the basis, which includes all q_i with i <= j,
    give Hq[:, j] above the diagonal, and Hermiticity gives the rest.
    Returns (store, count, Hq), the basis being the first `count` rows of
    the flattened store (the rows after `count` are zero) and Hq a host
    (count, count) array.

    With residual=True a fourth entry (j0, R) describes what the basis
    leaves out: R[j,k] = <r_j|r_k>, with r_j the part of P(H-E0)P|q_j>
    outside the basis for the vectors j0 <= j < count of the last block,
    the only ones whose image leaves it. A Ritz vector with coefficients u
    then has the residual norm ||(P(H-E0)P - w)|u>|| =
    sqrt(u[j0:]^dagger R u[j0:]), zero for an eigenstate; R is empty when
    the Krylov space was exhausted before nkry blocks."""
    xp = backend.xp()
    N = ts.flat_size()
    starts = real_starts(ts, starts)
    # preallocated and filled row by row: a vstack per accepted vector
    # copies the whole basis each time
    dtype = np.result_type(*[x.dtype for w in starts for x in w],
                           ts.As[0].dtype, *[W.dtype for W in ts.Ws])
    nmax = min(N, len(starts))
    store = xp.zeros((nmax, N), dtype=dtype)
    count = 0
    Hu = None  # Hu[i, j] = <q_i|P(H-E0)P|q_j> for i <= j, set below

    def add(block, cols=None):
        """Orthogonalize the candidates of `block` against the basis and
        keep the ones that are left. The starts are tangent vectors and
        are gauged first. The Krylov candidates (cols, their columns of
        Hu) are flat images P(H-E0)P|q_j>, already gauged by apply_many,
        and the first projection pass gives their overlaps with the n0
        vectors of the basis before this block, which are Hu[:n0, j]."""
        nonlocal count, store
        n0 = count
        for k, w in enumerate(block):
            if count == nmax:
                if cols is not None:  # the overlaps the rest still owes Hu
                    Hu[:n0, cols[k:]] = _overlap_columns(
                        store, xp.stack(block[k:]), n0)
                break
            if cols is None:
                w = ts.flat(ts.gauge(w))
            scale = float(xp.linalg.norm(w))
            for rep in range(2):
                c, w = backend.project_rows(store, w, count)
                if rep == 0 and cols is not None:
                    Hu[:n0, cols[k]] = backend.to_host(c)[:n0]
                w = ts.flat(ts.gauge(ts.unflat(w)))
            nw = float(xp.linalg.norm(w))
            if scale > 0.0 and nw > breakdown * scale:
                store = backend.setrow(store, count, w / nw)
                count += 1

    # The starts that survive their own orthogonalization size the store:
    # no later block is larger than the one before it, so count*(nkry+1)
    # rows hold the whole basis. For S, S^dagger pairs of Hermitian
    # operators (TaskKondoSpectrum's Sx, Sy, Sz) half the starts are
    # duplicates, and sizing by len(starts) doubled a store of tens of GB.
    add(starts)
    first = store
    nmax = min(N, count * (nkry + 1))
    store = xp.zeros((nmax, N), dtype=dtype)
    for k in range(count):
        store = backend.setrow(store, k, first[k])
    del first
    done = 0
    Hu = np.zeros((nmax, nmax), dtype=complex)  # Hu[i, j] for i <= j
    left = None
    for step in range(nkry + 1):
        n = count
        if n == done:
            break
        last = step == nkry
        Hb = []
        nb = ts.batch_size()
        for j0 in range(done, n, nb):  # the whole block at once
            j1 = min(n, j0 + nb)
            hs = ts.flat_many(ts.apply_h_many(ts.unflat_many(store[j0:j1])))
            if last and not residual:  # no projection pass to read them from
                Hu[:n, j0:j1] = _overlap_columns(store, hs, n)
            Hb.extend(hs[j - j0] for j in range(j0, j1))
            del hs
        if last:
            if residual:
                # two passes, as in add(): the remainder of a vector that
                # is already in the space is then rounding noise; the
                # first one gives the overlaps of Hu, as in add()
                out = []
                for j, h in zip(range(done, n), Hb):
                    for rep in range(2):
                        c, h = backend.project_rows(store, h, n)
                        if rep == 0:
                            Hu[:n, j] = backend.to_host(c)[:n]
                    out.append(backend.to_host(h))
                out = np.array(out)
                left = (done, out.conj() @ out.T)
            done = n
            break
        add(Hb, cols=np.arange(done, n))
        done = n
        del Hb
    upper = np.triu(Hu[:count, :count])
    Hq = upper + np.triu(upper, 1).conj().T
    Hq = 0.5 * (Hq + Hq.conj().T)  # real diagonal
    if not np.iscomplexobj(np.zeros(0, dtype=dtype)):
        Hq = Hq.real
    if residual:
        if left is None:  # nothing leaves an exhausted space
            left = (count, np.zeros((0, 0)))
        return store, count, Hq, left
    return store, count, Hq


def projected_spectrum(ts, ops, nkry, residual=False, column=False):
    """Energies and operator matrices in the basis {|GS>} + the Ritz states
    of P(H-E0)P on the block Krylov space started from P O|GS> and
    P O^dagger|GS> for every operator in `ops`.

    ops: list of (W, Wdag) MPO-array pairs, the operator and its adjoint.
    Returns (e, mats): e[0] = 0 for the ground state and e[1:] the Ritz
    energies relative to E0, and for each operator the matrix
    mats[k][a,b] = <a|O_k|b> in that basis, in which the O_k|GS> column is
    exact and the tangent-tangent block is the projection P O_k P.

    With column=True only row and column 0 are built, <a|O_k|GS> and
    <GS|O_k|a>, and the tangent-tangent block is left at zero: what a sum
    over transitions out of the ground state reads (the second-order dI/dV),
    without the n applications of every O_k and the n x n overlaps the
    block costs. Everything else is the same code, so the column is the
    one of the full build, bit for bit.

    With residual=True a third entry gives, for every state of the basis,
    the residual norm ||(P(H-E0)P - e)|a>|| of block_krylov: zero for the
    ground state and for a Ritz state that is an eigenstate of P H P, and
    of the order of the distance to its neighbours for one that only
    stands for a group of eigenstates the Krylov space has not resolved."""
    b = [ts.project_operator(W) for W, _ in ops]  # P O|GS>
    bd = [ts.project_operator(Wd) for _, Wd in ops]  # P O^dagger|GS>
    store, n, Hq, (jl, R) = block_krylov(ts, b + bd, nkry, residual=True)
    ws, U = np.linalg.eigh(Hq)
    res = np.zeros(n + 1)
    if len(R):
        Ul = U[jl:, :]
        res[1:] = np.sqrt(np.abs(np.einsum("jf,jk,kf->f", Ul.conj(), R, Ul)))
    Uh = U.conj().T
    e = np.concatenate([[0.0], ws])
    def overlaps(v):  # <q_i|v> for the n basis vectors, on the host
        return backend.to_host(backend.row_overlaps(store, v, n))[:n]

    xp = backend.xp()
    # A real basis (a real state and H, real_starts) has a real U as well,
    # and there an operator that is a phase times a real one (S_y = i times
    # a real operator) is applied as the real one, with the phase put back
    # on the matrix: real arithmetic in every application, in the overlaps
    # and in the Ritz transform, and S_y comes out imaginary with exact
    # zeros (_exact_phase_form), the form conductance.py reads to run its
    # own sums in real arithmetic. A complex basis keeps complex arithmetic.
    real_basis = not np.issubdtype(np.dtype(store.dtype), np.complexfloating)
    mats = []
    for (W, Wd), bk, bdk in zip(ops, b, bd):
        Wc = compress_mpo(W)
        phase, Wc = real_up_to_phase(Wc) if real_basis else (1.0, Wc)
        op = ts.operator(Wc)
        M = np.zeros((n + 1, n + 1), dtype=complex)
        M[0, 0] = phase * ts.expectation(op)
        M[1:, 0] = Uh @ overlaps(ts.flat(bk))
        M[0, 1:] = (Uh @ overlaps(ts.flat(bdk))).conj()
        if not column:
            real = real_basis and not any(xp.iscomplexobj(x) for x in Wc)
            G = np.zeros((n, n), dtype=float if real else complex)  # <q_i|O|q_j>
            nb = ts.batch_size(op)
            for j0 in range(0, n, nb):  # the basis in stacks, not vector by vector
                j1 = min(n, j0 + nb)
                ys = ts.flat_many(ts.apply_many(op, ts.unflat_many(store[j0:j1])))
                G[:, j0:j1] = _overlap_columns(store, ys, n)
                del ys
            M[1:, 1:] = phase * (Uh @ G @ U)
        if real_basis:
            M = _exact_phase_form(M, phase)
        mats.append(M)
    if residual:
        return e, mats, res
    return e, mats


def ritz_vector(basis, u):
    """The Ritz vector with Krylov coefficients u, as a tangent vector;
    basis is (ts, store, m)."""
    ts, store, m = basis
    xp = backend.xp()
    c = np.zeros(store.shape[0], dtype=np.result_type(u.dtype, store.dtype))
    c[:m] = u
    return ts.unflat(xp.asarray(c) @ store)


def correlator_poles(Ms, Ws, Wb, Wadag=None, nkry=100, variance=False,
                     pad_bonds=None, expand_bonds=None, expand_rounds=1):
    """The poles of <GS|A delta(w-H+E0) B|GS> from site arrays alone: Ms
    the ground-state MPS, Ws the Hamiltonian MPO, Wb the MPO of B and Wadag
    that of A^dagger (None when A = B^dagger). Shared by every backend that
    can export its site tensors. Returns (w, weights, variances, dim), the
    variances None unless variance=True. pad_bonds (None, an int or "full")
    pads the bonds of the MPS first, see the module docstring; dim is then
    the dimension of the padded tangent space. expand_bonds (None or an
    int) adds at most that many directions per bond from the two-site
    discarded space of B|GS> and A^dagger|GS> instead (module docstring,
    "Expanded bonds"), repeated up to expand_rounds times."""
    ops = [Wb] + ([] if Wadag is None else [Wadag])
    ts = TangentSpace(Ms, Ws, pad_bonds=pad_bonds, expand_bonds=expand_bonds,
                      expand_ops=ops, expand_rounds=expand_rounds)
    b = ts.project_operator(Wb)
    a = None if Wadag is None else ts.project_operator(Wadag)
    ws, weights, U, (store, m) = poles(ts, b, nkry, a=a)
    var = None
    if variance and m > 0:
        var = np.array([ts.two_site_variance(ritz_vector((ts, store, m), U[:, k]))
                        for k in range(len(ws))])
    return ws, weights, var, ts.dimension()


def spectrum_from_arrays(Ms, Ws, ops, nkry, pad_bonds=None,
                         expand_bonds=None, expand_rounds=1, residual=False,
                         column=False):
    """projected_spectrum from site arrays: ops is a list of (MPO of O,
    MPO of O^dagger). Returns (e, mats, dim), and with residual=True the
    residual norm of every state as a fourth entry. pad_bonds and
    expand_bonds as in correlator_poles; the expansion is fed by every O
    and O^dagger. column=True builds row and column 0 only
    (projected_spectrum)."""
    ts = TangentSpace(Ms, Ws, pad_bonds=pad_bonds, expand_bonds=expand_bonds,
                      expand_ops=[W for pair in ops for W in pair],
                      expand_rounds=expand_rounds)
    e, mats, res = projected_spectrum(ts, ops, nkry, residual=True,
                                      column=column)
    if residual:
        return e, mats, ts.dimension(), res
    return e, mats, ts.dimension()
