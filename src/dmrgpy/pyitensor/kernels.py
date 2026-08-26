"""Optional, OFF-BY-DEFAULT-ON-CPU JAX-accelerated contraction kernel for
the effective-Hamiltonian matvec inside DMRG's Lanczos eigensolve
(dmrg.py's _local_ground_state/_local_ground_state_penalized) and TDVP's
Krylov exponentiation (tdvp.py's _lanczos_expm_multiply).

Read this before flipping USE_JAX on: measured directly, on CPU, at
maxdim 60/150/300 on a 14-site chain (after dmrg.py's own contraction-
order fix -- see tensor.py's contract_many() -- made this matvec the
actual bottleneck rather than something else masking its cost), the JAX
path is 1.5x-2.5x *slower* than the plain NumPy path, not faster, and the
gap only narrows (not crosses over) as maxdim grows in that range. The
per-call numpy<->jax array conversion and XLA dispatch overhead
apparently exceeds the compute saved for tensors this small -- this
wasn't caught when the module was first written because at the time
dmrg.py's environment-building bug (now fixed) dominated total runtime by
such a wide margin that the matvec's own cost, JAX or not, barely
registered either way. USE_JAX therefore defaults to on only when JAX
itself reports a non-CPU device (see _detect_default_use_jax() below) --
UNTESTED on real GPU/TPU hardware (none was available while writing this),
so treat that heuristic as a reasonable guess, not a verified claim; on a
CPU-only install this now behaves exactly as if the module didn't exist.
Re-confirmed after the plain-path fix below (see its own section): on a
dynamical-METTS run JAX was still ~4x slower than the new plain default
(210s vs. 52s at n=10, nsamples=40, nt=10) -- the conclusion above holds
against the faster baseline too, not just the old one.

If JAX *is* enabled, the mechanism is: express the whole contraction as
one fused jax.numpy.einsum call (XLA fuses the operands itself, unlike
numpy.einsum -- confirmed separately that a single numpy.einsum call for
the same contraction, even with the path precomputed once via
numpy.einsum_path and reused rather than resolved every call, is over 4
orders of magnitude slower than plain pairwise ITensor.__mul__/tensordot
calls, since numpy's einsum evaluator doesn't dispatch >2-operand
contractions to BLAS the way tensordot does), JIT-compiled once at
*module level* (_einsum_jax below, subscripts marked as a static/hashable
argument) so the compiled-code cache is shared across every bond/sweep
with a matching contraction shape, not just within one bond's Lanczos
loop -- a naive jax.jit applied fresh inside each call would only ever
get reused within a single bond, since each such closure is a distinct
Python function object. compile_contraction() labels each distinct Index
deterministically (first-occurrence order, operands always listed the
same way) so structurally identical bonds produce byte-identical
subscript strings and hit the same cache entry.

== plain path (always on, the default when JAX/numba are unavailable or
disabled) ==

Without JAX or numba, dmrg.py's two_site_heff()/one_site_heff()/
zero_site_heff() get _make_matvec_planned() (below): the same precomputed
transpose+reshape+matmul plan the numba path uses (_plan_contraction_chain,
shared by both), just executed directly in NumPy instead of through a
compiled closure. This package originally used a bare pairwise
ITensor.__mul__ chain here (`t = t*L; t = t*H_i; ...`), which re-derives
that same plan (mul_plan()'s Index-matching, from tensor.py) from scratch
on *every* matvec call -- cheap enough per call to be invisible for DMRG's
ground-state sweep (a handful of calls per bond), but the dominant cost for
TDVP/METTS, which call the same matvec build thousands of times more per
run (many time steps x many Lanczos/Krylov iterations x many METTS
samples): confirmed via cProfile on a dynamical-METTS run, 57% of
cumulative time was in ITensor.__mul__ and its mul_plan()/Index-matching
machinery before this fix. Precomputing the plan once per bond/loop instead
of once per call closes most of that gap with no new dependency and no
compile tax: confirmed directly, 1.4x-2x end-to-end wall-time win on real
dynamical METTS runs (9.15s->6.68s at n=8; 152.0s->74.6s at n=10,
nsamples=40, nt=10), bit-identical results (this reorders *when* the plan
is computed, not *which* pairs get contracted first, so unlike tensor.py's
contract_many() fix elsewhere in this package it does not perturb
floating-point summation order).

== numba path (opt-in; see below for why it isn't the default) ==

Unlike JAX, a numba-compiled kernel operates directly on the same numpy
arrays already in memory -- no host<->device array conversion, no
separate dispatch runtime -- so it doesn't carry the per-call overhead
responsible for JAX's measured CPU regression. Measured directly: a
single representative pairwise contraction (a (60,5,60) environment
tensor against a (5,5,2,2) MPO tensor, matching maxdim=60/mpo-bond=5/
physical-dim=2 shapes from _extend_left/two_site_heff) was ~1.7x-2x
faster via a numba-compiled transpose+reshape+matmul than via
np.tensordot on the same arrays -- numpy's tensordot carries real Python-
level shape/axes-analysis overhead on every call that a compiled kernel
skips, even though both ultimately call the same underlying BLAS matmul
for the actual arithmetic (this is *not* numba's raw loops beating BLAS:
a naive 5-nested-loop numba kernel for the same contraction was 2x
*slower* than tensordot, confirmed directly -- the win only appears once
the kernel itself also reshapes to a 2D matmul and lets BLAS do the
arithmetic, same as tensordot does internally, just without tensordot's
own per-call Python overhead on top).

Getting that win in practice requires the transpose permutation and
reshape target to be baked into the compiled function as closure-captured
constants, not passed as ordinary runtime arguments -- confirmed
directly, passing the *same* permutation as a plain function argument
instead of a closure constant was consistently *slower* than tensordot,
not faster, apparently because numba's optimizer can't treat a
runtime-typed tuple argument's transpose as constant-foldable the way it
does for a captured closure value known to be fixed at compile time. So
_get_numba_contractor() below compiles and caches (in-process only, see
its own docstring for why not also on-disk via cache=True -- a real,
confirmed crash risk, not just a theoretical one) one specialized
function per distinct (permutation, reshape) pattern -- mirroring
compile_contraction()/_einsum_jax's per-pattern JIT-cache-reuse strategy
for the JAX path, for the same underlying reason (a bond's contraction
*shape* recurs every sweep once maxdim saturates, so the compile cost is
paid once and amortized over the whole run, within one process).

USE_NUMBA nonetheless defaults to False, and is worth even less now than
when this section was first written, because the baseline it has to beat
(the plain path above) is no longer the naive per-call ITensor.__mul__
chain -- it's the same precomputed-plan technique numba itself uses,
just without JIT. Re-measured directly against *that* baseline (not the
old naive one):

- Ground-state DMRG, independent chains (e.g. a parameter sweep: a fresh
  Spin_Chain per point, each starting from its own random MPS -- see
  mpscpp3's own randomMPS()-based default_mps(), which pyitensor's DMRG
  mirrors): the "steady state" this section used to describe *does not
  reliably occur*. Confirmed directly at n=14/maxdim=60: the in-process
  numba contractor cache kept growing across 4 independent gs_energy()
  calls (106 -> 144 -> 194 -> 210 entries, never leveling off) because
  each fresh random starting MPS produces its own bond-dimension growth
  trajectory through the sweeps before saturating at maxdim, so the exact
  edge-bond shapes hit differ call to call -- there is no fixed set of
  shapes to amortize a compile over. All 4 calls took 6.5s-42.6s each,
  vs. the plain path's steady ~0.7-0.8s; at n=24/maxdim=300 a single
  fresh-chain call took 133.5s with numba vs. 25.2s plain (5.3x slower,
  not the previously-reported "wash" -- that comparison was against the
  slower naive baseline).
- TDVP/METTS (a single long run that reuses the same matvec shapes many
  times once bond dimension saturates -- the one scenario where numba's
  in-process cache *does* stabilize): confirmed directly on a dynamical
  METTS run (n=10, nsamples=40, nt=10), numba's own compile tax is now a
  much larger fraction of a much smaller total: 328.8s cold vs. 51.8s
  plain (6.3x worse) on the very first process, and even fully warmed up
  (repeated in the same process after compilation) numba was still only
  ~11% faster than plain (46.1s vs 51.8s) -- a real but marginal edge,
  down from the 2x this section used to claim, because the plain path no
  longer leaves that much of numba's old advantage on the table.

So: still off by default, and now only worth opting into
(`kernels.USE_NUMBA = True`, set once at the start of a script) for a
narrower case than before -- a single long-running TDVP/METTS-like
session at small-to-moderate bond dimension where the same shapes recur
many times *within one run*, and even there the win is modest (~10%),
not the roughly-2x this section previously reported (that number was
real for numba-vs-naive-ITensor.__mul__, just no longer the comparison
that matters after the plain-path fix above). Independent-chain workloads
(parameter sweeps, or anything that builds a fresh random-start MPS each
time) should not enable it at all -- the per-call compile tax is paid
over and over with no steady state to amortize into.

Never a hard dependency either way: everything in this package works with
only NumPy/SciPy installed.
"""

import functools
import string

import numpy as np

try:
    import jax
    import jax.numpy as jnp
    # ITensor amplitudes and Hamiltonians are complex128 throughout (see
    # tensor.py's module docstring); JAX defaults every array to float32/
    # complex64 unless this is set, which would silently corrupt
    # cutoff-based truncation and every energy/overlap computed downstream.
    jax.config.update("jax_enable_x64", True)
    _HAVE_JAX = True
except ImportError:
    _HAVE_JAX = False

try:
    import numba
    _HAVE_NUMBA = True
except ImportError:
    _HAVE_NUMBA = False


def _detect_default_use_jax():
    """Always False. This used to auto-enable JAX whenever jax.devices()
    reported a non-CPU device, on the theory that a GPU's compute win might
    outrun the per-call dispatch overhead responsible for the measured CPU
    regression (see this module's docstring). That theory has now been
    tested on real hardware and it is wrong, so the heuristic is gone --
    following the instruction the docstring itself left for this outcome
    ("replaced with an explicit opt-in only, never automatic").

    Measured on an NVIDIA H200 (jax
    0.7.1 with its CUDA plugin, benchmarks/gpu/pyitensor_gpu_probe.py):
    ground-state DMRG on a Heisenberg chain was 5-11x *slower* on this
    path than on the plain one -- 0.26s vs 2.91s (n=16, maxm=60), 0.48s
    vs 2.25s (n=20, maxm=100), 0.53s vs 3.50s (n=24, maxm=200), same
    energies to 8e-14. So the auto-enable was not a harmless guess: on a
    GPU node it silently made every pyitensor run several times slower.

    The reason is the conversion, not the device. The same job's primitive
    microbenchmark (benchmarks/gpu/gpu_microbench.py) has the H200 doing
    the two-site matvec 5.5x faster than one Xeon 6248 core at chi=64 and
    688x faster at chi=1024 -- while a single host->device->host round trip
    for one theta-sized array costs *more* than the matvec itself at every
    chi above 64 (4.3ms vs 1.5ms at chi=512; 22.0ms vs 8.3ms at chi=1024).
    This path converts per call, so it pays that round trip on every
    Lanczos iteration and cannot win however fast the device is. Winning
    requires arrays that stay resident on the device, which is a port and
    not a flag -- see docs/pyitensor_gpu_port_plan.md.
    """
    return False


# Runtime on/off switch, independent of whether jax is importable. Always
# off by default now (see _detect_default_use_jax above: the automatic
# version made GPU runs 5-11x slower, and CPU runs 1.5-2.5x slower). Set it
# True yourself to measure this path; available() reflects what will
# actually be used, not just what's importable.
USE_JAX = _detect_default_use_jax()

# Off by default despite the real ~1.7x-2x per-contraction win (see this
# module's docstring): numba's own compile-time tax outweighs it for this
# library's actual one-shot-script usage pattern, and is only recovered
# in a sustained multi-call session. Opt in explicitly for that case.
USE_NUMBA = False

_LETTERS = string.ascii_letters  # 52 distinct single-character einsum labels


def available():
    """Whether make_matvec() will actually dispatch to the JAX path."""
    return _HAVE_JAX and USE_JAX


def available_numba():
    """Whether make_matvec() will actually dispatch to the numba path.
    Checked before available() -- see make_matvec()'s dispatch order."""
    return _HAVE_NUMBA and USE_NUMBA


def compile_contraction(operand_inds, output_inds):
    """operand_inds: list of Index-tuples, one per array operand (in a
    fixed, deterministic order -- callers must always list operands the
    same way for the same kind of contraction, see this module's
    docstring on why that's what makes JIT cache reuse work at all).
    output_inds: the desired final Index order. Returns a numpy/jax
    einsum subscript string: each distinct Index (by id+plev, i.e.
    Index.__eq__/__hash__) gets one letter, assigned in first-occurrence
    order while scanning operand_inds then output_inds."""
    label = {}

    def get_label(ind):
        if ind not in label:
            if len(label) >= len(_LETTERS):
                raise ValueError("compile_contraction: more than {} distinct "
                                  "indices, out of einsum labels".format(len(_LETTERS)))
            label[ind] = _LETTERS[len(label)]
        return label[ind]

    in_subs = [''.join(get_label(ind) for ind in inds) for inds in operand_inds]
    out_sub = ''.join(get_label(ind) for ind in output_inds)
    return ','.join(in_subs) + '->' + out_sub


if _HAVE_JAX:
    @functools.partial(jax.jit, static_argnames=("subscripts",))
    def _einsum_jax(subscripts, v, pieces):
        return jnp.einsum(subscripts, v, *pieces)


_numba_contractor_cache = {}


def _get_numba_contractor(perm_a, shape_a2, perm_b, shape_b2):
    """A numba-compiled `f(A, B) -> (A transposed+reshaped to shape_a2) @
    (B transposed+reshaped to shape_b2)`, specialized (permutation/shape
    baked in as closure constants, not arguments -- see this module's
    docstring for why that distinction is what makes this faster than
    tensordot rather than slower) and cached (in-process, via
    _numba_contractor_cache) per distinct pattern so a matvec chain that
    calls the same pairwise contraction shape on every Lanczos/Krylov
    iteration only pays the compile cost once per process.

    Deliberately NOT `cache=True` (numba's own on-disk cache) -- confirmed
    directly this is unsafe for this specific pattern: `_contract` is a
    distinct closure (same source text, different perm/shape constants
    baked in) every time this function compiles a new one, and numba's
    per-source-location on-disk overload index has no eviction, so it
    accumulates one entry per distinct closure *forever* across every
    process that ever ran this code (908 stale `.nbc` files were found
    for this one closure in a worktree used across several days of
    sessions). Past some accumulation, loading the disk cache started
    returning a corrupt overload and crashing with `RuntimeError: In
    'NRT_adapt_ndarray_to_python', 'descr' is NULL` deep inside a real
    dynamical-METTS TDVP run -- reproduced reliably in the affected
    worktree, and confirmed fixed both by deleting the accumulated cache
    files and by dropping `cache=True` entirely. The in-process dict cache
    above is what actually matters for this module's own "sustained
    multi-call session" use case (see the module docstring) -- cache=True
    only ever aimed to avoid recompilation on a *later, separate* process,
    a narrower benefit not worth this failure mode."""
    key = (perm_a, shape_a2, perm_b, shape_b2)
    f = _numba_contractor_cache.get(key)
    if f is not None:
        return f

    @numba.njit(fastmath=True)
    def _contract(A, B):
        At = A.transpose(perm_a).copy().reshape(shape_a2)
        Bt = B.transpose(perm_b).copy().reshape(shape_b2)
        return At @ Bt

    _numba_contractor_cache[key] = _contract
    return _contract


def _plan_contraction_chain(pieces, order_in, order_out):
    """Precompute the fixed sequence of pairwise transpose+reshape+matmul
    steps for v(order_in) * pieces[0] * pieces[1] * ... , using mul_plan()
    (tensor.py) on index *structure* alone -- called once per matvec build
    (once per bond/Lanczos-or-Krylov loop), not once per iteration, since
    pieces/order_in/order_out are fixed for the lifetime of one such loop
    and only the numeric contents of the vector being multiplied change
    from call to call. Shared by both the numba path (_make_matvec_numba)
    and the plain-NumPy path (_make_matvec_planned) below -- the only
    difference between them is whether the per-step transpose+reshape+
    matmul is executed via a numba-compiled closure or directly in NumPy.

    Returns (steps, final_perm) where steps is a list of (perm_a, shape_a2,
    perm_b, shape_b2, out_shape, piece_array) and final_perm reorders the
    chain's natural trailing free-index order into order_out."""
    from .tensor import mul_plan

    cur_inds = tuple(order_in)
    steps = []
    for p in pieces:
        a_axes, b_axes, a_free, b_free = mul_plan(cur_inds, p.inds)
        a_free_dims = tuple(cur_inds[i].dim for i in a_free)
        a_axes_dims = tuple(cur_inds[i].dim for i in a_axes)
        b_free_dims = tuple(p.inds[j].dim for j in b_free)
        a_free_size = int(np.prod(a_free_dims)) if a_free_dims else 1
        a_axes_size = int(np.prod(a_axes_dims)) if a_axes_dims else 1
        b_free_size = int(np.prod(b_free_dims)) if b_free_dims else 1
        perm_a = tuple(a_free + a_axes)
        perm_b = tuple(b_axes + b_free)
        out_shape = a_free_dims + b_free_dims
        steps.append((perm_a, (a_free_size, a_axes_size),
                      perm_b, (a_axes_size, b_free_size),
                      out_shape, p.array))
        cur_inds = tuple(cur_inds[i] for i in a_free) + tuple(p.inds[j] for j in b_free)

    final_perm = tuple(cur_inds.index(ind) for ind in order_out)
    return steps, final_perm


def _make_matvec_numba(pieces, order_in, shape_in, order_out):
    steps, final_perm = _plan_contraction_chain(pieces, order_in, order_out)
    out_shape = tuple(ind.dim for ind in order_out)
    contractors = [
        (_get_numba_contractor(perm_a, shape_a2, perm_b, shape_b2), step_shape, piece_arr)
        for perm_a, shape_a2, perm_b, shape_b2, step_shape, piece_arr in steps]

    def matvec(v_flat):
        cur = v_flat.reshape(shape_in)
        for contractor, step_shape, piece_arr in contractors:
            cur = contractor(cur, piece_arr).reshape(step_shape)
        return cur.transpose(final_perm).reshape(out_shape).reshape(-1)

    return matvec


def _make_matvec_planned(pieces, order_in, shape_in, order_out):
    """Same precomputed transpose+reshape+matmul plan as _make_matvec_numba
    (_plan_contraction_chain), executed directly in NumPy instead of through
    a numba-compiled closure -- no JIT compile tax, no host<->device
    transfer, always available. This is make_matvec()'s default fallback
    (see its docstring): it replaces what used to be a bare ITensor.__mul__
    chain that re-derived this same structural plan (mul_plan()'s Index
    matching, from tensor.py) from scratch on *every* call. That cost is
    negligible for DMRG's ground-state sweep (a handful of matvec calls per
    bond) but dominates for TDVP/METTS, which call the same matvec build
    thousands of times more per run (many time steps x many Lanczos/Krylov
    iterations x many METTS samples) -- confirmed directly via cProfile on a
    dynamical-METTS run: 57% of cumulative time was in ITensor.__mul__ and
    its mul_plan()/Index-matching machinery before this change.

    Confirmed directly: 1.4x-2x end-to-end wall-time win on real dynamical
    METTS runs (9.15s->6.68s at n=8; 152.0s->74.6s at n=10, nsamples=40,
    nt=10), bit-identical results against the ED cross-check in
    examples/finite_temperature/dynamical_metts_VS_ED. Unlike tensor.py's
    contract_many() fix elsewhere in this package, this does not change
    floating-point summation order (same tensordot-equivalent transpose/
    reshape/matmul per step, in the same left-to-right order as before) --
    it only removes redundant replanning, so it carries none of that
    change's excited-state-convergence-margin risk."""
    steps, final_perm = _plan_contraction_chain(pieces, order_in, order_out)

    # Hoist the *operand* side of every step out of the closure. `Bt`
    # depends only on `piece_arr` and `perm_b`, both fixed for the lifetime
    # of this matvec, yet it used to be recomputed on every call -- and
    # since a transposed array is generally not contiguous, that `.reshape`
    # forces a full copy of the operator tensor each time, not a view. For
    # an ordinary DMRG ground-state solve that is thousands of redundant
    # copies of the (largest) environment tensors: cProfile on a 30-site
    # fermionic chain at maxm=30 showed 5777 matvec calls carrying 4.37s of
    # self time and 122k reshape calls.
    #
    # Precomputing is exact, not an approximation: same arrays, same
    # left-to-right order, same matmuls, so results stay bit-identical --
    # this carries none of the summation-order risk a re-association would.
    # NOTE: every array operation below is a *method* (.transpose/.reshape/@),
    # never a numpy free function, so this path works unchanged whether the
    # arrays are numpy or JAX device arrays (backend.py). np.transpose() on a
    # device array would copy it to the host per call.
    prepared = []
    for perm_a, shape_a2, perm_b, shape_b2, step_shape, piece_arr in steps:
        Bt = piece_arr.transpose(perm_b).reshape(shape_b2)
        # skip a no-op transpose on the vector side too (common when the
        # planned axis order already matches)
        ident_a = (perm_a == tuple(range(len(perm_a))))
        prepared.append((perm_a, ident_a, shape_a2, Bt, step_shape))
    trivial_final = (final_perm == tuple(range(len(final_perm))))

    def matvec(v_flat):
        cur = v_flat.reshape(shape_in)
        for perm_a, ident_a, shape_a2, Bt, step_shape in prepared:
            At = cur if ident_a else cur.transpose(perm_a)
            cur = (At.reshape(shape_a2) @ Bt).reshape(step_shape)
        if not trivial_final:
            cur = cur.transpose(final_perm)
        return cur.reshape(-1)

    return matvec


def make_matvec(pieces, order_in, shape_in, order_out):
    """pieces: the fixed ITensor operators sandwiching the vector (e.g.
    [L, H_i, H_j, R], with None entries already filtered out by the
    caller). order_in/shape_in describe how a flat vector reshapes into
    the contraction's "v" operand; order_out is the result's Index order
    (see dmrg.py's two_site_heff()/one_site_heff() for why this can be a
    *different* Index identity than order_in even though it's the same
    shape -- L/R's own bra-side legs, not v's own link legs).

    Returns matvec(flat_vector) -> flat_vector: dispatches to the numba
    chain (_make_matvec_numba) if numba is available (checked first --
    see this module's docstring, numba wins on CPU where JAX regresses
    since it has no array-transfer overhead), else a single fused, JIT-
    compiled jax.numpy.einsum call if JAX is available and enabled, else
    the plain-NumPy precomputed-plan chain (_make_matvec_planned) -- always
    available, no opt-in required, see its own docstring for why this is
    make_matvec()'s default rather than a bare ITensor.__mul__ chain."""
    if available_numba():
        return _make_matvec_numba(pieces, order_in, shape_in, order_out)

    if available():
        operand_inds = [order_in] + [p.inds for p in pieces]
        subscripts = compile_contraction(operand_inds, order_out)
        piece_arrays = tuple(jnp.asarray(p.array) for p in pieces)

        def matvec(v_flat):
            v = jnp.asarray(v_flat).reshape(shape_in)
            out = _einsum_jax(subscripts, v, piece_arrays)
            return np.asarray(out).reshape(-1)

        return matvec

    return _make_matvec_planned(pieces, order_in, shape_in, order_out)
