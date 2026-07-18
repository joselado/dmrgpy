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

Without JAX or numba (see below), dmrg.py's two_site_heff()/one_site_heff()
get the exact pairwise ITensor.__mul__ chain (`t = t*L; t = t*H_i; ...`)
that existed before this module did -- never routed through numpy.einsum,
per the BLAS-dispatch gap above.

== numba path (default ON when importable) ==

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
_get_numba_contractor() below compiles and caches one specialized
function per distinct (permutation, reshape) pattern -- mirroring
compile_contraction()/_einsum_jax's per-pattern JIT-cache-reuse strategy
for the JAX path, for the same underlying reason (a bond's contraction
*shape* recurs every sweep once maxdim saturates, so the compile cost is
paid once and amortized over the whole run). Confirmed safe across
process restarts that a numba cache=True disk cache correctly
distinguishes different closures sharing the same source text by their
different captured constants, rather than colliding on the first compiled
variant seen (a real risk this pattern could otherwise hit silently).

USE_NUMBA nonetheless defaults to False. The per-contraction win above is
real but so is numba's own compile-time cost, and it is NOT amortized
away by cache=True the way the docstring above might suggest -- confirmed
directly: a first-ever @njit call in a fresh process pays a large
one-time LLVM/typing-system warm-up cost (~3s measured in isolation for
one specific closure) *on top of* per-signature compilation; cache=True
eliminates the *recompilation* of a previously-seen signature on a later
process (~3s -> ~0.2s measured), but not that first process-wide warm-up
tax, and DMRG naturally exercises several *distinct* contraction shapes
within a single run (edge bonds have smaller bond dimension than
saturated bulk bonds) each needing their own cache entry. Measured
end-to-end on this package's own regression suite -- exactly the
one-shot-script usage pattern CLAUDE.md describes as this library's real
verification workflow (`cd examples/... && python <script>.py`), not a
long-lived multi-calculation session -- selftest_dmrg.py (a handful of
small DMRG runs) took 1.83s with USE_NUMBA defaulted on vs 1.15s with it
off: a ~60% *regression*, not a win, because the fixed compile tax
dominates a script this short. A sustained, multi-call session (repeated
gs_energy() calls reusing the same shapes) does recover the steady-state
win measured above -- confirmed directly, a 14-site/maxdim=60 Heisenberg
chain's second and later gs_energy() calls in the same process ran at
~0.78-0.80s, at or slightly below the compiled ITensor v3 backend's own
~0.81s, after the first call's ~0.997s paid the compile tax -- so this is
a genuine case where "on by default" is the wrong call despite the
underlying kernel being a real, verified improvement, exactly the same
shape of finding as the JAX CPU regression above. Opt in explicitly
(`kernels.USE_NUMBA = True`, set once at the start of a script) for
workloads that call gs_energy()/tevol_method()/etc. repeatedly in one
process at small-to-moderate bond dimension -- parameter sweeps, or
anything long-running enough to amortize the fixed compile cost.

"Small-to-moderate bond dimension" matters: the ~1.7x-2x per-contraction
win itself is specific to tensors small enough that numpy's own
Python-level dispatch overhead is a meaningful fraction of the call cost.
Measured directly at a genuinely large scale (n=24, maxdim=300 -- closer
to where this library's real workloads live than the maxdim=60 benchmark
above) numba showed *no* win even with compile cost fully excluded from
the comparison: 105.7s without numba vs 108.1s with, a wash/slight loss,
not an improvement -- at that size the actual BLAS matmul FLOPs dominate
the per-call cost, leaving little Python-overhead fraction left for a
compiled kernel to shave off. So opting in is only worth it for small-
bond-dimension, sustained (many-call) sessions specifically -- it
recovers rough parity with (not a clear win over) the compiled ITensor v3
backend there, not a general-purpose speedup at any scale.

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
    """Best-effort default: only auto-enable JAX if it reports a non-CPU
    device, where the per-call dispatch overhead responsible for the
    measured CPU regression (see this module's docstring) may not apply.
    UNTESTED on real GPU/TPU hardware -- if it turns out JAX is a net loss
    there too, this heuristic should be replaced with an explicit opt-in
    only (set USE_JAX = True yourself, never automatic)."""
    if not _HAVE_JAX:
        return False
    try:
        return any(d.platform != "cpu" for d in jax.devices())
    except Exception:
        return False


# Runtime on/off switch, independent of whether jax is importable -- flip
# to True to force JAX on (e.g. on GPU hardware this heuristic doesn't
# recognize), or to False to force plain NumPy even where JAX would
# otherwise auto-enable. available() reflects what will actually be used,
# not just what's importable.
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
    tensordot rather than slower) and cached per distinct pattern so a
    matvec chain that calls the same pairwise contraction shape on every
    Lanczos iteration only pays the compile cost once."""
    key = (perm_a, shape_a2, perm_b, shape_b2)
    f = _numba_contractor_cache.get(key)
    if f is not None:
        return f

    @numba.njit(cache=True, fastmath=True)
    def _contract(A, B):
        At = np.transpose(A, perm_a).copy().reshape(shape_a2)
        Bt = np.transpose(B, perm_b).copy().reshape(shape_b2)
        return At @ Bt

    _numba_contractor_cache[key] = _contract
    return _contract


def _plan_numba_chain(pieces, order_in, order_out):
    """Precompute the fixed sequence of pairwise transpose+reshape+matmul
    steps for v(order_in) * pieces[0] * pieces[1] * ... , using mul_plan()
    (tensor.py) on index *structure* alone -- called once per matvec build
    (once per bond), not once per Lanczos iteration. Returns (steps,
    final_perm, out_shape) where steps is a list of (contractor,
    out_shape_2factors, piece_array) and final_perm/out_shape reorder the
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
        contractor = _get_numba_contractor(
            perm_a, (a_free_size, a_axes_size), perm_b, (a_axes_size, b_free_size))
        out_shape = a_free_dims + b_free_dims
        steps.append((contractor, out_shape, p.array))
        cur_inds = tuple(cur_inds[i] for i in a_free) + tuple(p.inds[j] for j in b_free)

    final_perm = tuple(cur_inds.index(ind) for ind in order_out)
    return steps, final_perm


def _make_matvec_numba(pieces, order_in, shape_in, order_out):
    steps, final_perm = _plan_numba_chain(pieces, order_in, order_out)
    out_shape = tuple(ind.dim for ind in order_out)

    def matvec(v_flat):
        cur = v_flat.reshape(shape_in)
        for contractor, step_shape, piece_arr in steps:
            cur = contractor(cur, piece_arr).reshape(step_shape)
        return np.transpose(cur, final_perm).reshape(out_shape).reshape(-1)

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
    compiled jax.numpy.einsum call if JAX is available and enabled,
    else the same pairwise ITensor.__mul__ chain this package always
    used."""
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

    from .tensor import ITensor

    def matvec(v_flat):
        t = ITensor(tuple(order_in), v_flat.reshape(shape_in))
        for p in pieces:
            t = t * p
        return t.transpose_to(order_out).reshape(-1)

    return matvec
