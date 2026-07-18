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

Without JAX (the CPU default), dmrg.py's two_site_heff()/one_site_heff()
get the exact pairwise ITensor.__mul__ chain (`t = t*L; t = t*H_i; ...`)
that existed before this module did -- never routed through numpy.einsum,
per the BLAS-dispatch gap above.

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

_LETTERS = string.ascii_letters  # 52 distinct single-character einsum labels


def available():
    """Whether make_matvec() will actually dispatch to the JAX path."""
    return _HAVE_JAX and USE_JAX


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


def make_matvec(pieces, order_in, shape_in, order_out):
    """pieces: the fixed ITensor operators sandwiching the vector (e.g.
    [L, H_i, H_j, R], with None entries already filtered out by the
    caller). order_in/shape_in describe how a flat vector reshapes into
    the contraction's "v" operand; order_out is the result's Index order
    (see dmrg.py's two_site_heff()/one_site_heff() for why this can be a
    *different* Index identity than order_in even though it's the same
    shape -- L/R's own bra-side legs, not v's own link legs).

    Returns matvec(flat_vector) -> flat_vector: a single fused, JIT-
    compiled jax.numpy.einsum call if JAX is available (see this module's
    docstring for the caching story), otherwise the same pairwise
    ITensor.__mul__ chain this package always used."""
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
