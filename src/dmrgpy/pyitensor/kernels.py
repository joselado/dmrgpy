"""Optional JAX-accelerated contraction kernel for the one hot loop in this
whole package: the effective-Hamiltonian matvec inside DMRG's Lanczos
eigensolve (dmrg.py's _local_ground_state/_local_ground_state_penalized)
and TDVP's Krylov exponentiation (tdvp.py's _lanczos_expm_multiply). Both
call the *same* matvec dozens of times per bond (niter ~ 30-50), so it is
by a wide margin the most-executed piece of arithmetic in the library --
everywhere else (AutoMPO compilation, sweep bookkeeping, SVD truncation)
runs once or a handful of times per bond, not dozens.

Without JAX, dmrg.py's two_site_heff()/one_site_heff() build this matvec
exactly as before this module existed: a chain of pairwise
ITensor.__mul__() calls (`t = t*L; t = t*H_i; ...`), each of which
dispatches to numpy.tensordot and therefore BLAS. That's deliberately
*not* routed through numpy.einsum here, despite this module's whole
purpose being "express the contraction as one fused call" for the JAX
path below -- confirmed directly, expressing the same contraction as a
single numpy.einsum call (even with the contraction path precomputed
once via numpy.einsum_path and reused, avoiding re-solving it on every
call) was over 4 orders of magnitude slower than the pairwise chain for
a representative bond (approx. 45ms/call vs 4us/call): numpy's einsum
evaluator doesn't dispatch >2-operand contractions to BLAS the way
tensordot does, no matter what path is fed to `optimize=`. So the numpy
path here is a pure pass-through to the original approach, not a
numpy.einsum reimplementation of it.

If JAX *is* available, the identical contraction *is* expressed as one
fused jax.numpy.einsum call, JIT-compiled -- there XLA's own compiler
does the operand-fusion decision itself rather than relying on einsum's
runtime path search, and does not have numpy's >2-operand-BLAS gap.

The JIT caching story is the other reason this is worth doing at all: a
naive `jax.jit` applied fresh inside each call to two_site_heff() would
only ever get reused across the niter Lanczos iterations *within one
bond*, since each such closure is a distinct Python function object and
JAX's cache is keyed on function identity as well as argument
shapes/dtypes. Compiling once at *module level* instead (_einsum_jax
below, with the einsum subscript string marked as a static/hashable
argument) means the cache is shared across every bond in every sweep
that happens to share the same contraction pattern -- which, once bond
dimensions saturate at maxdim after the first sweep or two, is most of
them. compile_contraction() labels each distinct Index deterministically
(by first-occurrence order, scanning operands in the same fixed order
every call), so two structurally identical bonds produce byte-identical
subscript strings and therefore hit the same compiled cache entry, not
just equivalent-but-distinct ones.

Never a hard dependency: everything in this package works with only
NumPy/SciPy installed (see USE_JAX/available() below) -- this module is
purely a performance option layered on top of already-correct code, and
the numpy path is exactly the code path that existed before it.
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

# Runtime on/off switch, independent of whether jax is importable -- set
# to False to force the plain-NumPy path even when JAX is available (e.g.
# to isolate whether a numerical discrepancy is JAX-related), or flip back
# to True after importing jax late. available() reflects what will
# actually be used, not just what's importable.
USE_JAX = _HAVE_JAX

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
