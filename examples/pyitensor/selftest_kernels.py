# Self-check for the optional JAX-accelerated kernel layer
# (dmrgpy.pyitensor.kernels). Runs regardless of whether JAX is actually
# installed: if it isn't, available() is False and the JAX-specific checks
# are skipped, but make_matvec()'s NumPy fallback path (identical to what
# ran before this module existed) is still exercised by every other
# selftest_*.py script's DMRG/TDVP runs.
import numpy as np

from dmrgpy.pyitensor import AutoMPO, randomMPS, to_mpo
from dmrgpy.pyitensor import kernels
from dmrgpy.pyitensor.dmrg import _all_right_environments, two_site_heff
from dmrgpy.pyitensor.sites import SiteX


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def _build_bond1_heff():
    np.random.seed(0)
    n = 6
    sites = SiteX([2] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    H = to_mpo(ampo, cutoff=1e-14, maxdim=5000)
    psi = randomMPS(sites, 12)
    right_env = _all_right_environments(H, psi)
    # bond i=1 (spanning sites 1,2) needs R covering sites i+2..N = 3..N
    R, Rbra = right_env[3]
    return H, psi, R, Rbra


def test_compile_contraction_labels_deterministically():
    from dmrgpy.pyitensor.index import Index
    a, b, c = Index(2, "Site"), Index(3, "Link"), Index(4, "Link")
    subs1 = kernels.compile_contraction([(a, b), (b, c)], (a, c))
    subs2 = kernels.compile_contraction([(a, b), (b, c)], (a, c))
    check("compile_contraction is deterministic for the same Index objects", subs1 == subs2)
    check("compile_contraction produces a valid einsum subscript shape", subs1 == "ab,bc->ac")


def test_jax_and_numpy_matvec_agree():
    if not kernels._HAVE_JAX:
        print("[SKIP] JAX not installed -- nothing to compare against NumPy")
        return
    H, psi, R, Rbra = _build_bond1_heff()

    kernels.USE_JAX = True
    matvec_jax, order_in, shape, x0 = two_site_heff(None, None, H, psi, 1, R, Rbra)
    kernels.USE_JAX = False
    matvec_np, *_ = two_site_heff(None, None, H, psi, 1, R, Rbra)
    kernels.USE_JAX = kernels._HAVE_JAX  # restore default

    out_jax = matvec_jax(x0)
    out_np = matvec_np(x0)
    rel_err = np.max(np.abs(out_jax - out_np)) / np.max(np.abs(out_np))
    check("JAX and NumPy matvec kernels agree to near machine precision "
          "(relative error={:.2e})".format(rel_err), rel_err < 1e-10)


def test_numpy_fallback_is_the_pairwise_chain_not_einsum():
    # Regression test: an earlier version of make_matvec() routed the
    # NumPy fallback through a single fused numpy.einsum call (mirroring
    # the JAX path), reasoning it would be a modest win via einsum's own
    # path optimizer. Confirmed directly instead to be catastrophically
    # slower (~4 orders of magnitude) than the pairwise ITensor.__mul__
    # chain for a representative bond, because numpy's einsum evaluator
    # doesn't dispatch >2-operand contractions to BLAS the way tensordot
    # does, regardless of the contraction path used. This just checks the
    # NumPy path still completes a modest number of calls quickly instead
    # of silently regressing back to that.
    import time
    kernels.USE_JAX = False
    H, psi, R, Rbra = _build_bond1_heff()
    matvec, order_in, shape, x0 = two_site_heff(None, None, H, psi, 1, R, Rbra)
    kernels.USE_JAX = kernels._HAVE_JAX
    t0 = time.time()
    for _ in range(50):
        matvec(x0)
    dt = time.time() - t0
    check("NumPy fallback matvec: 50 calls complete in well under a second "
          "(took {:.4f}s)".format(dt), dt < 1.0)


if __name__ == "__main__":
    test_compile_contraction_labels_deterministically()
    test_jax_and_numpy_matvec_agree()
    test_numpy_fallback_is_the_pairwise_chain_not_einsum()
    print("All pyitensor kernels checks passed.")
