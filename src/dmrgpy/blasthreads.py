"""Thread-count control for the BLAS/LAPACK underneath dmrgpy.

== Why this exists ==

DMRG is a great many *small* dense linear-algebra calls, not a few large
ones. A two-site tensor at bond dimension 30 with a spin-1/2 site is a
60x60 matrix; that is the size that hits `eigh`, `svd` and `matmul` tens of
thousands of times in one ground-state solve. At that size a multithreaded
BLAS can spend more time in thread barriers than in arithmetic, so the
library default -- one thread per core -- is at best not helping, and on a
busy or shared machine is a large loss.

Two separate effects, worth keeping apart.

The *mechanism* is that more threads is slower at this matrix size, full
stop. Measured on a 14-core host with MKL (numpy 2.4.6), going from one
thread to just two -- no oversubscription involved:

    60x60 complex svd     1.799 ms  ->  2.812 ms   (1.6x slower)
    60x60 complex eigh    0.774 ms  ->  1.792 ms   (2.3x slower)
    60x60 complex matmul  0.122 ms  ->  0.219 ms   (1.8x slower)

End to end, though, that mechanism alone is worth much less, because
dmrgpy's own Python/tensor-bookkeeping overhead partly hides it. On a
30-site fermionic chain at maxm=30, one thread against two:

    itensor_version="python"   12.9 s  ->  14.6 s   (1.13x)
    itensor_version=3           2.5 s  ->   2.5 s   (no difference)

The dramatic numbers show up only under *oversubscription*. On the same
host while another job held 10 of the 14 cores, letting MKL spin up its
default 12-14 threads cost:

    itensor_version="python"    9.4 s  ->  28.2 s
    itensor_version=3           1.6 s  ->  26.7 s

That is the case this module is really for: a shared cluster node, or
several dmrgpy runs launched in parallel, where every process's BLAS
believes it owns the machine. The compiled C++ backend degrades worse than
the pure-Python one, which is not a paradox -- it has less Python overhead
to hide the barrier cost behind.

None of this is a defect in either backend, and it is invisible unless you
look for it: nothing errors, the answers are correct, the run is just
several times slower than it should be. Note also that absolute timings
here come from a busy shared host, so treat the ratios as indicative rather
than exact, and measure on your own machine before tuning around them.

== Using it ==

    from dmrgpy import blasthreads

    with blasthreads.limit(1):
        e0 = chain.gs_energy()

    blasthreads.info()      # what the process is currently set to

`limit()` is a context manager and restores the previous setting on exit.
It is deliberately NOT applied automatically: threading does pay off once
the tensors get large enough, and dmrgpy has no business silently changing
a process-wide setting the caller may have chosen deliberately (a script
doing large-bond-dimension work, or its own outer parallelism, wants the
opposite of this).

Runtime control needs `threadpoolctl`, which is not a declared dependency;
without it `limit()` is a no-op that says so once. The environment-variable
route (`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`, or `OPENBLAS_NUM_THREADS=1`)
must be set *before* numpy is imported, so it belongs in the shell rather
than here -- but prefer it when you can: some BLAS libraries fix part of
their threading setup on first use, so a runtime limit applied later is
less effective than starting the process pinned. Measured here, the
env-var route reached ~1.8ms on the 60x60 svd where `limit(1)` from inside
an already-warm process reached ~2.9ms.
"""

import os
import warnings
from contextlib import contextmanager

_WARNED = [False]


def _controller():
    """threadpoolctl's ThreadpoolController, or None if unavailable."""
    try:
        import threadpoolctl
    except ImportError:
        return None
    try:
        return threadpoolctl.ThreadpoolController()
    except Exception:
        return None


@contextmanager
def limit(n=1):
    """Run the enclosed block with at most `n` BLAS/LAPACK/OpenMP threads.

    Restores the previous limits on exit, including if the block raises. A
    no-op (with a one-time warning) when `threadpoolctl` isn't installed --
    see this module's docstring for the environment-variable alternative.
    """
    ctrl = _controller()
    if ctrl is None:
        if not _WARNED[0]:
            _WARNED[0] = True
            warnings.warn(
                "dmrgpy.blasthreads.limit() needs threadpoolctl to change "
                "thread counts at runtime; leaving them alone. Install "
                "threadpoolctl, or set MKL_NUM_THREADS/OMP_NUM_THREADS/"
                "OPENBLAS_NUM_THREADS=1 in the environment before importing "
                "numpy. See dmrgpy/blasthreads.py for why this matters.",
                RuntimeWarning, stacklevel=3)
        yield
        return
    with ctrl.limit(limits=n):
        yield


def info():
    """A list of {user_api, internal_api, num_threads, filepath} dicts for
    every threadpool-backed library loaded in this process; an empty list if
    `threadpoolctl` isn't available."""
    ctrl = _controller()
    if ctrl is None:
        return []
    return ctrl.info()


def env_hint():
    """A one-line shell prefix that would pin every common BLAS to a single
    thread, for printing in a benchmark banner or an error message."""
    return ("MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
            "NUMEXPR_NUM_THREADS=1")


def current_hint():
    """A short human-readable summary of how threads are configured now --
    the relevant environment variables if set, otherwise what threadpoolctl
    reports."""
    envs = {k: v for k, v in os.environ.items()
            if k.endswith("_NUM_THREADS")}
    if envs:
        return ", ".join("{}={}".format(k, v) for k, v in sorted(envs.items()))
    got = info()
    if not got:
        return "unknown (threadpoolctl not installed, no *_NUM_THREADS set)"
    return ", ".join("{}={}".format(d.get("internal_api"), d.get("num_threads"))
                     for d in got)
