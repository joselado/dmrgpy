"""mpsalgebra.check_bicgstab: every BiCGSTAB solve of applyinverse on DMRG
measures its true relative residual ||A x - wf||/||wf|| once, warns above
BICGSTAB_WARN times the requested delta and raises above BICGSTAB_RAISE,
where x is worse than x = 0. Before, it returned whatever the recurrence
held: on itensor_version=3, the default, z - H with z = E0 + 1 + 0.1i on 10
sites at maxm=12 came back with ||x||/||b|| = 3e64 and no warning
(2026-11-21; the lead docs/development.md had left open since 2026-11-15)."""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, multioperator, spinchain

BACKENDS = ["python"] + ([3] if cppext.available(3) else [])


def heisenberg(n, iv, maxm):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=iv)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, 10
    sc.inverse_solver = "bicgstab"
    e0 = sc.gs_energy()
    b = sc.applyoperator(sc.Sz[n//2], sc.get_gs())
    return sc, h, e0, b


@pytest.mark.parametrize("iv", BACKENDS)
def test_diverged_solve_raises(iv):
    """z - H under truncation (maxm=12 on 10 sites): the recurrence leaves
    the solution, and after 20 iterations x is worse than x = 0"""
    sc, h, e0, b = heisenberg(10, iv, 12)
    sc.cvm_nit = 20
    A = (e0 + 1. + 0.1j)*multioperator.identity() - h
    with pytest.raises(RuntimeError, match="worse than x = 0"):
        sc.applyinverse(A, b)


@pytest.mark.parametrize("iv", BACKENDS)
def test_converged_solve_is_quiet(iv):
    """Every state held (maxm=30 on 8 sites): the true residual is at the
    tolerance (measured 1.1e-5 at delta = 1e-5), no warning"""
    sc, h, e0, b = heisenberg(8, iv, 30)
    A = (e0 + 1. + 0.1j)*multioperator.identity() - h
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        x = sc.applyinverse(A, b)
    assert np.isfinite(abs(x.dot(x)))


@pytest.mark.parametrize("iv", BACKENDS)
def test_floor_of_the_truncation_warns(iv):
    """H - e below the spectrum is positive definite, but at maxm=12 on 10
    sites the recurrence stops on its own truncated residual while the
    true one is 1e-2: a warning, not a raise"""
    sc, h, e0, b = heisenberg(10, iv, 12)
    A = h - (e0 - 0.3)*multioperator.identity()
    with pytest.warns(RuntimeWarning, match="BiCGSTAB"):
        sc.applyinverse(A, b)
