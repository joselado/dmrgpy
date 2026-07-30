"""Regression test for a real, pre-existing bug found and fixed while
building the KPM energy-truncation feature (see test_kpm_energy_
truncation*.py): mpscpp2/mpscpp3's chain_session.h check_kpm_moment()
(the divergence guard for a too-tight kpm_scale) used to call ITensor's
Error(...) macro, which prints a message and then calls abort()
unconditionally (itensor/util/error.h) -- it never actually throws,
despite ITError being a real std::runtime_error-derived exception type.
Confirmed directly: hitting this from Python (e.g. a too-narrow
kpm_scale) killed the whole interpreter with SIGABRT instead of raising
a catchable exception, even though check_kpm_moment's own message and
every existing except-RuntimeError-shaped test around it assumed
otherwise.

Fixed by changing just that one call site (not the many other
Error(...) calls in either file, which really do indicate broken
internal preconditions) to `throw ITError(...)` directly -- ITError
derives from std::runtime_error, which pybind11 auto-translates into
Python's RuntimeError, matching the pyitensor backend's own equivalent
check (chain.py's _check_kpm_moment, a plain `raise RuntimeError(...)`)
exactly.

This test exists specifically to catch a regression back to the old
Error(...)-based (uncatchable, process-killing) behavior -- it must
never be "fixed" by loosening it back to expecting a crash.
"""
import numpy as np
import pytest

from dmrgpy import cppext, spinchain


@pytest.mark.parametrize("itensor_version", [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2), reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    pytest.param("python", id="python"),
])
def test_kpm_divergence_raises_catchable_runtime_error(itensor_version):
    """A deliberately too-narrow kpm_scale must raise a plain, catchable
    RuntimeError with the expected message -- not crash the process --
    on every backend that has this guard."""
    n = 4
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    if itensor_version == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.kpm_scale = 0.3  # narrow enough to reliably trigger the guard
    name = (sc.Sz[0], sc.Sz[0])
    with pytest.raises(RuntimeError, match="KPM moments diverging"):
        sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                     es=np.linspace(0.3, 1.0, 71), delta=0.05)
    # If we get here (rather than the whole test process having been
    # killed by an uncaught abort), the fix is doing its job.
