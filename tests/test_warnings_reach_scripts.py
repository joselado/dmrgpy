"""Importing dmrgpy must not silence warnings in the caller's script.

pytest resets the warning filters for each test, so the check runs in a
fresh interpreter, the way a plain script imports dmrgpy."""
import os
import subprocess
import sys

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")

SCRIPT = r"""
import sys, warnings
sys.path.insert(0, sys.argv[1])
import dmrgpy.spinchain, dmrgpy.ysr
from dmrgpy.algebra import kpmextrapolate
catch_all = [f for f in warnings.filters if f[0] == "ignore" and f[1] is None
             and f[2] is Warning and f[3] is None]
assert not catch_all, catch_all
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("default")
    warnings.warn("reaches the script", UserWarning)
assert len(w) == 1
print("ok")
"""


def test_import_installs_no_catch_all_ignore_filter():
    out = subprocess.run([sys.executable, "-c", SCRIPT, SRC],
                         capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().endswith("ok")


def test_power_extrapolation_stays_quiet():
    """The fits inside the moment extrapolation keep their own silence."""
    import warnings
    import numpy as np
    from dmrgpy.algebra import kpmextrapolate
    mus = np.cos(0.3*np.arange(40))/(1.0 + np.arange(40))
    mus[5] = 0.0 # the log of a zero moment
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        f, finv = kpmextrapolate.fit_power_transform(mus)
        assert np.all(np.isfinite(finv(f(mus[:10]))))
