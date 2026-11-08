"""The four-point tensor <Cdag_i C_j Cdag_k C_l> with ctmode="full" on
julia_live: one exact MPO per tuple and inner(psi',M,psi), in Julia
(mpsjulialive/gs_extra.jl). It fell through to ctmode="explicit" there.
Checked element by element, repeated indices included, against julia_live's
"explicit", "python"'s "full" and ED."""

import contextlib
import io
import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy.entropytk import correlationentropy as ce

from _helpers import julia_live_param


N = 6


def _chain(**kw):
    fc = fermionchain.Fermionic_Chain(N, **kw)
    h = 0
    for i in range(N-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
        h = h + 0.7*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    for i in range(N):
        h = h + 0.05*(i-2.5)*fc.N[i] # breaks the reflection symmetry
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 64, 12
    return fc


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


@pytest.mark.parametrize("version", [julia_live_param()])
def test_julia_full_four_point_tensor(version):
    np.random.seed(0)
    fj = _chain(itensor_version=version)
    wj = _quiet(fj.get_gs)
    # "sweep" is the default on a julia_live spinless chain, so "full" is
    # asked for by name below
    assert ce._four_correlation_tensor_default_ctmode(wj) == "sweep"
    full = _quiet(lambda: wj.get_four_correlation_tensor(ctmode="full"))
    slow = _quiet(lambda: wj.get_four_correlation_tensor(ctmode="full",
                                                         accelerate=False))
    expl = _quiet(lambda: wj.get_four_correlation_tensor(ctmode="explicit"))
    np.random.seed(0)
    fp = _chain(itensor_version="python")
    wp = _quiet(fp.get_gs)
    pyfull = _quiet(lambda: wp.get_four_correlation_tensor(ctmode="full"))
    fe = _chain()
    we = _quiet(lambda: fe.get_gs(mode="ED"))
    ed = _quiet(lambda: we.get_four_correlation_tensor())
    assert full.shape == (N, N, N, N)
    assert np.max(np.abs(ed)) > 0.1
    for ref in (slow, expl, pyfull, ed):
        assert np.max(np.abs(full - ref)) < 1e-8
    # the repeated-index entries, by name: <N_i N_k> = <Cdag_i C_i Cdag_k C_k>
    for i in range(N):
        for k in range(N):
            nn = np.real(_quiet(lambda: we.dot(
                (fe.N[i]*fe.N[k])*we)))
            assert full[i, i, k, k] == pytest.approx(nn, abs=1e-8)
