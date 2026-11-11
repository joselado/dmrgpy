"""Regressions of the 2026-11 Kondo hole hunt, cluster dispatch: the mode is
validated before order="LL" qualifies it, n_gs is validated before anything
is built, and the implicit-default TypeError names every TASK backend."""
import numpy as np
import pytest
from dmrgpy import spinchain
from dmrgpy.taskdmrg import BACKENDS


def chain(n=3, itensor_version="python"):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


def test_mistyped_mode_with_LL_is_a_value_error_naming_the_modes():
    sc = chain()
    with pytest.raises(ValueError, match='mode must be "ED" or "DMRG"'):
        sc.get_kondo_spectrum(np.linspace(-1, 1, 3), mode="ed", order="LL")
    with pytest.raises(ValueError, match='mode must be "ED" or "DMRG"'):
        sc.get_kondo_spectrum(np.linspace(-1, 1, 3), mode="ed", order=3)


def test_LL_on_DMRG_is_still_not_implemented():
    sc = chain()
    with pytest.raises(NotImplementedError, match="mode=\"ED\" only"):
        sc.get_kondo_spectrum(np.linspace(-1, 1, 3), mode="DMRG", order="LL")


@pytest.mark.parametrize("n_gs", [0, -1, 1.5, True])
def test_n_gs_is_checked_before_anything_is_built(monkeypatch, n_gs):
    from dmrgpy.kondospectrumtk import stepfunctions, taskkondo, edkondo
    def boom(*a, **k): raise AssertionError("built before n_gs was checked")
    monkeypatch.setattr(stepfunctions, "FBuilder", boom)
    monkeypatch.setattr(taskkondo, "TaskKondoSpectrum", boom)
    monkeypatch.setattr(edkondo, "KondoSpectrum", boom)
    sc = chain()
    sc.get_gs = boom # a solve would also be building
    # T>0 builds the thermal FBuilder on the TASK route, so it is the case
    # the check has to come before
    for T in (0., 0.5):
        for kw in ({}, {"submode": "KPM", "delta": 1e-2,
                        "es": np.linspace(-3, 3, 9)}):
            with pytest.raises(ValueError,
                               match="n_gs must be a positive integer"):
                sc.get_kondo_spectrum(np.linspace(-1, 1, 3), mode="DMRG",
                                      T=T, order=3, n_gs=n_gs, **kw)


def test_implicit_default_error_names_every_task_backend():
    sc = chain()
    with pytest.raises(TypeError) as e:
        sc.get_kondo_spectrum(np.linspace(-1, 1, 3), mode="DMRG", T=0., delta=1e-2)
    msg = str(e.value)
    for b in BACKENDS:
        assert repr(b) in msg
    assert "julia_live" in msg
