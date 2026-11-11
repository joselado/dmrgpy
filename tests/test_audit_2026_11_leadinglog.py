"""Audit 2026-11 (leading-log cluster): the flow in real block arithmetic,
the tabulated exchange sum of the tip vertex, and an impurity with no spin.
kondospectrumtk/leadinglog.py, docs/new_methods/kondo_leading_log/."""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from dmrgpy.kondospectrumtk import leadinglog
from test_kondo_spectrum_tmatrix import MODELS, EVS, MUB, S12, _ks
from test_kondo_leading_log import _TWO_PASS, _spin_half_chain, W0, G0


def _dIdV(ks, **kw):
    return leadinglog.leading_log_dIdV(ks, EVS[::5], -0.05, U=0.25, **kw)


def test_the_block_form_is_real_for_a_real_hamiltonian_only():
    for name in MODELS:
        ks = _ks(*MODELS[name])
        assert leadinglog._block_is_real(leadinglog._spin_ops(ks)), name
    # a field along y makes the eigenbasis complex: Sx and Sz are not real
    ks = _ks(2*MUB*10*S12[1], S12)
    assert not leadinglog._block_is_real(leadinglog._spin_ops(ks))


@pytest.mark.parametrize("name", list(MODELS))
def test_real_and_complex_routes_agree(name, monkeypatch):
    """The real route is the complex one in other arithmetic: the same
    equations, so dI/dV agrees inside the integration tolerance (the step
    sequence differs). The complex route is the code from before."""
    ks = _ks(*MODELS[name])
    real = _dIdV(ks)
    monkeypatch.setattr(leadinglog, "_block_is_real", lambda S, tol=0.: False)
    cplx = _dIdV(ks)
    assert real == pytest.approx(cplx, rel=1e-9)
    assert real == pytest.approx(_TWO_PASS[name], rel=1e-9)


def test_a_field_along_y_takes_the_complex_route_and_gives_the_rotated_values():
    """B along y is B along x turned about z, and the sample couples
    isotropically, so dI/dV of the free S=1/2 is the one of Bx=10 T. The
    eigenbasis is complex here, so this goes through the complex route."""
    ks = _ks(2*MUB*10*S12[1], S12)
    assert not leadinglog._block_is_real(leadinglog._spin_ops(ks))
    assert _dIdV(ks) == pytest.approx(_TWO_PASS["S=1/2, Bx=10 T"], rel=1e-9)


def test_the_real_route_gives_the_flowed_vertex_of_the_complex_one(monkeypatch):
    ks = _ks(*MODELS["S=1, D=-1 meV, B=(1,0,2) T"])
    oms = np.array([2.2e-3, -0.7e-3, 0.05e-3])
    A, info = leadinglog.tip_vertex(ks, oms, 0.2, -0.05, omega0=W0, Gamma0=G0)
    monkeypatch.setattr(leadinglog, "_block_is_real", lambda S, tol=0.: False)
    Ac, infoc = leadinglog.tip_vertex(ks, oms, 0.2, -0.05, omega0=W0, Gamma0=G0)
    assert np.max(np.abs(A - Ac)) < 1e-9*np.max(np.abs(Ac))
    assert info["max_coupling"] == pytest.approx(infoc["max_coupling"], rel=1e-9)


def test_gamma0_must_be_positive_on_the_real_route():
    ks = _ks(*MODELS["S=1/2, Bx=10 T"])
    assert leadinglog._block_is_real(leadinglog._spin_ops(ks))
    with pytest.raises(ValueError):
        leadinglog.tip_vertex(ks, np.array([1e-3]), 0., -0.05, Gamma0=0.)
    with pytest.raises(ValueError):
        leadinglog.tip_vertex(ks, np.array([1e-3]), 0., -0.05, Gamma0=-1e-6)


def test_a_stopped_flow_on_the_real_route_is_the_complex_one(monkeypatch):
    """Jrho_s=-0.2 on the S=1 model stops the flow: the stop cutoff, the
    closed-form tail and the NaN mask follow the complex route."""
    ks = _ks(*MODELS["S=1, D=-1 meV, B=(1,0,2) T"])
    oms = np.array([2.2e-3, -0.7e-3, 0.05e-3])
    A, info = leadinglog.tip_vertex(ks, oms, 0., -0.2, omega0=W0, Gamma0=G0)
    monkeypatch.setattr(leadinglog, "_block_is_real", lambda S, tol=0.: False)
    Ac, infoc = leadinglog.tip_vertex(ks, oms, 0., -0.2, omega0=W0, Gamma0=G0)
    assert info["D_stop"] is not None
    assert info["D_stop"] == pytest.approx(infoc["D_stop"], rel=1e-9)
    assert np.max(np.abs(A - Ac)) < 1e-8*np.max(np.abs(Ac))


# ---------------------------------------------------------------------
# the tabulated exchange sum

def test_the_table_sum_is_the_sum_on_the_fly():
    rng = np.random.default_rng(1)
    d, n, nl = 7, 9, 12
    e = np.sort(rng.normal(size=d))*1e-3
    w = rng.normal(size=n)*2e-3
    Q = np.ascontiguousarray(rng.normal(size=(d, nl, d)))
    for D in (3e-2, 1e-4, 4e-6):
        ref = np.empty((n, d, nl)); got = np.empty((n, d, nl))
        leadinglog._hole_sum(Q, w, e, D, -1., W0, G0, ref)
        R, Wt = leadinglog._hole_tables(w, e, W0, G0)
        leadinglog._hole_sum_table(Q, R, Wt, D, got)
        assert np.max(np.abs(got - ref)) < 1e-14*np.max(np.abs(ref))


def test_the_table_is_off_above_the_cap_and_the_numbers_do_not_move(monkeypatch):
    ks = _spin_half_chain(4)
    w = np.array([1e-3, -2e-3, 3e-4]); e = np.asarray(ks.e)
    assert leadinglog._hole_tables(w, e, W0, G0) is not None
    monkeypatch.setattr(leadinglog, "HOLE_TABLE_BYTES", 16*len(w)*ks.dim**2 - 1)
    assert leadinglog._hole_tables(w, e, W0, G0) is None
    off = leadinglog.leading_log_dIdV(ks, EVS[::4], -0.05, U=0.25)
    monkeypatch.setattr(leadinglog, "HOLE_TABLE_BYTES", 256e6)
    on = leadinglog.leading_log_dIdV(ks, EVS[::4], -0.05, U=0.25)
    assert on == pytest.approx(off, rel=1e-12)


# ---------------------------------------------------------------------
# an impurity with no spin

def _no_spin_ks(occupied_only=False):
    d = 4
    Z = np.zeros((d, d), complex)
    Sx, Sy, Sz = Z.copy(), Z.copy(), Z.copy()
    if occupied_only:
        # spin on the two excited states, none on the occupied ground state
        Sx[1, 2] = Sx[2, 1] = .5
        Sy[1, 2], Sy[2, 1] = -.5j, .5j
        Sz[1, 1], Sz[2, 2] = .5, -.5
    return SimpleNamespace(e=np.array([0., 1e-3, 1.5e-3, 3e-3]),
                           p=np.array([1., 0, 0, 0]), dim=d, T=0., kB=8.617e-5,
                           Sx=Sx, Sy=Sy, Sz=Sz)


@pytest.mark.parametrize("occupied_only", [False, True])
def test_a_state_with_no_spin_gives_the_potential_scattering(occupied_only):
    """No spin on the occupied states: dI/dV is 2 pi 4 U^2, the second
    order, with no NaN and no warning, and the tip coupling is zero."""
    ks, U = _no_spin_ks(occupied_only), 0.25
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        dIdV, info = leadinglog.leading_log_dIdV(
            ks, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=U, return_info=True)
    assert dIdV == pytest.approx(2*np.pi*4*U**2, rel=1e-9)
    assert np.all(info["coupling"] == 0.) and info["D_stop"] is None


def test_an_empty_orbital_in_a_sector_gives_the_second_order():
    from dmrgpy import fermionchain
    fc = fermionchain.Spinful_Fermionic_Chain(2)
    fc.set_hamiltonian(-2e-3*sum(fc.Cdag[s]*fc.C[2+s] + fc.Cdag[2+s]*fc.C[s]
                                 for s in (0, 1)))
    fc.set_conserved_sector(Nf=0)
    eVs = np.array([-2e-3, 5e-4, 2e-3])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ll = fc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0.,
                                   order="LL", mode="ED")[1]
    assert ll == pytest.approx(2*np.pi*4*0.25**2, rel=1e-9)
