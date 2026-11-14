"""Item 15 of docs/kondo_speedup_plan.md: the leading-log flow on the
blocks of total S^z (leading_log_dIdV(sz_blocks=True), opt-in). The
labels come from the tip spin, the blocked flow gives the dense one to
roundoff; a degenerate level that mixes S^z sectors (B=0, a crossing in a
field) is rotated to the eigenbasis of S^z_tot and gives the dense flow
to its rtol; a transverse field warns and runs the dense flow bit for
bit."""

import warnings

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.multioperator import MO2matrix
from dmrgpy.kondospectrumtk import edkondo, leadinglog


def _chain(n=5, bz=0.3e-3, bx=0., spin="S=1/2"):
    sc = spinchain.Spin_Chain([spin]*n)
    h = 0
    for i in range(n - 1):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + bz*sc.Sz[i] + bx*sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc


EV = np.linspace(-3e-3, 3e-3, 21) + 1e-5


@pytest.mark.parametrize("n,spin", [(5, "S=1/2"), (3, "S=1")])
def test_the_labels_are_sz_tot_up_to_a_constant(n, spin):
    sc = _chain(n, spin=spin)
    ks = edkondo.KondoSpectrum(sc, 0, 0.)
    M, comp = leadinglog._sz_labels(leadinglog._spin_ops(ks))
    Szt = np.asarray(MO2matrix(sum(sc.Sz[i] for i in range(n)), sc.get_ED_obj()).todense())
    sz = np.real(np.einsum("ia,ij,ja->a", ks.vs.conj(), Szt, ks.vs))
    assert len(set(comp)) == 1
    assert np.ptp(M - sz) < 1e-9


def test_the_blocked_bracket_is_the_dense_one():
    """On random vertices with the block structure of labels in two
    components, which no product joins."""
    rng = np.random.default_rng(1)
    d = 12
    M = rng.integers(-2, 3, d)
    comp = np.where(np.arange(d) < 8, 0, 8)
    blocks = leadinglog._SzBlocks(M, comp)
    two_s = np.r_[np.ones(d), -np.ones(d)]
    L, C = np.r_[M, M], np.r_[comp, comp]
    keep = (((2*L + two_s)[:, None] == (2*L + two_s)[None, :])
            & (C[:, None] == C[None, :]))
    X, Y = rng.normal(size=(2, 2*d, 2*d))*keep
    dense = leadinglog._bracket(X, Y)
    assert np.max(np.abs(dense[~keep])) == 0.
    assert np.max(np.abs(blocks.bracket(X, Y) - dense)) < 1e-14


@pytest.mark.parametrize("T,tables", [(0., True), (0., False), (2.0, True)])
def test_sz_blocks_give_the_dense_flow(T, tables, monkeypatch):
    """At T=0 with the compact tables of the hole sum and with the kernel
    on the fly; at T>0 the products of G and A_on only."""
    if not tables:
        monkeypatch.setattr(leadinglog, "HOLE_TABLE_BYTES", 0)
    ks = edkondo.KondoSpectrum(_chain(), 0, T)
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert info["sz_blocks"] and not info["sz_rotated"]
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))


@pytest.mark.parametrize("n,T", [(5, 0.), (5, 2.0), (4, 0.)])
def test_b0_multiplets_run_rotated_to_sz_tot(n, T):
    """At B=0 every multiplet mixes S^z sectors in the eigenbasis, so the
    flow runs in the eigenbasis of S^z_tot and the tip vertex is rotated
    back: the dense flow to its rtol (measured 2.6e-12 to 4.1e-11 of the
    peak), as for any change of basis inside the degenerate levels."""
    ks = edkondo.KondoSpectrum(_chain(n, bz=0.), 0, T)
    assert leadinglog._sz_labels(leadinglog._spin_ops(ks)) is None
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert info["sz_blocks"] and info["sz_rotated"]
    assert np.max(np.abs(got - ref)) < 2e-10*np.max(np.abs(ref))


def test_a_level_crossing_in_a_field_is_rotated():
    """9 spins with J = 1 meV and 0.3 meV along z: a magnon level lies
    exactly 5B below the ferromagnetic multiplet, four pairs of levels of
    S^z 5 apart are degenerate, and the tip spin shows no blocks. The
    rotation to S^z_tot keeps the energies and finds them."""
    ks = edkondo.KondoSpectrum(_chain(9), 0, 0.)
    S = leadinglog._spin_ops(ks)
    assert leadinglog._sz_labels(S) is None
    cols = np.nonzero(ks.p > 0)[0]
    U, e, Sr, new, Wc = leadinglog._sz_rotation(ks, cols)
    assert np.max(np.abs(np.sort(e) - ks.e)) < 1e-10*np.ptp(ks.e)
    assert leadinglog._sz_labels(Sr) is not None
    assert len(new) == len(cols)
    assert np.max(np.abs(Wc.T @ Wc - np.eye(len(cols)))) < 1e-12


def test_a_flow_that_stops_on_the_blocks():
    """The quadrature below D_stop runs dense on the state of the blocks.
    dI/dV of a stopped flow is NaN at every bias on this spectrum, so the
    tip vertex itself is compared."""
    ks = edkondo.KondoSpectrum(_chain(), 0, 0.)
    cols = np.nonzero(ks.p > 0)[0]
    kw = dict(pole_threshold=0.06, cols=cols)
    ref, info = leadinglog.tip_vertex(ks, EV, 0.2, -0.05, **kw)
    got, info2 = leadinglog.tip_vertex(ks, EV, 0.2, -0.05, sz_blocks=True, **kw)
    assert info["D_stop"] is not None
    assert info2["D_stop"] == pytest.approx(info["D_stop"], rel=1e-9)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))


@pytest.mark.parametrize("bz,bx", [(0., 0.3e-3), (0.3e-3, 0.1e-3)])
def test_a_transverse_field_warns_and_runs_dense(bz, bx):
    """H does not conserve S^z_tot: no labels, the rotation is refused, a
    warning, and the dense flow bit for bit."""
    ks = edkondo.KondoSpectrum(_chain(bz=bz, bx=bx), 0, 0.)
    assert leadinglog._sz_labels(leadinglog._spin_ops(ks)) is None
    assert leadinglog._sz_rotation(ks, np.nonzero(ks.p > 0)[0]) is None
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2)
    with pytest.warns(RuntimeWarning, match="sz_blocks"):
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert not info["sz_blocks"]
    assert np.array_equal(got, ref)


def test_get_kondo_spectrum_takes_sz_blocks():
    sc = _chain()
    _, ref = sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, U=0.2, order="LL")
    _, got = sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, U=0.2, order="LL",
                                   sz_blocks=True)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))
    with pytest.raises(TypeError, match="sz_blocks"):
        sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, U=0.2, order=3,
                              sz_blocks=True)


def test_the_blocks_do_less_work():
    """On 6 spins the products of the flow are 7% of the dense ones, and
    the exchange sum of the tip reads the final states within one S^z of
    the ground state only."""
    ks = edkondo.KondoSpectrum(_chain(6), 0, 0.)
    S = leadinglog._spin_ops(ks)
    blocks = leadinglog._SzBlocks(*leadinglog._sz_labels(S))
    assert blocks.work() < 0.08
    cols = np.nonzero(ks.p > 0)[0]
    sums = blocks.hole_sum(cols, EV, ks.e, 20e-3, 5e-6)
    assert len(sums.rows) < 0.8*len(ks.e)
    assert sums.tables is not None


def test_sz_blocks_on_the_device_warn_and_run_the_device_flow():
    """The blocks are for the CPU; on the JAX backend the flow is the one of
    the device, with a warning."""
    pytest.importorskip("jax")
    from dmrgpy.pyitensor import backend
    ks = edkondo.KondoSpectrum(_chain(4), 0, 0.)
    name = backend.backend_name()
    backend.set_backend("jax")
    try:
        ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2)
        with pytest.warns(RuntimeWarning, match="for the CPU"):
            got = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True)
    finally:
        backend.set_backend(name)
    assert np.array_equal(got, ref)
