"""Item 15 of docs/kondo_speedup_plan.md: the leading-log flow on the
blocks of total S^z (leading_log_dIdV(sz_blocks=), automatic from
SZ_BLOCKS_FROM states by the author's pick of 2026-11-14). The
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
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
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
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert info["sz_blocks"] and info["sz_rotated"]
    assert np.max(np.abs(got - ref)) < 2e-10*np.max(np.abs(ref))


def test_a_level_crossing_in_a_field_is_rotated():
    """9 spins with J = 1 meV and 0.3 meV along z: a magnon level lies
    exactly 5B below the ferromagnetic multiplet, so four pairs of levels
    of S^z 5 apart are degenerate. Which basis eigh returns inside them is
    rounding luck (the first Triton run got a mixed one, a later local
    run a pure one), so the test turns each pair by 45 degrees, which is
    still an eigenbasis: the tip spin then shows no blocks, and the
    rotation to S^z_tot keeps the energies and finds them."""
    sc = _chain(9)
    ks = edkondo.KondoSpectrum(sc, 0, 0.)
    Szt = np.asarray(MO2matrix(sum(sc.Sz[i] for i in range(9)), sc.get_ED_obj()).todense())
    sz = np.real(np.einsum("ia,ij,ja->a", ks.vs.conj(), Szt, ks.vs))
    W = np.ptp(ks.e)
    pairs = [(a, a + 1) for a in range(len(ks.e) - 1)
             if ks.e[a + 1] - ks.e[a] < 1e-12*W and abs(sz[a] - sz[a + 1]) > 0.5]
    assert len(pairs) == 4
    for a, b in pairs:
        va, vb = ks.vs[:, a].copy(), ks.vs[:, b].copy()
        ks.vs[:, a], ks.vs[:, b] = (va + vb)/np.sqrt(2.), (va - vb)/np.sqrt(2.)
    ks.Sx, ks.Sy, ks.Sz = (ks._eigenbasis_spin_operator(x) for x in ("Sx", "Sy", "Sz"))
    S = leadinglog._spin_ops(ks)
    assert leadinglog._sz_labels(S) is None
    cols = np.nonzero(ks.p > 0)[0]
    U, e, Sr, new, Wc = leadinglog._sz_rotation(ks, cols)
    assert np.max(np.abs(np.sort(e) - ks.e)) < 1e-10*W
    assert leadinglog._sz_labels(Sr) is not None
    assert len(new) == len(cols)
    assert np.max(np.abs(Wc.T @ Wc - np.eye(len(cols)))) < 1e-12


@pytest.mark.parametrize("bz,bx", [(0., 0.3e-3), (0.3e-3, 0.1e-3)])
def test_a_transverse_field_warns_and_runs_dense(bz, bx):
    """H does not conserve S^z_tot: no labels, the rotation is refused, a
    warning, and the dense flow bit for bit."""
    ks = edkondo.KondoSpectrum(_chain(bz=bz, bx=bx), 0, 0.)
    assert leadinglog._sz_labels(leadinglog._spin_ops(ks)) is None
    assert leadinglog._sz_rotation(ks, np.nonzero(ks.p > 0)[0]) is None
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
    with pytest.warns(RuntimeWarning, match="sz_blocks"):
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert not info["sz_blocks"]
    assert np.array_equal(got, ref)


def test_get_kondo_spectrum_takes_sz_blocks():
    sc = _chain()
    _, ref = sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, U=0.2, order="LL",
                                   sz_blocks=False)
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
        ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
        with pytest.warns(RuntimeWarning, match="for the CPU"):
            got = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True)
    finally:
        backend.set_backend(name)
    assert np.array_equal(got, ref)


def _hubbard(cls, n, bz, nf):
    fc = cls(n)
    h = 0
    for i in range(n - 1):
        for s in ("up", "dn"):
            C, Cd = getattr(fc, "C" + s), getattr(fc, "Cdag" + s)
            h = h + 1e-3*(Cd[i]*C[i+1] + Cd[i+1]*C[i])
    for i in range(n):
        h = h + 4e-3*fc.Nup[i]*fc.Ndn[i] - 2e-3*(fc.Nup[i] + fc.Ndn[i]) + bz*fc.Sz[i]
    fc.set_hamiltonian(h)
    if nf is not None:
        fc.set_conserved_sector(Nf=nf)
    return fc


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("nf", [None, 3])
@pytest.mark.parametrize("bz", [0., 0.3e-3])
def test_hubbard_chains_run_on_the_blocks(native, nf, bz):
    """Both spinful classes, with and without an Nf sector (without one the
    tip spin conserves N too, so the blocks fall in several components),
    at B=0 (rotated) and in a field along z."""
    from dmrgpy import fermionchain
    cls = (fermionchain.Spinful_Fermionic_Chain_Native if native
           else fermionchain.Spinful_Fermionic_Chain)
    ks = edkondo.KondoSpectrum(_hubbard(cls, 3, bz, nf), 0, 0.)
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True,
                                                return_info=True)
    assert info["sz_blocks"] and info["sz_rotated"] == (bz == 0.)
    tol = 2e-10 if bz == 0. else 1e-12
    assert np.max(np.abs(got - ref)) < tol*np.max(np.abs(ref))


def test_a_stopped_flow_at_b0_runs_rotated():
    """Rotation, stop and rotation back: D_stop is read in the rotated basis
    (a maximum over entries) and moved by 6e-10 relative when measured."""
    ks = edkondo.KondoSpectrum(_chain(bz=0.), 0, 0.)
    cols = np.nonzero(ks.p > 0)[0]
    kw = dict(pole_threshold=0.06, cols=cols)
    ref, info = leadinglog.tip_vertex(ks, EV, 0.2, -0.05, sz_blocks=False, **kw)
    got, info2 = leadinglog.tip_vertex(ks, EV, 0.2, -0.05, sz_blocks=True, **kw)
    assert info2["sz_blocks"] and info2["sz_rotated"]
    assert info["D_stop"] is not None
    assert info2["D_stop"] == pytest.approx(info["D_stop"], rel=1e-8)
    assert np.max(np.abs(got - ref)) < 1e-10*np.max(np.abs(ref))


def test_a_field_along_y_warns_on_the_complex_route():
    sc = spinchain.Spin_Chain(["S=1/2"]*3)
    h = 0
    for i in range(2):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(3):
        h = h + 0.3e-3*sc.Sz[i] + 0.1e-3*sc.Sy[i]
    sc.set_hamiltonian(h)
    ks = edkondo.KondoSpectrum(sc, 0, 0.)
    assert not leadinglog._block_is_real(leadinglog._spin_ops(ks))
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
    with pytest.warns(RuntimeWarning, match="real route"):
        got = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=True)
    assert np.array_equal(got, ref)


def test_the_default_takes_the_blocks_from_32_states():
    """sz_blocks=None: the blocks from SZ_BLOCKS_FROM = 32 states on, the
    dense flow below, and a silent dense flow when the spectrum has no
    blocks, bit for bit the one of sz_blocks=False."""
    assert leadinglog.SZ_BLOCKS_FROM == 32
    for n, used in ((4, False), (5, True)):
        ks = edkondo.KondoSpectrum(_chain(n), 0, 0.)
        _, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, return_info=True)
        assert info["sz_blocks"] == used
    ks = edkondo.KondoSpectrum(_chain(5, bx=0.1e-3), 0, 0.)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got, info = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, return_info=True)
    assert not info["sz_blocks"]
    ref = leadinglog.leading_log_dIdV(ks, EV, -0.05, U=0.2, sz_blocks=False)
    assert np.array_equal(got, ref)
