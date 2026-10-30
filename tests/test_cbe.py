"""Controlled bond expansion (CBE) on itensor_version="python".

dmrg_method="CBE" is one-site DMRG whose bonds are expanded before each
update by the part of H|psi> in the two-site discarded space (Gleis, Li,
von Delft, PRL 130, 246402 (2023)); tevol_method="TDVP_CBE" is one-site
TDVP with the same expansion before each local step (Li, Gleis, von Delft,
arXiv:2208.10972). Both live in pyitensor/cbe.py. What is pinned here:

- the expansion kernel leaves the state unchanged and adds orthonormal
  directions orthogonal to the kept ones, for both selections;
- CBE-DMRG reaches ED and two-site DMRG on a Heisenberg chain, a spinful
  Hubbard chain (fermion signs) and inside a conserved sector;
- from a product state, one-site DMRG without expansion stays a product
  state while CBE reaches the ground state;
- CBE-TDVP follows ED after a quench from a product state, where TDVP_GSE
  stops growing its bonds after tdvp_gse_sweeps steps;
- every other backend refuses both settings, naming "python".

Every DMRG case pins its sweep schedule and its random seed.
"""

import types
import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain, groundstate, spinchain, timedependent
from dmrgpy.pyitensor import cbe
from dmrgpy.pyitensor import chain as chainmod
from dmrgpy.pyitensor.mpsalgebra import _link_at, inner
from dmrgpy.pyitensor.sweeps import Sweeps

warnings.filterwarnings("ignore", message=".*not compiled.*")


def _heisenberg(n, **kw):
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python", **kw)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h


def _bonds(psi):
    return [_link_at(psi, i, i + 1).dim for i in range(1, psi.length())]


def _neel(sc, n):
    """The Neel product state (all bonds 1), as the ground state of a
    staggered field, on the chain's own site indices."""
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    wfed = sc.get_gs(mode="ED")
    assert max(_bonds(wf.cpp_handle)) == 1
    return wf, wfed


# ---------------------------------------------------------------- kernel

@pytest.mark.parametrize("selection", cbe.SELECTIONS)
@pytest.mark.parametrize("side", ["left", "right"])
def test_expansion_keeps_the_state_and_adds_orthogonal_directions(selection, side):
    np.random.seed(3)
    n = 8
    sc, h = _heisenberg(n, maxm=6, nsweeps=3)
    sc.gs_energy()
    s = sc._session
    H, psi = s.H, s.wf0.copy()
    from dmrgpy.pyitensor.tdvp import _gauge_center
    from dmrgpy.pyitensor.dmrg import _all_left_environments, _all_right_environments
    i = 4
    if side == "right":   # center on i, bond (i, i+1) expanded by B rows
        _gauge_center(psi, i)
        L, Lb = _all_left_environments(H, psi)[i - 1]
        R, Rb = _all_right_environments(H, psi)[i + 2]
        before = psi.copy()
        added = cbe.expand_bond_right(psi, H, i, L, Lb, R, Rb, maxdim=6,
                                      delta=1.0, selection=selection)
        iso = psi.A(i + 1).transpose_to([_link_at(psi, i + 1, i)] + [
            ix for ix in psi.A(i + 1).inds if ix != _link_at(psi, i + 1, i)])
        m = np.asarray(iso).reshape(iso.shape[0], -1)
        gram = m.conj() @ m.T
    else:                 # center on i+1, bond (i, i+1) expanded by A columns
        _gauge_center(psi, i + 1)
        L, Lb = _all_left_environments(H, psi)[i - 1]
        R, Rb = _all_right_environments(H, psi)[i + 2]
        before = psi.copy()
        added = cbe.expand_bond_left(psi, H, i, L, Lb, R, Rb, maxdim=6,
                                     delta=1.0, selection=selection)
        lk = _link_at(psi, i, i + 1)
        iso = psi.A(i).transpose_to([ix for ix in psi.A(i).inds if ix != lk] + [lk])
        m = np.asarray(iso).reshape(-1, iso.shape[-1])
        gram = m.conj().T @ m
    assert added > 0
    assert _bonds(psi)[i - 1] == _bonds(before)[i - 1] + added
    # the expanded isometry is an isometry, so the old directions and the
    # new ones are orthonormal and mutually orthogonal
    assert np.max(np.abs(gram - np.eye(gram.shape[0]))) < 1e-12
    # and the state is the same state
    assert abs(inner(psi, before) - inner(before, before)) < 1e-12
    assert abs(inner(psi, psi) - inner(before, before)) < 1e-12


# ---------------------------------------------------------------- CBE-DMRG

def test_cbe_split_is_an_exact_isometry_where_the_gram_route_is_not():
    """The expansion projects on 1 - A A^dag, a projector only for an exact
    isometry A. svd.svd's Gram route returns columns of small singular
    values that are orthonormal only to eps*(s_max/s)^2; cbe.svd always
    takes the full SVD. Before, a 40-site Heisenberg sweep at maxm=50 went
    from -17.5413 up to -15.16."""
    from dmrgpy.pyitensor.index import Index
    from dmrgpy.pyitensor.svd import svd as gram_svd
    from dmrgpy.pyitensor.tensor import ITensor
    rng = np.random.default_rng(0)
    a, b = Index(60, tags="Link"), Index(40, tags="Link")
    q1, _ = np.linalg.qr(rng.standard_normal((60, 40)))
    q2, _ = np.linalg.qr(rng.standard_normal((40, 40)))
    s = np.logspace(0, -6.5, 40)
    T = ITensor((a, b), (q1 * s) @ q2.T)

    def iso_error(U):
        m = np.asarray(U.transpose_to([a] + [ix for ix in U.inds if ix != a]))
        return np.max(np.abs(m.conj().T @ m - np.eye(m.shape[1])))

    assert iso_error(cbe.svd(T, [a], cutoff=0.0)[0]) < 1e-13
    assert iso_error(gram_svd(T, [a], cutoff=0.0)[0]) > 1e-6



def test_cbe_heisenberg_matches_ed_and_two_site():
    np.random.seed(1)
    sc, h = _heisenberg(12, maxm=40, nsweeps=8, dmrg_method="CBE")
    e_cbe = sc.gs_energy()
    e_ed = sc.gs_energy(mode="ED")
    np.random.seed(1)
    sc2, _ = _heisenberg(12, maxm=40, nsweeps=8)
    e_2s = sc2.gs_energy()
    assert e_cbe == pytest.approx(e_ed, abs=1e-6)
    assert e_cbe == pytest.approx(e_2s, abs=1e-6)


def test_cbe_naive_and_shrewd_selection_agree():
    np.random.seed(2)
    sc, h = _heisenberg(10, maxm=24, nsweeps=6)
    groundstate.send_hamiltonian(sc)
    s = sc._session
    energies = {}
    for selection in cbe.SELECTIONS:
        np.random.seed(2)
        psi = s._default_mps(4)
        sweeps = s._make_sweeps(6, 24)
        energies[selection] = cbe.dmrg_cbe(psi, s.H, sweeps, selection=selection)
    e_ed = sc.gs_energy(mode="ED")
    for e in energies.values():
        assert e == pytest.approx(e_ed, abs=1e-6)


def test_cbe_spinful_hubbard_matches_ed_and_two_site():
    """Fermion signs: the Jordan-Wigner strings are in the MPO, and the
    expansion never touches them."""
    n = 4
    res = {}
    for method in ("CBE", "TWO_SITE"):
        np.random.seed(4)
        fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version="python",
                                                  maxm=64, nsweeps=10,
                                                  dmrg_method=method)
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagdn[i] * fc.Cdn[i + 1]
        h = h + h.get_dagger()
        for i in range(n):
            h = h + 4.0 * (fc.Nup[i] - 0.5) * (fc.Ndn[i] - 0.5) + 0.3 * i * fc.Nup[i]
        fc.set_hamiltonian(h)
        res[method] = fc.gs_energy()
    e_ed = fc.gs_energy(mode="ED")
    assert res["CBE"] == pytest.approx(e_ed, abs=1e-6)
    assert res["CBE"] == pytest.approx(res["TWO_SITE"], abs=1e-6)


@pytest.mark.parametrize("two_sz", [0, 2])
def test_cbe_in_a_conserved_sector_matches_ed(two_sz):
    np.random.seed(5)
    sc, h = _heisenberg(10, maxm=40, nsweeps=8, dmrg_method="CBE")
    for i in range(10):
        h = h + 0.3 * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    e = sc.gs_energy()
    e_ed = sc.gs_energy(mode="ED")
    assert e == pytest.approx(e_ed, abs=1e-6)
    assert sc.vev(sum(sc.Sz)).real == pytest.approx(two_sz / 2, abs=1e-8)


def test_one_site_dmrg_without_expansion_is_stuck_and_cbe_is_not():
    """From the Neel product state, one-site DMRG can only rotate each
    site: the bonds stay 1 and the energy stays the product-state one,
    -(n-1)/4. The expansion is what lets the bond grow."""
    n = 10
    np.random.seed(6)
    sc, h = _heisenberg(n, maxm=32, nsweeps=6)
    neel, _ = _neel(sc, n)
    sc.set_hamiltonian(h)
    groundstate.send_hamiltonian(sc)
    s = sc._session
    sweeps = s._make_sweeps(6, 32)
    e_ed = sc.gs_energy(mode="ED")
    stuck = neel.cpp_handle.copy()
    e_stuck = cbe.dmrg_cbe(stuck, s.H, sweeps, expand=False, noisy_start=False)
    assert max(_bonds(stuck)) == 1
    assert e_stuck == pytest.approx(-(n - 1) / 4, abs=1e-8)
    free = neel.cpp_handle.copy()
    e_cbe = cbe.dmrg_cbe(free, s.H, sweeps, noisy_start=False)
    assert max(_bonds(free)) > 1
    assert e_cbe == pytest.approx(e_ed, abs=1e-6)
    # the same through the public warm start
    sc.dmrg_method = "CBE"
    assert sc.gs_energy(wf0=neel) == pytest.approx(e_ed, abs=1e-6)


def test_changing_dmrg_method_gives_a_fresh_solve(monkeypatch):
    """dmrg_method is part of the send-cache key, so a solved chain
    switched to CBE solves again, with CBE."""
    calls = []
    real = chainmod._dmrg_cbe
    monkeypatch.setattr(chainmod, "_dmrg_cbe",
                        lambda *a, **k: calls.append(1) or real(*a, **k))
    np.random.seed(7)
    sc, h = _heisenberg(8, maxm=20, nsweeps=4)
    e1 = sc.gs_energy()
    assert calls == []
    sc.dmrg_method = "CBE"
    e2 = sc.gs_energy()
    assert calls == [1]
    assert e2 == pytest.approx(e1, abs=1e-6)


def test_dmrg_method_is_validated_and_refused_where_it_does_not_exist():
    sc, h = _heisenberg(6, maxm=10, nsweeps=2, dmrg_method="cbe")
    with pytest.raises(ValueError, match="dmrg_method"):
        sc.gs_energy()
    with pytest.raises(TypeError):
        spinchain.Spin_Chain([2] * 4, itensor_version="python", dmrg_methd="CBE")
    for version in (2, 3):
        sc = spinchain.Spin_Chain([2] * 6, itensor_version=version)
        h = sum(sc.Sz[i] * sc.Sz[i + 1] for i in range(5))
        sc.set_hamiltonian(h)
        sc.dmrg_method = "CBE"
        with pytest.raises(NotImplementedError, match='"python"'):
            sc.gs_energy()
        with pytest.raises(NotImplementedError, match='"python"'):
            sc.get_gs()
    # julia_live: the same check, without starting Julia
    stub = types.SimpleNamespace(dmrg_method="CBE", cbe_delta=0.1,
                                 itensor_version="julia_live")
    with pytest.raises(NotImplementedError, match='"python"'):
        groundstate.check_dmrg_method(stub)


def test_excited_states_refuse_cbe():
    sc, h = _heisenberg(6, maxm=10, nsweeps=2, dmrg_method="CBE")
    with pytest.raises(NotImplementedError, match="TWO_SITE"):
        sc.get_excited_states(n=2)


# ---------------------------------------------------------------- CBE-TDVP

def _quench(method, n=10, maxm=32, nt=40, dt=0.1):
    np.random.seed(8)
    sc, h = _heisenberg(n, maxm=maxm, nsweeps=6, tevol_method=method)
    wf, wfed = _neel(sc, n)
    h1 = h + 0.2 * sc.Sz[0]
    sc.set_hamiltonian(h1)
    op = sc.Sz[n // 2] + sc.Sz[n // 2 + 1] * 0.5
    _, z = timedependent.evolve_and_measure(sc, operator=op, wf=wf, nt=nt, dt=dt)
    _, zed = timedependent.evolve_and_measure(sc, mode="ED", operator=op,
                                              wf=wfed, nt=nt, dt=dt)
    return np.real(np.asarray(z)), np.real(np.asarray(zed))


def test_cbe_tdvp_quench_from_a_product_state_matches_ed():
    """Bond dimension 1 at t=0: one-site TDVP alone would never leave the
    product manifold. CBE-TDVP reaches the two-site TDVP accuracy, about
    1e-5 here, set by the time step; TDVP_GSE stops expanding after
    tdvp_gse_sweeps=3 steps and is two orders of magnitude worse."""
    z_cbe, zed = _quench("TDVP_CBE")
    z_2s, _ = _quench("TDVP")
    z_gse, _ = _quench("TDVP_GSE")
    err_cbe = np.max(np.abs(z_cbe - zed))
    err_2s = np.max(np.abs(z_2s - zed))
    err_gse = np.max(np.abs(z_gse - zed))
    assert err_cbe < 5e-5
    assert err_cbe < 2 * err_2s + 1e-6
    assert err_cbe < err_gse / 10


def test_cbe_tdvp_complex_time_step_matches_exact():
    """Chain.tdvp_cbe_step, the one-step route TDZ takes, with a complex
    step, against exp(-i dz H) applied exactly."""
    import scipy.linalg as sla
    np.random.seed(9)
    n = 6
    sc, h = _heisenberg(n, maxm=64, nsweeps=6)
    wf, wfed = _neel(sc, n)
    sc.set_hamiltonian(h)
    groundstate.send_hamiltonian(sc)
    s = sc._session
    dz = 0.05 * np.exp(-0.3j)
    psi = wf.cpp_handle
    for _ in range(4):
        psi = s.tdvp_cbe_step(s.H, psi, dz)
    from dmrgpy.multioperator import MO2matrix
    Hm = np.array(MO2matrix(h, sc.get_ED_obj()).todense())
    v = sla.expm(-1j * 4 * dz * Hm) @ np.asarray(wfed.v).reshape(-1)
    # compare a few observables rather than amplitudes (gauge-free)
    for i in range(n):
        Sz = np.array(MO2matrix(sc.Sz[i], sc.get_ED_obj()).todense())
        exact = (v.conj() @ Sz @ v).real / (v.conj() @ v).real
        mpo = s._mpo(sc.Sz[i].to_terms())
        got = (inner(psi, mpo, psi) / inner(psi, psi)).real
        assert got == pytest.approx(exact, abs=1e-6)


def test_tdvp_cbe_is_refused_on_other_backends():
    for version in (2, 3):
        sc = spinchain.Spin_Chain([2] * 6, itensor_version=version)
        h = sum(sc.Sz[i] * sc.Sz[i + 1] for i in range(5))
        sc.set_hamiltonian(h)
        sc.tevol_method = "TDVP_CBE"
        with pytest.raises(NotImplementedError, match='"python"'):
            timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=3, dt=0.1)
        with pytest.raises(NotImplementedError, match='"python"'):
            timedependent.evolution_DC(sc, name=(sc.Sz[0], sc.Sz[0]), nt=3, dt=0.1)
    stub = types.SimpleNamespace(tevol_method="TDVP_CBE",
                                 itensor_version="julia_live")
    with pytest.raises(NotImplementedError, match='"python"'):
        timedependent.check_tevol_method(stub)


# ---------------------------------------------------------------- GPU rule

def test_cbe_runs_on_the_jax_backend():
    """cbe.py takes its array namespace from pyitensor/backend.py, so both
    routes run under JAX and give the NumPy numbers."""
    pytest.importorskip("jax", reason="the JAX backend needs jax installed")
    from dmrgpy.pyitensor import backend as bk

    def run():
        np.random.seed(10)
        sc, h = _heisenberg(6, maxm=12, nsweeps=4, dmrg_method="CBE",
                            tevol_method="TDVP_CBE")
        e = sc.gs_energy()
        wf, _ = _neel(sc, 6)
        sc.set_hamiltonian(h)
        _, z = timedependent.evolve_and_measure(sc, operator=sc.Sz[2], wf=wf,
                                                nt=4, dt=0.1)
        return e, np.real(np.asarray(z))

    try:
        e_np, z_np = run()
        bk.set_backend("jax")
        e_jx, z_jx = run()
    finally:
        bk.set_backend("numpy")
    assert e_jx == pytest.approx(e_np, abs=1e-9)
    assert np.max(np.abs(z_jx - z_np)) < 1e-9
