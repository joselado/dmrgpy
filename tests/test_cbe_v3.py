"""Controlled bond expansion (CBE) on itensor_version=3.

mpscpp3/cbe.h is a port of pyitensor/cbe.py (tests/test_cbe.py pins the
kernel itself); these tests pin the port through the public API:

- dmrg_method="CBE" reaches ED and two-site DMRG on a Heisenberg chain and
  a spinful Hubbard chain (fermion signs), and agrees with "python";
- from the Neel product state, a warm start whose bonds are all 1, CBE
  grows the bonds and reaches the ground state;
- tevol_method="TDVP_CBE" follows ED after a quench from a product state
  to the accuracy of two-site TDVP, where TDVP_GSE stops growing its bonds;
- the one-step route TDZ takes, with a complex step, matches exp(-i dz H);
- a conserved sector is refused on v3, since cbe.h works on dense indices.
- for a non-Hermitian generator (a Lindbladian) the local steps are
  Arnoldi, as on "python", and the two backends agree step for step.

Every DMRG case pins its sweep schedule.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, groundstate, mps, spinchain, timedependent

warnings.filterwarnings("ignore", message=".*not compiled.*")

pytestmark = pytest.mark.skipif(not cppext.available(3),
                                reason="the itensor_version=3 extension is not compiled")


def _heisenberg(n, version=3, **kw):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version, **kw)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h


def _neel(sc, n):
    """The Neel product state, as the ground state of a staggered field."""
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    wfed = sc.get_gs(mode="ED")
    assert max(wf.get_entropy(b) for b in range(n - 1)) < 1e-10
    return wf, wfed


def test_cbe_heisenberg_matches_ed_two_site_and_python():
    e = {}
    for version, method in ((3, "CBE"), (3, "TWO_SITE"), ("python", "CBE")):
        np.random.seed(1)
        sc, _ = _heisenberg(12, version=version, maxm=40, nsweeps=8,
                            dmrg_method=method)
        e[(version, method)] = sc.gs_energy()
    e_ed = sc.gs_energy(mode="ED")
    for value in e.values():
        assert value == pytest.approx(e_ed, abs=1e-6)


def test_cbe_spinful_hubbard_matches_ed():
    """Fermion signs: the Jordan-Wigner strings are in the MPO, and the
    expansion never touches them."""
    n = 4
    res = {}
    for method in ("CBE", "TWO_SITE"):
        fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version=3,
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


def test_cbe_from_a_product_state_grows_the_bonds():
    """One-site DMRG alone would stay in the product manifold, at the
    Neel energy -(n-1)/4; the expansion is what lets the bonds grow."""
    n = 10
    sc, h = _heisenberg(n, maxm=32, nsweeps=6)
    neel, _ = _neel(sc, n)
    sc.set_hamiltonian(h)
    sc.dmrg_method = "CBE"
    e = sc.gs_energy(wf0=neel)
    e_ed = sc.gs_energy(mode="ED")
    assert e == pytest.approx(e_ed, abs=1e-6)
    assert sc.get_gs().get_entropy(n // 2 - 1) > 0.5


def test_cbe_small_units_are_unit_scaled():
    """The solve goes through solver_hamiltonian(), like two-site dmrg()."""
    sc, h = _heisenberg(8, maxm=20, nsweeps=6, dmrg_method="CBE")
    sc.set_hamiltonian(1e-6 * h)
    assert sc.gs_energy() == pytest.approx(sc.gs_energy(mode="ED"), rel=1e-7)


def test_cbe_refuses_a_conserved_sector_on_v3():
    sc, h = _heisenberg(8, maxm=20, nsweeps=4, dmrg_method="CBE")
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(ValueError, match="CBE"):
        sc.gs_energy()


def _quench(method, n=10, maxm=32, nt=40, dt=0.1):
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
    import scipy.linalg as sla
    from dmrgpy.multioperator import MO2matrix
    n = 6
    sc, h = _heisenberg(n, maxm=64, nsweeps=6)
    wf, wfed = _neel(sc, n)
    sc.set_hamiltonian(h)
    groundstate.send_hamiltonian(sc)
    H = sc.toMPO(h)
    dz = 0.05 * np.exp(-0.3j)
    psi = wf.cpp_handle
    for _ in range(4):
        psi = sc._session.tdvp_cbe_step(H.cpp_handle, psi, dz)
    psi = mps.MPS(sc, cpp_handle=psi)
    Hm = np.array(MO2matrix(h, sc.get_ED_obj()).todense())
    v = sla.expm(-1j * 4 * dz * Hm) @ np.asarray(wfed.v).reshape(-1)
    for i in range(n):
        Sz = np.array(MO2matrix(sc.Sz[i], sc.get_ED_obj()).todense())
        exact = (v.conj() @ Sz @ v).real / (v.conj() @ v).real
        got = (psi.dot(sc.Sz[i] * psi) / psi.dot(psi)).real
        assert got == pytest.approx(exact, abs=1e-6)


def test_mpo_hermiticity_test_on_v3():
    sc, h = _heisenberg(5, maxm=10, nsweeps=2)
    groundstate.send_hamiltonian(sc)
    s = sc._session
    assert s.mpo_is_hermitian(sc.toMPO(h).cpp_handle)
    assert s.mpo_is_hermitian(sc.toMPO(1e-6 * h).cpp_handle)
    assert not s.mpo_is_hermitian(sc.toMPO(h + 1e-3j * sc.Sz[2]).cpp_handle)
    assert not s.mpo_is_hermitian(sc.toMPO(1j * sc.Sx[0] * sc.Sy[1]).cpp_handle)
    assert s.mpo_is_hermitian(sc.toMPO(sc.Sx[0] * sc.Sy[1]).cpp_handle)


def _lindblad_cbe_steps(version, dt, nst):
    """nst CBE-TDVP steps exp(dt*L) of a driven XXZ Lindbladian on the
    doubled chain, from a product state; returns the doubled chain, the
    evolved MPS, the exact vector and L."""
    import scipy.linalg as sla
    from dmrgpy.lindblad import Lindblad_Chain
    from dmrgpy.multioperator import MO2matrix
    n = 3
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    sc.maxm = 64
    sc.set_hamiltonian(sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                           + 0.5 * sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1)))
    jumps = [0.6 * (sc.Sx[0] + 1j * sc.Sy[0]), 0.6 * (sc.Sx[n - 1] - 1j * sc.Sy[n - 1])]
    jumps += [0.4 * sc.Sz[i] for i in range(n)]
    lc = Lindblad_Chain(sc, jumps)
    d = lc.doubled
    d.maxm, d.nsweeps = 64, 6
    # a product state: the ground state of single-site fields only
    d.set_hamiltonian(sum((-1) ** (i // 2) * d.Sz[i] + 0.3 * d.Sx[i]
                          for i in range(2 * n)))
    wf, wfed = d.get_gs(), d.get_gs(mode="ED")
    d.set_hamiltonian(1j * lc.L)  # exp(-i dt (iL)) = exp(dt L)
    groundstate.send_hamiltonian(d)
    G = d.toMPO(1j * lc.L)
    assert not d._session.mpo_is_hermitian(G.cpp_handle) if version == 3 else True
    psi = wf.cpp_handle
    for _ in range(nst):
        psi = d._session.tdvp_cbe_step(G.cpp_handle, psi, dt)
    Lm = np.array(MO2matrix(lc.L, d.get_ED_obj()).todense())
    v = sla.expm(nst * dt * Lm) @ np.asarray(wfed.v).reshape(-1)
    return d, mps.MPS(d, cpp_handle=psi), v


def _sz_profile(d, psi):
    return np.array([psi.dot(d.Sz[i] * psi) / psi.dot(psi)
                     for i in range(d.ns)])


def test_cbe_tdvp_on_a_lindbladian_follows_python_and_exact():
    """Arnoldi local steps on v3 reproduce "python" step for step (the
    Lanczos applyExp route is 9e-10 away here), and both approach exp(t L)
    linearly in dt, 0.075 at dt=0.5 and 0.012 at dt=0.1 on <Sz_i>."""
    from dmrgpy.multioperator import MO2matrix
    d3, psi3, v = _lindblad_cbe_steps(3, 0.1, 20)
    dp, psip, _ = _lindblad_cbe_steps("python", 0.1, 20)
    z3, zp = _sz_profile(d3, psi3), _sz_profile(dp, psip)
    assert np.max(np.abs(z3 - zp)) < 1e-11
    assert abs(psi3.dot(psi3) - psip.dot(psip)) < 1e-11
    ed = d3.get_ED_obj()
    zed = np.array([(v.conj() @ np.array(MO2matrix(d3.Sz[i], ed).todense()) @ v)
                    / (v.conj() @ v) for i in range(d3.ns)])
    assert np.max(np.abs(z3 - zed)) < 2e-2


def test_an_extension_built_before_cbe_refuses_cbe_and_still_runs_two_site():
    """A v3 .so built before cbe.h has no set_dmrg_method: TWO_SITE runs
    as before, and CBE raises naming the rebuild."""
    import types
    sc, h = _heisenberg(6, maxm=10, nsweeps=2)
    real = sc._session
    sc._session = types.SimpleNamespace()
    sc.dmrg_method = "CBE"
    with pytest.raises(NotImplementedError, match="rebuild"):
        groundstate.check_dmrg_method(sc)
    sc.tevol_method = "TDVP_CBE"
    with pytest.raises(NotImplementedError, match="rebuild"):
        timedependent.check_tevol_method(sc)
    sc._session = real
