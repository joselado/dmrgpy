"""Regression tests of the 2026-11-21 hole hunt, cluster esr-krylov
(docs/audit_2026_11_21_hole_hunt.md, findings 1 and 2 and the
submode="KRYLOV" half of finding 3).

The basis of get_esr_spectrum(submode="KRYLOV") was the populated
eigenstates and the Krylov space of V|m>, closed under H and V but not
under the Korringa operators S^a_k, nor under a second application of V.
With the tip moment along the field (H and V both conserve Sz_tot) the T+
and T- of the S-T0 line of a dimer were missing, and the line came out
-5.0 times the eigenstate one (exactly 0 at zero tip field); with a tilted
tip an empty level two Sz steps away was missing on 4 spins (2.1e-2 of the
peak); and the |+1> of a saturated spin 1 was never in the basis (36% low).

What is pinned: on ED the Krylov route agrees with the eigenstate route
(submode=None) on the three cases, to the 1e-8 of tests/test_esr_stm.py;
its closure guard refuses, naming submode=None, when the closing stage is
taken away; and on DMRG, which has no closing stage, the route refuses the
chains whose partners lie in a sector it cannot reach, and agrees with ED
on the 4-spin chain whose partners it holds, to the 5e-5 of its own test.
"""
import numpy as np
import pytest

from dmrgpy import cppext, esrstm, spinchain
from dmrgpy.algebra import algebra
from dmrgpy.esrtk import spectrum

H = 4.135667696e-6 # eV/GHz
KW = dict(site=0, T=0.03, T1=100., T2=50., f_ref=20., rabi=1e-4)
TIP_Z = [0., 0., 1.]
TIP_TILTED = [0.5, 0., 0.866]
BACKENDS = ["python"] + ([3] if cppext.available(3) else [])


def chain(n, J, B, Bt, iv=None):
    """n spins 1/2, energies in GHz: J S_i.S_i+1 + B Sz_tot + Bt Sz_0, so
    that H conserves Sz_tot"""
    kw = {} if iv is None else dict(itensor_version=iv)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    h = 0
    for i in range(n-1):
        h = h + H*J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                     + sc.Sz[i]*sc.Sz[i+1])
    h = h + H*B*sum(sc.Sz) + H*Bt*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 1e-14
    return sc


def levels(sc, V):
    """ED energies in GHz above the ground state, <Sz_tot> and |<n|V|0>|^2"""
    ed = sc.get_ED_obj()
    Hm = np.asarray(algebra.todense(ed.get_hamiltonian()))
    Szt = np.asarray(algebra.todense(ed.get_operator(sum(sc.Sz))))
    Vm = np.asarray(algebra.todense(ed.get_operator(V)))
    w, v = np.linalg.eigh(Hm)
    mz = np.real(np.einsum("in,ij,jn->n", v.conj(), Szt, v))
    return (w - w[0])/H, mz, np.abs(v.conj().T @ Vm @ v[:, 0])**2


def dimer_window(Bt):
    """The S-T0 line of the dimer at J = 30 GHz, sqrt(J^2 + Bt^2) GHz, give
    or take 0.15 GHz"""
    f0 = np.sqrt(30.**2 + Bt**2)
    return np.linspace(f0 - 0.15, f0 + 0.15, 601)


def tilted_line():
    """The line of finding 2: 4 spins at B = 8.5 GHz, the strongest line
    of the tilted tip into Sz = -1 above the one Sz = -2 level, and that
    level, in GHz above the ground state"""
    sc = chain(4, 30., 8.5, 3.0)
    e, mz, vV = levels(sc, esrstm._tip_operator(sc, 0, TIP_TILTED))
    e2 = np.min(e[np.abs(mz + 2.) < 1e-6])
    cand = np.flatnonzero((np.abs(mz + 1.) < 1e-6) & (e > e2)
                          & (vV > 1e-6*np.max(vV[1:])))
    return e[cand[np.argmax(vV[cand])]], e2


def assert_close(a, b, tol):
    pk = np.max(np.abs(a))
    assert pk > 0.
    assert np.max(np.abs(a - b)) < tol*pk


@pytest.mark.parametrize("Bt", [3.0, 0.0])
def test_dimer_with_the_tip_along_the_field_on_ed(Bt):
    """Finding 1: the guide's dimer (J = 30 GHz, B = 5 GHz, tip field Bt
    along z, 30 mK, tip=[0,0,1]). T- lies below the window, T+ between the
    window and emax, and neither is reached by V = S^z_0 from the singlet;
    before, the route kept 2 of the 4 states"""
    sc = chain(2, 30., 5., Bt)
    f = dimer_window(Bt)
    _, a = sc.get_esr_spectrum(f, tip=TIP_Z, **KW)
    _, b = sc.get_esr_spectrum(f, tip=TIP_Z, submode="KRYLOV", **KW)
    assert_close(a, b, 1e-8)


@pytest.mark.parametrize("Bt", [3.0, 0.0])
def test_closure_guard_on_ed(Bt, monkeypatch):
    """Without the closing stage T+ is missing, which the record measured
    as 1.2e-4 of the peak at Bt = 3 GHz (about 6e-6 of the relaxation rate
    of T0, by hand, and about 5e-5 at Bt = 0). The guard refuses and names
    submode=None"""
    monkeypatch.setattr(spectrum, "_closing_krylov",
                        lambda Hs, Q1, S0, nkry, top, tol: Q1[:, :0])
    sc = chain(2, 30., 5., Bt)
    with pytest.raises(ValueError, match="submode=None"):
        sc.get_esr_spectrum(dimer_window(Bt), tip=TIP_Z, submode="KRYLOV",
                            **KW)


@pytest.mark.parametrize("Bt", [3.0, 0.0])
def test_dmrg_krylov_refuses_the_dimer(Bt):
    """Finding 1 on DMRG ("python"; v3 runs two sites on ED): T- comes
    from DMRG, since it lies below the window, but T+ lies between the
    window and emax, in a sector the tangent space of the singlet never
    reaches from V = S^z_0. The guard refuses (before, 1.2e-4 of the peak
    off, silently)"""
    np.random.seed(0)
    sc = chain(2, 30., 5., Bt, iv="python")
    with pytest.raises(ValueError, match="submode=None"):
        sc.get_esr_spectrum(dimer_window(Bt), tip=TIP_Z, mode="DMRG",
                            submode="KRYLOV", nstates=3, **KW)


@pytest.mark.parametrize("iv", BACKENDS)
def test_dmrg_krylov_refuses_a_partner_in_the_window(iv):
    """Finding 1, window W2 of its 4-spin chain: T- lies inside the window,
    so DMRG does not give it, and V = S^z_0 never leaves Sz_tot = 0; the
    route was 5.5 of the peak off. The guard refuses"""
    np.random.seed(0)
    sc = chain(4, 30., 5., 3.0)
    e, mz, v0 = levels(sc, sc.Sz[0])
    f0 = e[np.argmax(v0[1:]) + 1]
    etm = np.min(e[np.abs(mz - (mz[0] - 1.)) < 1e-6])
    assert etm < f0 # T- below the final state of the S-T0 line
    f = np.linspace(etm - 1.0, f0 + 0.5, 3001)
    scd = chain(4, 30., 5., 3.0, iv=iv)
    with pytest.raises(ValueError, match="submode=None"):
        scd.get_esr_spectrum(f, tip=TIP_Z, mode="DMRG", submode="KRYLOV",
                             nstates=6, **KW)


def test_tilted_tip_two_sz_steps_away_on_ed():
    """Finding 2: with a tilted tip V reaches Sz = 0 and +-1 from the
    ground state, and the empty Sz = -2 level lies below the final state
    of the line into Sz = -1, below the window; before, 2.1e-2 of the peak
    off"""
    f0, e2 = tilted_line()
    assert e2 < f0
    f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
    sc = chain(4, 30., 8.5, 3.0)
    _, a = sc.get_esr_spectrum(f, tip=TIP_TILTED, **KW)
    _, b = sc.get_esr_spectrum(f, tip=TIP_TILTED, submode="KRYLOV", **KW)
    assert_close(a, b, 1e-8)


@pytest.mark.parametrize("iv", BACKENDS)
def test_tilted_tip_on_dmrg_holds_its_partners(iv):
    """The chain of finding 2 on DMRG: the Sz = -2 level comes from DMRG
    (below the window), the tangent space of the ground state reaches
    every other state below emax, and the guard, which runs since H
    conserves Sz_tot, lets the call through: against ED"""
    np.random.seed(0)
    f0, _ = tilted_line()
    f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
    _, a = chain(4, 30., 8.5, 3.0).get_esr_spectrum(f, tip=TIP_TILTED, **KW)
    sc = chain(4, 30., 8.5, 3.0, iv=iv)
    _, b = sc.get_esr_spectrum(f, tip=TIP_TILTED, mode="DMRG",
                               submode="KRYLOV", nstates=10, **KW)
    assert_close(a, b, 5e-5)


@pytest.mark.parametrize("d", [0.0, 0.005])
@pytest.mark.parametrize("T", [0.0, 0.05])
def test_second_photon_of_a_spin_one_on_ed(d, T):
    """Finding 3, the submode="KRYLOV" half: a spin 1 in a field (rungs at
    20 GHz, or 10 MHz apart with d), saturated (rabi = 5e-3). The window
    reaches 2 F0 + 0.5 GHz, so both routes keep |+1>; before, the Krylov
    space of V|-1> was span{|-1>, |0>} and the second rung was missing
    (one pair against two, the line 36% low)"""
    F0, th = 20.0, 0.7
    sc = spinchain.Spin_Chain(["S=1"])
    Hm = H*F0*sc.Sz[0]
    if d != 0.:
        Hm = Hm + H*d*sc.Sz[0]*sc.Sz[0]
    sc.set_hamiltonian(Hm)
    f = np.append(np.linspace(19.9, 20.1, 41), 2*F0 + 0.5)
    kw = dict(tip=[np.sin(th), 0., np.cos(th)], T=T, T1=100., T2=50.,
              f_ref=F0, rabi=5e-3, return_poles=True)
    _, a, pa = sc.get_esr_spectrum(f, **kw)
    _, b, pb = sc.get_esr_spectrum(f, submode="KRYLOV", **kw)
    assert len(pa["frequency"]) == 2
    assert len(pb["frequency"]) == 2
    assert_close(a, b, 1e-8)
