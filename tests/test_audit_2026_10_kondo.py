"""Regressions of the 2026-10 hole hunt (docs/audit_2026_10_hole_hunt.md),
cluster kondo: the potential term's es coverage warning on fermion chains
(finding 2), the TASK warning for a state below the ground state (finding
3, and a state below it that the impurity spin does not reach), order=3
on a short chain and on v2 (finding 4), and nkry on the short-chain
fallback (finding 6)."""
import warnings
import numpy as np
import pytest
from dmrgpy import fermionchain, spinchain
from dmrgpy.kondospectrumtk import potentialdc, conductance
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum


def hubbard(n=3, U=2., b=0.3, itensor_version=3):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(
            n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(n):
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) + b*fc.Sz[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 64, 20, 0.
    return fc


def heisenberg(n, itensor_version=3, b=0.):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + b*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 64, 20, 0.
    return sc


def coverage_warnings(fc, top):
    eV = np.linspace(-1.2, 1.2, 7)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        potentialdc.third_order_potential_dIdV_dc(
                fc, 1, eV, -0.05, 0.3, omega0=4., Gamma0=1e-2,
                submode="KPM", delta=0.05, es=np.linspace(-0.3, top, 400))
    return [x for x in w if "es grid" in str(x.message)]


def test_sum_rule_on_a_hubbard_site_is_the_measured_spin_weight():
    fc = hubbard()
    ks = KondoSpectrum(fc, 1, T=0.)
    exact = sum((M@M)[0, 0].real for M in (ks.Sx, ks.Sy, ks.Sz))
    assert potentialdc._spin_casimir(fc, 1) == pytest.approx(exact, abs=1e-8)
    assert exact < 0.75 - 1e-3 # charge fluctuations, not S(S+1)


def test_narrow_es_warns_on_a_hubbard_chain_and_a_wide_one_does_not():
    fc = hubbard()
    assert len(coverage_warnings(fc, 0.8)) == 1 # spectrum reaches 6.8
    assert not coverage_warnings(fc, 12.)


def test_task_warns_when_handed_a_state_above_the_ground_state():
    sc = heisenberg(5, itensor_version="python", b=0.1)
    es, wfs = sc.get_excited_states(n=2)
    sc.set_gs(wfs[1]) # the tangent space of it reaches the ground state
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sc.get_kondo_spectrum(np.linspace(-0.3, 0.3, 5), T=0., mode="DMRG",
                              Jrho_s=-0.05, omega0=4., Gamma0=1e-2, kB=1.)
    assert any("not the ground state" in str(x.message) for x in w)


def test_task_on_the_ground_state_does_not_warn():
    sc = heisenberg(5, itensor_version="python", b=0.1)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sc.get_kondo_spectrum(np.linspace(-0.3, 0.3, 5), T=0., mode="DMRG",
                              Jrho_s=-0.05, omega0=4., Gamma0=1e-2, kB=1.)
    assert not any("not the ground state" in str(x.message) for x in w)


def _task_warnings(monkeypatch, coupled):
    # a Ritz state at -0.5 next to the ground state and one excitation;
    # the impurity spin reaches it from row 0 only when coupled
    from dmrgpy import taskdmrg
    from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
    e = np.array([0., -0.5, 0.3])
    M = np.diag([0.5, 0.5, -0.5]).astype(complex)
    M[0, 2] = M[2, 0] = 0.4
    if coupled:
        M[1, 2] = M[2, 1] = 0.1
    monkeypatch.setattr(taskdmrg, "session_spectrum",
                        lambda chain, ops, nkry: (e, [M, 0*M, M], 10))
    sc = heisenberg(3, itensor_version="python")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        TaskKondoSpectrum(sc, 0, nkry=5)
    return [str(x.message) for x in w]


def test_task_state_below_e0_that_the_spin_does_not_reach_is_named_so(
        monkeypatch):
    # a promoted sector state: rounding seeds a lower charge sector, which
    # carries no matrix element to the ground state (Ne=21 next to Ne=20)
    msgs = _task_warnings(monkeypatch, coupled=False)
    assert any("do not enter the spectrum" in m for m in msgs)
    assert not any("not the ground state" in m for m in msgs)


def test_task_state_below_e0_reached_through_an_excitation_still_warns(
        monkeypatch):
    msgs = _task_warnings(monkeypatch, coupled=True)
    assert any("not the ground state" in m for m in msgs)


def test_order3_with_a_correlator_submode_on_a_short_chain_uses_ed():
    # v3 below three sites falls back to ED: the Kondo term is then the
    # exact one, so order=3 minus order=2 is conductance's Kondo term
    sc = heisenberg(2, b=0.2)
    eV = np.linspace(-1., 1., 9)
    kw = dict(T=0., mode="DMRG", submode="KPM", delta=0.05,
              es=np.linspace(-1., 3., 400), omega0=4., Gamma0=1e-2, kB=1.,
              Jrho_s=-0.05)
    _, d3 = sc.get_kondo_spectrum(eV, order=3, **kw)
    _, d2 = sc.get_kondo_spectrum(eV, order=2, **kw)
    ref = conductance.third_order_kondo_dIdV(
            KondoSpectrum(sc, 0, T=0.), eV, -0.05, omega0=4., Gamma0=1e-2)
    assert np.max(np.abs((d3 - d2) - ref)) < 1e-10


def test_order3_with_a_correlator_submode_on_v2_raises():
    sc = heisenberg(4, itensor_version=2)
    with pytest.raises(NotImplementedError):
        sc.get_kondo_spectrum(np.linspace(-1., 1., 5), T=0., mode="DMRG",
                              order=3, delta=0.05, es=np.linspace(-1, 3, 50))


@pytest.mark.parametrize("nkry", [-5, "x", 2.5])
def test_nkry_is_checked_on_a_chain_too_short_for_the_session(nkry):
    sc = heisenberg(2)
    with pytest.raises(ValueError):
        sc.get_kondo_spectrum(np.linspace(-1., 1., 5), T=0., mode="DMRG",
                              nkry=nkry)
