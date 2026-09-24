"""Regression tests for the `kondo` cluster of the 2026-09-24 hole hunt
(docs/audit_2026_09_24_hole_hunt.md, findings 11 to 14).

  11. `get_kondo_spectrum(mode="ED")` accepted `**kwargs` and read none of
      them, so a misspelled physical parameter (`Jrho=` for `Jrho_s=`)
      silently fell back to its default. It now raises TypeError.
  12. `get_kondo_spectrum(mode="DMRG")` at an accidental ground-state
      degeneracy returned the spectrum of whichever member of the manifold
      the random start converged to, anywhere in [1.0, 2.0] against the
      T->0+ value 1.5 that mode="ED" defines. `n_gs=` now averages over
      the manifold; the default (n_gs=1) is the single converged state.
  13. The potential-interference DMRG term weighted every point of `es`
      with its first spacing, 101.7 times the exact value on a smooth grid
      dense at the line. It is a trapezoid rule on the actual grid now.
  14. The potential term needs `es` to reach past the top of the S_k
      spectrum, not just past max|eV| like the second-order term; the
      example that violated it is fixed, and the term now warns when the
      sum rule sum_k int S_kk = S(S+1) says the grid misses weight.

Findings 11, 13 and 14 run on ED and exact Lehmann correlators, so no
convergence enters; finding 12 runs on itensor_version="python", where the
single-state pick was reproduced by the hunter, at delta=1e-4 (the values
checked sit on plateaus 1 meV away from any line, so the broadening does
not move them).
"""

import warnings

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc

G = 2.0
MUB = 5.7883818066e-5 # eV/T
TP = 2*np.pi # get_kondo_spectrum returns 2*pi times dI/dV in e^2 T0^2/hbar


def zeeman_spin(B=10.0):
    """The single S=1/2 of tests/test_kondo_spectrum_potentialdc.py"""
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*B*sc.Sz[0])
    return sc


def crossing_chain(itensor_version="python", D=1e-3):
    """Finding 12's chain: an S=1 impurity with D*Sz^2 + D*Sz, at the
    crossing where |0> and |-1> are degenerate, next to two field-polarized
    S=1/2 spectators (three sites, so every DMRG backend runs it)"""
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"],
                              itensor_version=itensor_version)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc


# ------------------------------------------------------------ finding 11

@pytest.mark.parametrize("bad", [dict(Jrho=-0.05), dict(delta=1e-3),
                                 dict(submode="CVM"), dict(n_gs=2)])
def test_ed_mode_rejects_keywords_it_does_not_read(bad):
    sc = zeeman_spin(B=0.0)
    eVs = np.array([-4e-3, 0.0, 4e-3])
    # the keyword list itself, not the explanation after it (which names
    # the mode="DMRG" parameters anyway)
    with pytest.raises(TypeError, match=r"argument\(s\): %s\. " % list(bad)[0]):
        sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0,
                              omega0=20e-3, **bad)


def test_ed_mode_names_every_unknown_keyword_sorted():
    sc = zeeman_spin(B=0.0)
    with pytest.raises(TypeError, match=r"argument\(s\): Jrho, delta, u\. "):
        sc.get_kondo_spectrum(np.array([0.0]), site=0, u=0.25, delta=1e-3,
                              Jrho=-0.05)


def test_ed_mode_correct_spelling_is_unchanged():
    # the paper's Fig. 3b/7b zero-bias peak, pinned in
    # test_kondo_spectrum_paper_fig7.py at 1.137 in the figure's units
    sc = zeeman_spin(B=0.0)
    _, d = sc.get_kondo_spectrum(np.array([0.0]), site=0, Jrho_s=-0.05,
                                 T=1.0, omega0=20e-3)
    assert d[0]/TP == pytest.approx(1.1374, abs=1e-4)


# ------------------------------------------------------------ finding 12

def test_n_gs_averages_the_degenerate_manifold_like_ed():
    np.random.seed(0)
    sc = crossing_chain("python")
    eVs = np.array([-1e-3, 0.0, 1e-3])
    es = np.linspace(-30e-3, 30e-3, 3001)
    kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode="KPM",
              delta=1e-4, es=es)
    ks = KondoSpectrum(sc, 0, T=0.0)
    assert np.allclose(ks.p[:3], [0.5, 0.5, 0.0]) # an exact crossing
    _, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
    assert np.allclose(ed/TP, 1.5, atol=1e-10)
    # n_gs=1 (the default) is the single converged state: at this crossing
    # its zero-bias value is 1-<Sz_0> of that state, whichever member or
    # superposition of the two the random start reached
    _, one = sc.get_kondo_spectrum(eVs, **kw)
    sz = np.real(sc.vev(sc.Sz[0]))
    assert np.allclose(one/TP, 1.0 - sz, atol=2e-3)
    # n_gs=2 is the manifold average, mode="ED"'s T->0+ limit, whichever
    # basis of the manifold DMRG returns
    _, avg = sc.get_kondo_spectrum(eVs, n_gs=2, **kw)
    assert np.allclose(avg/TP, 1.5, atol=5e-3)
    # and the chain's own ground state is put back afterwards, on both the
    # Python side (vev) and the session (the KPM correlator reads the
    # session's state, so the restore has to reach it too)
    assert np.real(sc.vev(sc.Sz[0])) == pytest.approx(sz, abs=1e-12)
    _, again = sc.get_kondo_spectrum(eVs, **kw)
    assert np.allclose(again, one, rtol=0., atol=1e-9)


def test_n_gs_averages_the_two_time_kondo_term_too(monkeypatch):
    """order=3: the third-order Kondo term, a three-point function of the
    ground state, is linear in its density matrix as well, so n_gs=2 must
    give the average of the ED two-time term over the two degenerate
    eigenvectors, on the literal same (t2,tau) grid (the construction of
    test_kondo_spectrum_dmrgtwotime.py's grid-consistent reference; a
    deliberately coarse grid, since it is the same for both sides). The
    second-order term is stubbed to zero, so only the two-time term runs,
    and the bias points are few because the K_W kernel quadrature, not the
    time evolution, is what this costs."""
    import dmrgpy.kondospectrumtk.secondorder_dc as secondorder_dc
    from dmrgpy.kondospectrumtk.edtwotimeref import _levi_civita_coeff_G_chunk
    from dmrgpy.kondospectrumtk.twotime import kondo_term_from_two_time
    eVs = np.array([0.0, 1e-3, 2e-3]) # the term is even in eV
    monkeypatch.setattr(secondorder_dc, "second_order_dIdV_dc",
                        lambda *a, **kw: np.zeros(len(eVs)))
    np.random.seed(0)
    sc = crossing_chain("python")
    ks = KondoSpectrum(sc, 0, T=0.0)
    omega0, Gamma0, Jrho_s = 2e-3, 5e-6, 0.05
    n_t2_half, n_tau_half = 2, 3
    dt2, dtau = 25./Gamma0/10, (2*np.pi/2e-5)/15
    _, d = sc.get_kondo_spectrum(
            eVs, site=0, Jrho_s=Jrho_s, T=0.0, order=3, mode="DMRG",
            omega0=omega0, Gamma0=Gamma0, es=np.array([0.]), n_gs=2,
            dt2=dt2, n_t2_half=n_t2_half, dtau=dtau, n_tau_half=n_tau_half)
    term = d/(4*np.pi*Jrho_s)
    t2_grid = dt2*np.arange(-n_t2_half, n_t2_half+1)
    tau_grid = dtau*np.arange(-n_tau_half, n_tau_half+1)
    G = []
    for first in (0, 1): # each degenerate eigenvector as "the" ground state
        order = [first, 1-first] + list(range(2, ks.dim))
        k = KondoSpectrum.__new__(KondoSpectrum)
        k.e = ks.e[order]
        for a in ("Sx", "Sy", "Sz"):
            setattr(k, a, getattr(ks, a)[np.ix_(order, order)])
        G.append(_levi_civita_coeff_G_chunk(k, t2_grid, tau_grid))
    # the two members' three-point functions do differ (on a 9x13 grid their
    # zero-bias terms were -11.98 and -17.17), and the term is linear in G,
    # so the reference is the term of the averaged G
    assert np.max(np.abs(G[0] - G[1])) > 0.1*np.max(np.abs(G[0]))
    ref = kondo_term_from_two_time(t2_grid, tau_grid,
                                   iter([(t2_grid, 0.5*(G[0] + G[1]))]),
                                   eVs, omega0, Gamma0)
    assert np.allclose(term, ref, rtol=0., atol=1e-6)


def test_n_gs_1_never_touches_the_excited_states(monkeypatch):
    """The default path is the pre-fix one: no excited-state solve and no
    change of ground state (the three terms are stubbed out, as in
    test_kondo_spectrum.py's wiring test)"""
    import dmrgpy.kondospectrumtk.secondorder_dc as secondorder_dc
    eVs = np.linspace(-1e-3, 1e-3, 3)
    monkeypatch.setattr(secondorder_dc, "second_order_dIdV_dc",
                        lambda *a, **kw: np.full(len(eVs), 1.0))
    sc = zeeman_spin(B=5.0)
    def forbidden(*a, **kw):
        raise AssertionError("n_gs=1 must not solve for excited states")
    monkeypatch.setattr(sc, "get_excited_states", forbidden)
    monkeypatch.setattr(sc, "set_gs", forbidden)
    for extra in ({}, dict(n_gs=1)):
        _, d = sc._get_kondo_spectrum_dmrg(eVs, 0, 0.5, 0.0, 2.0, 20e-3,
                                           5e-6, order=2, es=np.array([0.]),
                                           **extra)
        assert np.allclose(d, 1.0)


@pytest.mark.parametrize("n_gs", [0, -1, 1.5, True, "2"])
def test_n_gs_must_be_a_positive_integer(n_gs):
    sc = zeeman_spin(B=5.0)
    with pytest.raises(ValueError, match="n_gs"):
        sc.get_kondo_spectrum(np.array([0.0]), site=0, T=0.0, order=2,
                              mode="DMRG", es=np.array([0.0]), n_gs=n_gs)


def test_n_gs_needs_an_mps_session():
    # a 1-site "python" chain runs its correlators on ED (two-site DMRG
    # has no update to make there), so there is no session state to set
    sc = spinchain.Spin_Chain(["1/2"], itensor_version="python")
    sc.set_hamiltonian(G*MUB*5.0*sc.Sz[0])
    with pytest.raises(NotImplementedError, match="n_gs"):
        sc.get_kondo_spectrum(np.array([0.0]), site=0, T=0.0, order=2,
                              mode="DMRG", es=np.array([0.0]), n_gs=2)


# ------------------------------------------------------------ finding 13

def test_potential_term_on_a_non_uniform_grid_matches_the_exact_sum():
    sc = zeeman_spin(B=10.0)
    ks = KondoSpectrum(sc, site=0, T=0.0)
    eVs = np.linspace(-1e-3, 2e-3, 21)
    Jrho_s, U = 0.1, 0.3
    ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
    peak = np.max(np.abs(ref))
    line = ks.e[1] # the 1.1577 meV Zeeman line
    # smooth sinh map concentrating points at the line: spacing 1.8e-8
    # there and 1.9e-6 at the ends, so the first spacing is 100 times the
    # local one at the line (the shipped rule returned 101.7x the exact
    # value at the peak on this grid)
    u = np.linspace(np.arcsinh((-1e-3 - line)/2e-5),
                    np.arcsinh((3e-3 - line)/2e-5), 12000)
    es = line + 2e-5*np.sinh(u)
    kw = dict(T0=1.0, mode="ED", submode="ED", delta=2e-6)
    mine = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, es=es, **kw)
    # 7.6e-4 of the peak, the Lorentzian broadening delta=2e-6, the same
    # residual as on the uniform grid of test_kondo_spectrum_potentialdc.py
    assert np.max(np.abs(mine - ref)) < 1e-3*peak
    uniform = third_order_potential_dIdV_dc(
            sc, 0, eVs, Jrho_s, U, es=np.linspace(-1e-3, 3e-3, 40_000), **kw)
    assert np.max(np.abs(mine - uniform)) < 1e-5*peak # measured 1.6e-7


# ------------------------------------------------------------ finding 14

def _potential_warnings(sc, eVs, es):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = third_order_potential_dIdV_dc(sc, 0, eVs, 0.1, 0.3, T0=1.0,
                                            mode="ED", submode="ED",
                                            delta=2e-6, es=es)
    return out, [w for w in caught
                 if "third_order_potential_dIdV_dc" in str(w.message)]


def test_potential_term_warns_when_es_misses_spectral_weight():
    sc = zeeman_spin(B=10.0)
    eVs = np.linspace(-1e-3, 2e-3, 7)
    # stops below the 1.16 meV line, which holds 0.5 of S(S+1)=0.75: the
    # result is then almost entirely missing, at biases the grid covers
    _, caught = _potential_warnings(sc, eVs, np.linspace(-1e-3, 1e-3, 8000))
    assert len(caught) == 1 and issubclass(caught[0].category, RuntimeWarning)
    assert "S(S+1) = 0.75" in str(caught[0].message)
    # the test suite's own grid covers everything, and stays silent
    _, caught = _potential_warnings(sc, eVs, np.linspace(-1e-3, 3e-3, 40_000))
    assert caught == []
