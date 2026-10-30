"""The Kondo screening cloud of an Anderson impurity at the end of a
tight-binding lead, Holzner et al., PRB 80, 205114 (2009),
arXiv:0906.2933 (kondocloud.py)."""

import warnings

import numpy as np
import pytest

from dmrgpy import kondocloud


def _lead(L=6, chain="native", **kw):
    p = dict(ed=-0.5, U=1.0, Gamma=0.5)
    p.update(kw)
    fc = kondocloud.anderson_lead(L, itensor_version="python", chain=chain, **p)
    fc.maxm = 64
    fc.nsweeps = 12
    fc.cutoff = 1e-12
    return fc


@pytest.mark.parametrize("chain", ["native", "interleaved"])
def test_cloud_dmrg_against_ed(chain):
    fc = _lead(chain=chain)
    ed = fc.get_kondo_cloud(mode="ED")
    dm = fc.get_kondo_cloud()
    assert np.max(np.abs(dm["SdSi"] - ed["SdSi"])) < 1e-6
    assert np.max(np.abs(dm["Sigma"] - ed["Sigma"])) < 1e-6
    assert dm["SdSd"] == pytest.approx(ed["SdSd"], abs=1e-6)
    assert dm["xi"] == ed["xi"]


@pytest.mark.parametrize("ed", [-0.5, -0.2, 0.3])
def test_singlet_sum_rule_is_exact(ed):
    # <S_tot^2> = 0 gives sum_i <S_d.S_i> = -<S_d^2>, so Sigma(L-1) = 0
    fc = _lead(L=6, ed=ed, half_filling=True)
    r = fc.get_kondo_cloud(mode="ED")
    assert r["Sigma"][0] == 1.0
    assert abs(r["sum_rule"]) < 1e-12
    assert np.sum(r["SdSi"]) == pytest.approx(-r["SdSd"], abs=1e-12)


def test_an_odd_chain_warns_and_breaks_the_sum_rule():
    with pytest.warns(UserWarning, match="odd"):
        fc = _lead(L=5)
    with pytest.warns(UserWarning, match="not a total-spin singlet"):
        r = fc.get_kondo_cloud(mode="ED")
    assert abs(r["sum_rule"]) > 1e-2


def test_norm_is_three_quarters_of_single_occupancy():
    fc = _lead()
    r = fc.get_kondo_cloud(mode="ED")
    m = fc.get_kondo_screening(mode="ED")
    assert r["SdSd"] == pytest.approx(0.75*m["M_imp"], abs=1e-12)
    assert r["mu"] == pytest.approx(0.25*m["M_imp"], abs=1e-12)
    assert r["nd"] == pytest.approx(1.0, abs=1e-12) # particle-hole symmetry


def test_connected_form_equals_the_full_form_in_zero_field():
    # SU(2): in a singlet <S_d.S_i>/<S_d.S_d> = <Sz_d Sz_i>/<Sz_d Sz_d>
    fc = _lead()
    a = fc.get_kondo_cloud(mode="ED")
    b = fc.get_kondo_cloud(mode="ED", connected=True)
    assert np.max(np.abs(a["Sigma"] - b["Sigma"])) < 1e-12


def test_field_dmrg_against_ed_and_the_moment_shrinks():
    fc = _lead(B=0.2)
    ed = fc.get_kondo_cloud(mode="ED", connected=True)
    dm = fc.get_kondo_cloud(connected=True)
    assert np.max(np.abs(dm["Sigma"] - ed["Sigma"])) < 1e-6
    r0 = _lead().get_kondo_cloud(mode="ED")
    assert ed["mu"] < r0["mu"]


@pytest.mark.parametrize("chain", ["native", "interleaved"])
def test_resonant_level_against_wick(chain):
    # U=0: <S_d.S_i> = -3/2 G_di^2 and <S_d.S_d> = 3/2 G_dd (1-G_dd).
    # A metallic lead needs a large bond dimension: at L=8 and maxm=100 the
    # far correlations are off by 5e-6, so we take maxm=256, which holds
    # the state of L=8 exactly.
    L, G, t = 8, 0.5, 1.0
    h = np.zeros((L, L))
    h[0, 1] = h[1, 0] = -np.sqrt(G*t)
    for i in range(1, L-1):
        h[i, i+1] = h[i+1, i] = -t
    ev, U = np.linalg.eigh(h)
    g = U[:, ev < 0] @ U[:, ev < 0].T
    fc = _lead(L=L, chain=chain, ed=0.0, U=0.0, Gamma=G)
    fc.maxm = 256
    r = fc.get_kondo_cloud()
    assert r["SdSd"] == pytest.approx(1.5*g[0, 0]*(1-g[0, 0]), abs=1e-6)
    assert np.max(np.abs(r["SdSi"] + 1.5*g[0, 1:]**2)) < 1e-6


def test_screening_length_definition():
    s = np.array([1.0, 0.5, 0.3, 0.12, 0.08, 0.0])
    assert kondocloud.screening_length(s, 0.9) == 4
    assert kondocloud.screening_length(s, 0.5) == 1
    assert np.isnan(kondocloud.screening_length(s[:4], 0.9))
    r = _lead().get_kondo_cloud(mode="ED", a=[0.5, 0.9])
    assert set(r["xi"]) == {0.5, 0.9}


def test_kondo_temperature_matches_the_paper():
    # the paper quotes xi_0.9 = 112 at U/Gamma = 5, L = 500, and p(0.9) = 6.7
    TK = kondocloud.kondo_temperature(1.0, 0.2)
    assert 6.7/TK == pytest.approx(112, rel=0.06)
    assert kondocloud.kondo_length(1.0, 0.2) == pytest.approx(2.0/TK)
    # away from ed = -U/2 the exponent shrinks, so T_K grows
    assert kondocloud.kondo_temperature(1.0, 0.2, ed=-0.25) > TK


def test_screening_length_scales_as_one_over_TK():
    # The paper's central result, xi_a = p(a)/T_K. At U/Gamma = 2 and 3,
    # xi_0.75 has saturated by L = 40 (the same at L = 60, maxm = 80), and
    # xi_0.75 T_K is the same for both, 1.61 and 1.60. The paper's own
    # numbers (xi_0.75 = 29 at U/Gamma = 5, L = 500) give p(0.75) = 1.8.
    p = []
    for ratio in [2.0, 3.0]:
        fc = _lead(L=40, Gamma=1.0/ratio)
        fc.nsweeps = 10
        fc.cutoff = 1e-10
        r = fc.get_kondo_cloud(a=0.75)
        p.append(r["xi"]*kondocloud.kondo_temperature(1.0, 1.0/ratio))
    assert p[1] == pytest.approx(p[0], rel=0.05)
    assert 1.4 < p[0] < 2.0
