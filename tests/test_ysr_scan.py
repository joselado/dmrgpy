"""Warm-started scans of the subgap energies (ysrtk/subgapscan.py): every
point of a Delta or Dz scan, and the bisection of Delta_c, against ED and
against a cold ysr.subgap_energies at the same point."""
import warnings

import numpy as np
import pytest

from dmrgpy import ysr
from dmrgpy.ysrtk import subgapscan

N, J, LAM = 4, 0.6, 4.0


def _reference(r, i, S, mode, itensor_version):
    fc = ysr.kondo_chain(N, J=J, Delta=r["Delta"][i], Dz=r["Dz"][i],
                         Lambda=LAM, S=S, Upen=r["Upen"], JH=r["JH"],
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    return ysr.subgap_energies(fc, mode=mode, sectors=list(r["sectors"][i]))


def _check_against_ed(r, S, itensor_version):
    for i in range(len(r["Delta"])):
        ed = _reference(r, i, S, "ED", itensor_version)
        assert set(r["sectors"][i]) == set(ed["sectors"])
        for s, e in ed["sectors"].items():
            assert r["sectors"][i][s] == pytest.approx(e, abs=1e-8)
        assert r["E_ysr"][i] == pytest.approx(ed["E_ysr"], abs=1e-8)
        assert r["ground_state"][i] == ed["ground_state"]


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("S,Dz", [("1/2", 0.0), ("1", 0.05)])
def test_delta_scan_against_ed_and_cold(S, Dz, itensor_version):
    """Every warm point equals ED to 1e-8 and the cold subgap_energies
    of the same point (cold_every=1) to 1e-9."""
    r = subgapscan.subgap_scan(N, Deltas=[0.05, 0.1, 0.2], J=J, Lambda=LAM,
                               S=S, Dz=Dz, maxm=64, nsweeps=15, cold_every=1,
                               itensor_version=itensor_version)
    assert [p["warm"] for p in r["points"]] == [False, True, True]
    assert sorted(r["cold"]) == [1, 2]
    assert r["nsweeps_warm"] == 8 # the default at nsweeps=15
    assert r["warm_cold"] < 1e-9
    _check_against_ed(r, S, itensor_version)
    first = _reference(r, 0, S, "DMRG", itensor_version)
    for s, e in first["sectors"].items():
        assert r["sectors"][0][s] == pytest.approx(e, abs=1e-9)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_dz_scan_against_ed(itensor_version):
    """A Dz scan through Dz=0 at fixed Upen and JH, every 2Sz from 0 to
    3 at every point; the ground state moves from Sz=+-1 to Sz=0."""
    Dzs = [-0.1, -0.05, 0.0, 0.05, 0.1]
    r = subgapscan.subgap_scan(N, Dzs=Dzs, Delta=0.1, J=J, Lambda=LAM,
                               S="1", maxm=64, nsweeps=15,
                               itensor_version=itensor_version)
    assert r["Upen"] == pytest.approx(2*(1 + J) + 2*2*0.1)
    assert all(set(s) == {0, 1, 2, 3} for s in r["sectors"])
    _check_against_ed(r, "1", itensor_version)
    assert r["Sz_ground"][0] == 2 and r["Sz_ground"][-1] == 0


def test_energy_differences_do_not_depend_on_the_penalties():
    """On the spin-S multiplet Upen is zero and JH a constant, so E_ysr
    of a scan at fixed Upen, JH equals subgap_energies with kondo_chain's
    defaults at that Dz."""
    r = subgapscan.subgap_scan(N, Dzs=[0.0, 0.08], Delta=0.1, J=J,
                               Lambda=LAM, S="1", mode="ED",
                               itensor_version="python")
    for i, Dz in enumerate([0.0, 0.08]):
        fc = ysr.kondo_chain(N, J=J, Delta=0.1, Dz=Dz, Lambda=LAM, S="1",
                             itensor_version="python")
        ed = ysr.subgap_energies(fc, mode="ED", sectors=[0, 1, 2, 3])
        assert r["E_ysr"][i] == pytest.approx(ed["E_ysr"], abs=1e-10)


@pytest.mark.parametrize("S", ["1/2", "1"])
def test_find_transition_against_ed(S):
    """The bisection on v3 with warm starts visits the same points and
    gives the same Delta_c as on ED, and the bracket ends solved cold
    agree with the warm ones."""
    kw = dict(ratio=1.02, Lambda=LAM, S=S)
    # 3 sites: the 4-site chain at J=0.6 has no screened phase
    ed = subgapscan.find_transition(3, J, 0.05, 0.5, mode="ED", **kw)
    dm = subgapscan.find_transition(3, J, 0.05, 0.5, maxm=64, nsweeps=15,
                                    check_cold=True, itensor_version=3, **kw)
    assert ed["x_lo"] > 0 > ed["x_hi"] and ed["hi"]/ed["lo"] <= 1.02
    assert ed["lo"] < ed["Delta_c"] < ed["hi"] and ed["monotonic"]
    assert np.allclose(dm["Delta"], ed["Delta"], rtol=1e-12)
    assert dm["Delta_c"] == pytest.approx(ed["Delta_c"], rel=1e-7)
    assert dm["Delta_c_cold"] == pytest.approx(dm["Delta_c"], rel=1e-8)
    assert dm["warm_cold"] < 1e-9


def test_scan_arguments():
    kw = dict(J=J, Lambda=LAM, mode="ED", itensor_version="python")
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, Deltas=[0.1], Dzs=[0.0], Delta=0.1, **kw)
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, **kw)
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, Dzs=[0.0, 0.1], **kw) # no Delta
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, Deltas=[0.1, 0.0], **kw)
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, Deltas=[0.1], sectors=[0, 2], **kw)
    with pytest.raises(ValueError): # Delta has no consumer in a Delta scan
        subgapscan.subgap_scan(N, Deltas=[0.1, 0.2], Delta=0.5, **kw)
    with pytest.raises(ValueError): # nor Dz in a Dz scan
        subgapscan.subgap_scan(N, Dzs=[0.0, 0.1], Delta=0.1, S="1", Dz=0.7,
                               **kw)
    with pytest.raises(ValueError):
        subgapscan.subgap_scan(N, Deltas=[0.1], nsweeps_warm=0, **kw)
    with pytest.raises(ValueError):
        subgapscan.find_transition(N, J, 0.1, 0.05, Lambda=LAM, mode="ED",
                                   itensor_version="python")


def test_metastable_warning_at_every_point():
    """The quasiparticle-gap check of subgap_energies runs at every point:
    a point whose |E_ysr| exceeds the bath gap warns, naming that point."""
    ens = {0: 0.0, 1: 0.5}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = subgapscan.summarize(ens, 1, qp_gap=0.1, reaches_gap=True,
                                   where=" (Delta = 0.1, Dz = 0)")
    assert out["ground_state"] == "singlet" and out["E_ysr"] == 0.5
    assert len(w) == 1 and "Delta = 0.1" in str(w[0].message)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        subgapscan.summarize(ens, 1, qp_gap=0.1, reaches_gap=False)


def test_default_warm_sweeps():
    """The default nsweeps_warm is the larger of 8 and nsweeps//3, at
    most nsweeps, on every backend; a given value is kept."""
    assert [subgapscan.warm_sweeps(ns) for ns in (20, 15, 30, 4)] == \
        [8, 8, 10, 4]
    assert subgapscan.warm_sweeps(20, 3) == 3


@pytest.mark.parametrize("mesh", ["mean", "zp"])
def test_quasiparticle_gap_on_the_mesh_of_the_chain(mesh):
    """The gap that the metastable warning compares against is that of
    kondo_chain on the same mesh (no solve)."""
    n = 20 # t_{n-1} is about 1e-6: below 0.1 Delta at 1e-4, above at 1e-6
    sc = subgapscan._SectorChains(n, J, LAM, 1.0, "1/2", 0.0, 0.0, mesh,
                                  None, None, 64, 20, None, [0, 1], "DMRG",
                                  {"itensor_version": "python"}, 1e-3, 0.0)
    for Delta in (1e-6, 1e-4, 1e-2): # 1e-6: the chain does not reach it
        fc = ysr.kondo_chain(n, J=J, Delta=Delta, Lambda=LAM, mesh=mesh,
                             itensor_version="python")
        qp, reaches = sc.gap(Delta)
        assert qp == pytest.approx(fc.qp_gap, rel=1e-12)
        assert reaches == fc.reaches_gap
    assert sc.gap(1e-4)[1] and not sc.gap(1e-6)[1]


def test_auto_warm_sweep_arguments():
    assert subgapscan.warm_sweeps(15, "auto") == "auto"
    kw = dict(J=J, Lambda=LAM, mode="ED", itensor_version="python")
    for bad in ("Auto", "fast", True):
        with pytest.raises(ValueError, match="nsweeps_warm"):
            subgapscan.subgap_scan(N, Deltas=[0.1], nsweeps_warm=bad, **kw)
    for bad in (0.0, -1e-6, np.nan):
        with pytest.raises(ValueError, match="warm_tol"):
            subgapscan.subgap_scan(N, Deltas=[0.1], nsweeps_warm="auto",
                                   warm_tol=bad, **kw)


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_auto_warm_sweeps_against_ed_and_cold(itensor_version):
    """nsweeps_warm="auto" on the S=1 Delta scan at Dz=0.05, the one where
    "python" needs the most warm sweeps. On d22 (run warm1) it used 4
    sweeps per point on v3 and 4 to 6 on "python", against 8 for the
    fixed default, and gave ED to 3.2e-13 (v3) and 2.3e-12 ("python") and
    the cold solve to 7.6e-15 and 2.3e-12."""
    r = subgapscan.subgap_scan(N, Deltas=[0.05, 0.1, 0.2], J=J, Lambda=LAM,
                               S="1", Dz=0.05, maxm=64, nsweeps=15,
                               cold_every=1, nsweeps_warm="auto",
                               itensor_version=itensor_version)
    used = r["nsweeps_warm"]
    assert used[0] is None and all(4 <= s <= 15 for s in used[1:])
    if itensor_version == 3:
        assert max(used[1:]) <= 8 # the fixed default at nsweeps=15
    for p in r["points"][1:]: # every block is a sweep from the one before
        assert all(len(b) == (s + 1)//2 for b, s in
                   zip(p["blocks"].values(), p["sweeps"].values()))
    assert r["warm_cold"] < 1e-9
    _check_against_ed(r, "1", itensor_version)
