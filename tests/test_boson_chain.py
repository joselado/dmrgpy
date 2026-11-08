"""Functional coverage for bosonchain.Bosonic_Chain, distilled from
examples/bosonic_chain and examples/boson_maxnb_v3_VS_ED down to small,
fast systems, comparing ED against DMRG on ITensor v2, v3, and the
pure-Python backend (see _helpers.py). Only Bosonic_Chain is covered
here, not SpinBoson_Chain (its ED backend is an unfinished stub, see
pyboson/boson.py).
"""
import numpy as np
import pytest

from dmrgpy import bosonchain
from _helpers import julia_live_param
from _julia_site_helpers import check_julia_live_against_ed

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3, setup_backend

# Looser than the 1e-6 spin/fermion tests get away with: mpscpp3's DMRG
# seeds from an actual random MPS (default_mps() = randomMPS(...), see
# CLAUDE.md/get_sites.h), not a fixed/seeded starting point. Spin/fermion
# tests dodge this because their dim-2 sites fit losslessly under the
# default maxm regardless of starting point; boson's dim-4+ sites don't,
# so occasional runs converge to a slightly different local optimum.
DMRG_TOL = 1e-4


def test_boson_dimer_hopping_energy():
    """Two-site boson hopping Hamiltonian H = t(Adag_0 A_1 + h.c.) at the
    default per-site dimension (maxnb=4, max 3 bosons/site). Unlike the
    spinless-fermion dimer (test_fermion_chain.py), bosons have no
    exclusion principle, so the ground state isn't the single-particle
    bonding orbital at E=-t: it piles multiple bosons into the bonding
    mode, up to what each site's occupation cutoff allows, giving the
    golden value -2*sqrt(3) checked here. Also exercises the
    itensor_version==3, ns<3 ED fallback (same mpscpp3 two-site-dmrg()
    limitation as spin/fermion dimers)."""
    n = 2
    bc = bosonchain.Bosonic_Chain(n)
    h = bc.Adag[0] * bc.A[1]
    h = h + h.get_dagger()

    golden = -2 * np.sqrt(3)
    e_ed, e_v2, e_v3 = energy_ed_v2_v3(bc, h)
    assert e_ed == pytest.approx(golden, abs=1e-8)
    assert e_v2 == pytest.approx(golden, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(golden, abs=DMRG_TOL)


def test_boson_chain_energy_and_density_default_maxnb():
    """4-site Bose-Hubbard-like chain (hopping + onsite repulsion) at the
    default per-site dimension (maxnb=4, i.e. site type code 104 on
    every backend): ground-state energy and total occupation must agree
    between ED and DMRG on both ITensor v2 and v3."""
    n = 4
    bc = bosonchain.Bosonic_Chain(n)
    bc.maxm = 120 # dim-4 sites need more bond dimension than the dim-2
    bc.nsweeps = 50 # sites spin/fermion tests get away with at the default
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 1.0 * bc.N[i] * (bc.N[i] - 1.0)

    Nop = 0
    for i in range(n):
        Nop = Nop + bc.N[i]

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(bc, h)
    assert e_v2 == pytest.approx(e_ed, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(e_ed, abs=DMRG_TOL)

    n_ed, n_v2, n_v3 = vev_ed_v2_v3(bc, h, Nop)
    assert n_v2.real == pytest.approx(n_ed.real, abs=DMRG_TOL)
    assert n_v3.real == pytest.approx(n_ed.real, abs=DMRG_TOL)


def test_boson_chain_get_density_and_fluctuation():
    """get_density/get_density_fluctuation (bosonchain.py) must agree
    between ED and DMRG (v3), and the total density must match the
    total-N vev computed the direct way."""
    n = 3
    bc = bosonchain.Bosonic_Chain(n)
    bc.maxm = 120
    bc.nsweeps = 50
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 0.5 * bc.N[i] * (bc.N[i] - 1.0)
    bc.set_hamiltonian(h)

    d_ed = bc.get_density(mode="ED")
    fluc_ed = bc.get_density_fluctuation(mode="ED")

    bc.set_hamiltonian(h)
    bc.setup_cpp(3)
    d_v3 = bc.get_density(mode="DMRG")
    fluc_v3 = bc.get_density_fluctuation(mode="DMRG")

    assert d_v3 == pytest.approx(d_ed, abs=DMRG_TOL)
    assert fluc_v3 == pytest.approx(fluc_ed, abs=DMRG_TOL)
    assert np.all(fluc_ed >= -1e-8) # fluctuations are non-negative


def test_boson_chain_nondefault_maxnb_matches_ed_on_v3_and_python():
    """Regression test for the maxnb-wiring fix (bosonchain.py's
    Bosonic_Chain now encodes each site's requested dimension into its
    DMRG type code, 100+dim, instead of always building the fixed
    4-level BosonFourSite regardless of maxnb -- see
    mpscpp3/get_sites.h and, on the pure-Python side,
    pyitensor/sites/boson.py's get_boson_site() factory). Uses a
    heterogeneous, non-default maxnb so that a hardcoded dim=4 site
    would silently disagree with ED. v2 is not checked: it still only
    understands the single fixed-dim-4 boson type code (104), see
    examples/boson_maxnb_v3_VS_ED."""
    n = 3
    maxnb = [3, 5, 3]
    bc = bosonchain.Bosonic_Chain(n, maxnb=maxnb)
    assert bc.sites == [103, 105, 103]

    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 0.7 * bc.N[i] * (bc.N[i] - 1.0)

    bc.set_hamiltonian(h)
    e_ed = bc.gs_energy(mode="ED")

    bc.set_hamiltonian(h)
    bc.maxm = 120
    bc.nsweeps = 50
    bc.setup_cpp(3)
    e_v3 = bc.gs_energy(mode="DMRG")
    assert e_v3 == pytest.approx(e_ed, abs=1e-4)

    bc.set_hamiltonian(h)
    setup_backend(bc, "python")
    e_py = bc.gs_energy(mode="DMRG")
    assert e_py == pytest.approx(e_ed, abs=1e-6)


@pytest.mark.parametrize("maxnb,t", [(2, 1.0), (3, 1.0), (5, 0.2)])
@pytest.mark.parametrize("v", [julia_live_param()])
def test_bose_hubbard_julia_live_matches_ed(maxnb, t, v):
    """julia_live builds an ITensors.jl Boson site of dimension maxnb
    (site-type code 100+maxnb), with the occupation projectors N0, N1, ...
    of sites_custom.jl: 4-site Bose-Hubbard chain against ED. The maxnb=5
    case has a weak hopping: at t=1 the julia_live solver, not the site,
    stays 1e-5 to 1e-2 above the ED energy from a random start."""
    n = 4
    bc = bosonchain.Bosonic_Chain(n, maxnb=[maxnb]*n)
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i]*bc.A[i + 1]
    h = -t*(h + h.get_dagger())
    mu = [0.5, 0.6, 0.4, 0.55]
    for i in range(n):
        h = h + 1.0*bc.N[i]*(bc.N[i] - 1.0) - mu[i]*bc.N[i]
    obs = [bc.N[0], bc.Adag[0]*bc.A[2], bc.D[1][0], bc.D[2][maxnb - 1],
           bc.N[1]*bc.N[3]]
    check_julia_live_against_ed(bc, h, obs, maxm=40, nsweeps=30)
