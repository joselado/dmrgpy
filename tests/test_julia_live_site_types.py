"""julia_live site types for the higher spins (S=3/2, 2, 5/2, defined in
mpsjulialive/sites_custom.jl) against mode="ED", and the error for a
site-type code julia_live does not know."""
import numpy as np
import pytest

from dmrgpy import spinchain

from _helpers import julia_live_param
from _julia_site_helpers import check_julia_live_against_ed


def _heisenberg(sc):
    n = sc.ns
    h = 0
    for i in range(n - 1):
        h = h + sc.SS(i, i + 1)
    return h + 0.3*sc.Sz[0] + 0.2*sc.Sx[n - 1] + 0.1*sc.Sy[1]


@pytest.mark.parametrize("spin", ["3/2", "2", "5/2"])
@pytest.mark.parametrize("v", [julia_live_param()])
def test_higher_spin_heisenberg_matches_ed(spin, v):
    """Heisenberg chain of 4 spins S with fields along x, y and z on
    different sites, so Sx, Sy and Sz all enter the Hamiltonian"""
    sc = spinchain.Spin_Chain([spin]*4)
    obs = [sc.Sz[0], sc.Sx[3], sc.Sy[1], sc.Sz[0]*sc.Sz[2],
           sc.Sy[1]*sc.Sy[2]]
    check_julia_live_against_ed(sc, _heisenberg(sc), obs, maxm=40)


@pytest.mark.parametrize("v", [julia_live_param()])
def test_ladder_names_equal_sx_sy_combinations(v):
    """Sp/Sm, which ITensors.jl does not name on S=1/2, and ISy are the
    combinations Sx+iSy, Sx-iSy and i*Sy of sites_custom.jl, measured on
    the same julia_live state"""
    sc = spinchain.Spin_Chain(["1/2", "1", "3/2", "1/2"])
    sc.setup_julia()
    sc.maxm = 30
    sc.set_hamiltonian(_heisenberg(sc))
    for i, j in ((0, 1), (1, 2), (2, 3)):
        spi = sc.get_operator("Sp", i)
        smj = sc.get_operator("Sm", j)
        ref = (sc.Sx[i] + 1j*sc.Sy[i])*(sc.Sx[j] - 1j*sc.Sy[j])
        assert abs(sc.vev(spi*smj) - sc.vev(ref)) < 1e-8, (i, j)
        isy = sc.get_operator("ISy", i)
        assert abs(sc.vev(isy*sc.Sx[j]) - 1j*sc.vev(sc.Sy[i]*sc.Sx[j])) < 1e-8


@pytest.mark.parametrize("v", [julia_live_param()])
def test_unknown_site_code_raises_naming_it(v):
    """S=3 is site-type code 7, which no DMRG backend builds"""
    with pytest.raises(Exception, match="code 7"):
        spinchain.Spin_Chain(["S=3"]*3, itensor_version="julia_live")
