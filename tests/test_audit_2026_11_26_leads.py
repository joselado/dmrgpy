"""Two leads left open by earlier audits, closed on 2026-11-26.

- edtk/edchain.py's EDchain.get_operator printed "Unrecognized operator in
  EDchain" and returned None on a type it cannot build, so the error came
  later and far away, in algebra.py (the 2026-09-25b record's companion
  tail at edchain.py:187). It now raises a TypeError naming the type, as
  the mixed-chain ED backend already did.
- On the lower-level get_dynamical_correlator_MB route a pair name= next
  to i=/j= was dropped silently, on ROOTN and CVM alike, where the public
  route raises (2026-09-25 open items, second pass finding 16).
"""
import numpy as np
import pytest

from dmrgpy import spinchain, bosonchain


def _heis(sc):
    h = 0
    for k in range(sc.ns - 1):
        h = h + sc.Sx[k]*sc.Sx[k+1] + sc.Sy[k]*sc.Sy[k+1] + sc.Sz[k]*sc.Sz[k+1]
    return h


def test_ed_get_operator_refuses_an_unknown_type():
    bc = bosonchain.Bosonic_Chain(3, maxnb=[3]*3)
    bc.set_hamiltonian(bc.A[0]*bc.Adag[1] + bc.A[1]*bc.Adag[0] + bc.N[2])
    ed = bc.get_ED_obj()
    with pytest.raises(TypeError, match="ndarray"):
        ed.get_operator(np.zeros(3))


@pytest.mark.parametrize("submode", ["ROOTN", "CVM"])
def test_lower_level_route_refuses_sites_next_to_a_pair(submode):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sc.set_hamiltonian(_heis(sc))
    with pytest.raises(TypeError, match="i=/j="):
        sc.get_dynamical_correlator_MB(submode=submode,
                                       name=(sc.Sz[1], sc.Sz[1]), i=1, j=1,
                                       es=np.linspace(0, 2, 3), delta=0.2)


def test_lower_level_route_keeps_sites_of_a_string_name():
    """a string name with i=/j= is still read, and agrees with the pair"""
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sc.set_hamiltonian(_heis(sc))
    kw = dict(submode="ROOTN", es=np.linspace(0, 2, 5), delta=0.2)
    y = np.asarray(sc.get_dynamical_correlator_MB(name="ZZ", i=1, j=1,
                                                  **kw)[1])
    y_11 = np.asarray(sc.get_dynamical_correlator_MB(
        name=(sc.Sz[1], sc.Sz[1]), **kw)[1])
    assert np.max(np.abs(y - y_11)) < 1e-10
