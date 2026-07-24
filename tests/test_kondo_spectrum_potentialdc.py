"""Coverage for the T=0 third-order potential-interference dI/dV computed
via the dynamical correlator (kondospectrumtk/potentialdc.py) instead of
the explicit excited-state sum (conductance.third_order_potential_dIdV),
following Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.

At T=0 the excited-state sum collapses to a convolution of the T=0
dynamical structure factor against the F0 kernel (see potentialdc.py's
module docstring) -- this is the DMRG-side route to this term
(mode="DMRG", submode="KPM"/"CVM") that needs no diagonalization beyond
the ground state, filling the gap Spin_Chain.get_kondo_spectrum's
mode="DMRG" path previously raised NotImplementedError for. Tested here
with mode="ED", submode="ED" (exact except for the numerical `delta`
broadening), since that is what's actually runnable without a compiled
C++ DMRG backend -- but the function itself is backend-agnostic.
"""
import numpy as np

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc

G = 2.0
MUB = 5.7883818066e-5 # eV/T


def test_third_order_potential_dc_matches_excited_state_sum():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    ks = KondoSpectrum(sc, site=0, T=0.0)

    eVs = np.linspace(-1e-3, 2e-3, 21)
    delta = 2e-6
    es = np.linspace(-1e-3, 3e-3, 40_000) # must cover the Zeeman gap (~1.16e-3)
    Jrho_s, U = 0.1, 0.3
    mine = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, T0=1.0,
                                          mode="ED", submode="ED",
                                          delta=delta, es=es)
    ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)

    assert np.max(np.abs(mine-ref)) < 0.01
    assert np.max(np.abs(mine-ref)) < 0.01*np.max(np.abs(ref))


def test_third_order_potential_dc_requires_es():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    import pytest
    with pytest.raises(ValueError):
        third_order_potential_dIdV_dc(sc, 0, np.array([0.0]), 0.1, 0.3,
                                       mode="ED", submode="ED")
