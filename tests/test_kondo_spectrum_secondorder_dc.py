"""Coverage for the T=0 second-order dI/dV computed via the dynamical
correlator (kondospectrumtk/secondorder_dc.py) instead of the explicit
excited-state sum (conductance.second_order_dIdV), following Ternes, New
J. Phys. 17 063016 (2015), arXiv:1505.04430.

At T=0 only the ground state is populated, so
sum_f |<f|S_alpha|GS>|^2 Theta0(eV-eps_f0) is exactly a Theta0-weighted
cumulative integral of the T=0 dynamical structure factor -- this exists
as the natural DMRG-side route to the second-order term (via
mode="DMRG", submode="KPM"/"CVM") that needs no diagonalization beyond
the ground state, unlike the excited-state-sum route. Tested here with
mode="ED", submode="ED" (exact except for the numerical `delta`
broadening), since that is what's actually runnable without a compiled
C++ DMRG backend -- but the function itself is backend-agnostic (any
mode/submode chain.get_dynamical_correlator accepts).
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import second_order_dIdV
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc

G = 2.0
MUB = 5.7883818066e-5 # eV/T


def test_second_order_dc_matches_excited_state_sum():
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    ks = KondoSpectrum(sc, site=0, T=0.0)

    eVs = np.linspace(-1e-3, 2e-3, 21)
    delta = 2e-6
    es = np.linspace(-1e-3, 3e-3, 40_000) # must cover the Zeeman gap (~1.16e-3)
    mine = second_order_dIdV_dc(sc, 0, eVs, T0=1.0, U=0.2, mode="ED",
                                 submode="ED", delta=delta, es=es)
    ref = second_order_dIdV(ks, eVs, T0=1.0, U=0.2)

    assert np.max(np.abs(mine-ref)) < 0.05
    assert np.max(np.abs(mine-ref)) < 0.01*np.max(np.abs(ref))


def test_second_order_dc_zero_at_U0_below_threshold():
    """Sanity check independent of the excited-state-sum reference: well
    below the Zeeman threshold and with no potential scattering, only the
    elastic channel is open, giving the same pi/2 plateau derived by hand
    in test_kondo_spectrum.py's analogous finite-T test."""
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    es = np.linspace(-1e-3, 3e-3, 40_000) # must cover the Zeeman gap (~1.16e-3)
    val = second_order_dIdV_dc(sc, 0, np.array([0.0]), T0=1.0, U=0.0,
                                mode="ED", submode="ED", delta=2e-6, es=es)[0]
    assert val == pytest.approx(np.pi/2, abs=0.03)
