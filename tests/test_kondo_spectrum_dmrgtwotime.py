"""Coverage for the DMRG (itensor_version=3) two-time-correlator
construction of the third-order Kondo term
(kondospectrumtk/dmrgtwotime.py), following Ternes, New J. Phys. 17
063016 (2015), arXiv:1505.04430.

IMPORTANT: dmrgtwotime.py was developed and written against the same
verified DMRG API already used elsewhere in this codebase (toMPO,
tdvp_step, MPS operator application/overlap -- see that module's own
docstring for the exact API citations) but could NOT be executed or
numerically validated in the environment it was developed in (no
compiled C++ backend: cppext.available(3) is False there). This test is
therefore the actual first real check of that code -- skipped entirely
when no compiled itensor_version=3 backend is available, rather than
silently claiming coverage that was never exercised.

Compares against the fully-validated ED two-time reference
(edtwotimeref.py, itself checked to ~0.01-1%% against the excited-state-
sum third_order_kondo_dIdV) on the smallest nontrivial system (a single
S=1/2 impurity), with deliberately modest, fast grid parameters.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy import cppext
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.edtwotimeref import two_time_kondo_term_ed

pytestmark = pytest.mark.skipif(
    not cppext.available(3),
    reason="the DMRG two-time construction needs a compiled itensor_version=3 backend",
)

G = 2.0
MUB = 5.7883818066e-5 # eV/T


def test_two_time_kondo_term_dmrg_matches_ed_reference():
    from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg

    sc = spinchain.Spin_Chain(["1/2"], itensor_version=3)
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    sc.get_gs() # ensure the ground state (and self.e0) are computed

    ks = KondoSpectrum(sc, site=0, T=0.0)

    Jrho_s, T0 = -0.05, 1.0
    omega0, Gamma0 = 2e-3, 5e-6
    eVs = np.linspace(-2e-3, 2e-3, 9)

    dmrg_term = two_time_kondo_term_dmrg(
            sc, 0, eVs, omega0=omega0, Gamma0=Gamma0,
            dt2=25./Gamma0/200, n_t2_half=200,
            dtau=(2*np.pi/2e-5)/500, n_tau_half=500)
    mine = 4*np.pi*T0**2*Jrho_s*dmrg_term

    ed_term = two_time_kondo_term_ed(
            ks, eVs, omega0=omega0, Gamma0=Gamma0,
            t2_width=25/Gamma0, t2_npts=40_000, t2_batch=10_000,
            tau_width=2*np.pi/2e-5, tau_npts=1_000)
    ref = 4*np.pi*T0**2*Jrho_s*ed_term

    assert np.max(np.abs(mine-ref)) < 0.1
