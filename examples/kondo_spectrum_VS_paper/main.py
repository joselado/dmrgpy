# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression example for Spin_Chain.get_kondo_spectrum (kondospectrumtk/),
# the third-order STM/Kondo perturbation-theory dI/dV spectrum of Ternes,
# New J. Phys. 17 063016 (2015), arXiv:1505.04430, for a single S=1/2
# impurity spin in a Zeeman field -- the paper's own simplest worked
# example (Figs. 2 and 3).
#
# Two of the paper's own closed-form equations for the numerical building
# blocks (eq. "step-fkt" for the temperature-broadened step Theta(x), and
# equ. "F_2" for the temperature-broadened Kondo log function F(eps,T))
# turned out not to reproduce the physics/figures the paper itself
# describes (Theta(x) as printed diverges as x->-inf instead of
# saturating at 0; F(eps,T) as printed drops the temperature-dependent
# broadening the surrounding text and Fig. F describe). kondospectrumtk/
# uses corrected forms instead, re-derived directly from the paper's own
# unambiguous defining equations (eq. "current" for Theta; equ. "F_1" for
# F) -- see kondospectrumtk/stepfunctions.py's module docstring for the
# full derivation notes. This example checks those corrected pieces
# quantitatively against the paper's own plotted values (digitized by eye
# from Fig. F(b)), and the assembled spectrum against the qualitative
# shapes described in the text and shown in Figs. 2 and 3.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.stepfunctions import F

G = 2.0
MUB = 5.7883818066e-5 # eV/T

# --- F(eps,T) vs Fig. F(b): T=1K, omega0=200meV, Gamma0=5ueV --------------
vals = F(np.array([0., 10e-3]), T=1.0, omega0=0.2, Gamma0=5e-6)
assert abs(vals[0] - 7.5) < 0.1, vals[0] # peak height read off the figure
assert abs(vals[1] - 3.1) < 0.2, vals[1] # value at eV-eps_m=10meV

# --- second order: S=1/2 Zeeman step, T=1K, g=2 (Fig. 2) ------------------
def make_chain(B):
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*B*sc.Sz[0])
    return sc

B = 10.0 # Tesla
zeeman = G*MUB*B
sc = make_chain(B)
eVs = np.linspace(-3*zeeman, 3*zeeman, 121)
eV, dIdV2 = sc.get_kondo_spectrum(eVs, site=0, T=1.0, order=2)

# symmetric temperature-broadened step centered at +-the Zeeman energy:
# low plateau (only the elastic channel open) well inside +-zeeman, high
# plateau (elastic + inelastic both open) well outside it
inside = np.abs(eV) < 0.3*zeeman
outside = np.abs(eV) > 2*zeeman
assert np.allclose(dIdV2, dIdV2[::-1], atol=1e-6) # symmetric in eV
assert np.mean(dIdV2[outside]) > 2*np.mean(dIdV2[inside])

# --- third order Kondo term: zero-field zero-bias resonance (Fig. 3a/b) ---
sc0 = make_chain(0.0)
span = 3e-3
eVs3 = np.linspace(-span, span, 61)
_, d2_only = sc0.get_kondo_spectrum(eVs3, site=0, T=1.0, order=2)
_, d3_total = sc0.get_kondo_spectrum(eVs3, site=0, Jrho_s=-0.05, T=1.0, order=3)
# with U=0 second order is exactly flat at B=0 (fully degenerate levels);
# the third-order Kondo term must add real structure (a zero-bias feature)
assert np.allclose(d2_only, d2_only[0])
assert np.max(np.abs(d3_total - d2_only)) > 0.05*d2_only[0]

# --- third order potential-interference term: bias asymmetry (Fig. 3c/d) --
# Isolate eq. "U-M"'s own contribution (d_U - d_noU): the third-order
# Kondo term alone (d_noU) is t->s current only (see
# third_order_kondo_dIdV's docstring) and so is not itself perfectly
# symmetric in eV either, which would confound a direct d_noU-vs-d_U
# comparison -- the potential-interference term specifically is what eq.
# "U-M" says should be bias-asymmetric.
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
scB = make_chain(10.0)
eVs4 = np.linspace(-3e-3, 3e-3, 41)
ksB = KondoSpectrum(scB, site=0, T=1.0)
dU_only = third_order_potential_dIdV(ksB, eVs4, Jrho_s=-0.05, U=0.25, T0=1.0)
assert np.allclose(third_order_potential_dIdV(ksB, eVs4, Jrho_s=-0.05, U=0.0, T0=1.0), 0.)
assert not np.allclose(dU_only, dU_only[::-1], atol=1e-8) # bias-asymmetric

print("All Kondo-spectrum checks against the paper's Figs. F/2/3 passed.")

# --- T=0: exact closed-form limit, consistent with small-but-finite T ------
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum as KS0
scB2 = make_chain(10.0)
eVs5 = np.linspace(-2e-3, 2e-3, 21)
_, dIdV_T0 = scB2.get_kondo_spectrum(eVs5, site=0, Jrho_s=-0.05, T=0.0, order=3)
_, dIdV_smallT = scB2.get_kondo_spectrum(eVs5, site=0, Jrho_s=-0.05, T=1e-3, order=3)
assert np.max(np.abs(dIdV_T0 - dIdV_smallT)) < 1e-2

# --- T=0 third-order Kondo term via the excited-state-free two-time
# construction (kondospectrumtk/twotime.py + edtwotimeref.py), the
# building block a DMRG (itensor_version=3) implementation reuses on top
# of real TDVP time evolution instead of ED's eigenbasis-exact evolution
# -- checked here against the ordinary excited-state-sum reference, using
# deliberately modest/fast grid parameters (see that module's own tests
# for the resolution/accuracy tradeoffs).
from dmrgpy.kondospectrumtk.conductance import third_order_kondo_dIdV as tokd
from dmrgpy.kondospectrumtk.edtwotimeref import two_time_kondo_term_ed

scC = make_chain(10.0)
ksC = KS0(scC, site=0, T=0.0)
eVs6 = np.linspace(-2e-3, 2e-3, 9)
omega0_tt, Gamma0_tt = 2e-3, 5e-6
twotime_term = 4*np.pi*1.0**2*(-0.05)*two_time_kondo_term_ed(
        ksC, eVs6, omega0=omega0_tt, Gamma0=Gamma0_tt,
        t2_width=25/Gamma0_tt, t2_npts=40_000, t2_batch=10_000,
        tau_width=2*np.pi/2e-5, tau_npts=1_000)
excited_state_sum_term = tokd(ksC, eVs6, -0.05, T0=1.0, omega0=omega0_tt,
                               Gamma0=Gamma0_tt)
assert np.max(np.abs(twotime_term - excited_state_sum_term)) < 0.02

print("All T=0 / two-time-correlator checks passed.")
