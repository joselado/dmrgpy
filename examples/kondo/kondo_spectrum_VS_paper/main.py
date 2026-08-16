# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

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
#
# Figs. 3b/3d are plotted in ABSOLUTE units (e^2 T0^2/h), which makes them
# a much stronger check than the shapes alone: their zero-bias peak values
# (1.13 at U=0, 1.39 at U=0.25) simultaneously pin down every prefactor in
# the third-order spectrum. Asserting on them here caught three real
# normalization/structure bugs at once -- see kondospectrumtk/
# conductance.py's module docstring for what they were.
import numpy as np
import matplotlib.pyplot as plt
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
# the third-order Kondo term must add a resonance centered exactly at
# eV=0 and symmetric in bias -- both tunneling directions contribute with
# the same sign (eq. "sym_z"), and every per-process curve of Fig. 3a is
# drawn symmetric about zero bias
assert np.allclose(d2_only, d2_only[0])
assert np.max(np.abs(d3_total - d2_only)) > 0.05*d2_only[0]
assert np.allclose(d3_total, d3_total[::-1], atol=1e-10) # even in eV
assert np.argmax(d3_total) == len(eVs3)//2 # peak exactly at zero bias

# --- third order potential-interference term: bias asymmetry (Fig. 3c/d) --
# Isolate eq. "U-M"'s own contribution: its two tunneling directions enter
# with opposite signs (eq. "asym_U"), so unlike every other term in the
# spectrum it is ODD in eV -- which is exactly the antisymmetric "sum"
# curve of Fig. 3c, and the sole source of bias asymmetry in Fig. 3d.
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
scB = make_chain(10.0)
eVs4 = np.linspace(-3e-3, 3e-3, 41)
ksB = KondoSpectrum(scB, site=0, T=1.0)
dU_only = third_order_potential_dIdV(ksB, eVs4, Jrho_s=-0.05, U=0.25, T0=1.0)
assert np.allclose(third_order_potential_dIdV(ksB, eVs4, Jrho_s=-0.05, U=0.0, T0=1.0), 0.)
assert not np.allclose(dU_only, dU_only[::-1], atol=1e-8) # bias-asymmetric
assert np.allclose(dU_only, -dU_only[::-1], atol=1e-10) # specifically odd

# --- absolute amplitudes vs Figs. 3b/3d ------------------------------------
# Figs. 3b/3d are plotted in absolute units (e^2 T0^2/h) for a single
# S=1/2 at B=0, T=1K, Jrho_s=-0.05, U=0 resp. 0.25 and the paper-default
# omega0=20meV; get_kondo_spectrum returns 2*pi times those units (its
# second-order plateau is 2*pi*0.75 where Fig. 3b's dashed curve reads
# 0.75). Reading the zero-bias peaks off the published figures gives 1.13
# and 1.39. These pin down the three prefactors nothing else here
# constrains: the Levi-Civita coefficient Im[X]/2, the sum over both
# tunneling directions, and the elastic potential channel's 4*|U|^2.
zero = np.array([0.0])
_, pk_b = sc0.get_kondo_spectrum(zero, site=0, Jrho_s=-0.05, U=0.0, T=1.0,
                                  order=3, omega0=20e-3)
_, pk_d = sc0.get_kondo_spectrum(zero, site=0, Jrho_s=-0.05, U=0.25, T=1.0,
                                  order=3, omega0=20e-3)
print("Fig. 3b zero-bias peak: %.3f (paper: 1.13)"%(pk_b[0]/(2*np.pi)))
print("Fig. 3d zero-bias peak: %.3f (paper: 1.39)"%(pk_d[0]/(2*np.pi)))
assert abs(pk_b[0]/(2*np.pi) - 1.13) < 0.02
assert abs(pk_d[0]/(2*np.pi) - 1.39) < 0.02

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

fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))

axes[0, 0].plot(eV, dIdV2, color="tab:blue")
axes[0, 0].axvline(zeeman, color="gray", ls=":")
axes[0, 0].axvline(-zeeman, color="gray", ls=":")
axes[0, 0].set_xlabel("eV")
axes[0, 0].set_ylabel("dI/dV (2nd order)")
axes[0, 0].set_title("Zeeman step (Fig. 2)")

# Fig. 3b itself: the zero-field Kondo peak splitting with magnetic field
eVs_fig3 = np.linspace(-4e-3, 4e-3, 161)
for Bfield in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]:
    _, dfig = make_chain(Bfield).get_kondo_spectrum(
            eVs_fig3, site=0, Jrho_s=-0.05, U=0.0, T=1.0, order=3, omega0=20e-3)
    axes[0, 1].plot(eVs_fig3*1e3, dfig/(2*np.pi), label="B=%g T"%Bfield)
_, dfig2 = make_chain(10.0).get_kondo_spectrum(eVs_fig3, site=0, T=1.0, order=2)
axes[0, 1].plot(eVs_fig3*1e3, dfig2/(2*np.pi), "k--", lw=1,
                 label="2nd order, B=10 T")
axes[0, 1].axhline(1.13, color="gray", lw=0.8, ls=":")
axes[0, 1].set_xlabel("eV (meV)")
axes[0, 1].set_ylabel(r"dI/dV ($e^2T_0^2/h$)")
axes[0, 1].set_title("Kondo peak splitting (Fig. 3b);\ndotted: paper's B=0 peak 1.13")
axes[0, 1].legend(fontsize=6)

axes[0, 2].plot(eVs4, dU_only, color="tab:red")
axes[0, 2].axhline(0, color="gray", lw=0.8, ls=":")
axes[0, 2].set_xlabel("eV")
axes[0, 2].set_ylabel("dI/dV (potential-interference term)")
axes[0, 2].set_title("bias asymmetry (Fig. 3c/d)")

axes[1, 0].plot(eVs5*1e3, dIdV_T0, "o-", label="T=0 (exact)")
axes[1, 0].plot(eVs5*1e3, dIdV_smallT, "x--", label="T=1e-3 (finite)")
axes[1, 0].set_xlabel("eV (meV)")
axes[1, 0].set_ylabel("dI/dV (3rd order)")
axes[1, 0].set_title("T=0 exact vs small-T limit")
axes[1, 0].legend(fontsize=8)

axes[1, 1].plot(eVs6*1e3, twotime_term, "o-", label="two-time (T=0, TDVP-style)")
axes[1, 1].plot(eVs6*1e3, excited_state_sum_term, "x--", label="excited-state sum")
axes[1, 1].set_xlabel("eV (meV)")
axes[1, 1].set_ylabel("3rd-order Kondo term")
axes[1, 1].set_title("two-time vs excited-state-sum cross-check")
axes[1, 1].legend(fontsize=8)

eps_grid = np.linspace(-10e-3, 10e-3, 200)
axes[1, 2].plot(eps_grid*1e3, F(eps_grid, T=1.0, omega0=0.2, Gamma0=5e-6), color="tab:purple")
axes[1, 2].plot([0, 10], [7.5, 3.1], "o", color="black", label="Fig. F(b) reference points")
axes[1, 2].set_xlabel("eV-eps_m (meV)")
axes[1, 2].set_ylabel("F(eps,T)")
axes[1, 2].set_title("temperature-broadened Kondo log F (Fig. F(b))")
axes[1, 2].legend(fontsize=8)

fig.suptitle("get_kondo_spectrum vs Ternes, New J. Phys. 17 063016 (2015)")
fig.tight_layout()
fig.savefig("kondo_spectrum_VS_paper.png", dpi=150)
print("saved plot to kondo_spectrum_VS_paper.png")
