# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#####################################################################
### Spin-1 chains: the dynamical structure factor S(k,w) of the    ###
### Haldane phase, and its exactly solvable AKLT point.            ###
#####################################################################
# Infinite_Many_Body_Chain.spectral_weights gives each quasiparticle
# branch's exact delta-peak residue |<k,a|O(k)|Psi>|^2. This example
# runs it on the two canonical S=1 chains, because between them they
# pin down every part of the machinery against something exact:
#
#   AKLT   H = sum_j [ S_j.S_{j+1} + (1/3)(S_j.S_{j+1})^2 ]
#          The valence-bond-solid ground state is EXACTLY a bond
#          dimension 2 MPS, so at D=2 there is no variational error at
#          all and three closed forms are available to check against:
#            * the static structure factor, from <S_0.S_r> = 4(-1/3)^r;
#            * Arovas, Auerbach & Haldane's single-mode dispersion
#              w(k) = (5/27)(5+3cos k) [PRL 60, 531 (1988)], which the
#              FIRST MOMENT sum_a E_a w_a / total must reproduce exactly,
#              since S^z_k|Psi> lies entirely inside the tangent space;
#            * the SU(2) content: the D^2(d-1) = 8 branches split into a
#              magnon TRIPLET and a quintuplet at every momentum, and a
#              spin-1 operator reaches only the triplet.
#
#   Haldane  H = sum_j S_j.S_{j+1}   (the pure S=1 Heisenberg chain)
#          No closed form, but sharp literature numbers: gap
#          Delta = 0.410479 at k=pi, e0/site = -1.401484, correlation
#          length xi ~ 6. The single-mode ansatz is known to be very good
#          near k=pi and to degrade away from it, which the sum-rule
#          fraction measures directly.
#
# Note the AKLT biquadratic term below is the FULL (S.S)^2 =
# sum_{ab} S^a_i S^b_i S^a_j S^b_j, not the diagonal-only
# sum_a (S^a_i S^a_j)^2 of examples/spin_models/bilinear_biquadratic --
# only the former has the exact valence-bond-solid ground state.

D_AKLT = 2      # exact for AKLT
D_HALDANE = 10  # xi ~ 6, so this is converged to ~1e-3 on the gap


def spin1_chain(biquadratic, D):
    ic = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
    Sc = [ic.SxC[0], ic.SyC[0], ic.SzC[0]]
    Sr = [ic.SxR[0], ic.SyR[0], ic.SzR[0]]
    h = Sc[0]*Sr[0] + Sc[1]*Sr[1] + Sc[2]*Sr[2]
    if biquadratic != 0.0:
        for a in range(3):
            for b in range(3):
                h = h + biquadratic*Sc[a]*Sc[b]*Sr[a]*Sr[b]
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"     # the excitation ansatz needs the mixed gauge
    ic.maxm = D
    ic.maxiter = 800
    ic.vumps_nrestarts = 6
    return ic


def aklt_szz(k):
    """Exact AKLT sum_r e^{ikr} <S^z_0 S^z_r>, from <S_0.S_r>=4(-1/3)^r."""
    x = -1./3.
    return 2./3. + (8./3.)*(x*np.cos(k) - x*x)/(1 - 2*x*np.cos(k) + x*x)


def aklt_sma(k):
    """Arovas-Auerbach-Haldane, w(k)=(5/27)(5+3cos k) for H=sum_j P^(2);
    the Hamiltonian above is 2*sum_j P^(2) + const, hence the factor 2."""
    return 2.0*(5./27.)*(5. + 3.*np.cos(k))


# =====================  AKLT: exact at D=2  ==========================
ic = spin1_chain(1./3., D_AKLT)
print("AKLT   e0 = %.12f   (exact -2/3 = %.12f)  converged=%s"
      % (ic.gs_energy(), -2./3., ic.converged))

ks = np.linspace(0.05, np.pi, 40)
dim = D_AKLT*D_AKLT*(3-1)
aklt_E, aklt_tot, aklt_mom, aklt_trip = [], [], [], []
for k in ks:
    e, w, tot = ic.spectral_weights("Sz", k, n=dim, return_total=True)
    # The triplet is identified by ENERGY, not by being lowest: the two
    # multiplets cross near k~0.9, below which the quintuplet is lower.
    trip = np.abs(e - aklt_sma(k)) < 1e-6
    aklt_E.append(e)
    aklt_tot.append(tot)
    aklt_mom.append(np.sum(e*w)/tot)
    aklt_trip.append(w[trip].sum()/tot)
aklt_E = np.array(aklt_E)
aklt_tot, aklt_mom = np.array(aklt_tot), np.array(aklt_mom)
aklt_trip = np.array(aklt_trip)

print("  max |S^zz(k) - exact|                = %.2e" % np.max(np.abs(aklt_tot - aklt_szz(ks))))
print("  max |first moment - AAH SMA|         = %.2e" % np.max(np.abs(aklt_mom - aklt_sma(ks))))
print("  min triplet-summed weight fraction   = %.12f" % aklt_trip.min())
assert np.max(np.abs(aklt_tot - aklt_szz(ks))) < 1e-8
assert np.max(np.abs(aklt_mom - aklt_sma(ks))) < 1e-6
assert aklt_trip.min() > 1.0 - 1e-9   # a spin-1 operator reaches only the triplet


# ==============  Haldane: the pure S=1 Heisenberg chain  =============
ich = spin1_chain(0.0, D_HALDANE)
print()
# `converged` often comes back False here, and that is expected rather than a
# failure: VUMPS's own D>1 convergence robustness is a documented open
# limitation (pyitensor/vumps.py's "Convergence robustness" docstring
# section), and the flag reports that its residual never dropped below the
# internal tolerance, not that the state is wrong -- e0 lands within ~3e-4 of
# the literature value and the sum rule below holds to ~1e-7 regardless. The
# sum rule in particular is an identity of the ANSATZ, not of the ground
# state's accuracy, so it holds even where the ground state is poor (checked
# down to D=4, where the gap is off by 50%).
print("Haldane e0 = %.10f   (literature -1.401484038971)  converged=%s"
      % (ich.gs_energy(), ich.converged))

ksh = np.linspace(0.0, 2*np.pi, 49)
dimh = min(12, D_HALDANE*D_HALDANE*2)
all_e, all_w, tot_h = [], [], []
for k in ksh:
    e, w, tot = ich.spectral_weights("Sz", k, n=dimh, return_total=True)
    all_e.append(e) ; all_w.append(w) ; tot_h.append(tot)
all_e, all_w, tot_h = np.array(all_e), np.array(all_w), np.array(tot_h)
band = all_e[:, 0]
# The lowest three branches ARE the magnon triplet here (unlike AKLT, which
# has a multiplet crossing). Their individual weights are basis-arbitrary
# within the near-degenerate multiplet, so only the SUM is physical -- see
# spectral_weights' own docstring, and note the splitting is set by the
# variational error, so do NOT group them by an energy tolerance.
frac_h = np.where(tot_h > 1e-12, all_w[:, :3].sum(axis=1)/np.maximum(tot_h, 1e-300), np.nan)

ipi = int(np.argmin(np.abs(ksh - np.pi)))
gap = band[ipi]
print("  gap E(pi)   = %.8f   (literature Haldane gap 0.4104789, err %+.1e)"
      % (gap, gap - 0.4104789))
print("  S^zz(pi)    = %.6f" % tot_h[ipi])
print("  triplet exhausts %.2f%% of the k=pi sum rule" % (100*frac_h[ipi]))

# Independent cross-check of the sum rule against a real-space correlator
# sum -- machinery that shares no code at all with the excitation ansatz.
# 120 sites, not 80: the Haldane correlation length is xi ~ 6, so the tail
# beyond r=80 is still ~1e-6 of the total -- borderline against the
# tolerance below, while r=120 leaves ~1e-9.
mean = ich.vev("Sz", 0).real
C = [ich.correlator("Sz", 0, "Sz", r).real - mean*mean for r in range(120)]
ref = C[0] + 2.0*sum(np.cos(np.pi*r)*C[r] for r in range(1, 120))
print("  S^zz(pi) from the real-space correlator sum = %.6f  (diff %.1e)"
      % (ref, abs(ref - tot_h[ipi])))
assert abs(ref - tot_h[ipi]) < 1e-5
assert abs(gap - 0.4104789) < 0.01      # D=10 with xi~6
assert frac_h[ipi] > 0.90               # single-mode is excellent at k=pi

# ------- S^zz(k,w): broadened from those same exact delta peaks ---------
delta = 0.08
es = np.linspace(0.0, 6.0, 400)
lor = (delta/np.pi)/((es[None, :, None] - all_e[:, None, :])**2 + delta**2)
S_map = np.einsum('kea,ka->ke', lor, all_w)

# dynamical_structure_factor is the convenience wrapper that does exactly
# this; run it on a coarse grid to demonstrate the public API and confirm
# the two agree, rather than paying for the whole fine grid twice.
kc = ksh[::8]
_, ec, Sc = ich.dynamical_structure_factor("Sz", ks=kc, energies=es,
                                            delta=delta, n=dimh)
assert np.max(np.abs(Sc - S_map[::8])) < 1e-10

# Where does the weight the magnon branch does NOT carry actually sit? Not
# at the two-magnon threshold 2*Delta, even though that is where the true
# continuum starts: a single-B tangent-space ansatz cannot represent two
# well-separated magnons at all (the multi-particle extension is Osborne &
# McCulloch, arXiv:2408.17117, not implemented here). Measured directly at
# k=pi with EVERY branch enumerated at D=8: the magnon triplet takes 97.6%,
# and the remaining 2.4% sits between w=2.7 and w=7.7 -- inside the
# two-magnon band [2*Delta, 2*max eps] = [0.82, ~5.0] but concentrated near
# its top, with nothing at all near its lower edge. So read the map as a
# single-mode result: quantitative on the magnon branch, and silent (not
# zero -- silent) about the continuum.
two_delta = 2*0.4104789
carried = all_w[:, :3].sum(axis=1)
print("  weight above the magnon branch, summed over k: %.4f of the total"
      % ((tot_h - carried).sum()/tot_h.sum()))
print("  S^zz(k) dynamic range across the zone: %.2e (k=0) .. %.3f (k=pi), factor %.1e"
      % (tot_h.min(), tot_h.max(), tot_h.max()/max(tot_h.min(), 1e-300)))
# k=0 must be an EXACT zero -- sum_j S^z_j is conserved. What is left is the
# ground state's own convergence error, not physics.
assert tot_h[0] < 1e-4

# ==========================  plots  ==================================
fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.4))

ax = axes[0, 0]
for a in range(dim):
    ax.plot(ks, aklt_E[:, a], "o", ms=3, color="C0",
            label="tangent-space branches" if a == 0 else None)
ax.plot(ks, aklt_sma(ks), "k--", lw=1.4,
        label=r"AAH single mode $\frac{10}{27}(5+3\cos k)$")
ax.set_xlabel("momentum $k$")
ax.set_ylabel("$E(k)$")
ax.set_title(r"AKLT at $D=2$ (exact): triplet + quintuplet, crossing at $k\approx0.6$")
ax.legend(fontsize=8)

ax = axes[0, 1]
ax.plot(ks, aklt_tot, "o", ms=4, label=r"AKLT, dmrgpy $\sum_a w_a(k)$")
ax.plot(ks, aklt_szz(ks), "k--", lw=1.4, label=r"AKLT, exact $S^{zz}(k)$")
ax.plot(ksh, tot_h, "s", ms=4, color="C3", label=r"Haldane, $D=%d$" % D_HALDANE)
ax.set_xlim(0, np.pi)
ax.set_xlabel("momentum $k$")
ax.set_ylabel(r"$S^{zz}(k)$ (sum rule total)")
ax.set_title("total spectral weight vs the exact structure factor")
ax.legend(fontsize=8)

ax = axes[1, 0]
# LOG colour scale, deliberately. S^zz(k) spans 5.3e5 between k=pi and k=0
# (measured: 3.31 vs 6.3e-6), so on a linear scale everything below ~1% of
# the k=pi peak is black and the magnon looks like it exists only near the
# zone boundary. It does not: the branch E_1(k) is well defined at every k
# (1.22 at k=0, 2.72 at k=pi/2), it is the matrix element that collapses --
# and at k=0 it collapses to an EXACT zero, since S^z(k=0) = sum_j S^z_j is
# a conserved charge and cannot create any excitation at all (the 6.3e-6 we
# see is that zero within this D=10 state's convergence). A linear scale is
# the honest picture of where scattering intensity actually is; a log scale
# is the honest picture of what the calculation produced, and it also makes
# the second branch near w~3 visible, which linear hides entirely.
im = ax.pcolormesh(ksh, es, np.maximum(S_map.T, 1e-12), shading="auto",
                   cmap="magma", norm=LogNorm(vmin=1e-4, vmax=S_map.max()))
ax.plot(ksh, band, "w--", lw=1.0, label="magnon branch (this calculation)")
ax.axhline(two_delta, color="cyan", lw=1.0, ls=":",
           label=r"$2\Delta$: true continuum edge, carries no ansatz weight")
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax.set_xlabel("momentum $k$")
ax.set_ylabel(r"energy $\omega$")
ax.set_title(r"$S^{zz}(k,\omega)$, Haldane chain, $D=%d$" % D_HALDANE)
ax.legend(loc="upper center", fontsize=7.5, framealpha=0.75)
fig.colorbar(im, ax=ax)

ax = axes[1, 1]
ax.plot(ksh, frac_h, "s-", ms=4, color="C3", label="Haldane (lowest triplet)")
ax.plot(ks, aklt_trip, "o-", ms=3, color="C0", label="AKLT (triplet, exact 1)")
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_xlim(0, np.pi)
ax.set_ylim(0.0, 1.08)
ax.set_xlabel("momentum $k$")
ax.set_ylabel("triplet weight / sum rule total")
ax.set_title(r"how much of $S^{zz}(k)$ the single mode carries")
ax.legend(fontsize=8, loc="lower right")

fig.tight_layout()
fig.savefig("haldane_structure_factor.png", dpi=150)
print("\nSaved plot to haldane_structure_factor.png")
