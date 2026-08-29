# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
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

ksh = np.linspace(0.1, np.pi, 33)
dimh = min(12, D_HALDANE*D_HALDANE*2)
band, frac_h, tot_h = [], [], []
for k in ksh:
    e, w, tot = ich.spectral_weights("Sz", k, n=dimh, return_total=True)
    band.append(e[0])
    # The lowest three branches ARE the magnon triplet here (unlike AKLT,
    # no crossing in this window); their weights are individually
    # basis-arbitrary within the near-degenerate multiplet, so only their
    # SUM is physical -- see spectral_weights' own docstring.
    frac_h.append(w[:3].sum()/tot)
    tot_h.append(tot)
band, frac_h, tot_h = np.array(band), np.array(frac_h), np.array(tot_h)

gap = band[-1]
print("  gap E(pi)   = %.8f   (literature Haldane gap 0.4104789, err %+.1e)"
      % (gap, gap - 0.4104789))
print("  S^zz(pi)    = %.6f" % tot_h[-1])
print("  triplet exhausts %.2f%% of the k=pi sum rule" % (100*frac_h[-1]))

# Independent cross-check of the sum rule against a real-space correlator
# sum -- machinery that shares no code at all with the excitation ansatz.
# 120 sites, not 80: the Haldane correlation length is xi ~ 6, so the tail
# beyond r=80 is still ~1e-6 of the total -- borderline against the
# tolerance below, while r=120 leaves ~1e-9.
mean = ich.vev("Sz", 0).real
C = [ich.correlator("Sz", 0, "Sz", r).real - mean*mean for r in range(120)]
ref = C[0] + 2.0*sum(np.cos(np.pi*r)*C[r] for r in range(1, 120))
print("  S^zz(pi) from the real-space correlator sum = %.6f  (diff %.1e)"
      % (ref, abs(ref - tot_h[-1])))
assert abs(ref - tot_h[-1]) < 1e-5
assert abs(gap - 0.4104789) < 0.01      # D=10 with xi~6
assert frac_h[-1] > 0.90                # single-mode is excellent at k=pi


# ==========================  plots  ==================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
for a in range(dim):
    ax.plot(ks, aklt_E[:, a], "o", ms=3, color="C0",
            label="tangent-space branches" if a == 0 else None)
ax.plot(ks, aklt_sma(ks), "k--", lw=1.4,
        label=r"AAH single mode $\frac{10}{27}(5+3\cos k)$")
ax.set_xlabel("momentum $k$")
ax.set_ylabel("$E(k)$")
ax.set_title(r"AKLT at $D=2$ (exact): triplet + quintuplet")
ax.legend(fontsize=8)

ax = axes[1]
ax.plot(ks, aklt_tot, "o", ms=4, label=r"dmrgpy $\sum_a w_a(k)$")
ax.plot(ks, aklt_szz(ks), "k--", lw=1.4, label=r"exact $S^{zz}(k)$, AKLT")
ax.plot(ksh, tot_h, "s", ms=4, color="C3", label=r"Haldane, $D=%d$" % D_HALDANE)
ax.set_xlabel("momentum $k$")
ax.set_ylabel(r"$S^{zz}(k)$ (sum rule total)")
ax.set_title("total spectral weight vs the exact structure factor")
ax.legend(fontsize=8)

ax = axes[2]
ax.plot(ksh, frac_h, "s-", ms=4, color="C3", label="Haldane (lowest triplet)")
ax.plot(ks, aklt_trip, "o-", ms=3, color="C0", label="AKLT (triplet, exact 1)")
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_ylim(0.0, 1.08)
ax.set_xlabel("momentum $k$")
ax.set_ylabel("triplet weight / sum rule total")
ax.set_title("how much of $S^{zz}(k)$ the single mode carries")
ax.legend(fontsize=8, loc="lower right")

fig.tight_layout()
fig.savefig("haldane_structure_factor.png", dpi=150)
print("\nSaved plot to haldane_structure_factor.png")
