# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#####################################################################
### Momentum-resolved dynamical structure factor S(k,w) directly  ###
### in the thermodynamic limit, from the quasiparticle branches'   ###
### own exact delta-peak spectral weights -- transverse-field      ###
### Ising model.                                                   ###
#####################################################################
# Infinite_Many_Body_Chain.spectral_weights returns, for each
# tangent-space excitation branch at momentum k, the residue
# |<k,a|O(k)|Psi>|^2 of the delta peak it contributes to S(k,w):
#
#     S(k,w) = sum_a weights[a] * delta(w - energies[a])
#
# Unlike kpm_finite/td_dynamical_correlator, which embed a finite window
# in the infinite chain, these peaks are exact in both momentum and
# energy -- no window, no broadening, no time truncation. What they miss
# is multi-particle continuum weight, and `return_total=True` measures
# exactly that: it gives the total over EVERY branch, which is the
# connected static structure factor sum_r e^{ikr}(<O_0 O_r> - <O>^2).
# The ratio weights.sum()/total is then the fraction of the response the
# single-mode picture actually accounts for.
#
# The physics that makes this a good demonstration: in the paramagnetic
# TFIM sigma^x is parity-ODD and creates exactly one quasiparticle, so
# its response IS the single coherent mode and the lowest branch carries
# ~100% of the sum rule. sigma^z is parity-EVEN and cannot reach that
# branch at all -- the selection rule is exact here, its weight on the
# lowest branch comes out at ~1e-21 -- so all of its response sits in
# the two-particle continuum, which the higher branches of the same
# H_eff(k) approximate. Both are computed below from the same ground
# state, and the continuum's exact free-fermion boundaries are overlaid
# on the branches that carry the sigma^z weight.
#
# Reference: Vanderstraeten, Haegeman & Verstraete, "Tangent-space
# methods for uniform matrix product states", arXiv:1810.07006,
# SciPost Phys. Lect. Notes 7 (2019), Sec. 6.

g = 2.5   # transverse field (g>1: gapped paramagnetic phase)
D = 4     # VUMPS bond dimension

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
# Sx,Sz have eigenvalues +-1/2, so -4*SxSx=-sigma^x sigma^x and
# -2*g*Sz=-g*sigma^z: the standard Pauli TFIM H=-sigma^x sigma^x - g sigma^z.
ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
ic.gs_method = "vumps"   # the excitation ansatz needs the mixed gauge
ic.maxm = D
ic.maxiter = 800
ic.vumps_nrestarts = 6
print("e0 = {:.10f}  (converged={})".format(ic.gs_energy(), ic.converged))

ks = np.linspace(-np.pi, np.pi, 61)
dim = D*D*(2-1)   # every branch of the ansatz at this bond dimension

# --- the band, its weight, and how much of the response it carries -----
bands, w_x, frac_x = [], [], []
all_e, all_wz, tot_z = [], [], []
for k in ks:
    e, wx, tx = ic.spectral_weights("Sx", k, n=1, return_total=True)
    ez, wz, tz = ic.spectral_weights("Sz", k, n=dim, return_total=True)
    bands.append(e[0])
    w_x.append(wx[0]) ; frac_x.append(wx[0]/tx)
    all_e.append(ez) ; all_wz.append(wz) ; tot_z.append(tz)
bands, w_x, frac_x = np.array(bands), np.array(w_x), np.array(frac_x)
all_e, all_wz, tot_z = np.array(all_e), np.array(all_wz), np.array(tot_z)
frac_z = all_wz/tot_z[:, None]

# Exact free-fermion band of the TFIM, eps(k) = 2*sqrt(1+g^2-2g*cos k),
# and the two-particle continuum eps(q)+eps(k-q) it generates.
def eps(k):
    return 2.0*np.sqrt(1.0 + g**2 - 2.0*g*np.cos(k))


exact = eps(ks)
qs = np.linspace(-np.pi, np.pi, 401)
two_particle = np.array([eps(qs) + eps(k - qs) for k in ks])
cont_lo, cont_hi = two_particle.min(axis=1), two_particle.max(axis=1)

print()
print("   k      E(k)       exact      w_x(k)     frac_x    w_z(lowest)")
for i in range(0, len(ks), 6):
    print("{:6.3f} {:10.6f} {:10.6f} {:10.6f} {:9.4f}   {:.3e}".format(
        ks[i], bands[i], exact[i], w_x[i], frac_x[i], all_wz[i, 0]))

# --- regression checks -------------------------------------------------
# 1) the dispersion is the exact free-fermion band
assert np.max(np.abs(bands - exact)) < 1e-4
# 2) the total weight is the connected static structure factor, summed in
#    real space from an entirely separate code path (vumps.py's own
#    mixed-gauge correlator, not the excitation ansatz)
mean = ic.vev("Sx", 0).real
Cr = [ic.correlator("Sx", 0, "Sx", r).real - mean*mean for r in range(60)]
for k, fx, wx in zip(ks, frac_x, w_x):
    ref = Cr[0] + 2.0*sum(np.cos(k*r)*Cr[r] for r in range(1, 60))
    assert abs(wx/fx - ref) < 1e-6
# 3) sigma^x is one coherent mode; sigma^z cannot reach it at all (an
#    exact parity selection rule, not merely a small number)
assert frac_x.min() > 0.99
assert frac_z[:, 0].max() < 1e-12
# 4) every branch that carries appreciable sigma^z weight lies inside the
#    exact two-particle continuum
heavy = frac_z > 1e-3
lo = np.repeat(cont_lo[:, None], dim, axis=1)
hi = np.repeat(cont_hi[:, None], dim, axis=1)
assert np.all(all_e[heavy] > lo[heavy] - 1e-6)
assert np.all(all_e[heavy] < hi[heavy] + 1e-6)

# --- the broadened S(k,w) map ------------------------------------------
delta = 0.08
ks_g, es_g, S = ic.dynamical_structure_factor(
    "Sx", ks=ks, energies=np.linspace(1.5, 8.5, 400), delta=delta, n=1)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
im = ax.pcolormesh(ks_g, es_g, S.T, shading="auto", cmap="magma")
ax.plot(ks, exact, "w--", lw=1.0, label="exact free-fermion band")
ax.set_xlabel("momentum $k$")
ax.set_ylabel(r"energy $\omega$")
ax.set_title(r"$S^{xx}(k,\omega)$, TFIM $g=%.1f$, $D=%d$" % (g, D))
ax.legend(loc="lower center", framealpha=0.7)
fig.colorbar(im, ax=ax)

ax = axes[1]
kk = np.repeat(ks[:, None], dim, axis=1)
mask = frac_z > 1e-4
ax.fill_between(ks, cont_lo, cont_hi, color="0.85", zorder=0,
                label="exact 2-particle continuum")
sc = ax.scatter(kk[mask], all_e[mask], c=frac_z[mask], cmap="viridis", s=26, zorder=3)
ax.plot(ks, exact, "r--", lw=1.2,
        label=r"1-particle band ($\sigma^z$ weight $\sim10^{-21}$)")
ax.set_xlabel("momentum $k$")
ax.set_ylabel(r"$E(k)$")
ax.set_ylim(0.0, cont_hi.max()*1.05)
ax.set_title(r"all branches, coloured by $\sigma^z$ weight fraction")
ax.legend(loc="lower center", fontsize=8, framealpha=0.8)
fig.colorbar(sc, ax=ax)

ax = axes[2]
i0 = int(np.argmin(np.abs(ks - 1.0)))
for opname, style, lab in (("Sx", "o-", r"$\sigma^x$"), ("Sz", "s--", r"$\sigma^z$")):
    _, w, tot = ic.spectral_weights(opname, ks[i0], n=dim, return_total=True)
    ax.plot(range(1, dim+1), np.cumsum(w)/tot, style, ms=4, label=lab)
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_ylim(-0.05, 1.1)
ax.set_xlabel("branches included (lowest first)")
ax.set_ylabel("fraction of the sum rule recovered")
ax.set_title(r"cumulative weight at $k=%.2f$" % ks[i0])
ax.legend()

fig.tight_layout()
fig.savefig("dynamical_structure_factor_tfim.png", dpi=150)
print("\nSaved plot to dynamical_structure_factor_tfim.png")
