# Infinite DMRG (iDMRG): VUMPS static correlators (vev/two_point_correlator),
# pure-Python (pyitensor) backend vs ITensor v3 (mpscpp3), on the
# transverse-field Ising model (TFIM).
#
# mpscpp3/chain_session.h's Chain::vumps_onsite_expectation/
# Chain::vumps_two_point_correlator are a line-for-line C++ port of
# pyitensor/vumps.py's own onsite_expectation/two_point_correlator --
# both compute <O>/<O_i O_j> directly from a converged VUMPS ground
# state's mixed-gauge {AC, AR} (Vanderstraeten, Haegeman, Verstraete,
# "Tangent-space methods for uniform matrix product states",
# arXiv:1810.07006, Eq.(34)/(37)-(39)), no eigenproblem needed, unlike
# pyitensor/idmrg.py's own dominant-right-fixed-point-based correlators.
# This example checks both backends land on the same ground-state energy
# AND the same two-point correlator as a function of both the bond
# dimension D and the separation r, and times each.
#
# H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] = -sigma^x_i sigma^x_{i+1} - g*sigma^z_i
# (Sx=sigma^x/2, Sz=sigma^z/2), the standard Pauli-matrix TFIM convention,
# g>1: gapped paramagnetic phase. Requires gs_method="vumps" on both
# backends (itensor_version=3's own gs_method="idmrg" has no correlator
# support at all, see ROADMAP.md).

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

G = 1.5  # transverse field (g>1: gapped paramagnetic phase)
D_VALUES = (1, 2, 3)
RS = list(range(0, 6))


def build_chain(itensor_version, D):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * G * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.maxm = D
    ic.maxiter = 800
    ic.vumps_nrestarts = 6
    return ic


e0 = {}
sz = {}
corr = {}
timing = {}
for version in ("python", 3):
    for D in D_VALUES:
        ic = build_chain(version, D)
        t0 = time.time()
        density = ic.gs_energy()
        t_gs = time.time() - t0
        t0 = time.time()
        sz[(version, D)] = ic.vev("Sz", 0).real
        corr[(version, D)] = np.array([ic.correlator("Sz", 0, "Sz", r).real for r in RS])
        t_corr = time.time() - t0
        e0[(version, D)] = density
        timing[(version, D)] = (t_gs, t_corr)
        print("itensor_version={:6} D={} e0={: .10f} converged={}  "
              "<Sz>={: .6f}  t_gs={:.2f}s t_correlator_scan={:.2f}s".format(
                  str(version), D, density, ic.converged, sz[(version, D)], t_gs, t_corr))

print()
print("-- cross-backend agreement (python vs v3) --")
for D in D_VALUES:
    de0 = abs(e0[("python", D)] - e0[(3, D)])
    dsz = abs(sz[("python", D)] - sz[(3, D)])
    dcorr = np.max(np.abs(corr[("python", D)] - corr[(3, D)]))
    print("D={}: |e0 diff|={:.2e}  |<Sz> diff|={:.2e}  max|correlator diff|={:.2e}".format(
        D, de0, dsz, dcorr))
    assert de0 < 1e-6, "python and v3 VUMPS energies disagree at D={}".format(D)
    assert dsz < 1e-6, "python and v3 VUMPS <Sz> disagree at D={}".format(D)
    assert dcorr < 1e-6, "python and v3 VUMPS correlators disagree at D={}".format(D)

print()
print("VUMPS static-correlator python-VS-v3 backend regression test PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

markers = {1: "o", 2: "s", 3: "^"}
colors = {"python": "tab:blue", 3: "tab:red"}
# v3 plotted as small filled markers, python as larger open (unfilled)
# markers on top -- the two agree to ~1e-14 (see the printed diffs above),
# so a naive same-style overlay would just hide one series entirely
# underneath the other.
for D in D_VALUES:
    ax1.plot(RS, corr[(3, D)], markers[D], color=colors[3],
              markersize=6, alpha=0.5 + 0.15 * D, label="v3 D={}".format(D))
for D in D_VALUES:
    ax1.plot(RS, corr[("python", D)], markers[D], color=colors["python"],
              markersize=11, markerfacecolor="none", markeredgewidth=1.5,
              linestyle="none", label="python D={}".format(D))
ax1.axhline(0, color="gray", lw=0.8, ls=":")
ax1.set_xlabel("distance $r$")
ax1.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax1.set_title("TFIM (g={}): static correlator vs $r$".format(G))
ax1.legend(fontsize=7)

for version, marker, label in (("python", "o-", "python"), (3, "s--", "v3")):
    ts = [timing[(version, D)][0] + timing[(version, D)][1] for D in D_VALUES]
    ax2.plot(D_VALUES, ts, marker, label=label)
ax2.set_xlabel("bond dimension $D$")
ax2.set_ylabel("wall time (gs_energy + correlator scan) [s]")
ax2.set_yscale("log")
ax2.set_title("Timing vs $D$")
ax2.legend()

fig.tight_layout()
fig.savefig("vumps_correlator_v3_VS_python.png", dpi=150)
print("\nSaved plot to vumps_correlator_v3_VS_python.png")
