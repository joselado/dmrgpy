# Infinite DMRG (iDMRG): VUMPS ground state + tangent-space/quasiparticle
# excitation ansatz, pure-Python (pyitensor) backend vs ITensor v3
# (mpscpp3), on the transverse-field Ising model (TFIM).
#
# mpscpp3/chain_session.h's Chain::vumps_ground_state/
# Chain::vumps_excitation_energies are a C++ port of pyitensor/vumps.py's
# and pyitensor/idmrg_excitations.py's own algorithms -- unlike
# Chain::idmrg_ground_state (built from genuine ITensor tensor-network
# objects), this port works entirely with dense row-major arrays closed
# over LAPACK (see VumpsResult's own doc comment in chain_session.h for
# why). This example checks both backends land on the same VUMPS ground-
# state energy AND the same momentum-resolved excitation dispersion, and
# both against the exact free-fermion single-magnon band (see
# examples/idmrg/excitation_gap_tfim/main.py, this example's own
# itensor_version="python"-only precursor).
#
# H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] = -sigma^x_i sigma^x_{i+1} - g*sigma^z_i
# (Sx=sigma^x/2, Sz=sigma^z/2), the standard Pauli-matrix TFIM convention.
# Exact free-fermion single-magnon dispersion: eps(k) = 2*sqrt(J^2 + g^2
# - 2*J*g*cos(k)), J=1. Requires gs_method="vumps" on both backends.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

J, G = 1.0, 2.5  # g > 1: gapped paramagnetic phase
D_VALUES = (1, 2, 3)
KS = np.linspace(0, np.pi, 13)
EXACT_DISPERSION = 2 * np.sqrt(J ** 2 + G ** 2 - 2 * J * G * np.cos(KS))


def build_chain(itensor_version, D):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "vumps"
    h = -4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * G * ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.maxm = D
    ic.etol = 1e-10
    ic.vumps_nrestarts = 6
    return ic


e0 = {}
dispersion = {}
timing = {}
for version in ("python", 3):
    for D in D_VALUES:
        ic = build_chain(version, D)
        t0 = time.time()
        density = ic.gs_energy()
        t_gs = time.time() - t0
        t0 = time.time()
        disp = np.array([ic.excitation_energies(k, n=1)[0].real for k in KS])
        t_exc = time.time() - t0
        e0[(version, D)] = density
        dispersion[(version, D)] = disp
        timing[(version, D)] = (t_gs, t_exc)
        print("itensor_version={:6} D={} e0={: .10f} converged={}  "
              "t_gs={:.2f}s t_excitation_scan={:.2f}s".format(
                  str(version), D, density, ic.converged, t_gs, t_exc))

print()
print("-- cross-backend agreement (python vs v3) --")
for D in D_VALUES:
    de0 = abs(e0[("python", D)] - e0[(3, D)])
    ddisp = np.max(np.abs(dispersion[("python", D)] - dispersion[(3, D)]))
    print("D={}: |e0 diff|={:.2e}  max|dispersion diff|={:.2e}".format(D, de0, ddisp))
    assert de0 < 1e-6, "python and v3 VUMPS energies disagree at D={}".format(D)
    assert ddisp < 1e-5, "python and v3 excitation dispersions disagree at D={}".format(D)

print()
print("-- agreement with the exact free-fermion dispersion (D=3) --")
err = np.max(np.abs(dispersion[(3, 3)] - EXACT_DISPERSION))
print("max|E(k) - exact| (v3, D=3) = {:.2e}".format(err))
assert err < 2e-2

print()
print("VUMPS + excitation ansatz python-VS-v3 backend regression test PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

ax1.plot(KS, EXACT_DISPERSION, "-", color="black", label="exact (free fermion)")
markers = {1: "o", 2: "s", 3: "^"}
colors = {"python": "tab:blue", 3: "tab:red"}
for version in ("python", 3):
    for D in D_VALUES:
        ax1.plot(KS, dispersion[(version, D)], markers[D], color=colors[version],
                  markersize=5, alpha=0.5 + 0.15 * D,
                  label="{} D={}".format(version if version == "python" else "v3", D))
ax1.set_xlabel("momentum $k$")
ax1.set_ylabel("$E(k)$")
ax1.set_title("TFIM (g={}): single-magnon dispersion".format(G))
ax1.legend(fontsize=7)

for version, marker, label in (("python", "o-", "python"), (3, "s--", "v3")):
    ts = [timing[(version, D)][0] + timing[(version, D)][1] for D in D_VALUES]
    ax2.plot(D_VALUES, ts, marker, label=label)
ax2.set_xlabel("bond dimension $D$")
ax2.set_ylabel("wall time (s): ground state + full $k$-scan")
ax2.set_xticks(D_VALUES)
ax2.set_title("wall time")
ax2.legend()

fig.suptitle("VUMPS + excitation ansatz: python (pyitensor) vs ITensor v3 (mpscpp3)")
fig.tight_layout()
fig.savefig("vumps_excitation_v3_VS_python.png", dpi=150)
print("saved plot to vumps_excitation_v3_VS_python.png")
