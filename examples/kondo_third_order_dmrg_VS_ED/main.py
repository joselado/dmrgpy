# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Third-order STM/Kondo tunneling spectrum (Ternes, New J. Phys. 17
# 063016 (2015), arXiv:1505.04430), T=0, computed two independent ways
# for the same chain and compared: mode="ED" (exact diagonalization,
# eigenstate-sum construction, kondospectrumtk/conductance.py) and
# mode="DMRG" (itensor_version=3, real TDVP time evolution, never
# diagonalizing beyond the ground state -- kondospectrumtk/dmrgtwotime.py).
# U=0 here, so this isolates the second-order term (computed via the
# dynamical correlator, kondospectrumtk/secondorder_dc.py) plus the
# third-order Kondo term specifically (the potential-interference term
# is exercised instead in examples/kondo_potential_term_dmrg_VS_ED).
#
# mode="DMRG"'s dynamical-correlator calls use a much coarser `delta`
# (2e-5) than this feature's test suite (2e-6): the underlying C++ KPM
# routine picks its own moment count from delta (not from an `n` kwarg,
# confirmed directly -- `n` is accepted by the Python wrapper but never
# forwarded to the compiled extension), so a very fine delta makes each
# correlator call take tens of seconds even for this tiny 3-site chain
# (confirmed directly: delta=2e-6 took ~40s per call: delta=2e-5 takes
# ~1s and is still fine enough to resolve this chain's ~1.16meV Zeeman
# splitting).
import numpy as np
from dmrgpy import spinchain

G = 2.0
MUB = 5.7883818066e-5 # eV/T

# 3-site S=1/2 chain: a Zeeman-split impurity at site 0 (the tip-coupled
# site) with a weak Heisenberg coupling to two spectator sites -- ITensor
# v3's two-site DMRG/TDVP cannot handle chains shorter than 3 sites, so
# this is the smallest chain this feature's DMRG path can run on
# (see tests/test_kondo_spectrum_dmrgtwotime.py's _build_chain).
sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=3)
h = G*MUB*10.0*sc.Sz[0]
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()

eVs = np.linspace(-2e-3, 2e-3, 41)
es = np.linspace(-3e-3, 3e-3, 800) # must cover the ~1.16meV Zeeman gap
Jrho_s = 0.05
omega0, Gamma0 = 2e-3, 5e-6
n_t2_half, n_tau_half = 10, 15
dt2 = 25./Gamma0/n_t2_half
dtau = (2*np.pi/2e-5)/n_tau_half

print("Computing the third-order Kondo spectrum via mode=\"ED\" ...")
_, dIdV_ed = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=Jrho_s, U=0.0, T=0.0,
                                    order=3, mode="ED")

print("Computing the third-order Kondo spectrum via mode=\"DMRG\" "
      "(real TDVP time evolution) ...")
_, dIdV_dmrg = sc.get_kondo_spectrum(
        eVs, site=0, Jrho_s=Jrho_s, U=0.0, T=0.0, order=3, mode="DMRG",
        submode="KPM", delta=2e-5, es=es, dt2=dt2, n_t2_half=n_t2_half,
        dtau=dtau, n_tau_half=n_tau_half)

diff = np.max(np.abs(dIdV_dmrg-dIdV_ed))
print("max |DMRG - ED| = %.3f (max |ED| = %.3f)"%(diff, np.max(np.abs(dIdV_ed))))
# a qualitative/order-of-magnitude check, not a tight precision test: the
# DMRG path carries KPM delta-broadening error (second-order term) and
# t2/tau-grid discretization error (third-order Kondo term) on top of
# what the ED path has
assert diff < 0.3*np.max(np.abs(dIdV_ed))

import matplotlib.pyplot as plt

plt.plot(eVs, dIdV_ed, marker="o", label="ED (exact)")
plt.plot(eVs, dIdV_dmrg, marker="s", label="DMRG (TDVP two-time)")
plt.xlabel("eV")
plt.ylabel("dI/dV")
plt.title("Third-order Kondo spectrum: ED vs DMRG")
plt.legend()
plt.show()
