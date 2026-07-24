# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Third-order potential-interference term of the STM/Kondo tunneling
# spectrum (eq. "U-M", Ternes, New J. Phys. 17 063016 (2015),
# arXiv:1505.04430), T=0, U!=0, computed two independent ways for the
# same chain and compared: the ED excited-state sum
# (conductance.third_order_potential_dIdV) and the DMRG dynamical-
# correlator construction (kondospectrumtk/potentialdc.py,
# itensor_version=3) added to fill the one remaining gap in this
# feature's mode="DMRG" path (it used to raise NotImplementedError for
# U!=0). Companion to examples/kondo_third_order_dmrg_VS_ED, which
# instead isolates the second-order term and the third-order *Kondo*
# term (U=0).
#
# See that other example's module docstring for why `delta` is set much
# coarser here (2e-5) than this feature's test suite default (2e-6): the
# compiled KPM extension picks its own moment count from delta, and a
# very fine delta makes each dynamical-correlator call take tens of
# seconds even for this tiny chain.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc

G = 2.0
MUB = 5.7883818066e-5 # eV/T

# same 3-site chain as examples/kondo_third_order_dmrg_VS_ED (ITensor
# v3's two-site DMRG/TDVP needs at least 3 sites)
sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=3)
h = G*MUB*10.0*sc.Sz[0]
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()
ks = KondoSpectrum(sc, site=0, T=0.0)

eVs = np.linspace(-2e-3, 2e-3, 41)
es = np.linspace(-3e-3, 3e-3, 800) # must cover the ~1.16meV Zeeman gap
Jrho_s, U = 0.05, 0.2

print("Computing the potential-interference term via the excited-state "
      "sum (ED, exact) ...")
dIdV_ed = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)

print("Computing the potential-interference term via the T=0 dynamical "
      "correlator (mode=\"DMRG\") ...")
dIdV_dmrg = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, T0=1.0,
                                           mode="DMRG", submode="KPM",
                                           delta=2e-5, es=es)

diff = np.max(np.abs(dIdV_dmrg-dIdV_ed))
print("max |DMRG - ED| = %.4f (max |ED| = %.4f)"%(diff, np.max(np.abs(dIdV_ed))))
# qualitative/order-of-magnitude check (KPM delta-broadening error), not
# a tight precision test -- see test_kondo_spectrum_potentialdc.py for
# the tight mode="ED" submode="ED" check of this same function
assert diff < 0.15*np.max(np.abs(dIdV_ed))

import matplotlib.pyplot as plt

plt.plot(eVs, dIdV_ed, marker="o", label="ED (exact)")
plt.plot(eVs, dIdV_dmrg, marker="s", label="DMRG (dynamical correlator)")
plt.xlabel("eV")
plt.ylabel("Potential-interference term of dI/dV")
plt.title("Third-order potential-interference term: ED vs DMRG")
plt.legend()
plt.show()
