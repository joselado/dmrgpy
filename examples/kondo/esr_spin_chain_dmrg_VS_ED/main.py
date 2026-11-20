# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# ESR-STM spectrum of an antiferromagnetic spin 1/2 chain in a field, with
# the tip on the first spin. Each line is a transition between many-body
# states, weighted by its matrix element, by the population difference and
# by the change of the spin under the tip. With DMRG the spectrum is read
# from the lowest eigenstates; we compare it with exact diagonalization.
import numpy as np
from dmrgpy import spinchain

h = 4.135667696e-6 # Planck constant in eV/GHz
n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
hh = 0
for i in range(n-1):
    J = 30.0*(1 + 0.2*np.cos(i)) # exchange in GHz
    hh = hh + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
hh = hh + 20.0*sum(sc.Sz) # Zeeman energy of 20 GHz
hh = hh + 2.0*sc.Sz[0] + 0.5*sc.Sx[0] # field of the tip on the first spin
sc.set_hamiltonian(h*hh) # energies in eV
sc.maxm = 64 # bond dimension
sc.nsweeps = 20 # number of DMRG sweeps

f = np.linspace(15., 30., 1501) # RF frequency in GHz
kw = dict(tip=[0.5, 0., 0.87], site=0, T=0.05, T1=100., T2=2., f_ref=20.,
          rabi=1e-4, vrf_over_vdc=0.1, p_cut=1e-4)
f, ied = sc.get_esr_spectrum(f, mode="ED", **kw)
f, idmrg = sc.get_esr_spectrum(f, mode="DMRG", nstates=24, **kw)
print("Max difference DMRG-ED =", np.max(np.abs(idmrg - ied)),
      "for a peak of", np.max(np.abs(ied)))

import matplotlib.pyplot as plt
plt.plot(f, ied, linewidth=4, label="ED")
plt.plot(f, idmrg, "--", linewidth=2, label="DMRG")
plt.xlabel("RF frequency (GHz)")
plt.ylabel("dI (units of G a V_DC)")
plt.legend()
plt.show()
