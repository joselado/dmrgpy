# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# ESR-STM spectrum of a 12-site antiferromagnetic spin chain at 50 mK. Only
# a few states are populated, but many final states lie in the window. With
# submode="KRYLOV", DMRG gives only the states below the window, and the
# final states of the lines come from a Krylov space started from V|m> in
# the tangent space of each populated state m. We compare it with ED.
import numpy as np
from dmrgpy import spinchain

h = 4.135667696e-6 # Planck constant in eV/GHz
n = 12 # number of sites
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

f = np.linspace(15., 30., 751) # RF frequency in GHz
kw = dict(tip=[0.5, 0., 0.87], site=0, T=0.05, T1=100., T2=2., f_ref=20.,
          rabi=1e-4, vrf_over_vdc=0.1, p_cut=1e-4)
f, ied = sc.get_esr_spectrum(f, mode="ED", **kw)
f, ikr = sc.get_esr_spectrum(f, mode="DMRG", submode="KRYLOV", nstates=6, **kw)
print("Max difference DMRG-ED =", np.max(np.abs(ikr - ied)),
      "for a peak of", np.max(np.abs(ied)))
assert np.max(np.abs(ikr - ied)) < 1e-4*np.max(np.abs(ied))

import matplotlib.pyplot as plt
plt.plot(f, ied, linewidth=4, label="ED")
plt.plot(f, ikr, "--", linewidth=2, label="DMRG, Krylov")
plt.xlabel("RF frequency (GHz)")
plt.ylabel("dI (units of G a V_DC)")
plt.legend()
plt.show()
