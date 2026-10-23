# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order Kondo term of the STM tunneling spectrum of a spin chain in a
# magnetic field (Ternes, New J. Phys. 17 063016). We compute it with
# exact diagonalization (ED) and with DMRG (two-time correlator with TDVP),
# on the compiled backend (v3) and on the pure Python one, and compare
# the time each one takes.
import time
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_kondo_dIdV
from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg

muB = 5.7883818066e-5 # Bohr magneton in eV/T
def chain(version):
    sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=version)
    h = 2.0*muB*10.0*sc.Sz[0] # Zeeman field of 10 T
    for i in range(2):
        h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 15 # bond dimension and sweeps
    return sc

eVs = np.linspace(-2e-3, 2e-3, 21) # bias voltages
omega0, Gamma0 = 2e-3, 2e-3 # cutoff and broadening
grid = dict(dt2=150.0, n_t2_half=20, dtau=150.0, n_tau_half=20) # time grid
curves = {}

t0 = time.time()
ks = KondoSpectrum(chain(2), site=0, T=0.0)
curves["ED"] = third_order_kondo_dIdV(ks, eVs, 1.0, T0=1.0, omega0=omega0, Gamma0=Gamma0)
print("ED time", time.time()-t0)
for version in [3, "python"]:
    t0 = time.time()
    sc = chain(version)
    sc.get_gs()
    d = two_time_kondo_term_dmrg(sc, 0, eVs, omega0=omega0, Gamma0=Gamma0, **grid)
    curves["DMRG %s" % version] = 4*np.pi*d
    print("DMRG", version, "time", time.time()-t0)

import matplotlib.pyplot as plt
plt.plot(eVs, curves["ED"], linewidth=4, label="ED")
plt.plot(eVs, curves["DMRG 3"], "--", linewidth=2, label="DMRG v3")
plt.plot(eVs, curves["DMRG python"], ":", linewidth=2, label="DMRG python")
plt.xlabel("eV")
plt.ylabel("third-order Kondo term of dI/dV")
plt.legend()
plt.show()
