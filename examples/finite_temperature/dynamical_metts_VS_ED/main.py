# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time correlator at finite temperature, C(t) = <Sz_0(t) Sz_0(0)>_T,
# of a Heisenberg chain in a field. It is computed with minimally
# entangled typical thermal states (METTS), averaged over random samples,
# and compared with the exact Boltzmann sum over all eigenstates (ED).
# At finite temperature this is what neutron scattering probes.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.edtk.edchain import EDOperator

n = 4 # number of sites
T = 1.0 # temperature
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_python() # python DMRG backend

# Heisenberg Hamiltonian in a field, H = sum_i S_i . S_{i+1} + 0.4 sum_i Sz_i
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + 0.4*sc.Sz[i]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension

name = (sc.Sz[0], sc.Sz[0])
ts, c1, err = sc.metts_dynamical_correlator(name, T, nt=8, dt=0.2,
        nsamples=150, nwarmup=30, dbeta_half_step=0.08, seed=42)

# exact result, sum over all eigenstates weighted by exp(-E/T)
ed = sc.get_ED_obj()
es, vs = ed.get_diagonalized_hamiltonian()
vs = np.array(vs)
a = vs.conj().T @ EDOperator(name[0], ed).SO @ vs
es = es - np.min(es)
w = np.exp(-es/T) ; w = w/w.sum() # Boltzmann weights
c2 = np.array([np.sum(w[:,None]*np.exp(1j*(es[:,None]-es[None,:])*t)*a*a.T) for t in ts])
print("Max difference METTS-ED =", np.max(np.abs(c1-c2)))

import matplotlib.pyplot as plt
plt.plot(ts, c2.real, linewidth=4, label="ED")
plt.errorbar(ts, c1.real, yerr=err, fmt="--", linewidth=2, label="METTS")
plt.xlabel("time t")
plt.ylabel("Re C(t)")
plt.legend()
plt.show()
