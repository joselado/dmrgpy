# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Finite temperature correlator C(t) = <Sz_0(t) Sz_0> of a Heisenberg chain
# in a field, with dynamical METTS (minimally entangled typical thermal
# states) on the Julia backend. It is averaged over random samples, so
# it carries a statistical error. We compare it with the exact sum
# over all the eigenstates, C(t) = sum_nm p_n |<n|Sz_0|m>|^2 e^{i(E_n-E_m)t}.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
T = 1.0 # temperature
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_julia()
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + 0.4*sc.Sz[i] # magnetic field
sc.set_hamiltonian(h)

name = (sc.Sz[0], sc.Sz[0])
ts, c, err = sc.metts_dynamical_correlator(name, T, nt=8, dt=0.2,
        nsamples=150, nwarmup=30, seed=42)

# exact reference with all the eigenstates
from dmrgpy.edtk.edchain import EDOperator
ed = sc.get_ED_obj()
es, vs = ed.get_diagonalized_hamiltonian()
es = es - np.min(es)
p = np.exp(-es/T) ; p = p/p.sum() # Boltzmann weights
U = np.array(vs)
A = U.conj().T @ EDOperator(name[0], ed).SO @ U
cex = [np.sum(p[:, None]*np.exp(1j*(es[:, None]-es[None, :])*t)*A*A.T) for t in ts]
print("Max difference =", np.max(np.abs(c-cex)))

import matplotlib.pyplot as plt
plt.plot(ts, np.real(cex), linewidth=4, label="exact")
plt.errorbar(ts, c.real, yerr=err, fmt="o", label="METTS")
plt.xlabel("time t")
plt.ylabel("Re C(t)")
plt.legend()
plt.show()
