# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Thermal energy <H>_T = Tr[H exp(-H/T)]/Tr[exp(-H/T)] of a Heisenberg
# chain. We compute it with a purified state (the chain is doubled with
# an ancilla, and the state is cooled down to temperature T), and compare
# with the exact Boltzmann average over all the eigenvalues.
import numpy as np
from dmrgpy import spinchain, thermal

n = 6 # number of sites
Ts = [0.1, 0.25, 0.5, 1.0, 2.0, 10.0] # temperatures

def hamiltonian(sc):
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.set_hamiltonian(hamiltonian(sc))
hm = sc.get_ED_obj().get_hamiltonian() # full matrix
hm = hm.toarray() if hasattr(hm,"toarray") else hm
ev = np.linalg.eigvalsh(hm)
exact = []
for T in Ts: # exact thermal average
    w = np.exp(-(ev-ev.min())/T)
    exact.append(np.sum(ev*w)/np.sum(w))

pur = []
for T in Ts:
    tc = thermal.Thermal_Spin_Chain(["S=1/2" for i in range(n)])
    h = hamiltonian(tc)
    tc.set_hamiltonian(h)
    tc.T = T # temperature
    tc.mode = "ED"
    wf = tc.get_gs() # purified thermal state
    pur.append(wf.dot(h*wf).real)
diff = np.max(np.abs(np.array(pur)-exact))
print("Max difference purification-exact =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(Ts, exact, linewidth=4, label="exact")
plt.plot(Ts, pur, "--", linewidth=2, label="purification")
plt.xscale("log")
plt.xlabel("temperature T")
plt.ylabel("thermal energy <H>")
plt.legend()
plt.show()
