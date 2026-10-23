# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Transverse spin spectrum S(w) = sum_n |<n|S+_0|GS>|^2 delta(w - E_n + E_0)
# of the S=1/2 Heisenberg chain. S+ changes Sz by one, so only the
# states of that one sector contribute. The SECTOR submode sums
# them directly, and here it must agree with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.operatornames import name2MO

n = 6 # number of sites
sc = spinchain.Spin_Chain([2 for i in range(n)])
sc.setup_python() # backend with quantum numbers

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 8 # number of DMRG sweeps

es = np.linspace(-0.2, 4.0, 400) # energies
delta = 0.1 # broadening of the peaks
name = (name2MO("Sm", sc)[0], name2MO("Sp", sc)[0]) # S-_0 and S+_0
(e1, s1) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name, es=es, delta=delta)
(e2, s2) = sc.get_dynamical_correlator(submode="SECTOR", name=name, nex=100, es=es, delta=delta)
print("Max difference SECTOR-ED =", np.max(np.abs(np.array(s2)-np.array(s1))))

import matplotlib.pyplot as plt
plt.plot(e1, np.real(s1), linewidth=4, label="ED")
plt.plot(e2, np.real(s2), "--", linewidth=2, label="SECTOR")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
