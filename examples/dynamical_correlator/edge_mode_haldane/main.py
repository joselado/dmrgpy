# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0)
# at the end of an open S=1 Heisenberg chain. The Haldane chain has
# spin-1/2 edge modes, so a low energy peak appears that becomes sharp
# when the chain is long. It is computed with the excited states (EX).
import numpy as np
from dmrgpy import spinchain

def get_dc(n):
    sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 20 # bond dimension
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(-0.2, 1, 300) # energies
    return sc.get_dynamical_correlator(name=name, nex=10, mode="DMRG",
            submode="EX", es=es, delta=0.02) # nex excited states

import matplotlib.pyplot as plt
for n in [4, 6, 8, 10]:
    x, y = get_dc(n)
    plt.plot(x, y.real + n/2, label="L = "+str(n)) # shifted for clarity
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
