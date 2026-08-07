# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy import mpsalgebra

np.random.seed(1) # reproducible random exchange couplings across sizes

sizes = [4, 6, 8, 10]
arnoldi_es, dmrg_es = [], []
for n in sizes:
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h +np.random.random()*sc.Sx[i]*sc.Sx[i+1]
        h = h +np.random.random()*sc.Sy[i]*sc.Sy[i+1]
        h = h +np.random.random()*sc.Sz[i]*sc.Sz[i+1]

    sc.set_hamiltonian(h)
    #print(sc.get_excited(n=40,mode="ED")) ; exit()
    e,wf = mpsalgebra.lowest_energy_arnoldi(sc,h,verbose=1)
    e_dmrg = sc.gs_energy()
    print("n =",n,"Arnoldi",e.real,"DMRG",e_dmrg)
    arnoldi_es.append(e.real)
    dmrg_es.append(e_dmrg)

plt.plot(sizes,arnoldi_es,"o-",label="Arnoldi")
plt.plot(sizes,dmrg_es,"s--",label="DMRG")
plt.xlabel("chain length n")
plt.ylabel("ground state energy")
plt.title("Arnoldi vs DMRG ground state energy")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()









