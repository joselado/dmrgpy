# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.environ["DMRGROOT"])
import os ; import sys ; sys.path.append(os.environ["PYQULAROOT"])

import numpy as np
from dmrgpy import spinchain

from pyqula import geometry

def get(n=2,bx=0.):
    g = geometry.single_square_lattice()
    m = g.get_supercell(n).get_hamiltonian(has_spin=False).get_hk_gen()([0.,0.,0.])
    
    n = m.shape[0] # number of sites in your chain
    spins = ["S=1/2" for i in range(n)] # create the sites
    sc = spinchain.Spin_Chain(spins) # create the chain
    
    h = 0
    # now define the Hamiltonian
    for i in range(n):
        for j in range(n):
            h = h + (-1)*m[i,j]*sc.Sz[i]*sc.Sz[j]
    
    h = 4*h + 2*bx*sum(sc.Sx)
    sc.set_hamiltonian(h)
    es = sc.get_excited(n=2,mode="ED")
    return es[1]-es[0]

import matplotlib.pyplot as plt

bxs = np.linspace(0.,4.,10)
for n in [2,3,4]:
    gs = [get(n=n,bx=bx) for bx in bxs]
    plt.plot(bxs,gs,marker="o",label=str(n)+"x"+str(n))


plt.xlabel("Bx")
plt.ylabel("Gap")
plt.legend()
plt.show()


