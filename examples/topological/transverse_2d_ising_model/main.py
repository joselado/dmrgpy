# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.environ["DMRGROOT"])
import os ; import sys ; sys.path.append(os.environ["PYQULAROOT"])

# Transverse-field Ising model on an L x L square lattice (periodic),
# H = -4 sum_ij t_ij Sz_i Sz_j + 2 Bx sum_i Sx_i, where t_ij is the
# neighbor matrix of the lattice. We compute the gap between the two
# lowest states with ED. It is small in the ordered phase, at small Bx.
import numpy as np
from dmrgpy import spinchain
from pyqula import geometry

def get_gap(L, bx):
    g = geometry.single_square_lattice()
    m = g.get_supercell(L).get_hamiltonian(has_spin=False).get_hk_gen()([0., 0., 0.])
    n = m.shape[0] # number of spins
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n):
        for j in range(n):
            h = h - 4*m[i, j]*sc.Sz[i]*sc.Sz[j]
    h = h + 2*bx*sum(sc.Sx)
    sc.set_hamiltonian(h)
    es = sc.get_excited(n=2, mode="ED")
    return es[1] - es[0]

import matplotlib.pyplot as plt
bxs = np.linspace(0, 4, 10)
for L in [2, 3]:
    plt.plot(bxs, [get_gap(L, bx) for bx in bxs], "o-", label=str(L)+"x"+str(L))
plt.xlabel("Bx")
plt.ylabel("gap")
plt.legend()
plt.show()
