import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
sys.path.append(os.environ["PYGRAROOT"]) # root for pygra
import matplotlib.pyplot as plt
import fermionchain
import geometry

n = 60
g = geometry.chain()
g = g.supercell(n)
h = g.get_hamiltonian()
h.add_rashba(0.5)
h.shift_fermi(2.0)
h.add_zeeman(0.3)
hk = h.get_hk_gen()

def fd(i,j):
    if i==j: return 0.1
    else: return 0.0

fc = fermionchain.Fermionic_Hamiltonian(len(g.r))
fc.set_spinful_hoppings(h.intra)
fc.set_pairing(fd)
pairs = [(0,i) for i in range(n)]
out = fc.correlator(pairs)
import matplotlib.pyplot as plt
plt.plot(range(n),out,c="blue",marker="o")
plt.show()
