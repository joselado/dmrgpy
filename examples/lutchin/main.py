# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
import fermionchain
import geometry
n = 40 
g = geometry.chain()
#g = g.supercell(n)
h = g.get_hamiltonian()
h.add_rashba(1.0)
h.shift_fermi(3.0)
h.add_zeeman(3.0)
h.add_swave(1.0)
h.get_bands(operator="sy")
exit()
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
