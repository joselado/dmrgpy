# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
ns = 10 # number of spinful fermionic sites
def get_fc(U):
    n = ns*2
    fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
    C = fc.C
    Cdag = fc.Cdag
    N = fc.N
    # add the Hamiltonian
    h = 0
    
    # add hopping
    for i in range(ns-1):
        for j in range(2):
          h = h + Cdag[2*i+j]*C[2*(i+1)+j]
    # add Hubbard
    h = h + h.get_dagger()
    for i in range(ns):
      d = 0
      for j in range(2): d = d + N[2*i+j]
      h = h + U*(d*d - 2*d) # Hubbard term
    fc.set_hamiltonian(h)
    print("Doing",U)
    return fc

print(get_fc(4.0).get_gap(mode="DMRG"))
#print(get_fc(4.0).get_gap(mode="DMRG"))
#exit()
us = np.linspace(0.0,10.0,20)
gs = [get_fc(u).get_gap() for u in us]
import matplotlib.pyplot as plt
plt.plot(us,gs,marker="o")
plt.show()


