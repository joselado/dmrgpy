# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

def get(JK):
    n = 10 # number of spinful fermionic sites
    fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
    h = 0
    mu = 0.
    for i in range(1,n-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1]
        h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(1,n): # hopping
        h = h + mu/2.*fc.Cdagup[i]*fc.Cup[i]
        h = h + mu/2.*fc.Cdagdn[i]*fc.Cdn[i]
    for i in [0]: # Hubbard, to enforce half filling
        h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + JK/2.*fc.Sx[0]*fc.Sx[1]
    h = h + JK/2.*fc.Sy[0]*fc.Sy[1]
    h = h + JK/2.*fc.Sz[0]*fc.Sz[1]
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.maxm = 40
    fc.set_hamiltonian(h) # set the hoppings
    
    wf = fc.get_gs()
    ss = wf.get_correlation_entropy_density()
    print(JK)
    return ss

Js = np.linspace(0.,2.,10) ; sss = [get(J) for J in Js]
fo = open("MAP.OUT","w")
for j in range(len(sss[0])):
    for i in range(len(Js)):
        fo.write(str(j)+"  ")
        fo.write(str(Js[i])+"  ")
        fo.write(str(sss[i][j])+"\n")
fo.close()
exit()
ss = get(0.4) # get the density
import matplotlib.pyplot as plt
plt.plot(range(len(ss)),ss)
plt.show()
