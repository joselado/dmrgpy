# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import spinchain
n = 4 # number of spin sites

# first let us create a Hubabrd model

fc = fermionchain.Spinful_Fermionic_Hamiltonian(n) # create the object
fc.maxm = 20
U = 6.0
def ft(i,j):
    if abs(i//2-j//2)==1 and i%2==j%2: return 1.0 # first neighbor coupling
    if i==j: return -2*U # set to half filling
    return 0.0

def fu(i,j):
    if i==j: return U
    return 0.0

fc.set_hoppings(ft) # add the term to the Hamiltonian
fc.set_hubbard_spinful(fu) # add the term to the Hamiltonian
pairs = [(0,i) for i in range(n)]
(esh,ch) = fc.get_dynamical_correlator(name="ZZ",mode="DMRG",i=0,j=1)
esh /= 2./U # scale by the estimated effective exchange
print(fc.get_density_spinful())

# now create the Heisenberg chain
sc = spinchain.Spin_Hamiltonian([2 for i in range(n)])
sc.set_exchange(lambda i,j: 1.0*(abs(i-j)==1))
(ess,cs) = sc.get_dynamical_correlator(name="ZZ",mode="DMRG",i=0,j=1)

import matplotlib.pyplot as plt

plt.scatter(esh,ch,c="red",label="Hubbard")
plt.plot(ess,cs,c="blue",label="Heisenberg")
#plt.plot(inds,mz2,c="green",label="Exact")
plt.legend()
plt.xlabel("Site")
plt.ylabel("Correlator")
plt.show()



