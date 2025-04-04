# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6

def get_dc(n):
    spins = ["S=1" for i in range(n)] # spin 1/2 heisenberg chain
    # create first neighbor exchange
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    
    sc.maxm = 20 # bond dimension
    
    name = (sc.Sz[0],sc.Sz[0]) # correlator to compute
    es = np.linspace(-0.2,1,300)
    delta = 2e-2
    (x,y) = sc.get_dynamical_correlator(name=name,
            nex = 10, # number of excited states
            mode="DMRG",submode="EX",es=es,delta=delta)
    return (x,y)


Ns = range(4,16,2)
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,5))

ii = 0
for n in Ns:
    c = ii/len(Ns) # colorcode
    print("Doing L=",n)
    (x,y) = get_dc(n) ; ii+= 1
    plt.plot(x,y.real+ii,c=(c,0.,1-c),label="L = "+str(n))

plt.legend()
plt.xlabel("Energy")
plt.ylabel("Dynamical correlator")
plt.tight_layout()
plt.show()












