# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

Mz = sum(sc.Sz) # total magnetization in Z
Mx = sum(sc.Sx) # total magnetization in X

#sc.itensor_version = "julia"


def get(mode="DMRG",z=1):
    print("Doing",z)
    sc.set_hamiltonian(Mz) # only magnetic field
    wf = sc.get_gs(mode=mode) # Z ferromagnetic wavefunction
    # The operator in the exponent MUST be Hermitian
    wf1 = sc.exponential(z*Mx,wf,mode=mode) # compute exp(Mx)*wf
    c = sc.overlap(wf,wf1,mode=mode) # compute <wf|exp(Mx)|wf>
    return c # return the number

# Compare the result between DMRG and ED
#print("ED",get("ED"))
#print("DMRG",get("DMRG"))
#exit()


zs = np.linspace(-1.0,1.0,20)
us0 = np.log([get(z=z,mode="DMRG").real for z in zs])
us1 = np.log([get(z=z,mode="ED").real for z in zs])

import matplotlib.pyplot as plt


plt.scatter(zs,us0,label="DMRG",marker="o",s=200,c="blue")
plt.scatter(zs,us1,label="ED",marker="o",s=60,c="red")
plt.legend()
plt.show()


