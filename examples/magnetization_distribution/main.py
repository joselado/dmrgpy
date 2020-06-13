# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 40 # number of sites in your chain
spins = [2 for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
#sc.itensor_version = "julia"

# now define the Hamiltonian
h = 0
for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1] # add exchange
for i in range(n): h = h - .5*sc.Sx[i] # add transverse field
#h = h - sc.Sz[n-1]*sc.Sz[0] # and apply periodic boundary conditions
sc.set_hamiltonian(sum(sc.Sz)) ; wf = sc.get_gs()
sc.set_hamiltonian(h) # and initialize the Hamiltonian
sc.get_gs(wf0=wf)
# setup some parameters
sc.maxm = 60
sc.kpmmaxm = sc.maxm


# now define the operator for which you want the distribution
M = (1/n)*sum(sc.Sz)
import time
t0 = time.time()
x,y = sc.get_distribution(X=M,delta=2e-2,kernel="lorentz") # compute a distribution
t1 = time.time()
print("Time",t1-t0)
x1,y1 = sc.get_distribution(X=M,delta=2e-2,kernel="jackson") # compute a distribution
#x1,y1 = sc.get_distribution(X=M,delta=1e-1,mode="ED") # compute a distribution

# plot the result and save it in a file
print("Integral DMRG=",np.trapz(y.real,dx=x[1]-x[0]))
#print("Integral ED=",np.trapz(y1.real,dx=x1[1]-x1[0]))
import matplotlib.pyplot as plt
np.savetxt("DISTRIBUTION.OUT",np.array([x,y.real]).T)
plt.plot(x,y.real,marker="o",label="DMRG Lorentz") # correlator using DMRG
plt.plot(x1,y1.real,marker="o",label="DMRG Jackson") # correlator using DMRG
plt.xlabel("magnetization")
plt.ylabel("distribution")
plt.legend()
plt.show()


