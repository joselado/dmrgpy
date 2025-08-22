# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
#sc.setup_julia()
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.get_gs()

sc.kpmmaxm = 20
sc.maxm = 20


from dmrgpy.kpmdmrg import kpm_wfa_wfb


#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
t1 = time.time()
name = (sc.Sz[i],sc.Sz[j])
es = np.linspace(-0.5,6,2000)
delta = 5e-2
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",name=name,es=es,delta=delta)
t2 = time.time()

# compute with the general KPM method

wf = sc.get_gs() # get the ground state
wfa = sc.Sz[i]*wf # first wavefunction
wfb = sc.Sz[j]*wf # second wavefunction
(x3,y3) = kpm_wfa_wfb(sc,wfa=wfa,wfb=wfb,delta=delta,X=h) # using the general KPM

x3 = x3 - sc.gs_energy() # shift by the ground state energy


# plot the results
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,y2.real,c="blue",label="Spin KPM")
plt.scatter(x3,y3.real,c="green",label="General KPM")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()












