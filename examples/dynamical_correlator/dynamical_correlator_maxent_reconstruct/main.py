# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 20
# create a random spin chain
spins = ["S=1/2" for i in range(n)] # spin 1 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

sc.kpmmaxm = 20 # KPM maxm
sc.maxm = 20 # KPM maxm
import time
i = 0
j = i

t1 = time.time()
name = (sc.Sx[i],sc.Sx[j])
es = np.linspace(-.5,5.,1000)
restart = False
if restart:
    (x0,y0) = sc.get_dynamical_correlator(name=name,delta=1e-1,es=es)
    np.savetxt("DIST.OUT",np.array([x0.real,y0.real]).T)
else:
    m = np.loadtxt("DIST.OUT").T
    x0,y0 = m[0],np.abs(m[1])
from dmrgpy.reconstruct import reconstruct_distribution
t0 = time.time()
(x1,y1) = reconstruct_distribution(x0,y0,n=8,bnds=None)
t1 = time.time()

print("Time in the reconstruction",t1-t0)





# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x0,y0.real,c="blue",label="Plain")
plt.scatter(x1,y1.real,c="red",label="Reconstructed")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()












