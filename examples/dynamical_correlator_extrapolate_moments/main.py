# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain


n = 10
# create a random spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
#sc.itensor_version = "julia"
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]*np.random.random()
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
#sc.get_gs()

sc.kpmmaxm = 20 # KPM maxm
sc.maxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
name = (sc.Sz[i],sc.Sz[j])
es = np.linspace(-0.5,6,2000)
delta = 1e-1
sc.kpm_extrapolate = False
(x1,y1) = sc.get_dynamical_correlator(name=name,es=es,delta=delta*2)
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)
t2 = time.time()
print("Time without extrapolation",t2-t1)

sc.kpm_extrapolate = True
(x3,y3) = sc.get_dynamical_correlator(name=name,es=es,delta=delta)
t3 = time.time()
print("Time with extrapolation",t3-t2)



# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.scatter(x1,y1.real,c="red",label="Reduced quality")
plt.plot(x2,y2.real,c="blue",label="No extrapolation")
plt.scatter(x3,y3.real,c="green",label="Extrapolation")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



