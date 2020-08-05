# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create a random spin chain
#spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
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
print("Time with DMRG",t2-t1)
#exit()


(x3,y3) = sc.get_dynamical_correlator(mode="ED",submode="ED",name=name,es=es,
        delta=delta)
t3 = time.time()
print("Time with ED",t3-t2)


(x4,y4) = sc.get_dynamical_correlator(mode="ED",submode="KPM",name=name,es=es,
        delta=delta)
t4 = time.time()
print("Time with KPM-ED",t4-t3)


# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,y2.real,c="blue",label="DMRG")
plt.scatter(x3,y3.real,c="green",label="ED")
plt.scatter(x4,y4.real,c="red",label="ED KPM")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



