# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain,spinchain
n = 4

# create first neighbor exchange
sc = fermionchain.Spinon_Chain(n) # create the spin chain
#sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)]) # create the spin chain
#sc.setup_julia()
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.get_gs()


#print(sc.vev(sc.Cdagup[0]*sc.Cup[1]))

#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
i,j = 0,0
t1 = time.time()
name = (sc.Sz[i],sc.Sz[j])
name = (sc.Cdagup[i],sc.Cup[j])
es = np.linspace(-0.5,6,2000)
delta = 5e-2
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",name=name,es=es,delta=delta)
t2 = time.time()
print("Time with DMRG",t2-t1)




# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,y2.real,c="blue",label="DMRG")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



