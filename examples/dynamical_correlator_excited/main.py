# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create a random spin chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
h = 0
for i in range(n-1):
      j = i+1
      p = np.random.random()
      h = h + p*(sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j])


#h = h + h.get_dagger()
sc.set_hamiltonian(h)
sc.nsweeps = 30
sc.maxm = 20
sc.kpmmaxm = 20

#print(sc.gs_energy(mode="ED"))
#print(sc.gs_energy(mode="DMRG"))
#exit()

#sc.kpmmaxm = 20 # KPM maxm

es0 = sc.get_excited(n=14,mode="DMRG")
es1 = sc.get_excited(n=14,mode="ED")
print("Difference in excited states",np.mean(np.abs(es0-es1)))

import time
i = np.random.randint(n)
j = np.random.randint(n)
name = (sc.Sz[i],sc.Sz[j])
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",submode="EX",name=name)
t2 = time.time()
print("Time with DMRG excited states",t2-t1)


(x3,y3) = sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name)
t3 = time.time()
print("Time with DMRG KPM",t3-t2)


(x4,y4) = sc.get_dynamical_correlator(mode="ED",submode="ED",name=name)
t4 = time.time()
print("Time with ED",t4-t3)



# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,np.abs(y2),c="blue",label="DMRG EX")
plt.scatter(x3,np.abs(y3),c="green",label="DMRG KPM")
plt.scatter(x4,np.abs(y4),c="red",label="ED")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



