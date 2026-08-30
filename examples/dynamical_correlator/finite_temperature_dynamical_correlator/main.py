# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create a random spin chain
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
def fj(i,j):
  if abs(i-j)==1: return 1.0
  else: return 0.0
# set_exchange() was removed; this is the Hamiltonian it built --
# the coupling summed over BOTH orderings of every pair, so a
# nearest-neighbour fj=1 gives 2*sum_<ij> S_i.S_j, not sum_<ij>.
h = 0
for i in range(n):
  for j in range(n):
    if fj(i,j)!=0.0: h = h + fj(i,j)*sc.SS(i,j)
sc.set_hamiltonian(h)

#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(mode="ED",
        i=i,j=j,name="ZZ")
t2 = time.time()
print("Time with T=0",t2-t1)


# submode="ED" (the exact Lehmann sum) is the only submode with a finite-T
# implementation -- the default "KPM" raises here. This script never got
# far enough to find that out: before name="ZZ" was made to work on the ED
# path it died one call earlier, inside EDOperator, with "RuntimeError: No
# active exception to reraise".
(x3,y3) = sc.get_dynamical_correlator(mode="ED",T=0.0005,submode="ED",
        i=i,j=j,name="ZZ")
t3 = time.time()
print("Time with finite T",t3-t2)


#(x4,y4) = sc.get_dynamical_correlator(mode="ED",submode="KPM",
#        i=i,j=j,name="ZZ")
#t4 = time.time()
#print("Time with KPM-ED",t4-t3)


# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,np.abs(y2),c="blue",label="T=0.0")
plt.scatter(x3,np.abs(y3),c="green",label="T=0.1")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()












