# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
# create a random spin chain
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain

def getsc(iv):
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
  sc.itensor_version = iv
  sc.get_gs()
  return sc

#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
j = i
t1 = time.time()
sc = getsc(2)
name = (sc.Sz[i],sc.Sz[j]) # name= takes a pair of MultiOperators, not a
                            # bare "ZZ" string -- see kpmdmrg.py's
                            # get_dynamical_correlator and
                            # tests/test_dynamical_correlator.py for the
                            # same convention
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",name=name)
t2 = time.time()
print("Time with ITensor2",t2-t1)


sc = getsc(3)
name = (sc.Sz[i],sc.Sz[j])
(x3,y3) = sc.get_dynamical_correlator(mode="DMRG",name=name)
t3 = time.time()
print("Time with ITensor3",t3-t2)




# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,np.abs(y2),c="blue",label="ITensor 2")
plt.scatter(x3,np.abs(y3),c="green",label="ITensor 3")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()












