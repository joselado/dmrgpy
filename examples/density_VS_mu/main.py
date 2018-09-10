import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain

n = 30 # number of spinful fermionic sites

def get_fc(mu):
  fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
  
  ####### Input matrices #######
  
  # Array with the hoppings and with hubbard couplings
  # These are the matrices that you have to modify
  hopping = np.zeros((n,n))
  hubbard = np.zeros((n,n))
  for i in range(n-1):  hopping[i,i+1] = 1. ; hopping[i+1,i] = 1.
  for i in range(n-3):  hopping[i,i+3] = 1./3. ; hopping[i+3,i] = 1./3.
  for i in range(n): U = 4.0 ; hubbard[i,i] = U/2. ; hopping[i,i] = -U + mu
  
  # The implemented Hamiltonian is
  # H = \sum_ij hopping[i,j] c^dagger_i c_j + hubbard[i,j] n_i n_j
  # with n_i = c^\dagger_{i,up} c_{i,up} + c^\dagger_{i,dn} c_{i,dn}
  
  # the previous matrices are for a half filled Hubbard chain
  
  ##############################
  
  
  # Setup the Many Body Hamiltonian
  fc.set_hoppings(lambda i,j: hopping[i,j]) # set the hoppings
  fc.set_hubbard(lambda i,j: hubbard[i,j]) # set the hubbard constants
  #fc.set_fields(lambda i: [0.,0.,0.2]) # set the hubbard constants
  
  #fc.nsweeps = 7
  
  
  # Compute the dynamical correlator defined by
  # <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>
  
  i = 0 # first index of the dynamical correlator
  j = 0 # second index of the dynamical correlator
  delta = 0.1 # energy resolution (approximate)
  fc.nsweeps = 14
  fc.kpmmaxm = 20 # maximum bond dimension in KPM
  return fc






mus = np.linspace(0.5,2.5,30)
ds = []
d2s = []
deltas = []
for mu in mus:
  fc = get_fc(mu)
  d = np.mean(fc.get_density())
  d2 = np.mean(fc.get_density_fluctuation())
  delta = np.abs(np.mean(fc.get_delta()))
  deltas.append(delta)
  ds.append(d)
  d2s.append(d2)
  print(mu,d,delta,d2)


ds = np.array(ds)
deltas = np.array(deltas)


np.savetxt("D_VS_MU.OUT",np.matrix([mus,deltas,ds,d2s]).T)

import matplotlib.pyplot as plt
plt.subplot(121)
plt.plot(mus,ds,marker="o")
plt.plot(mus,d2s,marker="o")
plt.subplot(122)
plt.plot(mus,deltas,marker="o")
plt.show()
