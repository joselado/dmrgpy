from __future__ import print_function
import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg


import spinchain


alphas = np.linspace(0.0,2.0,100)

def getf(alpha):
  def fun(i,j): # function for the exchange
    if i==(j+1): # first neighbors
      return np.identity(3) # AF
    if i==0 and j==(n-1): # first and last
      m = np.zeros((3,3)) # initialize
      m[2,2] = 1.0
      m[0,0] = np.cos(alpha)
      m[1,1] = np.cos(alpha)
      m[0,1] = np.sin(alpha)
      m[1,0] = -np.sin(alpha)
      return m # return matrix
    else: return np.zeros((3,3))
  return fun


es = []

for alpha in alphas:
  n = 20
  spins = [2 for i in range(n)] # spins
  sc = spinchain.Spin_Hamiltonian(spins) # create the chain
#  alpha = np.random.random()
  sc.set_exchange(getf(alpha*np.pi))
#  e = sc.gs_energy(mode="full")
  e = sc.gs_energy(mode="DMRG")
#  print(e0,e1)
#  exit()
  es.append(e)
  print(alpha,e)
#  exit()

es = np.array(es) # array
  
import matplotlib.pyplot as plt

plt.scatter(alphas,es)

plt.show()
