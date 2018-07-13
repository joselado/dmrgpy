from __future__ import print_function
import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg


import spinchain



def getf(alpha,g=1.0):
  def fun(i,j): # function for the exchange
    if i==(j+1): # first neighbors
      return g*np.identity(3) # AF
    if i==0 and j==(n-1): # first and last
      m = np.zeros((3,3)) # initialize
      m[2,2] = 1.0
      m[0,0] = np.cos(alpha)
      m[1,1] = np.cos(alpha)
      m[0,1] = np.sin(alpha)
      m[1,0] = -np.sin(alpha)
      return g*m # return matrix
    else: return np.zeros((3,3))
  return fun


es = []

ns = np.array(range(4,40,2))

for n in ns:
  alpha = 1./20.
  spins = [2 for i in range(n)] # spins
  sc0 = spinchain.Spin_Hamiltonian(spins) # create the chain
  sc1 = spinchain.Spin_Hamiltonian(spins) # create the chain
  # Add the exchange
  g=1.0
  sc0.set_exchange(getf(0.0,g=g))
  sc1.set_exchange(getf(alpha*np.pi,g=g))
  # compute energy
  e0 = sc0.gs_energy(mode="DMRG")
  e1 = sc1.gs_energy(mode="DMRG")
#  print(e0,e1)
#  exit()
  es.append(e1-e0)
  print(alpha,e0,e1,(e1-e0)*n/alpha)
#  exit()

es = np.array(es) # array
  
import matplotlib.pyplot as plt

plt.scatter(ns,es/((alpha/ns)))

plt.show()
