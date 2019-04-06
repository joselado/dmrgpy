# Add the root path of the dmrgpy library
# replace PATH_TO_DMRGPY by your actual path to the library
# PATH_TO_DMRGPY = /home/jose/programs/dmrgpy/src
#import os ; import sys ; sys.path.append(PATH_TO_DMRGPY)

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains


def get_energy(ind):
  """Given a certain ordering of the indexes, return the energy"""
  n = len(ind)
  spins = [2 for i in range(n)] # list with the different spins of your system
  sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object
  sc.maxm = 10 # bond dimension, controls the accuracy
  sc.nsweeps = 3
  def fj(i,j):
      if abs(ind[i]-ind[j])==1: return 1.0
      return 0.0
  sc.set_exchange(fj) # add the exchange couplings
  return sc.gs_energy()




# compute the energy with a random sorting

import random
n = 10 # length of your list (number of spins of a many body system)
ind0 = [i for i in range(n)] # initial sequence
ntries = 20 # do a certain number of random trials
emin = 1e10 # initalize as a big number
indmin = None # initialize as None


# Now do several random resorting of the indexes and compute their energy
for i in range(ntries):
  ind = random.sample(ind0,n) # randomize the sequence
  energy = get_energy(ind) # get the energy for this sequence
  print("Computing",energy)
  if energy<emin: # if better energy, store it
      emin = energy
      indmin = ind
print("\n\n\n")
# print the stochastically optimized result
print("Minimum energy found",emin)
print("Best sorting found",indmin)
print("\n")
print("Now the best energy that could be found (by construction)")
print("Energy for optimal sorting",get_energy(ind0))
print("Optimal sorting",ind0)



