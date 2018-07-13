import numpy as np

def select_state(sc,indexes):
  """Return the wavefunction which has certain indexes"""
  v = np.array(indexes) # vector
  print("Select state",v)
  vout = np.array([0. for i in range(sc.size)]) # zero vector
  for i in range(sc.size): # loop over basis
    d = v - sc.basis[i] # difference
    if np.sum(np.abs(d))<0.1: # if found
      vout[i] = 1.0
      return vout 


