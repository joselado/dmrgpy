import numpy as np

def linear_chain(n=5,s=.5):
  """Linear chain"""
  spins = [.5 for i in range(n)]
  return spins


def kitaev_coupling():
  """Return a list of matrices with the kitaev coupling"""
  m0 = np.matrix([[0. for i in range(3)] for j in range(3)]) # zero matrix
  ts = [m0.copy() for i in range(3)] # couplings
  ts[0][0,0] = 1.0
  ts[1][1,1] = 1.0
  ts[2][2,2] = 1.0
  return ts

