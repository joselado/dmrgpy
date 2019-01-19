from scipy.sparse import csc_matrix as csc
import numpy as np

def Av(A,v,T=False):
  """Compute matrix times vector"""
  if T: return csc(A)*csc(v).transpose() # return A*v
  else: return csc(A)*csc(v) # return A*v


def braket(v,w):
  """This is to compute the braket between two vectors"""
  return (csc(v).H*csc(w)).todense()[0,0] # return v*w


def vAw_braket(v,A,w):
  """Compute the braket between a vector and an operator"""
  return (csc(v)*csc(A)*csc(w).H).todense()[0,0]
