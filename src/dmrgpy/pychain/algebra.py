from scipy.sparse import csc_matrix as csc
import numpy as np

def Av(A,v):
  return csc(A)*csc(v) # return A*v


def braket(v,w):
  return (csc(v).H*csc(w)).todense()[0,0] # return v*w


def vAw_braket(v,A,w):
  return (csc(v)*csc(A)*csc(w).H).todense()[0,0]
