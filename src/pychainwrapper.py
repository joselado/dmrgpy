import numpy as np


# wrapper function for pychain

def get_full_hamiltonian(self):
  sc = get_pychain(self) # get pychain object
  def get_coupling(i,j):
    """Return the coupling between two sites"""
    for c in self.couplings:
      if i==c.i and j==c.j: return c.g
    return np.zeros((3,3))
  h = sc.add_tensor_interaction(get_coupling) # add interaction
  h = h + sc.add_exchange(self.fields) # add magnetic fields
  return h


def get_pychain(self):
  import pychain.build
  sc = pychain.build.Spin_chain()
  # the pychain library assumes that s=1/2 is for spin one-half
  # whereas in DMRG s = 2 is for S=1/2
  sc.build((np.array(self.sites)-1.)/2.,use_lib=False) 
  return sc

