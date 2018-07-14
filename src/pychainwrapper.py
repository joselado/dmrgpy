import numpy as np


# wrapper function for pychain

def get_full_hamiltonian(self):
  sc = get_pychain(self) # get pychain object
  h = sc.add_tensor_interaction(self.get_coupling) # add interaction
#  print(self.fields)
  h = h + sc.add_exchange(self.fields) # add magnetic fields
  return h


def get_pychain(self):
  import pychain.build
  sc = pychain.build.Spin_chain()
  # the pychain library assumes that s=1/2 is for spin one-half
  # whereas in DMRG s = 2 is for S=1/2
  sc.build((np.array(self.spins)-1.)/2.,use_lib=False) 
  return sc

