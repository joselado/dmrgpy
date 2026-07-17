
from __future__ import print_function
import scipy.sparse as sparse
import numpy as np
from .. import multioperator
from ..edtk import edchain

maxsize = 300000 # maximum matrix size


def get_dimension(spins):
  size = 1
  for s in spins:
    size *= 2*s+1
  return int(size)


class Spin_chain(edchain.EDchain):
  size = 0 # size of the Hamiltonian
  def build(self,spins):
    """Creates a spin chain, using as input the spins list"""
    if get_dimension(spins)>maxsize:
      print("Surpased maximum allowed dimension for ED",maxsize)
      print("Dimension of the requested Hilbert space",get_dimension(spins))
      raise
    from . import chain
    ###############################
    # now build the spin operators #
    ###############################
    sobj = chain.get_chain(spins) # return a list of classes with the spins
    self.nspins = len(spins)
    self.spins = spins
    self.basis = sobj.basis # store the basis
    self.sxi = sobj.sxi  # store different sx
    self.syi = sobj.syi  # store different sy
    self.szi = sobj.szi  # store different sz
    self.sx = sobj.sx  # store sx
    self.sy = sobj.sy  # store sy
    self.sz = sobj.sz  # store sz
    self.wf0 = None # ground state
    self.e0 = None # ground state energy
    self.hamiltonian = None # Hamiltonian, as a multioperator
    self.size = self.sxi[0].shape[0] # store size of the Hamiltonian
  def get_identity(self):
      return sparse.identity(self.sx.shape[0],dtype=np.complex128)
  def get_operator(self,name,i=0):
      """Return an operator"""
      if type(name)==multioperator.MultiOperator:
          return multioperator.MO2matrix(name,self) # return operator
      else:
        if name=="X" or name=="Sx": return self.sxi[i]
        elif name=="Y" or name=="Sy": return self.syi[i]
        elif name=="Z" or name=="Sz": return self.szi[i]
        elif name=="Id": return self.get_identity()
        else:
            print(name)
            raise
