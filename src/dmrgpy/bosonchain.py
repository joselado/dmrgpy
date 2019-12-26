import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .manybodychain import Many_Body_Hamiltonian
from .pyboson import boson

class Bosonic_Hamiltonian(Many_Body_Hamiltonian):
    """Bosonic Hamiltonian"""
    def __init__(self,n,maxnb=None):
        if maxnb is None: maxnb = [6 for i in range(n)] # maximum # of bosons
        self.maxnb = maxnb # maximum number of bosons
        Many_Body_Hamiltonian.__init__(self,[-1 for i in range(n)])
    def get_ED_obj(self):
        """Return the associated ED object"""
        if not self.has_ED_obj: # not computed
          if np.exp(np.sum(np.log(self.maxnb)))>10000: raise
          out = boson.bosonchain(self.maxnb)
          out.hamiltonian = self.hamiltonian
          self.ed_obj = out # store object
          return out
        else: return self.ed_obj


