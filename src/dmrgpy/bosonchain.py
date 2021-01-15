import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .manybodychain import Many_Body_Chain
from .pyboson import boson

class Bosonic_Chain(Many_Body_Chain):
    """Bosonic Hamiltonian"""
    def __init__(self,n,maxnb=None):
        if maxnb is None: maxnb = [4 for i in range(n)] # maximum # of bosons
        self.maxnb = maxnb # maximum number of bosons
        Many_Body_Chain.__init__(self,[104 for i in range(n)])
        self.use_ampo_hamiltonian = True # use ampo
        self.N = [self.get_operator("N",i) for i in range(self.ns)]
        self.D0 = [self.get_operator("N0",i) for i in range(self.ns)]
        self.D1 = [self.get_operator("N1",i) for i in range(self.ns)]
        self.D2 = [self.get_operator("N2",i) for i in range(self.ns)]
        self.D3 = [self.get_operator("N3",i) for i in range(self.ns)]
        self.A = [self.get_operator("A",i) for i in range(self.ns)]
        self.Adag = [self.get_operator("Adag",i) for i in range(self.ns)]
    def get_ED_obj(self):
        """Return the associated ED object"""
        if not self.has_ED_obj: # not computed
          if np.exp(np.sum(np.log(self.maxnb)))>10000: raise
          out = boson.bosonchain(self.maxnb)
          out.hamiltonian = self.hamiltonian
          self.ed_obj = out # store object
          return out
        else: return self.ed_obj


