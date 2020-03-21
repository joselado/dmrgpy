from ..algebra import algebra
from .. import multioperator
import scipy.sparse.linalg as slg
from .one2many import one2many
import numpy as np


class EDchain():
    """Generic class for an ED chain"""
    def __init__(self):
        self.operators = dict() # empty dictionary
        self.localdim = [] # empty list
        self.computed_gs = False
    def get_operator(self,name,i=0):
        """Return an operator"""
        if type(name)==multioperator.MultiOperator: # input is a MO
            return multioperator.MO2matrix(name,self)
        return self.operators[(name,i)] # return the operator
    def get_hamiltonian(self):
        """Return the Hamiltonian"""
        return self.get_operator(self.hamiltonian) # return operator
    def gs_energy(self):
        """Return ground state energy"""
        self.h = self.get_hamiltonian()
        return algebra.ground_state(self.h)[0]
    def get_gs(self):
        """Get ground state wavefunction"""
        if self.computed_gs: return self.wf0
        else: 
          e0,wf0 = algebra.ground_state(self.get_hamiltonian())
          self.wf0 = wf0
          self.e0 = e0
          self.computed_gs = True
          return self.wf0
    def vev(self,op):
        """Return a vacuum expectation value"""
        wf0 = self.get_gs()
        op = multioperator.MO2matrix(op,self) # return operator
        return algebra.braket_wAw(wf0,op)
    def get_excited(self,**kwargs):
        """Excited states"""
        h = self.get_hamiltonian()
        return algebra.lowest_eigenvalues(h,**kwargs)
    def create_operator(self,a,i=None,name=None):
        """Create the different operators"""
        ns = self.localdim # local dimensions
        if name is None: raise
        if i is None: raise
        if a.shape[0]!=ns[i]: raise
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex) for n in ns] # identities
        op = one2many(ids,a,i) # one to many body
        self.operators[(name,i)] = op # store in the dictionary
    def get_identity(self):
        ids = [np.identity(n,dtype=np.complex) for n in self.localdim] 
        op = one2many(ids,ids[0],0) # one to many body
        return op
    def get_dynamical_correlator(self,**kwargs):
        from . import dynamics
        return dynamics.get_dynamical_correlator(self,**kwargs)



