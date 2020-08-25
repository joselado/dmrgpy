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
    def get_excited_states(self,**kwargs):
        """Excited states"""
        h = self.get_hamiltonian()
        (es,ws) = algebra.lowest_states(h,**kwargs)
        ws = [State(w,self) for w in ws] # transform to states
        return (es,ws)
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
    def get_distribution(self,**kwargs):
        """Return a certain distribution"""
        from . import distribution
        return distribution.get_distribution(self,**kwargs)
    def exponential(self,h,wf):
        """Exponential of a wavefunction"""
        h = self.MO2matrix(h) # convert to matrix 
        return algebra.expm(h)@wf # return
    def MO2matrix(self,m): return multioperator.MO2matrix(m,self)
    def overlap(self,wf1,wf2): return np.dot(np.conjugate(wf1),wf2)
    def applyoperator(self,A,wf): return self.MO2matrix(A)@wf
    def random_state(self):
        """Return a random state"""
        n = self.get_hamiltonian().shape[0] # convert to matrix
        v = np.random.random(n)-.5 + 1j*(np.random.random(n)-.5)
        return State(v,self) # return the state




class State():
    """This is a dummy class to contain states"""
    def __init__(self,v,MBO):
        self.v = v # store the vector
        self.MBO = MBO # store the many-body object
    def __rmul__(self,a):
        """Multiply by something"""
        if type(a)==multioperator.MultiOperator: # multioperator
            A = self.MBO.MO2matrix(a)  # get the matrix
            w = A@self.v # multiply
            return State(w,self.MBO) # create a new object
        elif multioperator.isnumber(a):
            w = a*self.v # multiply
            return State(w,self.MBO) # create a new object
        else: raise # not implemented
    def __add__(self,a):
        if type(a)==State: return State(self.v + a.v,self.MBO)
        else: raise
    def __mul__(self,x):
        if multioperator.isnumber(x): 
            return State(x*self.v,self.MBO)
        else: raise
    def __truediv__(self,a):
        if multioperator.isnumber(a): # number
          return (1./a)*self
        else: raise
    def __sub__(self,a):
        return self + (-1)*a
    def __neg__(self):
        return (-1)*self
    def overlap(self,a):
        if type(a)==State: # state object
            return np.dot(np.conjugate(self.v),a.v)
        else: raise # not implemented
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def dot(self,a):
        return np.sum(np.conjugate(self.v)*a.v)





