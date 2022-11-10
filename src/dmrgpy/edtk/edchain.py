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
        elif type(name)==str: # string
            return self.operators[(name,i)] # return the operator
        else: # unrecognized type
            print("Unrecognized operator in EDchain",type(name))
    def get_hamiltonian(self):
        """Return the Hamiltonian"""
        return self.get_operator(self.hamiltonian) # return operator
    def gs_energy(self):
        """Return ground state energy"""
        return self.get_excited(n=1)[0]
    def get_gs(self,array_mode=True):
        """Get the ground state"""
#        if array_mode: 
#            print("This will be deprecated")
#            return self.get_gs_array()
#        else: 
        return State(self.get_gs_array(),self)
    def get_gs_array(self):
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
        wf0 = self.get_gs_array()
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
        wf = State(wf,self) # convert to State
        h = self.MO2matrix(h) # convert to matrix 
        return State(algebra.expm(h)@wf.v,self) # return
    def MO2matrix(self,m): return multioperator.MO2matrix(m,self)
    def overlap(self,wf1,wf2):
        return wf1.dot(wf2) 
    def applyoperator(self,A,wf): 
        wf = State(wf,self)
        return A*wf #return self.MO2matrix(A)@wf
    def random_state(self):
        """Return a random state"""
        n = self.get_operator("Id").shape[0] # dimension
        v = np.random.random(n)-.5 + 1j*(np.random.random(n)-.5)
        return State(v,self).normalize() # return the state




class State():
    """This is a dummy class to contain states"""
    def __init__(self,v,MBO):
        if type(v)==State:
            self.v = v.v.copy()
            self.MBO = v.MBO
        else:
            self.v = v # store the vector
            self.MBO = MBO # store the many-body object
    def __rmul__(self,a):
        """Multiply by something"""
        from ..algebra.algebra import ismatrix
        if type(a)==multioperator.MultiOperator: # multioperator
            A = self.MBO.MO2matrix(a)  # get the matrix
            w = A@self.v # multiply
            return State(w,self.MBO) # create a new object
        elif multioperator.isnumber(a):
            w = a*self.v # multiply
            return State(w,self.MBO) # create a new object
        elif ismatrix(a):
            w = a@self.v
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
    def aMb(self,M,b): return self.dot(M*b)
    def normalize(self):
        norm = np.sqrt(self.dot(self).real) # norm
        if norm>1e-8: return self/norm
        else: return None
    def get_correlation_entropy(self,**kwargs):
        from .. import entanglement
        return entanglement.get_correlation_entropy_from_wf(self,**kwargs)
    def applyinverse(self,a,**kwargs):
        if type(a)==multioperator.MultiOperator: # multioperator
            A = self.MBO.MO2matrix(a)  # get the matrix
        elif type(a)==EDOperator: # multioperator
            A = a.SO  # get the matrix
        else: raise
        w = algebra.applyinverse(A,self.v)
        return State(w,self.MBO) # create a new object





class EDOperator():
    """This is a dummy class for operators, so that it resembles the
    tensor network Static Operator. Useful for testing purposes"""
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        if type(MO)==multioperator.MultiOperator: # multioperator
            self.SO = MBO.get_operator(MO) # generate the static operator
        elif type(MO)==EDOperator:
            self.SO = MO.SO.copy() # dummy copy
        else: 
            print("Unrecognized type in EDOperator",type(MO))
            raise
    def __mul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==State: # input is an MPS
            return State(self.SO@v.v,self.MBO)
        elif type(v)==EDOperator: # input is an MPO
            out = self.copy() # copy
            out.SO = self.SO@v.SO # matrix multiplication
            return out
        elif type(v)==MultiOperator: # input is a multioperator
            return self*EDOperator(v,self.MBO)
        else: 
            print("Incompatible object",type(v))
            raise
    def __rmul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MultiOperator: # input is a multioperator
            return EDOperator(v,self.MBO)*self
        else: raise
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def trace(self):
        from ..algebra.algebra import trace
        return trace(self.SO)
    def get_dagger(self):
        out = self.copy()
        from ..algebra.algebra import dagger
        out.SO = dagger(self.SO)
        return out









