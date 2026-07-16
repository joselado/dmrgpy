from __future__ import print_function
from copy import copy as _shallow_copy
import os
import numpy as np
from . import entropy
from . import multioperator
import subprocess

class MPS():
    """Object for an MPS, backed by an opaque in-process extension handle
    (mpscpp2/chain_session.h's MPS, see mpscpp2/bindings.cc) -- nothing is
    read from or written to disk."""
    def __init__(self,MBO=None,name="psi_GS.mps",cpp_handle=None):
        if MBO is None:
            self.path = os.getcwd() # current directory
            self.MBO = None
        else:
            self.path = MBO.path # path to the many body object folder
            self.MBO = MBO
        self.name = name # initial name
        self.cpp_handle = cpp_handle # opaque in-process extension handle
        self.mode = "DMRG" # mode of the object
    def set_MBO(self,MBO):
        """Set the MBO"""
        self.path = MBO.path # path to the many body object folder
        self.MBO = MBO # set the object
    def dot(self,x):
        if self.MBO is not None: return self.MBO.overlap(self,x)
        else: raise
    def overlap(self,x):
        if self.MBO is not None: return self.MBO.overlap(self,x)
        else: raise
    def aMb(self,M,b):
        if self.MBO is not None: return self.MBO.aMb(self,M,b)
        return self.dot(M*b) # workaround
    def __radd__(self,x): return self + x
    def __add__(self,x):
        if x==0: return self # do nothing
        if self.MBO is not None: return self.MBO.summps(self,x)
        else: raise
    def __sub__(self,x):
        return self + (-1)*x
    def __neg__(self,x):
        return (-1)*x
    def __truediv__(self,x): return self*(1./x)
    def copy(self,name=None):
        """Copy this wavefunction"""
        if name is None:
          name = id_generator()+".mps" # create a new name
        # A shallow copy is enough, and deliberately used instead of a deep
        # one: the opaque cpp_handle is never mutated in place (every Chain
        # method takes wf by const reference and returns a brand-new MPS,
        # see chain_session.h), and self.MBO should stay shared -- there is
        # no reason to duplicate the whole chain object just to copy one
        # wavefunction. A real deepcopy would also choke on cpp_handle,
        # which has no pickle/deepcopy support.
        out = _shallow_copy(self)
        out.name = name
        return out
    def __deepcopy__(self,memo):
        """Defer to copy(): see the note there for why a shallow copy is
        already correct and a real deep copy would fail on cpp_handle (or
        recurse needlessly into self.MBO)."""
        return self.copy()
    def write(self,name=None,path=None):
        """No-op: the wavefunction already lives in self.cpp_handle,
        nothing to write to disk. Kept as a callable so existing call
        sites (e.g. self.execute(wf.write)) don't need special-casing."""
        pass
    def get_fermionic_parity(self,**kwargs):
        from .fermionicparity import get_fermionic_parity
        return get_fermionic_parity(self,**kwargs) # parity of the state
    def get_entropy(self,b=None):
        """Compute entanglement entropy in a bond"""
        if b is None: # compute all 
            return np.mean([self.get_entropy(i) for i in range(1,self.MBO.ns)])
        if self.MBO is not None: return self.MBO.get_entropy(self,b=b)
        else: raise
    def get_site_entropy(self,i):
        if self.MBO is not None: return self.MBO.get_site_entropy(self,i)
        else: raise # not implemented
    def get_bond_entropy(self,i,j=None):
        if j is None: j = i + 1
        if self.MBO is not None: return self.MBO.get_bond_entropy(self,i,j)
        else: raise # not implemented
    def get_correlation_matrix(self,**kwargs):
        from . import entanglement
        return entanglement.get_correlation_matrix(self.MBO,
                wf=self,**kwargs)
    def get_four_correlation_tensor(self,**kwargs):
        from . import entanglement
        return entanglement.get_four_correlation_tensor(self,**kwargs)
    def get_correlation_entropy(self,**kwargs):
        from . import entanglement
        return entanglement.get_correlation_entropy_from_wf(self,**kwargs)
    def get_correlation_entropy_density(self,**kwargs):
        from . import entanglement
        return entanglement.get_correlation_entropy_density(self.MBO,
                 wf=self,**kwargs)
    def get_CFT_central_charge(self):
        return entropy.central_charge(self)
    def get_pair_entropy(self,i,j):
        if self.MBO is not None: return self.MBO.get_pair_entropy(self,i,j)
        else: raise # not implemented
    def get_mutual_information(self,i,j):
        if self.MBO is not None: 
            return self.MBO.get_mutual_information(self,i,j)
        else: raise # not implemented
    def rename(self,name):
        self.execute(lambda: subprocess.run(["mv",self.name,name]))
        self.name = name
    def execute(self,f):
        pwd = os.getcwd() # path
        os.chdir(self.path) # go
        f()
        os.chdir(pwd) # go back
    def clean(self):
        self.execute(lambda: subprocess.run(["rm",self.name]))
        del self
    def norm(self):
        return np.sqrt(self.dot(self).real) # norm
    def get_conjugate(self):
        """Return the conjugate wavefunction"""
        from .mpsalgebra import conjugate_mps
        if self.MBO is None: raise
        return conjugate_mps(self.MBO,self)
    def normalize(self,tol=1e-8):
        """Normalize a wavefunction"""
        # This function has problems for too many sites due to
        # the orthogonality catastrophe. Probably MPS should be 
        # orthogonalized at the C++ level
        norm = np.sqrt(self.dot(self).real) # norm
        if norm>tol: return self*(1./norm)
        else: 
            print("WARNING, state is not normalizable. Returning None")
            return None
    def __rmul__(self,A):
        """Multiply by an operator"""
        if self.MBO is not None:
            if type(A)==multioperator.MultiOperator: # MO type
                return self.MBO.applyoperator(A,self) # apply the operator
            elif multioperator.isnumber(A):
                return A*multioperator.identity()*self
            else: raise
        else: raise
    def __mul__(self,x):
        if multioperator.isnumber(x):
            return x*multioperator.identity()*self
        else: raise
    def get_dm(self,**kwargs):
        """COmpute the density matrix"""
        from .dmtk.densitymatrix import dm
        return dm(self,**kwargs)





import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))



from .randommps import random_mps
from .randommps import orthogonal_random_mps
from .randommps import random_product_state









