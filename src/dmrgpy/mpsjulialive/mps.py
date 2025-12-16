from copy import deepcopy
import os
import numpy as np
from .. import entropy
from .. import multioperator
from .juliasession import Main as Mainjl

class MPS():
    """Object for an MPS"""
    # this is an MPS object using the
    # Julia live environment
    def __init__(self,jlmps,MBO=None):
        self.MBO = light_MBO(MBO) # clean and store
        self.jlmps = jlmps
    def set_MBO(self,MBO):
        """Set the MBO"""
        self.MBO = MBO # set the object
    def dot(self,x):
        return self.overlap(x)
    def overlap(self,x):
        return Mainjl.overlap(self.jlmps,x.jlmps)
    def aMb(self,M,b):
        return self.dot(M*b) # workaround
#        if self.MBO is not None: return self.MBO.aMb(self,M,b)
#        return self.dot(M*b) # workaround
    def __radd__(self,x): return self + x
    def __add__(self,x):
        if x==0: return self # do nothing
        if self.MBO is not None:
            jlmps = Mainjl.summps(self.jlmps,x.jlmps,self.MBO.maxm)
            wf3 = MPS(jlmps,MBO=self.MBO)
            return wf3
        else: raise
    def __sub__(self,x):
        return self + (-1)*x
    def __neg__(self,x):
        return (-1)*x
    def __truediv__(self,x): return self*(1./x)
    def copy(self,name=None):
        """Copy this wavefunction"""
        return MPS(self.jlmps,MBO=self.MBO)
#        out = deepcopy(self) # copy everything
#        return out
    def write(self,name=None,path=None):
        return # dummy method
    def get_entropy(self,b=None):
        """Compute entanglement entropy in a bond"""
        raise # not implemented
        if b is None: # compute all 
            return np.mean([self.get_entropy(i) for i in range(1,self.MBO.ns)])
        if self.MBO is not None: return self.MBO.get_entropy(self,b=b)
        else: raise
    def get_site_entropy(self,i):
        raise # not implemented
    def get_bond_entropy(self,i,j=None):
        raise # not implemented
    def get_correlation_entropy(self,**kwargs):
        from .. import entanglement
        return entanglement.get_correlation_entropy_from_wf(self,**kwargs)
    def get_correlation_entropy_density(self,**kwargs):
        from .. import entanglement
        return entanglement.get_correlation_entropy_density(self.MBO,
                 wf=self,**kwargs)
    def get_CFT_central_charge(self):
        raise # not implemented
    def get_pair_entropy(self,i,j):
        raise # not implemented
    def get_mutual_information(self,i,j):
        raise # not implemented
    def rename(self,name):
        return # dummy method
    def execute(self,f):
        return # dummy method
    def clean(self):
        return # dummy method
    def norm(self):
        return np.sqrt(self.dot(self).real) # norm
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
                from .mpo import MPO
                return MPO(A,MBO=self.MBO)*self
#                return self.MBO.applyoperator(A,self) # apply the operator
            elif multioperator.isnumber(A):
                jlmps = Mainjl.mpstimesscalar(A,self.jlmps)
                return MPS(jlmps,MBO=self.MBO)
#                return A*multioperator.identity()*self
            else: raise
        else: raise
    def __mul__(self,x):
        if multioperator.isnumber(x):
            return x*self
        else: raise
    def get_dm(self,**kwargs):
        """Compute the density matrix"""
        from ..dmtk.densitymatrix import dm
        return dm(self,**kwargs)


def random_mps(self):
    jlmps = Mainjl.random_state(self.jlsites)
    return MPS(jlmps,MBO=self)







def light_MBO(MBO):
    """Create an MBO object with cleaned up data"""
    return MBO
    out = MBO.copy()
    del out.wf0 # remove this
    del out.hamiltonian # remove this
    return out



class LMBO():
    """Light many-body object, with the key methods required"""
    def __init__(self,MBO):
        """Light many body-object"""
        self.nsweeps = MBO.nsweeps
        self.maxm = MBO.maxm
        self.cutoff = MBO.cutoff
        self.sites = MBO.sites
        self.jlsites = MBO.jlsites
    def toMPO(self,H):
        from .mpo import MPO
        return MPO(H,MBO=self)




