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
        """Create the different operator"""
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
    def get_dynamical_correlator(self,
            es=np.linspace(-1.0,10,500),delta=1e-1,
            name=None):
        """
        Compute the dynamical correlator
        """
        from ..algebra import kpm
        self.get_gs() # compute ground state
        if name is None: raise
        if type(name[0])==multioperator.MultiOperator: # multioperator
          A = name[0].get_dagger() # dagger
          A = self.get_operator(A)
          B = self.get_operator(name[1])
        else:
          raise # this is no longer used
        vi = A@self.wf0 # first wavefunction
        vj = B@self.wf0 # second wavefunction
        h = self.get_hamiltonian()
        m = -np.identity(h.shape[0])*self.e0+h # matrix to use
        emax = slg.eigsh(h,k=1,ncv=20,which="LA")[0] # upper energy
        scale = np.max([np.abs(self.e0),np.abs(emax)])*3.0
        n = 10*int(scale/delta) # number of polynomials
        (xs,ys) = kpm.dm_vivj_energy(m,vi,vj,scale=scale,
                                    npol=n*4,ne=n*10,x=es)
        return xs,np.conjugate(ys) # return correlator





