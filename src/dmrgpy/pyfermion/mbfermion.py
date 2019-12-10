from . import states
from ..pychain.spectrum import ground_state
import numpy as np
from scipy.sparse import csc_matrix,identity
import scipy.sparse.linalg as slg
from ..algebra import algebra
from .. import operatornames
from .. import multioperator
from .. import funtk


nmax = 15 # maximum number of levels

def get_spinless_hamiltonian(m0,hubbard=None):
    """Compute ground state energy"""
    MBf = MBFermion(len(m0)) # create many body fermion object
    MBf.add_hopping(m0)
#    h = MBf.one2many(m0) # get the single body matrix
#    h = states.one2many(m0) # single to many body Hamiltonian
    if hubbard is not None: # if hubbard given
        MBf.add_hubbard(hubbard)
    return MBf.h

def gs_energy(m0,spinless=True,hubbard=None):
    if spinless:
      h = get_spinless_hamiltonian(m0,hubbard=hubbard)
    else: raise
    return ground_state(h)[0] # return GS energy



class MBFermion():
    """
    Class for a many body fermionic Hamiltonian
    """
    def __init__(self,n,fconf = lambda x: True):
        """
        Initialize the object
        """
        self.n = n # number of different levels
        if n>nmax: raise # too big system
        self.c_dict = dict() # dictionary with the annhilation operators
        self.basis = states.generate_basis(self.n,fconf) # basis
        self.nMB = len(self.basis) # dimension of many body hamiltonian
        self.basis_dict = states.get_dictionary(self.basis) # dictionary
        self.h = csc_matrix(([],([],[])),shape=(self.nMB,self.nMB)) # Hamil
    def get_c(self,i):
        """
        Return the annhilation operator for site i in the many body basis
        """
        if i in self.c_dict: # if already computed
            return self.c_dict[i] # return matrix
        else: # not computed yet
            m = states.destroy(self.basis,self.basis_dict,self.n,i)
            self.c_dict[i] = m # store matrix
            return m
    def clean(self):
        """
        Initialize the Hamiltonian
        """
        self.h = self.get_zero()
    def get_zero(self):
        """
        Return the zero matrix
        """
        return csc_matrix(([],([],[])),shape=(self.nMB,self.nMB))
    def get_identity(self):
        """
        Return the identity
        """
        return identity(self.h.shape[0],dtype=np.complex)
    def add_hopping(self,m):
        """
        Add a single particle term to the Hamiltonian
        """
        self.h = self.h + self.one2many(m) # add contribution
    def add_multioperator(self,m):
        """
        Add a multioperator Hamiltonian
        """
        self.h = self.h + self.get_operator(m) # add the operator
    def get_hopping(self,m):
        """
        Return Hopping matrix
        """
        return self.one2many(m) # return the matrix
    def add_hubbard(self,hubbard):
        """
        Add a Hubbard term to the hamiltonian
        """
        if hubbard is None: return
        self.h = self.h + self.hubbard(hubbard) # add Hubbard term
    def add_vijkl(self,f):
        """
        Add a generalized interaction
        """
        self.h = self.h + self.get_vijkl(f)
    def get_vijkl(self,f):
        """
        Return the generalized interaction
        """
        return get_vijkl(self,f)
    def get_gs(self):
        """
        Return the ground state
        """
        e,wf = ground_state(self.h) # return GS
        self.energy = e # store energy
        self.wf0 = wf # store wavefunction
        return self.energy
    def get_excited(self,**kwargs):
        """Excited states"""
        return algebra.lowest_eigenvalues(MBF.h,**kwargs)
    def get_cd(self,i):
        """
        Return the creation operator for site i in the many body basis
        """
        return self.get_c(i).H # return the dagger
    def get_density(self,i):
        """
        Return the density operator
        """
        return self.get_cd(i)@self.get_c(i)
    def one2many(self,m0):
        """
        Convert a single body Hamiltonian into a many body one
        """
        m = csc_matrix(([],([],[])),shape=(self.nMB,self.nMB)) # initialize
        for i in range(len(m0)):
          for j in range(len(m0)):
              if abs(m0[i,j])>1e-7:
                m = m + self.get_cd(i)@self.get_c(j)*m0[i,j] # add contribution
        return m # return many body hamiltonian
    def get_pairing(self,fun):
        """
        Return a pairing term
        """
        m = self.get_zero() # get zero matrix
        out = funtk.fun2list(fun,self.n) # get list of pairings
        for o in out:
            i,j,delta = o[0],o[1],o[2] # get the parameters
            mt = self.get_c(i)@self.get_c(j)*delta
            m = m + mt + mt.H # add contributions
        return m # return matrix
    def add_pairing(self,fun):
        self.h = self.h + self.get_pairing(fun) # add contribution
    def hubbard(self,m0):
        """
        Return the many body matrix for certain hubbard couplings
        """
        m = csc_matrix(([],([],[])),shape=(self.nMB,self.nMB)) # initialize
        for i in range(len(m0)):
          for j in range(len(m0)):
              if abs(m0[i,j])>1e-7:
                ci = self.get_c(i)
                ni = ci.H*ci
                cj = self.get_c(j)
                nj = cj.H*cj
                m = m + ni*nj*m0[i,j] # add contribution
        return m # return many body hamiltonian
    def get_correlator(self,name="",pairs=[]):
        """
        Compute a set of correlators for a wavefunction
        """
        namei,namej = operatornames.recognize(name) # get the operator
        out = [] # empty list
        self.get_gs() # get ground state
        for p in pairs: # loop over pairs
            A = self.get_operator(namei,p[0]) # get matrix
            A = A@self.get_operator(namej,p[1]) # get matrix
            out.append(algebra.braket_wAw(self.wf0,A))
        return np.array(out) # return array
    def vev(self,A):
        """Return the ground state expectation value"""
        m = self.get_operator(A) # return the operator
        self.get_gs() # get ground state
        return algebra.braket_wAw(self.wf0,m) # return the overlap
    def excited_vev(self,A,**kwargs):
        m = self.get_operator(A) # return the operator
        wfs = algebra.lowest_eigenvectors(self.h,**kwargs)
        return np.array([algebra.braket_wAw(wf,m) for wf in wfs])
    def get_operator(self,name,i=0):
        """
        Return a certain operator
        """
        if name is None: return self.get_zero() # return zero operator
        from .. import multioperator
        if type(name)==multioperator.MultiOperator:
            return multioperator.MO2matrix(name,self) # return operator

#        elif type(name)==multioperator.MultiOperator: # Multioperator
#            out = self.get_identity()*0. # initialize
#            for o in name.op: # loop over operators
#                m = self.get_identity()
#                m = m*o[0] # get coefficient
#                for j in range(1,len(o)):
#                  m = m@self.get_operator(o[j][0],int(o[j][1]))
#                out = out + m
#            return out # return operator
        ### conventional procedure ###
        elif name=="density" or name=="N": return self.get_density(i)
        elif name=="C": return self.get_c(i)
        elif name=="Cdag": return self.get_cd(i)
        else: raise
    def get_dynamical_correlator(self,i=0,j=0,
            es=np.linspace(-1.0,10,500),delta=1e-1,
            name="densitydensity"):
        """
        Compute the dynamical correlator
        """
        from ..algebra import kpm
        from .. import operatornames
        self.get_gs() # compute ground state
        if type(name[0])==multioperator.MultiOperator: # multioperator
          A = self.get_operator(name[0])
          B = self.get_operator(name[1])
        else:
          namei,namej = operatornames.recognize(name) # get the operator
          namei = operatornames.hermitian(namei) # get the dagger
          A = self.get_operator(namei,i)
          B = self.get_operator(namej,j)
#        A = self.get_cd(i) # first operator
#        B = self.get_cd(j) # second operator
        vi = A@self.wf0 # first wavefunction
        vj = B@self.wf0 # second wavefunction
        m = -identity(self.h.shape[0])*self.energy+self.h # matrix to use
        emax = slg.eigsh(self.h,k=1,ncv=20,which="LA")[0] # upper energy
        scale = np.max([np.abs(self.energy),np.abs(emax)])*3.0
        n = int(scale/delta) # number of polynomials
        (xs,ys) = kpm.dm_vivj_energy(m,vi,vj,scale=scale,
                                    npol=n*4,ne=n*10,x=es)
        return xs,np.conjugate(ys)/scale*np.pi*2 # return correlator






def get_vijkl(self,f):
    """
    Return a generalized interaction in the many body basis
    """
    m = self.get_zero()
    if f is None: return m
    for i in range(self.n):
      for j in range(self.n):
        for k in range(self.n):
          for l in range(self.n):
              c = f(i,j,k,l) # get the value
              if np.abs(c)>1e-8: # non zero
                  m = m + c*self.get_cd(i)@self.get_c(j)@self.get_cd(k)@self.get_c(l)
    return m








