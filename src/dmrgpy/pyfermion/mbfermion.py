from . import states
from ..pychain.spectrum import ground_state
import numpy as np
from scipy.sparse import csc_matrix,identity
import scipy.sparse.linalg as slg
from ..algebra import algebra
from .. import operatornames
from .. import multioperator
from .. import funtk
from ..edtk import edchain

nmax = 20 # maximum number of levels

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



class MBFermion(edchain.EDchain):
    """
    Class for a many body fermionic Hamiltonian
    """
    def __init__(self,n,nf=None):
        """
        Initialize the object
        """
        super().__init__()
        if nf is not None: raise
        fconf = lambda x: True
        self.n = n # number of different levels
        if n>nmax: raise # too big system
        self.c_dict = dict() # dictionary with the annhilation operators
        self.basis = states.generate_basis(self.n,fconf) # basis
        self.nMB = len(self.basis) # dimension of many body hamiltonian
        self.basis_dict = states.get_dictionary(self.basis) # dictionary
        self.h = csc_matrix(([],([],[])),shape=(self.nMB,self.nMB)) # Hamil
        self.C = [self.get_c(i) for i in range(n)]
        self.Cdag = [self.get_cd(i) for i in range(n)]
    def get_c(self,i):
        """
        Return the annhilation operator for site i in the many body basis
        """
        if i in self.c_dict: # if already computed
            return self.c_dict[i] # return matrix
        else: # not computed yet
            m = states.destroy(self.basis,self.basis_dict,self.n,i)
#            m = m.H
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
        self.hamiltonian = m
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
    def get_excited(self,**kwargs):
        """Excited states"""
        return algebra.lowest_eigenvalues(self.h,**kwargs)
    def get_cd(self,i):
        """
        Return the creation operator for site i in the many body basis
        """
        return np.conjugate(self.get_c(i)).T # return the dagger
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
    def excited_vev(self,A,**kwargs):
        m = self.get_operator(A) # return the operator
        wfs = algebra.lowest_eigenvectors(self.h,**kwargs)
        return np.array([algebra.braket_wAw(wf,m) for wf in wfs])
    def test(self):
        test_commutation(self)
    def get_gs(self,**kwargs):
        out = super().get_gs(**kwargs)
        return Fermionic_State(out,out.MBO)
    def get_operator(self,name,i=0):
        """
        Return a certain operator
        """
        if name is None: return self.get_zero() # return zero operator
        from .. import multioperator
        if type(name)==multioperator.MultiOperator:
            return multioperator.MO2matrix(name,self) # return operator

        elif name=="density" or name=="N": return self.get_density(i)
        elif name=="C": return self.get_c(i)
        elif name=="Id": return self.get_identity()
        elif name=="Cdag": return self.get_cd(i)
        else: 
            print("Unrecognised operator",name)
            raise





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



def test_commutation(self):
    """Perform a test of the commutation relations"""
    C = [self.get_operator("C",i) for i in range(self.n)]
    Cdag = [self.get_operator("Cdag",i) for i in range(self.n)]
    Id = self.get_operator("Id")
    n = len(C) # number of sites
    ntries = 8 # number of tries
    for ii in range(ntries):
        i = np.random.randint(n)
        j = np.random.randint(n)
        d = Cdag[i]@Cdag[j] + Cdag[j]@Cdag[i]
        if np.max(np.abs(d))>1e-7: 
            print("Cdag,Cdag failed",i,j)
            raise
        d = C[i]@C[j] + C[j]@C[i]
        if np.max(np.abs(d))>1e-7: 
            print("C,C failed",i,j)
            raise
        if i==j: d = Cdag[i]@C[j] + C[j]@Cdag[i] - Id
        else: d = Cdag[i]@C[j] + C[j]@Cdag[i] 
        if np.max(np.abs(d))>1e-7:
            print("Cdag,C failed",i,j)
            raise
    print("Commutation test passed")



class Fermionic_State(edchain.State):
    """Special object for many-body fermions"""
    def get_dm(self,**kwargs):
        """Special fucntion to compute density matrix"""
        from ..dmtk.densitymatrix import dm_fermionic
        return dm_fermionic(self,**kwargs)


