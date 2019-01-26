from . import states
from ..pychain.spectrum import ground_state
import numpy as np
from scipy.sparse import csc_matrix


nmax = 15 # maximum number of levels

def get_spinless_hamiltonian(m0,hubbard=None):
    """Compute ground state energy"""
    MBf = MBFermion(len(m0)) # create many body fermion object
    h = MBf.one2many(m0) # get the single body matrix
#    h = states.one2many(m0) # single to many body Hamiltonian
    if hubbard is not None: # if hubbard given
        h = h + MBf.hubbard(hubbard) # add Hubbard term
    return h

def gs_energy(m0,spinless=True,hubbard=None):
    if spinless:
      h = get_spinless_hamiltonian(m0,hubbard=hubbard)
    else: raise
    return ground_state(h)[0] # return GS energy



class MBFermion():
    """
    Class for a many body fermionic Hamiltonian
    """
    def __init__(self,n):
        """
        Initialize the object
        """
        self.n = n # number of different levels
        if n>nmax: raise # too big system
        self.c_dict = dict() # dictionary with the annhilation operators
        self.basis = states.generate_basis(self.n,lambda x: True) # basis
        self.nMB = len(self.basis) # dimension of many body hamiltonian
        self.basis_dict = states.get_dictionary(self.basis) # dictionary
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
    def get_cd(self,i):
        """
        Return the creation operator for site i in the many body basis
        """
        return self.get_c(i).H # return the dagger
    def one2many(self,m0):
        """
        Convert a single body Hamiltonian into a many body one
        """
        m = csc_matrix(([],([],[])),shape=(self.nMB,self.nMB)) # initialize
        for i in range(len(m0)):
          for j in range(len(m0)):
              if abs(m0[i,j])>1e-7:
                m = m + self.get_cd(i)*self.get_c(j)*m0[i,j] # add contribution
        return m # return many body hamiltonian
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
    def correlator(self,pairs,wf):
        """
        Compute a set of correlators for a wavefunction
        """
        out = [] # output list
        for p in pairs: # loop over pairs
            m = self.get_cd(p[0])*self.get_c(p[1]) # get matrix
            raise # not finished yet










