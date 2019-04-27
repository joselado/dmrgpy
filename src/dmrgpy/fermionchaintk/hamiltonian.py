import numpy as np
from . import mop



def set_swave_pairing_spinful(self,fun):
    """
    Add onsite swave pairing to a spinful Hamiltonian
    The pairing term is of the form
    Delta_i c_{i,up} c_{i,down} + h.c.
    """
    def fp(i,j):
        if i//2==j//2 and i!=j: # same site, different spins
            if i<j: return fun(i//2) # pairing in that site
            else: return -fun(i//2) # pairing in that site
        else: return 0.0
    self.set_pairings_MB(fp) # set pairing



def set_hubbard_spinful(self,fun):
    """
    Add Hubbard interation in a spinful manner
    The Hubbard term will be defined as
    n_i n_j, with n_i = n_{i,up} + n_{i,,down}
    """
    def fh(i,j):
        """Return Hubbard"""
        ii = i//2 # index of the site without spin
        jj = j//2 # index of the site without spin
        return fun(ii,jj) # return the hubbard term
    self.set_hubbard_spinless(fh) # set hubbard



def set_hubbard_spinless(self,fun):
    """
    Hubbard term for spinless fermions
    """
    self.set_hubbard_MB(fun) # set hubbard


def set_exchange_spinful(self,fun):
    """
    Add exchange coupling to a spinful fermionic Hamiltonian
    """
    out = mop.get_zero("hamiltonian_multioperator")
    for i in range(self.ns//2): # loop over sites
      for j in range(self.ns//2): # loop over sites
          if np.max(np.abs(obj2matrix(fun(i,j)-fun(j,i))))>1e-8: raise
          m = obj2matrix(fun(i,j)) # get the matrix
          for k in range(3): # loop 
            for l in range(3): # loop
              out = out + m[k,l]*mop.get_si(j=k,i=i)*mop.get_si(j=l,i=j)
#    self.hamiltonian_multioperator = out + self.hamiltonian_multioperator
    from .. import multioperator
    f,mo = multioperator.MO2vijkl(out) # transform into vijkl
    self.set_vijkl(f) # set the function
    if mo is not None:
      self.hamiltonian_multioperator = mo + self.hamiltonian_multioperator







def obj2matrix(a):
    m = np.zeros((3,3))
    m[0][0] = 1.0
    m[1][1] = 1.0
    m[2][2] = 1.0
    if np.array(a).shape==(3,3): return a # return matrix
    else: return a*m  # return identity times input

