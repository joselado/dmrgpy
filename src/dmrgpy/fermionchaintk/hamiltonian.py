import numpy as np



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







