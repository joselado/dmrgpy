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
