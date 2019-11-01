# library to solve zn models using ED
import numpy as np
from ..edtk.one2many import one2many


class bosonchain():
    def __init__(self,ns):
        """Initialize"""
        self.nsites = len(ns) # number of sites
        self.ns = ns # list with the maximum number of bosons in each site
        self.create_operators() # initialize operators
    def create_operators(self):
        """Create the different operators"""
        dop = dict() # dictionary
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex) for n in ns] # identities
        ds = [] # empty list
        for n in ns:
            d = np.zeros((n,n),dtype=np.complex)
            for i in range(n-1): d[i,i+1] = 1.0 # one more boson
            ds.append(d) # store
        # now create the many body basis
        for i in range(self.nsites): 
            op = one2many(ids,ds[i])
            dop[("ad",i)] = op # creation
            dop[("a",i)] = op.H # annhilation
        self.operators = dop # store dictionary
    def get_operator(self,name,i=0):
        """Return an operator"""
        return self.operators[(name,i)] # return the operator


