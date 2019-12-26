# library to solve zn models using ED
import numpy as np
from ..edtk.one2many import one2many
from ..edtk.edchain import EDchain


class bosonchain(EDchain):
    def __init__(self,maxnb):
        """Initialize"""
        self.nsites = len(maxnb) # number of sites
        self.maxnb = maxnb # list with the maximum number of bosons in each site
        EDchain.__init__(self)
        self.create_operators() # initialize operators
    def create_operators(self):
        """Create the different operators"""
        dop = dict() # dictionary
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex) for n in self.maxnb] # identities
        ds = [] # empty list
        for n in self.maxnb:
            d = np.zeros((n,n),dtype=np.complex)
            for i in range(n-1): d[i,i+1] = np.sqrt(i+1) # one more boson
            ds.append(d) # store
        # now create the many body basis
        for i in range(self.nsites): 
            op = one2many(ids,ds[i],i)
            dop[("A",i)] = op # creation
            dop[("Adag",i)] = np.transpose(np.conjugate(op)) # annhilation
            dop[("density",i)] = np.transpose(np.conjugate(op))@op
        # create density operators
        self.operators = dop # store dictionary
    def get_identity(self):
        ids = [np.identity(n,dtype=np.complex) for n in self.maxnb] # identities
        return one2many(ids)




