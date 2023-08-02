# library to solve zn models using ED
import numpy as np
from ..edtk.one2many import one2many
from ..edtk.edchain import EDchain


class BosonChain(EDchain):
    def __init__(self,maxnb):
        """Initialize"""
        self.nsites = len(maxnb) # number of sites
        self.maxnb = maxnb # list with the maximum number of bosons in each site
        EDchain.__init__(self) # parent initialization
        self.create_operators() # initialize operators
    def create_operators(self):
        """Create the different operators"""
        dop = dict() # dictionary
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex_) for n in self.maxnb] # identities
        ds = [] # empty list
        for n in self.maxnb:
            d = np.zeros((n,n),dtype=np.complex_)
            for i in range(n-1): d[i,i+1] = np.sqrt(i+1) # one more boson
            ds.append(d) # store
        # now create the many body basis
        for i in range(self.nsites): 
            op = one2many(ids,ds[i],i)
            dop[("A",i)] = op # creation
            dop[("Adag",i)] = np.transpose(np.conjugate(op)) # annhilation
            dop[("density",i)] = np.transpose(np.conjugate(op))@op
            dop[("N",i)] = np.transpose(np.conjugate(op))@op
        for i in range(self.nsites): 
            op = one2many(ids,ids[i],i)
            dop[("Id",i)] = op
        # now create the occupation operators
        for i in range(self.nsites): # loop over sites
            ops = ids[i]*0.0 # initialize
            for n in range(ops.shape[0]): # loop over occupations
                name = "N"+str(n) # name
                op = ops*0.0 # initialize
                op[n] = 1.0
                op = one2many(ids,op,i) # to many-body
                dop[(name,i)] = op
        # create density operators
        self.operators = dop # store dictionary
    def get_identity(self):
        """Identity operator"""
        ids = [np.identity(n,dtype=np.complex_) for n in self.maxnb] # identities
        return one2many(ids)


class SpinBosonChain(EDchain):
    """Object that simultanously incorporates bosons and spins"""
    def __init__(self,sites):
        """Initialize the object"""
        from ..bosonchain import is_boson,get_site
        nsites = len(sites) # number of sites
        idsites = [get_site(s) for s in sites] # transform the name to a number
        isboson = [is_boson(s) for s in idsites] # if this site is a boson 
        spins = [] # empty list
        bosons = [] # empty list
        for i in range(nsites): # loop over sites
            if isboson[i]: bosons.append(sites[i]) # add to bosons
            else: spins.append(sites[i]) # add to spins
        self.nsites = len(nsites) # number of sites
        print("Spin Boson Chain is not functional yet")
        raise # this is not finished
        EDchain.__init__(self) # parent initialization
        ### Now create all the operators




bosonchain = BosonChain # for backwards compatibility

