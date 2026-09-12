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
        ids = [np.identity(n,dtype=np.complex128) for n in self.maxnb] # identities
        ds = [] # empty list
        for n in self.maxnb:
            d = np.zeros((n,n),dtype=np.complex128)
            for i in range(n-1): d[i,i+1] = np.sqrt(i+1) # one more boson
            ds.append(d) # store
        # now create the many body basis
        for i in range(self.nsites): 
            op = one2many(ids,ds[i],i)
            dop[("A",i)] = op # annihilation (d[i,i+1]=sqrt(i+1) lowers n)
            dop[("Adag",i)] = np.transpose(np.conjugate(op)) # creation
            dop[("density",i)] = np.transpose(np.conjugate(op))@op
            dop[("N",i)] = np.transpose(np.conjugate(op))@op
        for i in range(self.nsites): 
            op = one2many(ids,ids[i],i)
            dop[("Id",i)] = op
        # now create the occupation operators: "N<k>" is the projector
        # |k><k| onto the state with exactly k bosons on this site, the
        # same operator the DMRG backends build under that name
        # (pyitensor/sites/boson.py's build_matrix(dim,[(k+1,k+1,1.0)]),
        # ITensor's own BosonFourSite). Note op[n,n], not op[n]: op is a
        # (d,d) array, so op[n] would assign the whole ROW n and give
        # sum_m |n><m| -- non-Hermitian, not idempotent, and boson-number
        # changing. That is invisible on a number-conserving Hamiltonian
        # (every off-diagonal piece has exactly zero expectation value in
        # a state of definite total N) and gives negative "probabilities"
        # as soon as the Hamiltonian breaks boson-number conservation.
        for i in range(self.nsites): # loop over sites
            ops = ids[i]*0.0 # initialize
            for n in range(ops.shape[0]): # loop over occupations
                name = "N"+str(n) # name
                op = ops*0.0 # initialize
                op[n,n] = 1.0 # the projector |n><n| on this site
                op = one2many(ids,op,i) # to many-body
                dop[(name,i)] = op
        # create density operators
        self.operators = dop # store dictionary
    def get_identity(self):
        """Identity operator"""
        ids = [np.identity(n,dtype=np.complex128) for n in self.maxnb] # identities
        return one2many(ids)


class SpinBosonChain(EDchain):
    """Object that simultanously incorporates bosons and spins"""
    def __init__(self,sites):
        """Initialize the object"""
        # This class was never finished (the operator construction below
        # the raise does not exist), so say so up front rather than
        # falling over on the way there -- the line that used to raise
        # was `self.nsites = len(nsites)`, a TypeError on an int that
        # hid the intended message behind an unrelated traceback.
        raise NotImplementedError(
            "the ED backend for SpinBoson_Chain is not implemented yet; "
            "use mode=\"DMRG\" (itensor_version=3 or \"python\") on this "
            "chain, or Bosonic_Chain for a purely bosonic one")
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

