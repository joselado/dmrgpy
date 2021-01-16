import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .manybodychain import Many_Body_Chain
from .pyboson import boson

class Bosonic_Chain(Many_Body_Chain):
    """Bosonic Hamiltonian"""
    def __init__(self,n,maxnb=None):
        if maxnb is None: maxnb = [4 for i in range(n)] # maximum # of bosons
        self.maxnb = maxnb # maximum number of bosons
        Many_Body_Chain.__init__(self,[104 for i in range(n)])
        self.use_ampo_hamiltonian = True # use ampo
        self.N = [self.get_operator("N",i) for i in range(self.ns)]
        self.D0 = [self.get_operator("N0",i) for i in range(self.ns)]
        self.D1 = [self.get_operator("N1",i) for i in range(self.ns)]
        self.D2 = [self.get_operator("N2",i) for i in range(self.ns)]
        self.D3 = [self.get_operator("N3",i) for i in range(self.ns)]
        self.A = [self.get_operator("A",i) for i in range(self.ns)]
        self.Adag = [self.get_operator("Adag",i) for i in range(self.ns)]
    def get_ED_obj(self):
        """Return the associated ED object"""
        if not self.has_ED_obj: # not computed
          if np.exp(np.sum(np.log(self.maxnb)))>10000: raise
          out = boson.bosonchain(self.maxnb)
          out.hamiltonian = self.hamiltonian
          self.ed_obj = out # store object
          return out
        else: return self.ed_obj


from .spinchain import get_site as get_site_spin

def get_site(label):
    out = get_site_spin(label) # get the spin site
    if out is None: # this is not a spin, try a boson
        if label=="B": return 104
        elif label=="B4": return 104
        else: raise
    return out # otherwise


def is_boson(n):
    """Check if a certain site is a boson"""
    if 1<n<10: return False # nope, this is a spin
    elif 101<n<200: return True # yes
    else: 
        print("Unrecognized site",n)
        raise # unkown

class SpinBoson_Chain(Many_Body_Chain):
    """Bosonic Hamiltonian"""
    def __init__(self,sitesin,n=None,maxnb=None):
#        if maxnb is None: maxnb = [4 for i in range(n)] # maximum # of bosons
#        self.maxnb = maxnb # maximum number of bosons
        sites = [get_site(s) for s in sitesin] # get the labels
        Many_Body_Chain.__init__(self,sites) # initialize
        self.use_ampo_hamiltonian = True # use ampo
        self.N = [self.get_operator("N",i) for i in range(self.ns)]
        self.D0 = [self.get_operator("N0",i) for i in range(self.ns)]
        self.D1 = [self.get_operator("N1",i) for i in range(self.ns)]
        self.D2 = [self.get_operator("N2",i) for i in range(self.ns)]
        self.D3 = [self.get_operator("N3",i) for i in range(self.ns)]
        self.A = [self.get_operator("A",i) for i in range(self.ns)]
        self.Adag = [self.get_operator("Adag",i) for i in range(self.ns)]
        self.Sx = [self.get_operator("Sx",i) for i in range(self.ns)]
        self.Sy = [self.get_operator("Sy",i) for i in range(self.ns)]
        self.Sz = [self.get_operator("Sz",i) for i in range(self.ns)]
        # now depurate the operators
        for i in range(self.ns): # loop over sites
            if is_boson(sites[i]): # for bosonic sites
                self.Sx[i] = 0 # set to zero
                self.Sy[i] = 0 # set to zero
                self.Sz[i] = 0 # set to zero
            else: # for spin sites
                self.N[i] = 0 # set to zero
                self.D0[i] = 0 # set to zero
                self.D1[i] = 0 # set to zero
                self.D2[i] = 0 # set to zero
                self.D3[i] = 0 # set to zero
                self.A[i] = 0 # set to zero
                self.Adag[i] = 0 # set to zero
    def get_ED_obj(self):
        """Return the associated ED object"""
        raise # not implemented yet!
        if not self.has_ED_obj: # not computed
          if np.exp(np.sum(np.log(self.maxnb)))>10000: raise
          out = boson.bosonchain(self.maxnb)
          out.hamiltonian = self.hamiltonian
          self.ed_obj = out # store object
          return out
        else: return self.ed_obj

