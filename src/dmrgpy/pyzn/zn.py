# library to solve zn models using ED
import numpy as np
from ..edtk.one2many import one2many


class znchain():
    def __init__(self,ns):
        """Initialize"""
        self.nsites = len(ns) # number of sites
        self.ns = ns # list with the integers of the Zn model
        self.create_operators() # initialize operators
    def create_operators(self):
        """Create the different operators"""
        dop = dict() # dictionary
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex) for n in ns] # identities
        taus = [] # empty list
        for n in ns:
            tau = np.zeros((n,n),dtype=np.complex)
            omega = 1j*2.*np.pi/n
            for i in range(n): tau[i,i] = omega**i
            taus.append(tau) # store
        # now create the many body basis
        for i in range(self.nsites): dop[("tau",i)] = one2many(ids,taus[i])
        self.operators = dop # store dictionary
    def get_operator(self,name,i=0):
        """Return an operator"""
        return self.operators[(name,i)]


