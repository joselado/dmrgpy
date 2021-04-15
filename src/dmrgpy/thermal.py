# this library contains routines to perform thermal calculations
import numpy as np

from .spinchain import Spin_Chain

class Thermal_Spin_Chain():
    def __init__(self,sites,T=0.1,**kwargs):
        sitesT = []
        for s in sites: 
            sitesT += [s]
            sitesT += [s]
        n = len(sites) # number of sites
        self.MBChain = Spin_Chain(sitesT) # get the chain
        self.computed_gs = False
        self.T = T # temperature
        Sx = [self.MBChain.Sx[2*i] for i in range(n)] 
        Sy = [self.MBChain.Sy[2*i] for i in range(n)] 
        Sz = [self.MBChain.Sz[2*i] for i in range(n)] 
        self.all_Sx = self.MBChain.Sx
        self.all_Sy = self.MBChain.Sy
        self.all_Sz = self.MBChain.Sz
        self.Sx = Sx
        self.Sy = Sy
        self.Sz = Sz
        self.wf0 = None
        self.mode = "DMRG"
    def get_gs(self):
        """Compute the ground state"""
        if self.computed_gs: 
            return self.wf0
        else:
            h = 0. # initialize
            for i in range(len(self.Sx)):
                h = h + self.all_Sx[2*i]*self.all_Sx[2*i+1]
                h = h + self.all_Sy[2*i]*self.all_Sy[2*i+1]
                h = h + self.all_Sz[2*i]*self.all_Sz[2*i+1]
            self.MBChain.mode = self.mode # overwrite the mode
#            wf = self.MBChain.random_mps()
            if self.T>1e-5: # non-zero temperature
                self.MBChain.set_hamiltonian(h) # singlet Hamiltonian
                wf = self.MBChain.get_gs() # get the fully entangled WF
                wf0 = anneal(self.MBChain,self.hamiltonian,wf,self.T)
            elif self.T<1e-5:
                self.MBChain.set_hamiltonian(self.hamiltonian)
                wf0 = self.MBChain.get_gs() # get the fully entangled WF
            wf0 = wf0.normalize()
            self.computed_gs = True
            self.wf0 = wf0 # store
            self.MBChain.wf0 = wf0 # overwrite in the full object
            self.MBChain.hamiltonian = self.hamiltonian # overwrite
            return wf0
    def set_hamiltonian(self,h):
        self.hamiltonian = (h + h.get_dagger())/2.
        self.computed_gs = False

                

                
def anneal(sc,h,wf,T,dbeta=0.1):
    """Anneal a certain wavefunction using an exponential"""
    beta = 1./T # beta part
    n = int(beta/dbeta) # number of steps
    h0 = h*(beta/n) # this temperature step
    wf = wf.normalize() # normalize
    for i in range(n): # do as many steps
        print("Annealing, energy",wf.dot(h*wf).real)
        wf1 = (1-h0)*wf # apply the operator
        wf1 = wf1.normalize() # normalize
        if np.abs(1.-np.abs(wf1.dot(wf)))<1e-7: return wf
        wf = wf1 # redefine
    return wf







