# this library contains routines to perform thermal calculations
import numpy as np

from .spinchain import Spin_Chain
from . import multioperator

class Thermal_Spin_Chain():
    def __init__(self,sites,T=0.1,mode="DMRG",**kwargs):
        if T<0.: # a negative temperature used to be silently treated as T=0
            raise ValueError("Thermal_Spin_Chain: T must be >= 0, got "
                             +repr(T))
        # mode is the wrapper's own, and get_gs() writes it onto MBChain;
        # that used to be a hardcoded "DMRG", which overwrote a forwarded
        # mode="ED" so that Thermal_Spin_Chain(...,mode="ED") ran DMRG
        if mode is not None: # None leaves the chain on its default, as on any chain
            from .mode import _check_mode
            _check_mode(mode,"Thermal_Spin_Chain mode (mode=)")
        sitesT = []
        for s in sites: 
            sitesT += [s]
            sitesT += [s]
        n = len(sites) # number of sites
        # **kwargs was accepted and then dropped on the floor here, so
        # Thermal_Spin_Chain(...,itensor_version="python") silently built
        # the default (compiled) backend instead
        self.MBChain = Spin_Chain(sitesT,mode=mode,**kwargs) # get the chain
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
        self.mode = mode
    def get_gs(self):
        """Compute the ground state"""
        if self.computed_gs: 
            return self.wf0
        else:
            def terms():
                for i in range(len(self.Sx)):
                    yield self.all_Sx[2*i]*self.all_Sx[2*i+1]
                    yield self.all_Sy[2*i]*self.all_Sy[2*i+1]
                    yield self.all_Sz[2*i]*self.all_Sz[2*i+1]
            h = multioperator.msum(terms()) # initialize
            self.MBChain.mode = self.mode # overwrite the mode
#            wf = self.MBChain.random_mps()
            if self.T>1e-5: # non-zero temperature
                self.MBChain.set_hamiltonian(h) # singlet Hamiltonian
                wf = self.MBChain.get_gs() # get the fully entangled WF
                wf0 = anneal(self.MBChain,self.hamiltonian,wf,self.T)
                wf0 = wf0.normalize()
                # Through the public setters, so MBChain holds the
                # physical Hamiltonian and the annealed state as a state
                # set by hand: the next read takes it unswept, and every
                # correlator measures it. Assigning MBChain.wf0 and
                # MBChain.hamiltonian directly left the singlet
                # Hamiltonian on the session and in the send-cache (and,
                # on mode="ED", the singlet ED object), so gs_energy()
                # returned the singlet energy (-2.25 against <wf|H|wf> =
                # -0.4164 on a 3-site chain at T=1), a correlator re-solved
                # over the annealed state and vev() on ED read the singlet
                # state (<Sz0 Sz1> -0.1667 and 0.0 against -0.0694)
                self.MBChain.set_hamiltonian(self.hamiltonian)
                self.MBChain.set_gs(wf0)
            else: # T<=1e-5: plain ground state. `else`, not another elif:
                # the two branches used to leave T exactly 1e-5 (and, before
                # the check in __init__, any negative T) falling through
                # both, so wf0 was never assigned and the next line raised
                # UnboundLocalError. MBChain's own solved state is this
                # one already, so nothing is written back to it
                self.MBChain.set_hamiltonian(self.hamiltonian)
                wf0 = self.MBChain.get_gs() # get the fully entangled WF
                wf0 = wf0.normalize()
            self.computed_gs = True
            self.wf0 = wf0 # store
            return wf0
    def set_hamiltonian(self,h):
        self.hamiltonian = (h + h.get_dagger())/2.
        self.computed_gs = False

                

                
def anneal(sc,h,wf,T,dbeta=0.1):
    """Anneal a certain wavefunction using an exponential"""
    # Purification convention: |Psi(beta)> = e^{-beta*H/2}|Psi(0)>, so that
    # tracing out the ancilla out of |Psi(beta)><Psi(beta)| gives
    # rho ~ e^{-beta*H}, i.e. the physical thermal state at the requested
    # T = 1/beta. Applying the full beta here (as this used to) doses the
    # wavefunction with e^{-beta*H} instead of e^{-beta*H/2}, which after
    # tracing out the ancilla is the thermal state of T/2, not T -- confirmed
    # numerically against exact ED thermal averages (see docs/user_guide).
    beta_half = 1./(2.*T) # half-beta part
    n = int(beta_half/dbeta) # number of steps
    h0 = h*(beta_half/n) # this temperature step
    wf = wf.normalize() # normalize
    for i in range(n): # do as many steps
        print("Annealing, energy",wf.dot(h*wf).real)
        wf1 = (1-h0)*wf # apply the operator
        wf1 = wf1.normalize() # normalize
        if np.abs(1.-np.abs(wf1.dot(wf)))<1e-7: return wf
        wf = wf1 # redefine
    return wf







