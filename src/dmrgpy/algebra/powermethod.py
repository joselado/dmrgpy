import numpy as np


def power_method_several(self,H,verbose=0,n0=10,ntries=10,shift=0.0,
        orthogonal=None,unfiltered=False,error=1e-6,**kwargs):
    """Simple implementation of the power method"""
    from .arnolditk import random_state
    from .krylov import gram_smith,gram_smith_single
    def f(): # function to perform a single power method iteration
        wf = random_state(self,orthogonal=orthogonal) # take a random MPS
        eold = -1e20
        for i in range(n0): # number of tries
            if shift!=0.: wf = shift*wf + H*wf
            else: wf = H*wf # perform one iteration
            wf = wf.normalize() # normalize wavefunction
            if orthogonal is not None: # reorthogonalize
                wf = gram_smith_single(wf,orthogonal)
            ei = wf.aMb(H,wf)
            if verbose>2: print("Energy in this iteration",ei)
            if np.abs(eold-ei)<error:
                if verbose>2: print("Stopping PM")
                break
            eold = ei
        if verbose>1: print("PM energy",ei)
        return ei,wf # return energy and wavefunction
    es,wfs = [],[]
    for i in range(ntries): # make several tries
        if verbose>0: print("Power method #",i)
        e,wf = f() # use the power method
        es.append(e) # store
        wfs.append(wf.copy()) # store
    if unfiltered: return es,wfs # return without purifying
    wfs = gram_smith(wfs) # orthogonalize
    from .krylov import krylov2states
    wfs = krylov2states(H,wfs,**kwargs) # sort according to the criteria
    es = [wfi.aMb(H,wfi) for wfi in wfs] # energies
    if verbose>0: 
        print("All energies in PM",np.round(es,3))
    return es,wfs



def power_method(self,H,**kwargs):
    """Return a single wavefunction of the power method"""
    es,wfs = power_method_several(self,H,**kwargs) # get all
    es = np.abs(es)
    ind = np.where(es==np.max(es))[0][0] # index of the max energy
    return es[ind],wfs[ind].copy() # return this one



def multi_power_method(self,H,orthogonal=True,nwf=1,**kwargs):
    """Return several vectors using the power method"""
    if orthogonal: # orthogonal method
      return power_method_orthogonal(self,H,nwf=nwf,**kwargs)
    else: 
      return power_method_several(self,H,ntries=nwf,**kwargs)






def power_method_orthogonal(self,H,nwf=1,**kwargs):
    """Return N power method wavefucntions with orthogonalization"""
    es = [] # storage
    wfs = [] # storage
    for i in range(nwf):
        e,wf = power_method(self,H,orthogonal=wfs,ntries=1,**kwargs)
        es.append(e)
        wfs.append(wf.copy())
    return es,wfs



def estimate_radius(self,H,n0=10):
    """Given a Hamiltonian, make an estimate of the radius of the spectra"""
    from .arnolditk import random_state
    wf = random_state(self) # random state
    for i in range(n0): # warmup number of tries
        wf = wf.normalize() # normalize
        wf = H*wf # multiply
    return np.sqrt(np.abs(wf.dot(wf)))




