import numpy as np

from ..fermionchain import Fermionic_Chain

def dm(self,**kwargs):
    """Return the density matrix of a wavefunction"""
    if type(self.MBO)==Fermionic_Chain: # fermionic object
        return dm_fermionic(self,**kwargs)
    else: raise # not implemented


def dm_fermionic(self,inds=[0]):
    """Compute the density matrix of a system"""
    n = len(inds) # number of sites
    out = np.zeros((n,n),dtype=np.complex_) # complex output
    for i in range(n):
        for j in range(n):
            wfi = self.MBO.C[i]*self # first wavefunction
            wfj = self.MBO.C[j]*self # first wavefunction
            out[i,j] = wfj.dot(wfi) # store
    return out


