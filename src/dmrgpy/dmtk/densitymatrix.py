import numpy as np

from ..fermionchain import Fermionic_Chain

def dm(self,**kwargs):
    """Return the density matrix of a wavefunction"""
    if isinstance(self.MBO,Fermionic_Chain): # fermionic object
        return dm_fermionic(self,**kwargs)
    raise NotImplementedError("get_dm is only implemented for a "
            "Fermionic_Chain, got "+type(self.MBO).__name__)


def dm_fermionic(self,inds=[0]):
    """Compute the density matrix of a system"""
    inds = list(inds) # sites of the density matrix
    n = len(inds) # number of sites
    out = np.zeros((n,n),dtype=np.complex128) # complex output
    wfs = [self.MBO.C[k]*self for k in inds] # C_k |psi>, one per site
    for i in range(n):
        for j in range(n):
            out[i,j] = wfs[j].dot(wfs[i]) # <psi|Cdag_{inds[j]} C_{inds[i]}|psi>
    return out


