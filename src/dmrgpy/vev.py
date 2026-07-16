# compute vacuum expectation values using multioperators

from . import multioperator

import numpy as np

def power_vev(self,wf=None,n=4,X=None,**kwargs):
    """Compute the moments of an operator"""
    if wf is None: wf = self.get_gs(**kwargs) # ground state
    wfs = [wf.copy()] # wavefunctions
    wfi = wf.copy()
    for i in range(n): # get X**n |WF>
        wfi = X*wfi
        wfs.append(wfi.copy()) # apply the operator
    out = [0.0j for i in range(2*(n//2))] # output
    for i in range(n//2): # loop over powers wanted
        out[2*i] = wfs[i].dot(wfs[i])
        out[2*i+1] = wfs[i].dot(wfs[i+1])
    return np.array(out)



def multi_vev(self,MO,wf=None,npow=1,**kwargs):
    """
    Compute a VEV using multioperators
    """
    MO = multioperator.obj2MO(MO,name="vev_multioperator")
    if MO.name!="vev_multioperator": raise
    if npow==0: return 1.0
    if wf is None: wf = self.get_gs() # get the ground state
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    c = self._session.vev(MO.to_terms(),wf.cpp_handle,npow=int(npow))
    return c


def vev(*args,**kwargs):
    return multi_vev(*args,**kwargs)


def excited_vev(*args,**kwargs):
    return multi_vev(*args,excited=True,**kwargs)












