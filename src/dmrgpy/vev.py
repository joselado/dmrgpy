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

    For npow>1 the operator is applied npow-1 times, each application
    truncated to the chain's maxm and cutoff, so <MO^npow> carries that
    truncation: below the state's own bond dimension times the MPO's, the
    number depends on maxm (2026-09-24c audit, finding 7). The energy
    variance, where that error is the whole answer, goes through
    energy_variance() below instead.
    """
    MO = multioperator.obj2MO(MO,name="vev_multioperator")
    if MO.name!="vev_multioperator": raise
    if npow==0: return 1.0
    if wf is None: wf = self.get_gs() # get the ground state
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    c = self._session.vev(MO.to_terms(),wf.cpp_handle,npow=int(npow))
    return c


# A bond-dimension cap for an MPO application that is meant not to truncate
# at all: the result's bond dimension is bounded by the MPO's times the
# state's, far below this, so only the cutoff acts.
_UNCAPPED = 1000000


def energy_variance(self,h,wf=None):
    """<wf|(H-<H>)^2|wf> on a DMRG session backend, as ||(H-<H>)|wf>||^2.

    vev(H,npow=2) applies H to the state truncated to the chain's maxm and
    returns <wf|H|trunc(H|wf>)>, and <H^2>-<H>^2 is then a difference of two
    numbers of order E0^2 whose error is set by that truncation: a state
    solved and measured at the same maxm reported its variance 10 to 51
    times too low, one wider than maxm reported it at order one (2026-09-24c
    audit, finding 7). Subtracting <H> first makes the applied vector the
    small residual itself, so the relative cutoff acts on its own norm, and
    the application runs uncapped, so maxm does not act at all. Building
    H^2 as one MPO instead is exact too but O(L^3) to construct on
    "python" (302 s against 1.8 s at 80 sites)."""
    if wf is None: wf = self.get_gs()
    e = multi_vev(self,h,wf=wf)
    dh = h - e*multioperator.identity()
    old = self.maxm
    self.maxm = _UNCAPPED
    try:
        out = multi_vev(self,dh,wf=wf,npow=2)
    finally:
        self.maxm = old
        self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
        self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    return out


def vev(*args,**kwargs):
    return multi_vev(*args,**kwargs)


def excited_vev(*args,**kwargs):
    return multi_vev(*args,excited=True,**kwargs)












