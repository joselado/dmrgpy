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
    if getattr(self,"use_cpp_extension",False) and self._session is not None:
        return multi_vev_cpp_ext(self,MO,wf,npow)
    wf.write(name="wf_vev.mps") # write wavefunction
    taskd = MO.get_dict() # get the dictionary
    self.task["vev"] = "true" # do a VEV
    self.task["wf_vev"] = "wf_vev.mps" # WF to use VEV
    self.task["pow_vev"] = int(npow) # power
    self.write_task() # write the tasks in a file
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.execute(lambda: MO.write()) # write multioperator
    self.run() # perform the calculation
    m = self.execute(lambda: np.genfromtxt("VEV.OUT"))
    return m[0]+1j*m[1] # return result


def multi_vev_cpp_ext(self,MO,wf,npow):
    """
    Compute a VEV via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::vev), mirroring multi_vev()'s
    DMRG path exactly but with no file I/O.
    """
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    c = self._session.vev(MO.to_terms(),wf.cpp_handle,npow=int(npow))
    return c


def vev(*args,**kwargs):
    return multi_vev(*args,**kwargs)


def excited_vev(*args,**kwargs):
    return multi_vev(*args,excited=True,**kwargs)












