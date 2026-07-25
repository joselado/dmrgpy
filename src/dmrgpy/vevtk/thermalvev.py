import numpy as np
from ..algebra.algebra import maxsize


def thermal_vev(MB,Op,submode="EX",**kwargs):
    """Compute the thermal vev"""
    if submode=="EX": # excited states method
        return thermal_vev_ex(MB,Op,**kwargs) # thermal VEV
    else: raise




def thermal_vev_ex(MB,Op,mode="ED",T=0.,**kwargs):
    """Thermal VEV with excited states"""
    dim = MB.get_hamiltonian().shape[0] # full Hilbert space dimension
    n_requested = kwargs.get("n",None)
    if n_requested is None: # caller did not ask for a specific truncation
        if dim<=maxsize: # small enough to diagonalize exactly anyway
            kwargs["n"] = dim # use the full spectrum, no truncation
            n_used = dim
        else: n_used = 10 # matches algebra.lowest_states' own default
    else: n_used = n_requested
    truncated = n_used<dim # did we actually leave out any eigenstates?
    (es,wfs) = MB.get_excited_states(mode=mode,**kwargs) # get excited states
    es = es - np.min(es) # with respect to the GS
    if T>0.: # finite temperature
        beta = 1./T # beta
        P = np.exp(-beta*es) # Boltzman probabilities
        P = P/np.sum(P) # normalize probabilities
    else:
        print("Only finite temperature")
        raise # not implemented
    if truncated and P[-1]>1e-6*np.max(P):
        # the highest computed excited state still carries non-negligible
        # Boltzmann weight, so the states left out by the truncation are
        # not safe to ignore: summing only over "n_used" would silently
        # return a wrong thermal average instead of the true Tr(O exp(-beta H))/Z
        raise RuntimeError(
            "thermal_vev_ex: the "+str(n_used)+" excited states computed do not "
            "capture the full Boltzmann weight at T="+str(T)+" (highest included "
            "state has relative weight "+str(P[-1]/np.max(P))+" of the largest). "
            "Pass a larger n= to edchain.vev()/thermal_vev_ex() to include more "
            "excited states.")
    vals = np.array([wf.aMb(Op,wf) for wf in wfs]) # expectation values
    return np.sum(P*vals) # return the expectation value
