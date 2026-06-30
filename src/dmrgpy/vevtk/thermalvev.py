import numpy as np


def thermal_vev(MB,Op,submode="EX",**kwargs):
    """Compute the thermal vev"""
    if submode=="EX": # excited states method
        return thermal_vev_ex(MB,Op,**kwargs) # thermal VEV
    else: raise




def thermal_vev_ex(MB,Op,mode="ED",T=0.,**kwargs):
    """Thermal VEV with excited states"""
    (es,wfs) = MB.get_excited_states(mode=mode,**kwargs) # get excited states
    es = es - np.min(es) # with respect to the GS
    if T>0.: # finite temperature
        beta = 1./T # beta
        P = np.exp(-beta*es) # Boltzman probabilities
        P = P/np.sum(P) # normalize probabilities
    else: 
        print("Only finite temperature")
        raise # not implemented
    if T>np.max(es):
        print("Warning, highest excited state comparable to thermal energy, increase n=",T,np.max(es))
    vals = np.array([wf.aMb(Op,wf) for wf in wfs]) # expectation values
    return np.sum(P*vals) # return the expectation value
