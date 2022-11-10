import numpy as np
from ..algebra import algebra


def trace(self,A,**kwargs):
    """Return the full trace of an operator"""
    return self.toMPO(A,**kwargs).trace() # compute the trace


def inverse_trace(self,A,mode="ED",**kwargs):
    """Return trace of the inverse of an operator"""
    if mode=="ED": 
        return inverse_trace_ED(self,A)
    if mode=="DMRG":
#        return stochastic_inverse_trace(self,A,mode="ED",**kwargs)
        return stochastic_inverse_trace(self,A,mode="DMRG",**kwargs)



def inverse_trace_ED(self,A,**kwargs):
    """Return trace of the inverse of an operator"""
    M = self.toMPO(A,mode="ED").SO # extract the matrix
    L = M.shape[0] # dimension
    return algebra.trace(algebra.inv(M))/L # return the (average) trace


def stochastic_inverse_trace(self,A,n=10,delta=1e-2,mode="DMRG",**kwargs):
    """Compute the inverse trace using the stochastic method"""
    from ..randommps import random_states
    wfs = random_states(self,mode=mode,n=n,orthogonal=True) # states
    n = len(wfs) # overwrite
    wBw = [wf.dot(self.applyinverse(A,wf,delta=1e-4)) for wf in wfs] # inverses
    from scipy.stats import tstd
    print("Fluctuation",tstd(wBw)/np.sqrt(n))
    return np.sum(wBw)/len(wBw) # return the sum



