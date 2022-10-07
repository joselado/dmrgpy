import numpy as np
import scipy.linalg as lg
from .entropytk.correlationentropy import get_correlation_matrix


def get_correlation_entropy_density(self,**kwargs):
    """Compute the density of spinless correlation entropy"""
    from . import fermionchain
    cm = get_correlation_matrix(self,**kwargs)
    ce,ws = lg.eigh(cm) # return eigenvalues
    ws = ws.T # transpose
    ce[ce<1e-6] = 1. # set to 1, so that they have zero entropy
    ss = -ce*np.log(ce) # entropies 
    out = 0.*ws[0].real # initialize
    for i in range(len(ss)): out = out + ss[i]*np.abs(ws[i])**2
    if type(self)==fermionchain.Fermionic_Chain: pass
    elif type(self)==fermionchain.Spinful_Fermionic_Chain:
      out = [out[2*i]+out[2*i+1] for i in range(len(out)//2)]
    return np.array(out) # return the density

 

def get_correlation_eigenvalues(self,**kwargs):
    """Return the correlation eigenvalues"""
    cm = get_correlation_matrix(self,**kwargs)
    ce = lg.eigvalsh(cm) # return eigenvalues
    # a sanity check
    if np.max(ce)>1.01: raise # they should be at most 1.
    if np.min(ce)<-0.01: raise # they should be zero the smallest
    return ce


def get_correlation_entropy(self,**kwargs):
    """Return the correlation entropy"""
    vs = get_correlation_eigenvalues(self,**kwargs)
    out = 0.0
    for v in vs:
        if v>1e-6: out = out + v*np.log(v)
    return -out



def get_correlation_entropy_from_wf(self,**kwargs):
    """Return the correlation entropy of a wavefunction"""
    if self.MBO is not None: # if it has an MBO
        return get_correlation_entropy(self.MBO,wf=self,**kwargs)
    else: raise





def get_correlated_orbitals(self,ordered=True,**kwargs):
    """Return the most correlated orbitals"""
    cm = get_correlation_matrix(self,**kwargs)
    es,vs = lg.eigh(cm) # diagonalize
    vs = np.abs(vs.T)**2 # transpose
    esa = np.abs(es-.5)
    if ordered: vs = [y for (x,y) in sorted(zip(esa,vs),key=lambda x: x[0])]
    return vs


def get_correlated_density(self,n=4,**kwargs):
    """Return the most correlated orbitals"""
    vs = get_correlated_orbitals(self,**kwargs)
    return np.sum(vs[0:n],axis=0)



def get_second_order_correlation_entropy(self,**kwargs):
    """Return the second_order correlation entropy"""
    from .entropytk.correlationentropy import get_highorder_correlation_matrix
    cm = get_highorder_correlation_matrix(self,**kwargs)
#    print(np.sum(np.abs(cm-np.conjugate(cm.T))))
    vs = lg.eigvalsh(cm) # return eigenvalues
    out = 0.0
    for v in vs:
        if v>1e-6: out = out + v*np.log(v)
    out = -out
#    print(np.round(vs,3),np.round(out,2))
    return out
