import numpy as np
import scipy.linalg as lg
from .entropytk.correlationentropy import get_correlation_matrix

 

def get_correlation_eigenvalues(self,**kwargs):
    """Return the correlation eigenvalues"""
    cm = get_correlation_matrix(self,**kwargs)
    return lg.eigvalsh(cm) # return eigenvalues


def get_correlation_entropy(self,**kwargs):
    """Return the correlation entropy"""
    vs = get_correlation_eigenvalues(self,**kwargs)
    out = 0.0
    for v in vs:
        if v>1e-6: out = out + v*np.log(v)
    return -out


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
    print(np.sum(np.abs(cm-np.conjugate(cm.T))))
    vs = lg.eigvalsh(cm) # return eigenvalues
    print(np.round(vs,3))
    out = 0.0
    for v in vs:
        if v>1e-6: out = out + v*np.log(v)
    return -out
