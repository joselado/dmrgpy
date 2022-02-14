import numpy as np
import scipy.linalg as lg
from .. import fermionchain



def get_correlation_matrix(self,T=0.,**kwargs):
    """Compute the correlation matrix of a finite temperature state"""
    if T==0.: return get_correlation_matrix_zeroT(self,**kwargs)
    else: return get_correlation_matrix_finiteT(self,T=T,**kwargs)



def get_correlation_matrix_finiteT(self,T=1.,**kwargs):
    """Wrapper for finite temperature"""
    ### This is currently only implemented with ED
    raise # not finished yet
    n = len(self.C) # number of sites
    if n>14: raise # not implemented
    (es,wfs) = fc.get_excited_states(mode="ED",n=2**n)
    dms = [get_correlation_matrix_zeroT(self,wf=wf) for wf in wfs]
    es = es - np.min(es) # minus ground state energy
    return dm




def get_correlation_matrix_zeroT(self,operators=None,wf=None,**kwargs):
    """Compute the correlation matrix of a ground state"""
    if wf is None: wf = self.get_gs(**kwargs) # compute ground state
    if operators is None:
        if type(self)==fermionchain.Fermionic_Chain:
            operators = self.C # fermionic operators
        elif type(self)==fermionchain.Spinful_Fermionic_Chain:
            operators = self.C # fermionic operators
        else: raise
    # create the matrix
    n = len(operators)
    cm = np.zeros((n,n),dtype=np.complex)
    for i in range(n):
        A = operators[i].get_dagger()
        for j in range(i,n):
            B = operators[j]
            C = A*B
            print(i,j)
    #        out = self.vev(C,**kwargs)
            out = wf.dot(C*wf) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix

