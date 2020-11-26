import numpy as np
import scipy.linalg as lg
from . import fermionchain

def get_correlation_matrix(self,operators=None,wf=None,**kwargs):
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
        for j in range(n):
            B = operators[j]
            C = A*B
            print(i,j)
            cm[i,j] = self.vev(C,**kwargs)
    return cm # return matrix
  
