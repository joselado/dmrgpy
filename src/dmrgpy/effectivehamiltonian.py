import numpy as np
import scipy.linalg as lg

def get_effective_hamiltonian(self,**kwargs):
    """Return an effective Hamiltonian"""
    n = self.ns # number of sites 
    pairs = []
    for i in range(n): # loop
        for j in range(n): # loop
            pairs.append([i,j]) # store
    cs = self.get_correlator(pairs=pairs,apply_hamiltonian=True,**kwargs)
    norm = self.get_correlator(pairs=pairs,apply_hamiltonian=False,**kwargs)
    h = np.zeros((n,n),dtype=np.complex) # create matrix with coefficients
    b = np.zeros((n,n),dtype=np.complex) # create matrix with overlaps
    for k in range(len(cs)): # loop
        i,j = pairs[k] # get index
        h[i,j] = cs[k] # store
        b[i,j] = norm[k] # store
    return h,b # return effective hamiltonian
