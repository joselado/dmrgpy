import numpy as np
import scipy.linalg as lg



def get_correlation_matrix(self,T=0.,**kwargs):
    """Compute the correlation matrix of a finite temperature state"""
    if T==0.: return get_correlation_matrix_zeroT(self,**kwargs)
    else: return get_correlation_matrix_finiteT(self,T=T,**kwargs)



def get_correlation_matrix_finiteT(self,T=1.,**kwargs):
    """Wrapper for finite temperature"""
    ### This is currently only implemented with ED
    n = len(self.C) # number of sites
    if n>14: raise # not implemented
    (es,wfs) = self.get_excited_states(mode="ED",n=2**n)
    dms = [get_correlation_matrix_zeroT(self,wf=wf) for wf in wfs]
    es = es - np.min(es) # minus ground state energy
    dm = 0j # initialize
    Z = 0. # initialize
    for i in range(len(es)): # loop over energies
        dm = dm + dms[i]*np.exp(es[i]/T)
        Z = Z + np.exp(es[i]/T)
    return dm/Z # return matrix




def get_correlation_matrix_zeroT(self,operators=None,
                              basis ="electron", 
                              dmmode="fast",
                              wf=None,**kwargs):
    """Compute the correlation matrix of a ground state"""
    from .. import fermionchain
    if wf is None: wf = self.get_gs(**kwargs) # compute ground state
    if operators is None: # no operators provided
        if fermionchain.isfermion(self):
            if basis=="Nambu": 
              operators = [o for o in self.C] 
              operators += [o for o in self.Cdag] 
            else: # just normal basis
              operators = self.C # fermionic operators
        else: 
            print("Unrecognized type",type(self))
            raise
    # create the matrix
    if dmmode=="simple":
        return correlation_matrix_clean(operators,wf,self)
    elif dmmode=="fast":
        return correlation_matrix_fast(operators,wf)
    elif: raise # not implemented
    n = len(operators)
    cm = np.zeros((n,n),dtype=np.complex)
    for i in range(n):
        A = self.get_dagger(operators[i])
        for j in range(i,n):
            B = operators[j]
            wf1 = B*wf # first operator
            wf2 = A*wf1 # second operator
            out = wf.dot(wf2) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix



def correlation_matrix_clean(operators,wf,self):
    """Compute the correlation matrix of a wavefunction with the fastest 
    algorithm"""
    # create the matrix
    n = len(operators)
    cm = np.zeros((n,n),dtype=np.complex)
    for i in range(n):
        A = self.get_dagger(operators[i])
        for j in range(i,n):
            B = operators[j]
            wf1 = B*wf # first operator
            wf2 = A*wf1 # second operator
            out = wf.dot(wf2) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix

def correlation_matrix_fast(operators,wf):
    """Compute the correlation matrix of a wavefunction with the fastest 
    algorithm"""
    # create the matrix
    n = len(operators)
    wfs = [o*wf for o in operators] # compute all the wavefunctions
    cm = np.zeros((n,n),dtype=np.complex)
    for i in range(n):
        for j in range(i,n):
            out = wfs[i].dot(wfs[j]) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix



## High order correlation entropy


def get_highorder_correlation_matrix(self,operators=None,wf=None,**kwargs):
    """Compute the correlation matrix of a ground state"""
    from .. import fermionchain
    if wf is None: wf = self.get_gs(**kwargs) # compute ground state
    if operators is None:
        if type(self)==fermionchain.Fermionic_Chain:
            operators = self.C # fermionic operators
        elif type(self)==fermionchain.Spinful_Fermionic_Chain:
            operators = self.C # fermionic operators
        else: raise
    # create the matrix
    n = len(operators)
    cm = np.zeros((n,n,n,n),dtype=np.complex)
    for i in range(n):
        A = operators[i].get_dagger()
        for j in range(i+1,n):
            B = operators[j].get_dagger()
            for k in range(n):
                C = operators[k]
                for l in range(k+1,n):
                    D = operators[l]
                    Op = A*B*D*C
                    out = wf.dot(Op*wf) # overlap
                    cm[i,j,k,l] = out
#                    if np.abs(out)>1e-4:
#                       print(i,j,k,l,np.round(out,2))
#    print("Trace",np.sum([cm[i,i,i,i] for i in range(n)]))
    cm = four2two(cm) 
    if np.sum(np.abs(cm-np.conjugate(cm.T)))>1e-4: raise
    return cm # return matrix









def four2two(m):
    n = m.shape[0] # dimension
    out = np.zeros((n*n,n*n),dtype=np.complex) # output
    ii = 0
    for i in range(n-1):
      for j in range(i+1,n):
        jj = 0
        for k in range(n-1):
          for l in range(k+1,n):
              out[ii,jj] = m[i,j,k,l]
              jj += 1
        ii += 1
    return out[0:ii,0:ii]












