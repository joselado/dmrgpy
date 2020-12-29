import numpy as np
import scipy.linalg as lg

# routines to perform iterative diagonalization

def mpsarnoldi(self,H,wf=None,e=0.0,delta=1e-1,
        mode="SI",**kwargs):
    """Compute an eigenvector using the Arnoldi algorithm"""
    if mode=="SI": # target a specific energy
        M = H - (e+delta*1j) # shift
        Op = lambda x: self.applyinverse(M,x) # operator to apply
        def fe(es): # function return the right WF
            es = np.abs(es-e) # minimum energy
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="GS": # target the ground state
        shift = self.ns/2. # number of sites
        M = H-shift # start with the Hamiltonian
        Op = lambda x: M*x # operator to apply
        def fe(es): # function return the right WF
            return np.where(es==np.min(es))[0][0] # return the index
    else: raise
    return mpsarnoldi_iteration(self,Op,H,wf,fe,**kwargs)

def mpsarnoldi_iteration(self,Op,H,wf,fe,maxit=30,maxde=1e-4,
        verbose=1,n=6):
    """Single iteration of the restarted arnoldi algorithm"""
    wfs = []
    if wf is None:
        wf = self.random_mps() # random initial guess
        wf = wf.normalize()
    else: wfs.append(wf)
    for i in range(n-len(wfs)): # loop over Krylov vectors
        wf = Op(wf) # apply operator
        wf = gram_smith_single(wf,wfs) # orthogonalize
        if wf is None: break # stop the loop
        wfs.append(wf.copy()) # store
    nw = len(wfs) # number of wavefunction
    mh = np.zeros((nw,nw),dtype=np.complex) # output matrix
    for i in range(nw):
      for j in range(nw):
          mh[i,j] = wfs[j].dot(H*wfs[i]) # compute representation
    (es,vs) = lg.eigh(mh)
    inde = fe(es) # get the index
    v0 = vs.T[inde] # get the wavefunction
    wf = 0
    for i in range(nw): wf = wf + np.conjugate(v0[i])*wfs[i] # add
    wf = wf.normalize()
    ef = wf.dot(H*wf).real
    ef2 = wf.dot(H*(H*wf)).real
    error = ef2-ef**2
    if verbose>0: print("Error",error,"Energy",ef)
    if error<maxde: return wf
    if maxit<0: return wf
    if nw==1: return wf
    else: return mpsarnoldi_iteration(self,Op,H,wf,fe,maxit=maxit-1,
            maxde=maxde,verbose=verbose,n=n)



def gram_smith_single(w,ws):
    """Gram smith orthogonalization for a single wavefunction"""
    out = []
    n = len(ws)
    w = w.normalize()
    for wj in ws: # loop over stored wavefunctions
        w = w - wj.dot(w)*wj # remove the overlap with each WF
    return w.normalize()

