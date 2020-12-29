import numpy as np
import scipy.linalg as lg

# routines to perform iterative diagonalization

def mpsarnoldi(self,H,wf=None,e=0.0,delta=1e-1,
        mode="SI",P=None,**kwargs):
    """Compute an eigenvector using the Arnoldi algorithm"""
    if mode=="SI": # target a specific energy
        M = H - (e+delta*1j) # shift
        Op = lambda x: self.applyinverse(M,x) # operator to apply
        def fe(es): # function return the right WF
            es = np.abs(es-e) # minimum energy
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="GS": # target the ground state
        # for non-Hermitian matrices, this targets the eigenvalues
        # with most negative real part
        shift = self.ns/2. # number of sites
        M = H-shift # start with the Hamiltonian
        if P is None: Op = lambda x: M*x # operator to apply
        else: Op = lambda x: M*(P*x) # operator to apply
        def fe(es): # function return the right WF
            es = es.real # take the real part
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="MRGS": # target the most real ground state
        shift = self.ns/2. # number of sites
        M = H-shift # start with the Hamiltonian
        Op = lambda x: M*x # operator to apply
        def fe(es): # function to return the right WF
            esr = es.real # take the real part
            esr = esr-np.min(esr) # shift the minimum to zero
            esi = np.abs(es.imag) # take the imaginary part
            d = 1./(esr**2 +esi**2/delta + delta) # weight for each eigenvalue
            return np.where(d==np.max(d))[0][0] # return the index
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
    (es,vs) = diagonalize(mh) # diagonalize
    inde = fe(es) # get the index
    v0 = vs.T[inde] # get the wavefunction
    wf = 0
    for i in range(nw): wf = wf + np.conjugate(v0[i])*wfs[i] # add
    wf = wf.normalize()
    ef = wf.dot(H*wf)
    if np.abs(ef.imag)<1e-8: ef= ef.real
    ef2 = wf.dot(H*(H*wf))
    error = np.abs(ef2-ef**2) # compute the error
    if verbose>0: print("Error",error,"Energy",np.round(ef,4))
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



def diagonalize(mh):
    if np.max(np.abs(mh-np.conjugate(mh.T)))<1e-6:
        return lg.eigh(mh)
    else:
        return lg.eig(mh)
