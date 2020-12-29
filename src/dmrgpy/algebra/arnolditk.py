import numpy as np
import scipy.linalg as lg

# routines to perform iterative diagonalization

def mpsarnoldi(self,H,wf=None,e=0.0,delta=1e-1,
        n=6,maxde=1e-4,mode="SI",maxit=30,**kwargs):
    """Compute an eigenvector using the Arnoldi algorithm"""
    wfs = []
    if wf is None:
        wf = self.random_mps() # random initial guess
        wf = wf.normalize()
    else: wfs.append(wf)
    M = H - (e+delta*1j) # shift
    for i in range(n-len(wfs)): # loop over Krylov vectors
        wf = self.applyinverse(M,wf,**kwargs) # apply the inverse
        wf = gram_smith_single(wf,wfs) # orthogonalize
        if wf is None: break # stop the loop
        wfs.append(wf.copy()) # store
    nw = len(wfs) # number of wavefunction
    mh = np.zeros((nw,nw),dtype=np.complex) # output matrix
    for i in range(nw):
      for j in range(nw):
          mh[i,j] = wfs[j].dot(H*wfs[i]) # compute representation
    (es,vs) = lg.eigh(mh)
    es = es-e # difference
    v0 = [y for (x,y) in sorted(zip(np.abs(es),vs.T))][0]
    wf = 0
    for i in range(nw): wf = wf + np.conjugate(v0[i])*wfs[i] # add
    wf = wf.normalize()
    ef = wf.dot(H*wf).real
    ef2 = wf.dot(H*(H*wf)).real
    error = ef2-ef**2
    print("Error",error,"Energy",ef)
    if error<maxde: return wf
    if maxit<0: return wf
    if nw==1: return wf
    else: return mpsarnoldi(self,H,wf=wf,e=e,mode=mode,
            delta=delta,n=n,maxit=maxit-1,**kwargs)



def gram_smith_single(w,ws):
    """Gram smith orthogonalization for a single wavefunction"""
    out = []
    n = len(ws)
    w = w.normalize()
    for wj in ws: # loop over stored wavefunctions
        w = w - wj.dot(w)*wj # remove the overlap with each WF
    return w.normalize()

