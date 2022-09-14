import numpy as np
import scipy.linalg as lg

def krylov_matrix_representation(H,wfs):
    """Given a Krylov subspace, return the matrix representation"""
    accelerate=False
    nw = len(wfs) # number of wavefunction
    if nw==0: raise # something wrong
    mh = np.zeros((nw,nw),dtype=np.complex) # output matrix
    if accelerate:
        for i in range(nw):
            mh[i,i] = wfs[i].aMb(H,wfs[i]) # compute representation
        for i in range(nw-1):
            mh[i+1,i] = wfs[i+1].aMb(H,wfs[i]) # compute representation
            mh[i,i+1] = wfs[i].aMb(H,wfs[i+1]) # compute representation
    else:
        for i in range(nw):
          for j in range(nw):
              mh[i,j] = wfs[j].aMb(H,wfs[i]) # compute representation
    return mh



def recompute_energies(H,vs,wfs):
    """Given certain eigenvectors, recompute the energies"""
    eout = [] # empty list
    for v0 in vs: # loop over WF
        wf = 0
        for i in range(len(wfs)):
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wf = wf.normalize()
        eout.append(wf.aMb(H,wf)) # compute expectation value
    return np.array(eout) # return energies




def gram_smith_single(w,ws):
    """Gram smith orthogonalization for a single wavefunction"""
    if len(ws)==0: return w
    out = []
    n = len(ws)
    w = w.normalize()
    for wj in ws: # loop over stored wavefunctions
        w = w - wj.dot(w)*wj # remove the overlap with each WF
    return w.normalize()


def gram_smith(ws):
    """Gram smith orthogonalization"""
    out = []
    n = len(ws)
    for i in range(n):
        w = ws[i].copy() # copy wavefunction
        w = gram_smith_single(w,out) # orthogonalize
        out.append(w) # store
    return out


def diagonalize(mh):
    if np.max(np.abs(mh-np.conjugate(mh.T)))<1e-6:
        return lg.eigh(mh)
    else: # non Hermitian
        return lg.eig(np.conjugate(mh).T)



def rediagonalize(H,wfs):
    """Given certain eigenfunctions, rediagonalize a Hamiltonian"""
    n = len(wfs) # number of wavefunctions
    mh = np.zeros((n,n),dtype=np.complex) # empty Hamiltonian
    for i in range(n):
      for j in range(n):
          mh[i,j] = wfs[j].aMb(H,wfs[i]) # compute representation
    (es,vs) = diagonalize(mh) # diagonalize
    wfout = [] # empty list
    for j in range(n):
        v0 = vs.T[j] # get the wavefunction
        wf = 0
        for i in range(n): # loop over components
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wfout.append(wf.copy()) # store wavefunction
    return wfout


def most_mixed_wf(H,wfs,info=False):
    """Return the most mixed wavefunction"""
    if len(wfs)==1: return wfs[0].copy() # return wavefunction
    ef,wfk = krylov_eigenstates(H,wfs) # compute eigenstates
#    wfk = wfs
#    if info:
#        iden = krylov_matrix_representation(1.,wfk)
#        print("Krylov orthogonality") # get the representation
#        print(np.round(iden,1)) # get the representation
    ef = np.array([wfi.aMb(H,wfi) for wfi in wfk]) # compute energies
    ef2 = np.array([wfi.aMb(H,H*wfi) for wfi in wfk]) # compute energies square
    error = np.sqrt(np.abs(ef2-ef**2)) # compute the error
    ind = np.where(error==np.max(error))[0][0] # index of the maximum error
    if info: print("Index of the most mixed",ind)
    if info: print("Energy of the most mixed",np.round(ef[ind],2))
    if info: print("Eigenenergies",np.round(ef,2))
    if info: print("Fluctuations",np.round(error,2))
    return wfk[ind].copy() # return the most mixed WF


def krylov_eigenstates(H,wfs):
    """Return the eigenstates and eigenenergies of an operator"""
    mh = krylov_matrix_representation(H,wfs) # get the representation
    (es,vs) = diagonalize(mh) # diagonalize
    wfout = unitary_transformation(vs.T,wfs) # output eigenfunctions
    return es,wfout # return eigenstates



def unitary_transformation(vs,wfs):
    """Perform a unitary transformation"""
    wfout = [] # storage
    for v0 in vs: # loop over WF
        wf = 0
        for i in range(len(wfs)):
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wf = wf.normalize()
        wfout.append(wf.copy()) # store wavefunction
    return wfout # return transformed wavefunctions

