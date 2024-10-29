import numpy as np
import scipy.linalg as lg

def krylov_matrix_representation(H,wfs):
    """Given a Krylov subspace, return the matrix representation"""
    accelerate=False
    nw = len(wfs) # number of wavefunction
    if nw==0: raise # something wrong
    mh = np.zeros((nw,nw),dtype=np.complex128) # output matrix
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
        es,ws = lg.eig(np.conjugate(mh).T)
        return np.conjugate(es),ws



def rediagonalize(H,wfs):
    """Given certain eigenfunctions, rediagonalize a Hamiltonian"""
    n = len(wfs) # number of wavefunctions
    mh = np.zeros((n,n),dtype=np.complex128) # empty Hamiltonian
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
    weight = error + 1e-6 # normalize
    weight = weight/np.max(weight) # weight
    wfo = wfs[0]*0. # initialize
    for i in range(len(wfs)):
        phi = np.exp(1j*np.random.random()*np.pi*2) # random phase
        wfo = wfo + phi*weight[i]*wfs[i] # output
    wfo = wfo.normalize() # normalize wavefunction
    if info: print("Eigenenergies",np.round(ef,2))
    if info: print("Fluctuations",np.round(error,2))
    return wfo


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
        wf = 0*wfs[0] # initialize
        for i in range(len(wfs)):
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wf = wf.normalize()
        wfout.append(wf.copy()) # store wavefunction
    return wfout # return transformed wavefunctions




def generalized_diagonalize(H,wfs):
    """Return the eigenvalues and eigenvectors of a generalized
    eigenvalue problem"""
    mh = krylov_matrix_representation(H,wfs) # matrix representation
    b = krylov_matrix_representation(1.,wfs) # matrix representation
#    print(np.round(b,2))
    if np.max(np.abs(mh-np.conjugate(mh.T)))<1e-6:
        return lg.eigh(mh,b=b)
    else: # non Hermitian
        es,ws = lg.eig(np.conjugate(mh).T,b=b)
        return np.conjugate(es),ws 



def select_states(es,vs,fe,ne=1):
    """Select states according to a criteria"""
    elist = [e for e in es] # list with the energies
    vlist = [v for v in vs] # list with the states
    vstore = [] # empty list
    estore = [] # empty list
    einds = [] # indexes
    for i in range(ne): # loop over desired energies
        ie = fe(np.array(elist)) # get the desired index
        estore.append(elist[ie]) # store this energy
        vstore.append(vlist[ie]) # store this WF
        del elist[ie] # ignore in the next iteration
        del vlist[ie] # ignore in the next iteration
    return estore,vstore


def selectwf(es,vs,wfs,fe,ne=1):
    """Select the wavefunctions that should be returned"""
    estore,vstore = select_states(es,vs,fe,ne=ne)
    wfout = [] # output wavefunctions
    for v0 in vstore: # loop over WF
        wf = 0*wfs[0]
        for i in range(len(wfs)):
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wf = wf.normalize() # normalize
        wfout.append(wf.copy()) # store wavefunction
    eout = np.array(estore) # convert to array
    return eout,wfout




def krylov2states(H,wfs,criteria=None,**kwargs):
    """Given a Hamiltonian and a set of states, return the eigenstates
    that fufill a certain criteria"""
    if criteria is None: return wfs # return all
    mh = krylov_matrix_representation(H,wfs) # get the representation
    (es,vs) = diagonalize(mh) # diagonalize
    # select the wavefunctions
    ef,wf = selectwf(es,vs.T,wfs,criteria,**kwargs,ne=len(wfs)) 
    return wf




