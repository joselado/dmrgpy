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




def normalize_image(wf):
    """wf/||wf|| for wf = A|u>, the image of a unit vector under an
    operator, or None when that image is exactly zero.

    Such a state's norm carries A's units, so wf.normalize()'s absolute
    1e-8 floor is a floor on the operator's units, not a test for a
    vanishing state: with A written in small units (1e-9*H), a perfectly
    good H|u> came back as None and the next product raised TypeError
    (2026-09-25b audit, finding 22, through powermethod.estimate_radius).
    Dividing by the norm itself is right here and only here: an image of
    a unit vector is not a difference of two states, so there is no
    cancellation to normalize into noise. Anything that subtracts goes
    through normalize(tol=...) with a reference of its own. The norm is
    computed as normalize() computes it, so where normalize() would have
    divided, this divides by the same number."""
    nrm = np.sqrt(max(np.real(wf.dot(wf)),0.))
    if nrm==0.: return None
    return wf*(1./nrm)


def is_hermitian_matrix(mh,atol=1e-6):
    """Whether a small Krylov-space matrix is Hermitian:
    max|mh-mh^dagger| < atol*min(1,max|mh|), the zero matrix counting as
    Hermitian. It was an absolute max|mh-mh^dagger| < 1e-6, which every
    matrix of an operator written in small units passes, sending a
    non-Hermitian one to eigh (2026-09-25b audit, lead
    scale-arnolditk-absolute-stop, by reading). Below entries of order one
    the test is now relative to the matrix; at and above it, it is the old
    absolute one bit for bit, so no decision of the IRAM and Arnoldi
    routes on a Hamiltonian in ordinary units moves -- the one-sided rule
    e7b1196 set for v2/v3 (mo_terms.h's unit_scale_up)."""
    scale = np.max(np.abs(mh))
    if scale==0.: return True
    return bool(np.max(np.abs(mh-np.conjugate(mh.T))) < atol*min(1.,scale))


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
        if w is not None:
            out.append(w) # store
    return out


def diagonalize(mh):
    if is_hermitian_matrix(mh): # relative to mh's own size, see there
        return lg.eigh(mh)
    else: # non Hermitian
        es,ws = lg.eig(np.conjugate(mh).T)
        return np.conjugate(es),ws



def rediagonalize(H,wfs):
    """Given certain eigenfunctions, rediagonalize a Hamiltonian, returning
    both the refined eigenvalues and eigenvectors (as MPS linear
    combinations of the input wfs). Returning the eigenvalues here too
    lets callers avoid building the same O(n^2) representation matrix a
    second time just to get them (see excited.py's purify=True path)."""
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
    return es,wfout


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
    if is_hermitian_matrix(mh): # relative to mh's own size, see there
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




