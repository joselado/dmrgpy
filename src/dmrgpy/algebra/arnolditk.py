import numpy as np
import scipy.linalg as lg

# routines to perform iterative diagonalization

def mpsarnoldi(self,H,wf=None,e=0.0,delta=1e-1,
        mode="GS",P=None,
        recursive_arnoldi=False,
        nwf=1, 
        **kwargs):
    """Compute an eigenvector using the Arnoldi algorithm"""
    if mode=="ShiftInv": # target a specific energy with shift and invert
        M = H - (e+delta*1j) # shift
        Op = lambda x: self.applyinverse(M,x) # operator to apply
        def fe(es): # function return the right WF
            es = np.abs(es-e) # minimum energy
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="GS": # target the ground state
        # for non-Hermitian matrices, this targets the eigenvalues
        # with most negative real part
        shift = self.ns # number of sites
        M = H-shift # start with the Hamiltonian
        if P is None: Op = lambda x: M*x # operator to apply
        else: Op = lambda x: M*(P*x) # operator to apply
        def fe(es): # function return the right WF
            es = es.real # take the real part
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="SM": # smallest magnitude
        M = H+e # start with the Hamiltonian
        Op = lambda x: M*x # operator to apply
        def fe(es): # function return the right WF
            es = np.abs(es-e) # distance to the wanted eigenvalue
            return np.where(es==np.min(es))[0][0] # return the index
    elif mode=="SI": # smallest imaginary part
        M = H # start with the Hamiltonian
        Op = lambda x: M*x # operator to apply
        def fe(es): # function return the right WF
            esi = np.abs(es.imag) # imaginary part
            return np.where(esi==np.min(esi))[0][0] # return the index
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
    if nwf==1: 
        return mpsarnoldi_iteration(self,Op,H,fe,ne=1,**kwargs)
    else: 
        if recursive_arnoldi:
            wfout = [] # empty list
            eout = [] # empty list
            for i in range(nwf): # loop over desired wavefunctions
                ei,wfi = mpsarnoldi_iteration(self,
                        Op,H,fe,ne=1,
                        wfs=[],
                        wfskip=wfout,**kwargs)
                wfout.append(wfi.copy()) # store wavefunction
                eout.append(ei) # store wavefunction
            eout,wfout = sortwf(eout,wfout,fe) # resort the result
            return np.array(eout),wfout # return wavefunctions
        else:
          return mpsarnoldi_iteration(self,Op,H,fe,ne=nwf,**kwargs)


def mpsarnoldi_iteration(self,Op,H,fe,maxit=1,
        maxde=1e-4,
        ne=1, # number of energies to return
        maxdwf = 1e-3, # maximum change in the wavefunction
        wfskip=[], # wavefunctions to skip
        shift = 0.0, # shift for the operator
        verbose=0,
        wfs = [], # initial Krylov vectors
        n=10):
    """Single iteration of the restarted arnoldi algorithm"""
    if verbose>1:
        print("Eigenvalue shift",shift)
    if len(wfs)==0: # no vectors given
        wf = self.random_mps() # random initial guess
        wf = wf.normalize()
        wf0 = None
    else: 
        wf = wfs[0].copy() # take the "best" one
        wf0 = wf.copy() # store initial wavefunction
    for i in range(n-len(wfs)): # loop over Krylov vectors
        wf = Op(wf) + shift*wf # apply operator
        wf = gram_smith_single(wf,wfs+wfskip) # orthogonalize
        if verbose>1: print("Krylov vector #",i)
        if wf is None: 
            if verbose>1: print("Zero vector found, use a random one")
            wf = self.random_mps(orthogonal=wfs+wfskip) # random MPS
        wf = wf.normalize()
        wfs.append(wf.copy()) # store
    nw = len(wfs) # number of wavefunction
    if nw==0: raise # something wrong
    mh = krylov_matrix_representation(H,wfs) # get the representation
    if verbose>2:
        iden = krylov_matrix_representation(1.,wfs)
        print("Krylov orthogonality") # get the representation
        print(np.round(iden,1)) # get the representation
    (es0,vs) = diagonalize(mh) # diagonalize
    es = recompute_energies(H,vs.T,wfs) # recompute the energies
#    if verbose>1:
#        print("Bare energies",es0[0:ne])
#        print("Purified energies",es[0:ne])
    if maxit<2: # if last iteration has been reached
        eout,wfout = selectwf(es,vs.T,wfs,fe,ne=ne) # select the wavefunction
        return eout,wfout # return the results
    #########################################
    # if the algorithm is recursive, continue
    #########################################
    ef,wf = selectwf(es,vs.T,wfs,fe,ne=1) # select the wavefunction
    if np.abs(ef.imag)<1e-8: ef= ef.real
    wf2 = H.get_dagger()*wf
    ef2 = wf2.dot(H*wf)
    error = np.abs(ef2-ef**2) # compute the error
    dwf = 1.0 # difference with respect to the initial wf
    if wf0 is not None:
        dwf = 1.0 - np.abs(wf0.dot(wf)) # difference between WF
    if verbose>0: 
        print("Error in energy",error,"Energy",np.round(ef,4))
        if wf0 is not None: print("Error in WF",dwf)
    # stop according to several criteria
    if error<maxde or dwf<maxdwf: 
        eout,wfout = selectwf(es,vs.T,wfs,fe,ne=ne) # select the wavefunctions
        return eout,wfout
    # if this point is reached, recall
    nk = max([2,n//2]) # number of vector to keep
    eout,wfout = selectwf(es,vs.T,wfs,fe,ne=nk) # select nk "best" WF
    if nk==1: 
        wfout = [wfout]
        eout = [eout]
    if verbose>1: 
        print("Restarting with",nk,"wavefunctions")
        print("Restarting energies",eout)
        wfout = gram_smith(wfout) # orthogonalize again
    if verbose>2: 
        iden = krylov_matrix_representation(1.,wfout)
        print("Orthogonality") # get the representation
        print(np.round(iden,1)) # get the representation
    return mpsarnoldi_iteration(self,Op,H,fe,
            maxit=maxit-1,
            maxde=maxde,
            verbose=verbose,
            shift = eout[0], # shift operator
            n=n,
            wfs = wfout, # initial wavefunctions
            maxdwf=maxdwf,
            ne=ne,
            wfskip=wfskip)


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



def selectwf(es,vs,wfs,fe,ne=1):
    """Select the wavefunctions that should be returned"""
    elist = [e for e in es] # list with the energies
    vlist = [v for v in vs] # list with the eigenvectors
    vstore = [] # empty list
    estore = [] # empty list
    einds = [] # indexes
    for i in range(ne): # loop over desired energies
        ie = fe(np.array(elist)) # get the desired index
        estore.append(elist[ie]) # store this energy
        vstore.append(vlist[ie]) # store this WF
        del elist[ie] # ignore in the next iteration
        del vlist[ie] # ignore in the next iteration
    wfout = [] # output wavefunctions
    for v0 in vstore: # loop over WF
        wf = 0
        for i in range(len(wfs)): 
            wf = wf + np.conjugate(v0[i])*wfs[i] # add
        wf = wf.normalize()
        wfout.append(wf.copy()) # store wavefunction
    eout = np.array(estore) # convert to array
    if ne==1: return eout[0],wfout[0]
    else: return eout,wfout


def sortwf(es,wfs,fe):
    """Sort wavefunction according to a criteria"""
    elist = [e for e in es] # list with the energies
    listinds = [i for i in range(len(es))] # list with indexes
    inds = [] # indexes
    for i in range(len(es)): # loop over desired energies
        ie = fe(np.array(elist)) # get the desired index
        inds.append(listinds[ie]) # store this index
        del listinds[ie] # ignore in the next iteration
        del elist[ie] # ignore in the next iteration
    eout = [es[i] for i in inds] # pick these energies
    wfout = [wfs[i] for i in inds] # pick these energies
    return np.array(eout),wfout # return energies and wavefunctions








def gram_smith_single(w,ws):
    """Gram smith orthogonalization for a single wavefunction"""
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
        return lg.eig(mh)



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





