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
    elif mode=="LM": # largest magnitude
        Op = lambda x: H*x # operator to apply
        def fe(es): # function return the right WF
            es = np.abs(es) # maximum energy
            return np.where(es==np.max(es))[0][0] # return the index
    elif mode=="GS": # target the ground state
        # for non-Hermitian matrices, this targets the eigenvalues
        # with most negative real part
        M = H # start with the Hamiltonian
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
        wfskip=None, # wavefunctions to skip
        shift = 0.0, # shift for the operator
        verbose=0,
        n0=100, # number of warm up iterations
        wfs = None, # initial Krylov vectors
        n=10):
    """Single iteration of the restarted arnoldi algorithm"""
    if wfs is None: wfs=[]
    if wfskip is None: wfskip=[]
    if verbose>1:
        print("Eigenvalue shift",shift)
    if len(wfs)==0: # no vectors given
        wf = self.random_mps() # random initial guess
        wf = wf.normalize()
        for i in range(n0):
            if verbose>1: print("Pre-iteration",i)
            wfold = wf.copy() # copy
            wf = Op(wf) + shift*wf # apply operator
            wf = wf.normalize()
#            mix = 0.0 # quench oscillations
#            wf = mix*wf + (1-mix)*wfold # mix
#            wf = wf.normalize()
            if verbose>1: print("Eigen",wf.aMb(H,wf))
            wf = gram_smith_single(wf,wfskip) # orthogonalize
        wf0 = None
        if n==1 and ne==1: # power method
            return wf.dot(Op(wf)),wf # return WF and energy
    else: 
        wf = wfs[-1].copy() # take the "worst" one
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
    ef,wf = selectwf(es,vs.T,wfs,fe,ne=ne) # select the wavefunctions
    if verbose>1: # print the orthogonality matrix 
        iden = krylov_matrix_representation(1.,wf)
        print("Orthogonality") # get the representation
        print(np.round(iden,1)) # get the representation
    if maxit<2: return ef,wf # if last iteration has been reached
    #########################################
    # if the algorithm is recursive, continue
    #########################################
    ef,wf = selectwf(es,vs.T,wfs,fe,ne=ne) # select the wavefunctions
    if ne==1: wf = [wf] # if just one
    ef = np.array([wfi.aMb(H,wfi) for wfi in wf]) # compute energies
    ef2 = np.array([wfi.aMb(H,H*wfi) for wfi in wf]) # compute energies square
    error = np.sqrt(np.abs(ef2-ef**2)) # compute the error
#    dwf = 1.0 # difference with respect to the initial wf
#    if wf0 is not None:
#        dwf = 1.0 - np.abs(wf0.dot(wf)) # difference between WF
    if verbose>0: 
        print("Energies",np.round(ef,2))
        print("Error in energies",np.round(error,2))
#        if wf0 is not None: print("Error in WF",dwf)
    # stop according to several criteria
    if np.sum(error)<maxde: return ef,wf # return wavefunctions 
    # if this point is reached, recall
    nk = max([2,n//2]) # number of vector to keep
    eout,wfout = selectwf(es,vs.T,wfs,fe,ne=nk) # select nk "best" WF
    if nk==1: 
        wfout = [wfout]
        eout = [eout]
    if verbose>1: 
        print("Restarting with",nk,"wavefunctions")
        print("Restarting energies",eout)
    wfout = gram_smith(wfout) # orthogonalize the basis again
    return mpsarnoldi_iteration(self,Op,H,fe,
            maxit=maxit-1,
            maxde=maxde,
            verbose=verbose,
            shift = shift, # shift operator
            n=n,
            wfs = wfout, # initial wavefunctions
            maxdwf=maxdwf,
            ne=ne,
            wfskip=wfskip)

from .krylov import krylov_matrix_representation
from .krylov import recompute_energies
from .krylov import gram_smith_single
from .krylov import gram_smith
from .krylov import diagonalize
from .krylov import rediagonalize



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



def lowest_energy(self,H,**kwargs):
    """Compute the most negative energy of a Hamiltonian,
    assuming a Hermitian Hamiltonian"""
    emax,wf = mpsarnoldi(self,H,mode="LM",n0=20,n=4,
            **kwargs) # warm up
    # most negative
    return mpsarnoldi(self,H,mode="GS",shift=-np.abs(emax),**kwargs) 


def lowest_energy_non_hermitian(self,H,n=1,npm=30,**kwargs):
    """Compute the most negative energy of a Hamiltonian,
    assuming a non Hermitian Hamiltonian"""
    emax,wf = mpsarnoldi(self,H,mode="LM",n0=npm,n=1,nwf=1,
            **kwargs) # warm up
    # most negative real part
    return mpsarnoldi(self,H,mode="GS",n0=npm,n=3*n+1,nwf=n,maxit=6,
            shift=-np.abs(emax),**kwargs)



