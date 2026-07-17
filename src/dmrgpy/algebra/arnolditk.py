import numpy as np
import scipy.linalg as lg
from . import algebra
from . import powermethod

# routines to perform iterative diagonalization




arnoldimode = "DMRG" # mode of the arnoldi method


def mpsarnoldi(self,H,wf=None,e=0.0,delta=1e-3,
        mode="GS",P=None,
        recursive_arnoldi=False,
        nwf=1, 
        **kwargs):
    """Compute an eigenvector using the Arnoldi algorithm"""
    shift = 0. # zero shift unless stated otherwise
    if mode=="ShiftInv": # target a specific energy with shift and invert
        Mi = H - e # shift
        M = self.toMPO(H,mode=arnoldimode) # accelerate
        Op = lambda x: self.applyinverse(Mi,x,delta=delta,maxn=20) # operator to apply
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
        M = self.toMPO(H,mode=arnoldimode) # accelerate
        shift = -1.5*powermethod.estimate_radius(self,H,n0=10)
#        print("Estimated radius",shift)
        if P is None: 
            from ..multioperatortk.staticoperator import StaticOperator
            Op = lambda x: M*x + shift*x # operator to apply
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
    # Op is an affine function of M (Op=a*M+b) for every mode above except
    # "ShiftInv" (Op=(M-e)^-1): affine Op shares M's eigenvectors, so its
    # own cheap Hessenberg matrix can be diagonalized directly (only the
    # eigenvalues need an -- equally cheap -- correction back to M's
    # scale). ShiftInv needs the exact (but O(n^2)) matrix representation
    # of M itself -- see mpsarnoldi_iteration_single.
    op_is_affine = mode!="ShiftInv"
    if nwf==1: # just the ground state
        return mpsarnoldi_iteration(self,Op,M,fe,ne=1,
                shift=shift,maxde=delta,op_is_affine=op_is_affine,
                **kwargs)
    else:
        if recursive_arnoldi:
            wfout = [] # empty list
            eout = [] # empty list
            for i in range(nwf): # loop over desired wavefunctions
                ei,wfi = mpsarnoldi_iteration(self,
                        Op,M,fe,ne=1,
                        wfs=[],maxde=delta,op_is_affine=op_is_affine,
                        wfskip=wfout,**kwargs)
                wfout.append(wfi[0].copy()) # store wavefunction
                eout.append(ei[0]) # store wavefunction
            eout,wfout = sortwf(eout,wfout,fe) # resort the result
            return np.array(eout),wfout # return wavefunctions
        else:
          return mpsarnoldi_iteration(self,Op,M,fe,
                  maxde=delta,op_is_affine=op_is_affine,
                  ne=nwf,**kwargs)


def mpsarnoldi_iteration(self,Op,H,fe,
        verbose=0, # verbosity
        maxde=1e-3, # maximum error in the energies
        maxit=10, # maximum number of recursive iterations
        wfs = None, # initial Krylov vectors (only wfs[0], if any, is used
                     # as the restart seed -- see mpsarnoldi_iteration_single)
        nkry_min = None, # minimum number of krylov vectors
        nkry_max = None, # maximum number of krylov vectors
        ne=1, # number of energies to return
        **kwargs # other arguments
        ):
        """Restarted Arnoldi ("thick restart"): build a Krylov subspace of
        growing size, each outer iteration reseeded directly from the ne
        Ritz vectors found so far (not collapsed into a single combined
        vector) -- keeping all ne of them as independent seeds is what
        lets a simultaneous multi-eigenvalue search (ne>1) resolve
        (near-)degenerate eigenvalues, which a single restart vector
        cannot: a lone vector's Krylov chain only ever explores one
        direction inside a degenerate eigenspace. Optional deflation
        against wfskip covers excited states found in earlier calls."""
        if nkry_min is None: nkry_min = ne + 2 # default value
        if nkry_max is None: nkry_max = 2*ne + 4 # default value
        nkry_max = max(nkry_max,nkry_min) # keep the range well ordered
        nkry = nkry_min # initialize
        seeds = wfs # None on the first call, else the previous ne Ritz vectors
        for i in range(maxit): # loop over iterations
            if verbose>0:
                print("Arnoldi iteration #",i)
                print("Number of Krylov vectors",nkry)
            es,wfs,error = mpsarnoldi_iteration_single(self,Op,H,fe,
                           ne=ne,n=ne+nkry,seeds=seeds,verbose=verbose,**kwargs)
            seeds = wfs # reseed the next restart with all ne Ritz vectors
            dnk = np.abs(np.log(np.mean(error))/np.log(maxde)) # rescaled error
            dnk = np.min([dnk,1.]) # upper cutoff
            dnk = np.max([0,dnk]) # lower cutoff
            nkry = int(np.round(dnk*nkry_max + (1.-dnk)*nkry_min)) # new number of krylov vectors
            if verbose>0:
                print(maxde)
                print("Error in Arnoldi iteration",np.round(error,3))
                print("Mean error",np.round(np.mean(error),3))
                print("Krylov update",dnk)
            if np.max(error)<maxde: break # stop if the error is smaller than the threshold
        return es,wfs # return energies and wavefunctions



def mpsarnoldi_iteration_single(self,Op,H,fe,
        ne=1, # number of energies to return
        maxdwf = 1e-3, # maximum change in the wavefunction
        maxde=1e-3, # maximum error in energy
        wfskip=None, # wavefunctions to skip (locked/already found states)
        shift = 0.0, # shift for the operator (informational only)
        verbose=0, # verbosity
        mix = 0., # mixing for the new wavefunction
        n0=10, # number of warm up iterations
        seeds=None, # seed vector(s) for the Krylov chain, up to ne of them
        n=10, # dimension of krylov space
        op_is_affine=True, # is Op an affine function of H (Op=a*H+b)?
        ntries_pm=3 # unused, kept for backwards-compatible call sites
        ):
    """Build an Arnoldi (Krylov) chain of dimension n starting from seeds
    (or a fresh warm start if seeds is None/empty), and return the ne
    Ritz pairs selected by fe together with their residual error.

    The residual is obtained from the Arnoldi/Hessenberg relation
    H Q = Q Hn + beta * q_{n+1} e_n^T, i.e. for a Ritz pair (theta,y) of
    Hn, ||H y - theta y|| = |beta * y[-1]| -- a byproduct of the last
    Krylov coefficient and the last component of the Ritz eigenvector,
    computed with a single extra Op(x) application (to complete the last
    Hessenberg column) instead of one extra Op(x) application per
    requested eigenvalue (the previous approach recomputed H*wf for every
    selected Ritz vector every outer iteration).

    When Op is only an approximate/non-affine proxy for H (op_is_affine=
    False, currently just mode="ShiftInv", where Op=(H-e)^-1), the basis
    is still built by extending with Op (it converges towards the right
    invariant subspace fastest), but selection/eigenvalues must come from
    H's own O(n^2) matrix representation on that basis instead of Op's
    (cheap O(n) Hessenberg one) -- diagonalizing Op's own matrix would
    pick Ritz pairs by Op's eigenvalues, which fe cannot interpret."""
    if wfskip is None: wfskip=[]
    if verbose>1:
        print("Eigenvalue shift",shift)
    if not seeds: # no seed given, get ne of them with a warm-up power method
        seeds = arnoldi_warmup_multi(self,Op,H,n0,ne,wfskip,error=maxde*10,
                verbose=verbose) if n0>0 else \
                [random_state(self,orthogonal=wfskip if wfskip else None)
                        for _ in range(ne)]
    else:
        seeds = [gram_smith_single(wf,wfskip) if wfskip else wf.normalize()
                for wf in seeds]
    if mix!=0.: # finite mixing, to escape stagnation on restart
        seeds = [gram_smith_single((1.-mix)*wf + mix*random_state(self),
                    wfskip) if wfskip else
                 ((1.-mix)*wf + mix*random_state(self)).normalize()
                 for wf in seeds]
    wfs,hmat,beta = build_arnoldi_chain(self,Op,seeds,n,wfskip=wfskip,
            verbose=verbose)
    if op_is_affine: # cheap: hmat was already built for free from Op
        (es,vs) = diagonalize(hmat)
    else: # exact but O(n^2): H's own representation on the Op-built basis
        (es,vs) = diagonalize(krylov_matrix_representation(H,wfs))
    if verbose>2:
        iden = krylov_matrix_representation(1.,wfs)
        print("Krylov orthogonality") # get the representation
        print(np.round(iden,1)) # get the representation
    estore,vstore = krylov.select_states(es,vs.T,fe,ne=ne) # select Ritz pairs
    # residual from Op's Hessenberg relation -- exact for the common
    # affine modes (Op's residual is H's, up to the known affine scale),
    # an approximate proxy for mode="ShiftInv" (only used to drive the
    # adaptive Krylov-size/stopping heuristic, never the returned answer)
    error = np.array([np.abs(beta*v0[-1]) for v0 in vstore])
    wf = krylov.unitary_transformation(vstore,wfs) # build the wavefunctions
    # for the affine modes, es are Op's (possibly shifted/scaled)
    # eigenvalues -- recompute the reported energies against the caller's
    # true H (cheap: ne aMb calls, not the O(n^2) needed for selection)
    ef = np.array(estore) if not op_is_affine else \
            np.array([wfi.aMb(H,wfi) for wfi in wf])
    if verbose>0:
        print("Energies",np.round(ef,3))
        print("Residuals",np.round(error,6))
    return ef,wf,error


def arnoldi_warmup_multi(self,Op,H,n0,ne,wfskip,error=1e-2,verbose=0):
    """ne independent, mutually deflated warm-start vectors: restores the
    old algorithm's block start, needed for a simultaneous (ne>1)
    multi-eigenvalue search to reliably resolve (near-)degenerate
    eigenvalues -- a single shared warm-start vector's Krylov chain only
    ever explores one direction inside a degenerate eigenspace, however
    long it is grown (confirmed to occasionally miss a degenerate
    partner in mode="GS" non-Hermitian searches)."""
    seeds = []
    for i in range(ne):
        seeds.append(arnoldi_warmup(self,Op,H,n0,wfskip+seeds,
                error=error,verbose=verbose))
    return seeds


def arnoldi_warmup(self,Op,H,n0,wfskip,error=1e-2,verbose=0):
    """Cheap warm start: n0 power-iteration steps with Op (the same
    operator the Arnoldi chain itself uses -- correct for every mode,
    including shift-invert/projected ones where Op is not simply H*x),
    deflated against wfskip, to get a decent seed direction. Replaces
    running one independent power-method chain per desired eigenvalue
    (as the old block warm-start did) with a single shared one, and
    breaks out as soon as the Rayleigh quotient stops moving."""
    wf = random_state(self,orthogonal=wfskip if wfskip else None)
    eold = -1e20
    for i in range(n0):
        wf = Op(wf)
        if wfskip: wf = gram_smith_single(wf,wfskip)
        wf = wf.normalize()
        ei = wf.aMb(H,wf).real
        if verbose>2: print("Warmup energy",ei)
        if np.abs(eold-ei)<error: break
        eold = ei
    return wf


def build_arnoldi_chain(self,Op,seeds,n,wfskip=None,verbose=0):
    """Build an orthonormal Krylov (Arnoldi-like) basis of dimension n,
    starting from the given seeds (1 or more, mutually orthonormalized
    here) and extending with Op as needed to reach n vectors, deflating
    against wfskip throughout. Returns:
      wfs  -- the n orthonormal basis vectors
      hmat -- the n x n matrix representation, built incrementally from
              the (modified, twice-repeated) Gram-Schmidt coefficients
              instead of the O(n^2) all-pairs sandwich contractions
              krylov_matrix_representation would need. With a single seed
              this is exact in exact arithmetic (Op(wfs[k]) then has no
              component outside span(wfs[0..k+1]), the standard Krylov
              argument); with multiple independent seeds the block
              formed by their mutual couplings is generally dense rather
              than banded, which the row-at-a-time construction below
              already produces correctly (each row's Gram-Schmidt runs
              against the *entire* current basis, not just wfs[0..k]).
      beta -- the residual Hessenberg coefficient h_{n+1,n} for the last
              vector, i.e. the norm of Op(wfs[n-1]) after removing its
              component along wfs -- combined with a Ritz eigenvector's
              last component this gives the exact residual with a single
              extra Op(x) application (to fill in the last Hessenberg
              column), rather than one extra application per selected
              eigenvalue.
    Uses the same (transposed) convention as krylov_matrix_representation
    (hmat[k,i] = <wfs[i]|Op|wfs[k]>) so the existing diagonalize()/
    select_states()/unitary_transformation() helpers can be reused as-is.
    """
    if wfskip is None: wfskip = []
    wfs = gram_smith(seeds)[:n] # mutually orthonormalize the given seeds
    hmat = np.zeros((n,n),dtype=np.complex128)
    def norm(v): return np.sqrt(v.dot(v).real) # works for MPS and ED State alike
    def orthogonalize(v,basis):
        # modified Gram-Schmidt, repeated once ("twice is enough") --
        # restores orthogonality lost to roundoff/MPS truncation, cheaper
        # and more robust than solving a generalized eigenvalue problem
        # against a near-singular overlap matrix
        coeffs = np.zeros(len(basis),dtype=np.complex128)
        for _ in range(2):
            for i,b in enumerate(basis):
                c = b.dot(v)
                coeffs[i] += c
                v = v - c*b
        return v,coeffs
    for k in range(n-1):
        v = Op(wfs[k])
        if wfskip: v,_ = orthogonalize(v,wfskip) # deflate found states
        # orthogonalize against the *entire* current basis (not just
        # wfs[0..k]): with multiple independent seeds, Op(wfs[k]) can
        # have components along seeds with index >k too
        v,coeffs = orthogonalize(v,wfs)
        hmat[k,:len(wfs)] = coeffs
        if len(wfs)<n: # still need to grow the basis to reach size n
            beta = norm(v)
            if beta<1e-8: # invariant subspace hit, use a random vector
                if verbose>1: print("Invariant subspace found, using random vector")
                v = random_state(self,orthogonal=wfs+wfskip)
                beta = norm(v)
            # hmat follows krylov_matrix_representation's convention,
            # hmat[a,b] = <wfs[b]|Op|wfs[a]> -- so the coupling from
            # wfs[k] to the newly appended vector belongs at column
            # len(wfs) of row k (not the k+1,k transpose -- that slip
            # drops the coupling entirely and silently produces a wrong
            # Hessenberg matrix/wrong energies)
            hmat[k,len(wfs)] = beta
            wfs.append(v.normalize())
    # complete the last row: this single extra Op(x) call gives beta_last,
    # which yields the residual of ALL ne selected Ritz pairs at once
    # (see mpsarnoldi_iteration_single)
    v = Op(wfs[n-1])
    if wfskip: v,_ = orthogonalize(v,wfskip)
    v,coeffs = orthogonalize(v,wfs)
    hmat[n-1,:] = coeffs
    beta_last = norm(v)
    return wfs,hmat,beta_last


from .krylov import krylov_matrix_representation
from .krylov import gram_smith_single
from .krylov import gram_smith
from .krylov import diagonalize
from .krylov import rediagonalize
from . import krylov

selectwf = krylov.selectwf

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



def lowest_energy(self,H,n=1,**kwargs):
    """Compute the most negative energy of a Hamiltonian,
    assuming a Hermitian Hamiltonian"""
    return mpsarnoldi(self,H,mode="GS",nwf=n,**kwargs) 


def lowest_energy_non_hermitian(self,H,n=1,**kwargs):
    """Compute the most negative energy of a Hamiltonian,
    assuming a non Hermitian Hamiltonian"""
#    emax,wf = mpsarnoldi(self,H,mode="LM",n0=npm,n=1,nwf=1,maxit=1,
#            **kwargs) # warm up
    # most negative real part
    return mpsarnoldi(self,H,mode="GS",nwf=n,**kwargs)





def most_positive_energy(self,H,npm=10,verbose=0,dt=0.1,n=1,**kwargs):
    """Compute the most negative energy of a Hamiltonian,
    assuming a Hermitian Hamiltonian"""
    emax,wf = power_method(self,H,n0=npm,verbose=verbose) # PM estimate
    Hs = H -np.abs(emax) # shift the energy
    Hs = Hs/emax # scale the Hamiltonian to approx [0,1]
    Hexp = 1. + dt*Hs # approximation to the exponential
    # most negative
    (es,wfs) = mpsarnoldi(self,Hexp,mode="LM",n0=npm,
            n=2*n+1,nwf=n,maxit=3,
            verbose=verbose,**kwargs)
    es = [wfi.aMb(H,wfi) for wfi in wfs] # compute energies
    return np.array(es),wfs


def most_negative_energy(self,H,**kwargs):
    """Compute the most negative energy"""
    es,wfs = most_positive_energy(self,-H,**kwargs)
    return -es,wfs

#lowest_energy_non_hermitian = most_negative_energy



def random_state(self,orthogonal=None):
    """Generate a random state"""
    if orthogonal is None:
        if arnoldimode=="DMRG": return self.random_mps()
        elif arnoldimode=="ED": return self.get_ED_obj().random_state()
    else:
        while True: # infinite loop
            wf = random_state(self) # generate a random state
            wf = gram_smith_single(wf,orthogonal)
            if wf is not None: return wf # return this wavefunction


