from . import mps
import numpy as np


def best_gs(sc,n=1):
    """Compute many ground states, and retain only the best one"""
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        sc.computed_gs = False # initialize
        # force a fresh session solve: with the Hamiltonian unchanged the
        # send-cache in gs_energy_single would otherwise return the
        # previous iteration's cached energy instead of re-running DMRG
        sc._session_ham_cache = None
        e0 = sc.gs_energy() # ground state energy
        if e0<emin: # only keep this one if it actually improves on the best
            wf0 = sc.wf0 # copy wavefunction
            emin = e0
    sc.set_initial_wf(wf0) # set the wavefunction

def gs_energy_bestof(self,**kwargs):
    maxm = self.maxm
    maxm2 = maxm
    self.maxm = maxm2
    gs_energy_many(self,**kwargs)
    self.maxm = maxm
    return gs_energy_single(self,reconverge=True)

def gs_energy_many(self,n=20,**kwargs):
    """Compute many ground states, and retain only the best one"""
    wfs = [] # list with GS
    es = [] # list with energies
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        self.computed_gs = False # initialize
        self.gs_from_file = False
        self.skip_dmrg_gs = False
        self._session_ham_cache = None # force a fresh session solve (see best_gs)
        e0 = gs_energy_single(self,reconverge=False,
                **kwargs) # ground state energy
        if e0<emin: 
            wf0 = self.wf0.copy() # copy wavefunction
            emin = e0
        else: self.wf0.clean() # remove wavefunction
        print("Best of",n,i,e0,emin)
    wf0.rename("psi_GS.mps")
    self.set_initial_wf(wf0) # set the wavefunction
    print("Final energy",self.vev(self.hamiltonian).real)
    return emin

def gs_energy_single(self,wf0=None,reconverge=None,maxde=None,maxdepth=5):
    """
    Return the ground state energy via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain): the Hamiltonian/sweep params/
    wavefunction are passed as in-memory arguments to self._session.
    """
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    # Only re-send the Hamiltonian when it (or the MPO bond dimension it
    # is built with) actually changed since the last send to this same
    # session: the session's set_hamiltonian() invalidates its energy
    # and band-edge caches unconditionally, so an unconditional re-send
    # here turned every get_dynamical_correlator() call's internal
    # ground-state re-verification into a real warm re-sweep and forced
    # KPM to redo its band-edge DMRG on every call even with an
    # unchanged Hamiltonian. Keying on the session object itself (by
    # reference) makes the cache self-invalidating whenever a fresh
    # Chain is created (setup_cpp/setup_python/__deepcopy__); comparing
    # the to_terms() output (not the MultiOperator identity) catches
    # in-place mutation of self.hamiltonian. Every solver parameter that
    # a re-run would pick up (maxm, nsweeps, cutoff, noise, and the MPO
    # bond dimension the Hamiltonian is built with) is part of the key:
    # the session's energy cache survives a skipped re-send, so a user
    # bumping any of these between bare gs_energy() calls must get a
    # fresh solve, not the cached energy computed under the old params.
    terms = self.hamiltonian.to_terms()
    key = (self.maxm,self.nsweeps,self.cutoff,self.noise,
           max(self.maxm,self.mpomaxm),terms)
    cache = getattr(self,'_session_ham_cache',None)
    if cache is None or cache[0] is not self._session or cache[1]!=key:
        self._session.set_hamiltonian(terms)
        self._session_ham_cache = (self._session,key)
    if wf0 is not None:
        self._session.set_wavefunction(wf0.cpp_handle)
    if reconverge is not None: # overwrite skip_dmrg_gs
        self.skip_dmrg_gs = not reconverge # if the computation should be rerun
    out = self._session.gs_energy(skip_dmrg=self.skip_dmrg_gs)
    self.e0 = out # store ground state energy
    self.computed_gs = True
    self.sites_from_file = True
    self.gs_from_file = True
    self.skip_dmrg_gs = True
    wf0 = mps.MPS(MBO=self,cpp_handle=self._session.gs_wavefunction()).copy()
    self.set_initial_wf(wf0) # set the initial wavefunction
    if maxde is not None: # enforce a maximum fluctuation in the energy
      e = self.vev(self.hamiltonian)
      e2 = self.vev(self.hamiltonian,npow=2)
      de = np.sqrt(abs(e2-e**2)) # fluctuation in the energy
      de = de/self.ns # normalize by the number of sites
      if de>maxde and maxdepth>0: # if a maximum energy fluctuation
          maxm,nsweeps = self.maxm,self.nsweeps
          noise = self.noise
          print("Energy fluctuation = ",de,maxm)
          self.maxm = maxm*2
          self.nsweeps = 2 # just two sweeps
          self.noise = 0.0
          gs_energy_single(self,maxde=maxde,reconverge=True,
                  maxdepth=maxdepth-1) # execute again
          self.maxm = maxm
          self.nsweeps = nsweeps # restore
          self.noise = noise
    self.computed_gs = True # ground state has been computed
    return out # return energy


def gs_energy(self,**kwargs):
    if self.is_hermitian(self.hamiltonian): # put a check for Hermitian
        if self.itensor_version in (2,3,"python"): # C++ or pure-Python version
            return gs_energy_cpp(self,**kwargs)
        elif self.itensor_version=="julia_live":
            from .mpsjulialive.groundstate import get_gs_dmrg
            e0,wf0 = get_gs_dmrg(self,**kwargs)
            self.wf0 = wf0
            return e0
        else: raise
    else: # non-Hermitian excited states
        # Julia version has its own function
        if self.itensor_version=="julia_live":
            from .mpsjulialive.groundstate import get_gs_dmrg
            e0,wf0 = get_gs_dmrg(self,ishermitian=False,**kwargs)
            self.wf0 = wf0
            return e0
        elif self.itensor_version in (2,3,"python"): # real non-Hermitian DMRG
            from .nhdmrg import gs_energy_nhdmrg
            return gs_energy_nhdmrg(self,**kwargs)
        else: # any other backend falls back to Krylov
            es,ws = self.get_excited_states(n=1,**kwargs)
            self.computed_gs = True
            self.e0 = es[0]
            self.wf0 = ws[0].copy() # copy wavefunction
            return self.e0




def gs_energy_generalized(self,A,lam0=None):
    """Smallest generalized eigenvalue lambda solving H|psi>=lambda*A|psi>
    for this chain's own Hamiltonian (self.hamiltonian, already set via
    set_hamiltonian()) and a Hermitian positive-definite metric operator
    A (a MultiOperator, same calling convention as vev()). Stores the
    resulting wavefunction as this chain's ground state, mirroring
    gs_energy_single()'s own wf0 handling.

    Implemented for itensor_version="python" (pyitensor/dmrg.py's
    dmrg_generalized()) and itensor_version=3 (mpscpp3/chain_session.h's
    Chain::gs_energy_generalized, a line-for-line port of the same
    self-consistent Lagrange-multiplier algorithm against ITensor v3's
    own dmrg()/Sweeps/sum()) -- see either docstring for the derivation.
    mpscpp2 (itensor_version=2) and julia_live don't have this session
    method yet. There is also no ED implementation of this method at all
    (unlike vev()/gs_energy()/..., which all honor self.mode="ED" for
    cross-validation) -- self.mode="ED" is rejected explicitly below
    rather than silently ignored.

    Non-Hermitian self.hamiltonian is dispatched to a separate solver,
    nhdmrg.py's nhdmrg_generalized() (the non-Hermitian, complex-lambda,
    biorthogonal-quotient generalization of NH-DMRG, mirroring how this
    function itself generalizes plain gs_energy()) -- implemented on the
    same itensor_version="python"/3 pair as the Hermitian path above (no
    mpscpp2 support either way).

    CAVEAT (found via code review, not fixed -- see this codebase's usual
    "document the quirk" convention rather than adding a state-tracking
    flag threaded through every consumer): self.wf0/self.e0/computed_gs
    afterward hold the eigenvector/eigenvalue of the *shifted* problem
    H-lambda*A (or its biorthogonal NH counterpart), not a plain
    eigenstate of self.hamiltonian alone. Every other method that treats
    self.wf0 as an ordinary ground state -- get_excited_states() (its
    overlap-penalty anchor), any dynamical/KPM correlator, NH-KPM
    (nonhermitian/kpm.py) -- has no way to detect this and will silently
    build on the wrong reference state if called afterward without first
    recomputing a genuine ground state (gs_energy()/nhdmrg(), which reset
    wf0 to a real eigenstate of self.hamiltonian). Call
    gs_energy_generalized() as the last step of a calculation, or
    explicitly recompute the plain ground state before using any other
    method that reads self.wf0."""
    if self.mode=="ED":
        raise NotImplementedError(
            "gs_energy_generalized has no ED implementation -- unset "
            "self.mode (or set it to \"DMRG\") to use the DMRG solver")
    if self.itensor_version not in (3,"python"):
        raise NotImplementedError(
            "gs_energy_generalized is only implemented for "
            "itensor_version=3 or 'python' so far -- call "
            "chain.setup_cpp(version=3) or chain.setup_python() first")
    if self._session is None:
        # itensor_version==3 but no compiled extension for it (sites.py's
        # initialize() leaves self._session as None in that case, the
        # same "extension not compiled" state mode.py's own get_mode()
        # falls back to ED for elsewhere) -- there is no ED fallback for
        # this method, so fail with an actionable message instead of an
        # AttributeError on the calls below.
        raise RuntimeError(
            "gs_energy_generalized needs a compiled ITensor v3 extension "
            "(itensor_version=3) but none is available for this chain -- "
            "run `python install.py --itensor-version=3`, or call "
            "chain.setup_python() to use the pure-Python backend instead")
    if self.hamiltonian is None:
        raise RuntimeError("gs_energy_generalized called before set_hamiltonian")
    if not self.is_hermitian(self.hamiltonian):
        # Non-Hermitian H: dispatch to the NH-DMRG generalized solver
        # (nhdmrg.py's nhdmrg_generalized()/gs_energy_generalized_nhdmrg())
        # -- mirrors gs_energy()'s own non-Hermitian dispatch to NH-DMRG.
        # Without this branch, the local two-site solver in both backends
        # would silently Hermitize its effective-Hamiltonian matrix before
        # diagonalizing, producing a well-defined but physically
        # meaningless "eigenvalue" with no warning at all.
        #
        # Deliberately checked *before* the itensor_version==3-and-ns<3
        # guard below: that guard exists only because the Hermitian path
        # calls real ITensor v3 dmrg() (whose two-site sweep aborts the
        # whole process for short chains), but nhdmrg_generalized()'s own
        # two-site sweep is hand-rolled directly against
        # arnoldi_smallest_real/manual ITensor contractions and never
        # calls dmrg() at all -- confirmed directly, a 2-site chain runs
        # it without aborting -- so the non-Hermitian path must not be
        # rejected for a short chain it can actually handle fine.
        from .nhdmrg import gs_energy_generalized_nhdmrg
        return gs_energy_generalized_nhdmrg(self,A,lam0=lam0)
    if self.itensor_version==3 and self.ns<3:
        # ITensor v3's two-site dmrg() aborts the whole process (SIGABRT,
        # "LocalOp is default constructed") for chains shorter than 3
        # sites -- see mode.py's own itensor_version==3 guard, which
        # exists specifically to route every *other* DMRG entry point
        # around this by falling back to ED. There is no ED fallback
        # here, so this must be rejected explicitly (mirrored again,
        # defense in depth, by Chain::gs_energy_generalized itself on
        # the C++ side) rather than silently crashing the interpreter.
        # Hermitian-path-only (see the non-Hermitian branch above for why
        # NH-DMRG doesn't need this guard).
        raise RuntimeError(
            "gs_energy_generalized: ITensor v3's two-site DMRG can't "
            "handle a chain this short (n=%d < 3 sites) -- use "
            "itensor_version=\"python\" instead"%self.ns)
    from . import multioperator
    A = multioperator.obj2MO(A)
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    self._session.set_hamiltonian(self.hamiltonian.to_terms())
    if self.itensor_version=="python": # pyitensor accepts lam0=None directly
        session_lam0 = lam0
    else: # the compiled v3 binding takes a plain float, NaN meaning "unset"
        session_lam0 = float('nan') if lam0 is None else lam0
    lam = self._session.gs_energy_generalized(A.to_terms(),lam0=session_lam0)
    self.e0 = lam
    wf0 = mps.MPS(MBO=self,cpp_handle=self._session.gs_wavefunction()).copy()
    self.set_initial_wf(wf0) # set the initial wavefunction (resets computed_gs)
    self.computed_gs = True # ...so this must come after set_initial_wf, not
                             # before -- mirrors gs_energy_single's own
                             # ordering (its trailing re-assertion after the
                             # same set_initial_wf() reset), which an earlier
                             # version of this function got backwards: setting
                             # computed_gs=True *before* set_initial_wf() left
                             # it silently False on return, so the next
                             # ordinary gs_energy()/get_gs() call would see
                             # computed_gs==False and quietly re-run a plain
                             # ground-state DMRG solve (warm-started from this
                             # method's own wf0), overwriting self.wf0/self.e0
                             # with a different quantity than the caller asked
                             # for -- confirmed directly via execution.
    return lam


def gs_energy_cpp(self,policy="single",**kwargs):
    """Compute ground state with C++ DMRG"""
    if policy=="single":
        return gs_energy_single(self,**kwargs)
    if policy=="many":
        return gs_energy_many(self,**kwargs)
    else: raise






def lowest_energy(self,h):
    """Return the lowest energy of the Hamiltonian"""
    raise # not finished yet
    self.execute(lambda: h.write("hamiltonian.in")) # write Hamiltonian
    task = {"GS":"true",
            }
    self.task = task
    self.write_task()
    self.run() # perform the calculation
    out = self.execute(lambda: np.genfromtxt("GS_ENERGY.OUT"))
    return out



def get_gs_manifold(MBO,n=2,tol=1e-3,**kwargs):
    """Return the ground state manifold, i.e. all the states with the
    lowest energy"""
    (es,wfs) = MBO.get_excited_states(n=n,**kwargs)
    e0 = es[0] # ground state
    ngs = len(es[np.abs(es-e0)<tol]) # number of ground states
    if ngs<n: # all the GS found
        wfo = []
        for (e,w) in zip(es,wfs):
            if np.abs(e-e0)<tol: wfo.append(w)
        return wfo
    else: 
        print("Recalling with ",n+1,"states")
        return get_gs_manifold(MBO,n=n+1,tol=tol,**kwargs)




def set_gs(MBO,wf):
    """Set the ground state int he object"""
    mode = wf.mode # get the mode
    if mode=="DMRG": # DMRG mode
        MBO.computed_gs = True 
        MBO.wf0 = wf.copy() # set the wavefunction
    elif mode=="ED": # ED mode
        MBO.get_ED_obj() # generate the ED object
        MBO.ED_obj.computed_gs = True # comptued GS
        MBO.ED_obj.wf0 = wf.v.copy() # copy the array

