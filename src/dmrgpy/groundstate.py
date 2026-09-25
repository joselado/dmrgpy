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

def ramp_key(self):
    """The bond-ramp settings, as part of gs_energy_single()'s send-cache
    key: the session caches its own energy across a skipped Hamiltonian
    re-send, so a user changing the ramp between two bare gs_energy()
    calls must get a fresh solve rather than the energy computed under the
    old schedule -- exactly the same reasoning as for
    maxm/nsweeps/cutoff/noise there."""
    return (self.bond_ramp,self.bond_ramp_start,self.bond_ramp_fraction,
            self.bond_ramp_noise_decay)


def sector_key(self):
    """The conserved sector, as part of gs_energy_single()'s send-cache key:
    the session caches its own energy across a skipped Hamiltonian re-send,
    so a user changing sector between two bare gs_energy() calls must get a
    fresh solve rather than the previous sector's energy. Same reasoning as
    ramp_key() above (set_conserved_sector() also calls restart(), which
    clears the cache outright -- this keeps the key honest regardless)."""
    sector = getattr(self,"conserved_sector",None)
    return None if not sector else tuple(sorted(sector.items()))


def solver_key(self):
    """The solver parameters a stored ground state was computed under.

    Same fields as gs_energy_single()'s own send-cache key, minus the
    Hamiltonian itself: changing that goes through set_hamiltonian(),
    which resets computed_gs outright, restart=False included (it used to
    keep computed_gs, so every read but the public correlator answered for
    the old Hamiltonian, 2026-09-24c audit, finding 5). Used by gs_is_current() below to
    decide whether a cached self.e0/self.wf0 still answers the question
    being asked."""
    return (self.maxm,self.nsweeps,self.cutoff,self.noise,
            max(self.maxm,self.mpomaxm),ramp_key(self),sector_key(self))


def gs_is_current(self):
    """True if the stored ground state may be returned as-is.

    computed_gs alone is not enough: gs_energy()/get_gs() used to short-
    circuit on it, which returned before gs_energy_single() -- and hence
    before its send-cache, which does key on maxm/nsweeps/cutoff/noise/
    ramp/sector -- was ever entered. The textbook convergence check

        for m in [10,20,40]: sc.maxm = m; print(sc.gs_energy())

    therefore printed the m=10 energy three times, silently: no warning,
    a perfectly plausible flat curve, and (measured on a 6-site
    Heisenberg chain) an answer 5% above the true ground state. The
    fluctuation-based retry inside gs_energy_single() bumps maxm the same
    way and was equally affected. Comparing the solver parameters against
    the ones in force when the state was stored makes the short circuit
    honest while keeping a repeated call with unchanged parameters as
    cheap as it was before.

    A state the caller injected on a session backend (set_gs(),
    set_initial_wf(), set_initial_wf_guess(), see mark_injected() below)
    is never current until gs_energy_single() has taken it: that is the
    step that hands it to the DMRG session, which KPM, TD and the
    excited-state search read instead of the Python-side wf0, and that
    gives it its own energy <wf|H|wf>. Returning it before that step is
    what made every correlator measure the session's own solved state
    after a set_gs() (2026-09-24b hole hunt, findings 11 and 12).

    A state with no recorded key and no injection mark was put there by
    a backend with no session, a restored julia_live snapshot or its
    solve, and is returned as-is, since there is nothing to hand it to."""
    if not self.computed_gs: return False
    if pending_injection(self) is not None: return False # see above
    key = getattr(self,"_gs_solver_key",None)
    if key is None: return True # no session to hand it to, see above
    return key==solver_key(self)


def mark_injected(self,wf,reconverge=False,supplied=True):
    """Store `wf` (a copy of it) as this chain's state, marked as injected
    by the caller rather than produced by a solve.

    This is the contract the file-based backend had and the pybind port
    lost (2026-09-24b hole hunt, findings 11 and 12): the next ground-state
    read, gs_energy() or get_gs(), reaches gs_energy_single(), which hands
    a copy of the state to the session and either takes it unswept, with
    e0 = <wf|H|wf> (reconverge=False: set_gs(), set_initial_wf()), or
    sweeps from it (reconverge=True: set_initial_wf_guess()). Only the
    public setters call this; a solver storing its own result assigns
    self.wf0 directly, since its state is already the session's.

    The mark holds the injected object itself and counts only while
    self.wf0 is that object, so anything that replaces or clears wf0 --
    restart(), a solve, a backend switch -- retires it without having to
    know it exists. On a chain with no session (julia_live, or a C++
    backend that fell back to ED) there is nothing to hand the state to,
    and the old behaviour stands: set_gs() stores it as current and
    set_initial_wf() leaves the next solve to that backend. Returns
    whether the state was marked.

    supplied=False is for the library's own re-marks of a state it
    computed (promote_to_dense): the state is taken unswept like any
    other, but does not count as the caller's, which is what
    submode="SECTOR" asks (see state_supplied())."""
    self.wf0 = wf.copy()
    if getattr(self,"_session",None) is None:
        self._gs_injected = None
        return False
    self._gs_injected = ("reconverge" if reconverge else "skip",self.wf0,
                         bool(supplied))
    self.computed_gs = False
    self._gs_solver_key = None
    return True


def pending_injection(self):
    """"skip" or "reconverge" when the Python-side state was injected by
    the caller and has not been handed to the session yet, None
    otherwise (see mark_injected())."""
    mark = getattr(self,"_gs_injected",None)
    if mark is None or mark[1] is not self.wf0: return None
    return mark[0]


def state_supplied(self):
    """True when the chain's current state was set by the caller (set_gs(),
    set_initial_wf(), gs_energy(wf0=x, reconverge=False)) rather than
    solved, pending or already taken; retired by the next solve, by
    restart() and by set_hamiltonian(). submode="SECTOR", which measures
    the reference sector's own ground state, reads it to refuse a state it
    cannot see (2026-09-24c audit, finding 8)."""
    mark = getattr(self,"_gs_injected",None)
    if pending_injection(self) is not None:
        return mark[2] if len(mark)>2 else True
    return bool(getattr(self,"_gs_supplied",False))


def detached_copy(wf):
    """A copy of an MPS that no session call can mutate behind the
    caller's back. On itensor_version="python" set_wavefunction() stores
    the handle it is given and the next sweep rewrites that MPS's tensor
    list in place, so handing it the caller's own handle made
    gs_energy(wf0=x) move x itself onto the ground state (<H> of x went
    from -0.987 to -1.000 on a 3-site Heisenberg chain). MPS.copy()
    duplicates that list on "python" and is free on the C++ backends,
    whose MPS has value semantics."""
    return wf.copy()


def _session_parameters(self):
    """Hand the chain's solver parameters to the session (every entry
    point that may sweep needs them)."""
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    # Bond-dimension ramp for the ground-state sweep schedule, see
    # Many_Body_Chain.__init__ (manybodychain.py) for what it does and
    # Chain::make_sweeps_ramped() / pyitensor's _make_sweeps_ramped() for
    # the schedule itself. hasattr-guarded so an out-of-date compiled
    # extension (one built before this method existed) keeps working: it
    # then simply uses the C++-side defaults, which are the same as the
    # Python-side ones.
    if hasattr(self._session,"set_bond_ramp"):
        self._session.set_bond_ramp(self.bond_ramp,self.bond_ramp_start,
                                    self.bond_ramp_fraction,
                                    self.bond_ramp_noise_decay)


def _take_injected_state(self,wf,supplied=True):
    """Make `wf` this chain's ground state, unswept, on both sides.

    The session gets a detached copy with set_wavefunction(), and e0 is
    <wf|H|wf>, which is what the file-based backend's get_gs_energy()
    returned for a state read back under skip_dmrg_gs. The session's own
    energy is NOT used: set_wavefunction() drops it, so the session's
    gs_energy(skip_dmrg=True) would then run a sweep from the injected
    state, which is finding 13 of the 2026-09-24b hole hunt.

    The same drop used to leave one more sweep in the session: the lower
    band edge KPM rescales with (and the excited-state search sets its
    penalty weight from) was filled lazily by gs_energy(skip_dmrg=True),
    so the first KPM call after a push swept the pushed state in place
    (|<x|session>|^2 = 0.989 on "python", 0.987 on v3, for x not an
    eigenstate). 867e2b4 pre-filled both edges with excited_states(1)
    before every push, which cost an upper-edge -H solve and an energy
    fluctuation it threw away on the first read after an injection, 7.1 s
    on a 24-site "python" chain for a 0.04 s <x|H|x> (2026-09-24c audit,
    finding 9). Each session's minimum_energy() now solves from a fresh
    start when the state it holds has no energy, and puts the state back,
    so nothing needs filling here and the edges are paid by the first
    call that reads them. The edges are the Hamiltonian's, from a solve,
    not the injected state's energy; for a member of a degenerate ground
    manifold the two coincide."""
    _session_parameters(self)
    send_hamiltonian(self) # precondition: H on the session
    hermitian = self.is_hermitian(self.hamiltonian)
    self._session.set_wavefunction(detached_copy(wf).cpp_handle)
    e = self.aMb(wf,self.hamiltonian,wf)/self.overlap(wf,wf)
    if hermitian: e = float(np.real(e))
    self.e0 = e
    self.wf0 = wf
    self._gs_injected = None
    self._gs_supplied = supplied
    self.computed_gs = True
    self.sites_from_file = True
    self.gs_from_file = True
    self.skip_dmrg_gs = True
    self._gs_solver_key = solver_key(self)
    return e


def send_hamiltonian(self):
    """Put self.hamiltonian on the session, if it is not already there.

    Split out of gs_energy_single() so that the paths which drive the
    session DIRECTLY -- the real-time quench/evolve entry points in
    timedependent.py -- can establish the same precondition. They pass
    the Hamiltonian's terms to quench_tdvp()/quench_tebd()/... as an
    argument, but those C++/pyitensor methods start the trajectory from
    the session's OWN ground state (get_gs()), which needs set_hamiltonian
    to have been called. On a chain whose ground state had not been solved
    yet, that meant `Chain::gs_energy called before set_hamiltonian` --
    ITensor's Error(), i.e. abort(), taking the interpreter with it, from
    a plain `evolution_DC(mode="DMRG")` on a freshly built chain. It
    "worked" for every caller that happened to touch gs_energy() first,
    which is every test and example in this repo, and is why it survived
    both the 2026-08 and 2026-09 audits. Found writing the regression test
    for the 2026-09 audit's finding #7.

    The caching rationale below is the original from gs_energy_single, and
    is unchanged -- re-sending invalidates the session's energy and
    band-edge caches, so it must stay conditional.
    """
    key,terms = _send_key(self)
    if _sent(self,key): return
    if terms is None: # an already-built MPO
        if not hasattr(self._session,"set_hamiltonian_mpo"):
            raise NotImplementedError(
                "set_hamiltonian was given an already-built MPO "
                "(StaticOperator), which this backend cannot accept -- "
                "only itensor_version=3 implements set_hamiltonian_mpo. "
                "Pass a MultiOperator instead, or switch backend.")
        self._session.set_hamiltonian_mpo(self.hamiltonian.cpp_handle)
    else:
        self._session.set_hamiltonian(terms)
    self._session_ham_cache = (self._session,key)


def _send_key(self):
    """(key, terms): send_hamiltonian()'s cache key, and the term list to
    send (None for a Hamiltonian that is already an MPO, a StaticOperator,
    keyed on its handle instead)."""
    from .multioperatortk.staticoperator import StaticOperator
    base = (self.maxm,self.nsweeps,self.cutoff,self.noise,
            max(self.maxm,self.mpomaxm),ramp_key(self),sector_key(self))
    if isinstance(self.hamiltonian,StaticOperator):
        return base+(id(self.hamiltonian.cpp_handle),),None
    terms = self.hamiltonian.to_terms()
    return base+(terms,),terms


def hamiltonian_on_session(self):
    """True when self.hamiltonian, under the current solver parameters, is
    what this chain last sent to its current session, i.e. when
    send_hamiltonian() would not re-send it."""
    return _sent(self,_send_key(self)[0])


def _sent(self,key):
    cache = getattr(self,'_session_ham_cache',None)
    return (cache is not None and cache[0] is self._session
            and cache[1]==key)


def ground_state_on_session(self):
    """Make sure the ground state every session correlator reads is the
    chain's current one, on both sides, and that the session has the
    Hamiltonian; call it before driving the session directly.

    Three cases. A state that is current and whose Hamiltonian is the
    session's costs nothing, no session call at all, so a repeated
    correlator on a solved chain does not re-sweep. A state the caller
    injected (set_gs(), set_initial_wf()) is handed to the session by
    get_gs() through gs_energy_single(), unswept. And a state that is
    stored but was computed for a Hamiltonian the session no longer
    has (set_hamiltonian(restart=False)) is solved again, warm-started
    from the session's previous state on v2/v3 and from a random one on
    "python", whose session drops its state when the terms change (2026-09
    audit, finding 2), since the session answers with the Hamiltonian it
    holds.

    This replaces the set_initial_wf(self.wf0) that every correlator used
    to open with, the trigger of the file-based backend's hand-off of the
    stored state to the C++ program. The pybind port kept the trigger and
    dropped the hand-off, so the line only reset computed_gs and the
    next get_gs() put the session's own solved state back over whatever
    set_gs() had set (2026-09-24b hole hunt, finding 11), or, after a
    set_wavefunction() had dropped the session's energy, ran a real
    sweep from it (finding 13)."""
    current = gs_is_current(self)
    if current and hamiltonian_on_session(self): return # the cache hit
    if current: self.computed_gs = False # the stored state is not for this H
    self.get_gs()
    send_hamiltonian(self)


def gs_energy_single(self,wf0=None,reconverge=None,maxde=None,maxdepth=5):
    """
    Return the ground state energy via the in-process session
    (mpscpp2/mpscpp3's chain_session.h Chain, or pyitensor's): the
    Hamiltonian, sweep parameters and wavefunction are passed as in-memory
    arguments to self._session.

    Where the answer comes from, in order of precedence:

    - wf0=, an explicit start: a detached copy (detached_copy()) is handed
      to the session and swept from, or taken unswept with reconverge=False;
    - a state the caller injected (mark_injected()): taken unswept, with
      e0 = <wf|H|wf>, after set_gs()/set_initial_wf(), and swept from after
      set_initial_wf_guess();
    - otherwise the session's own state: its cached energy when it has one
      under the current Hamiltonian and parameters (skip_dmrg_gs, which
      reconverge=True overrides), a sweep from it when it does not.
    """
    supplied = True
    if wf0 is not None:
        mode = "skip" if reconverge is False else "reconverge"
        start = wf0.copy()
    else:
        mode = pending_injection(self)
        start = self.wf0
        if mode is not None and len(self._gs_injected)>2:
            supplied = self._gs_injected[2]
    if mode=="skip":
        out = _take_injected_state(self,start,supplied=supplied)
    else:
        _session_parameters(self)
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
        # A Hamiltonian that is already an MPO (a StaticOperator, e.g. from
        # toMPO() and MPO algebra) is handed to the session directly -- it
        # has no symbolic term list to key a cache on, and building one
        # would defeat the point of having assembled it as an MPO. Identity
        # of the handle plus the solver parameters is the cache key instead.
        send_hamiltonian(self) # precondition: H on the session (see above)
        if mode=="reconverge":
            # after send_hamiltonian, which on "python" drops the session's
            # state when the terms changed (2026-09 audit, finding 2)
            self._session.set_wavefunction(detached_copy(start).cpp_handle)
            skip = False
        elif reconverge is not None: skip = not reconverge
        else: skip = self.skip_dmrg_gs
        out = self._session.gs_energy(skip_dmrg=skip)
        self.e0 = out # store ground state energy
        self.sites_from_file = True
        self.gs_from_file = True
        self.skip_dmrg_gs = True
        # the solve's own result, the session's state already: assigned,
        # not injected (mark_injected() is for the public setters only)
        self.wf0 = mps.MPS(MBO=self,cpp_handle=self._session.gs_wavefunction()).copy()
        self._gs_injected = None
        self._gs_supplied = False
    self.computed_gs = True # ground state has been computed
    self._gs_solver_key = solver_key(self) # ...under these parameters
    if maxde is not None: # enforce a maximum fluctuation in the energy
      # the variance of the state itself, not of its truncation to maxm
      # (2026-09-24c audit, finding 7), and per site: maxde is a
      # fluctuation per site, while gs_energy_fluctuation() is the total
      from . import vev as _vev
      de = np.sqrt(abs(_vev.energy_variance(self,self.hamiltonian,wf=self.wf0)))
      de = de/self.ns # normalize by the number of sites
      if de>maxde and maxdepth>0: # if a maximum energy fluctuation
          maxm,nsweeps = self.maxm,self.nsweeps
          noise = self.noise
          ramp = self.bond_ramp
          print("Energy fluctuation = ",de,maxm)
          self.maxm = maxm*2
          self.nsweeps = 2 # just two sweeps
          self.noise = 0.0
          # No bond ramp here: this is a deliberate, already-warm
          # refinement of a converged state at doubled maxm over just two
          # sweeps, so there is no cheap-early-sweep phase to win and
          # spending one of the two sweeps below the target bond dimension
          # would simply halve what the retry does.
          self.bond_ramp = False
          gs_energy_single(self,maxde=maxde,reconverge=True,
                  maxdepth=maxdepth-1) # execute again
          mpo_cap = max(self.maxm,self.mpomaxm) # what the refined solve built H with
          self.maxm = maxm
          self.nsweeps = nsweeps # restore
          self.noise = noise
          self.bond_ramp = ramp
          self.computed_gs = True # ground state has been computed
          self._gs_solver_key = solver_key(self) # ...under these parameters
          # The session holds the refined state and its energy, and H built
          # at the same MPO cap, so record it as sent under the restored
          # parameters: left keyed on the doubled maxm, the next correlator
          # saw a Hamiltonian that was not on the session and re-solved at
          # the original maxm, discarding the refinement (2026-09-24c audit,
          # finding 6). Re-sending instead would drop the session's energy
          # cache and sweep the refined state at the original maxm anyway.
          if mpo_cap==max(self.maxm,self.mpomaxm):
              self._session_ham_cache = (self._session,_send_key(self)[0])
          out = self.e0 # the refined energy, the one left on the chain
    return out # return energy


def gs_energy(self,**kwargs):
    if self.is_hermitian(self.hamiltonian): # put a check for Hermitian
        if self.itensor_version in (2,3,"python"): # C++ or pure-Python version
            return gs_energy_single(self,**kwargs)
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
            if pending_injection(self)=="skip" and not kwargs:
                # set_gs()/set_initial_wf() on a non-Hermitian chain: the
                # state as given, as on the Hermitian route; NH-DMRG takes
                # no start state, so set_initial_wf_guess() still re-solves
                return _take_injected_state(self,self.wf0)
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
    dmrg_generalized()), itensor_version=3 (mpscpp3/chain_session.h's
    Chain::gs_energy_generalized, a line-for-line port of the same
    self-consistent Lagrange-multiplier algorithm against ITensor v3's
    own dmrg()/Sweeps/sum()) and itensor_version="julia_live"
    (mpsjulialive/generalized.jl's get_gs_generalized, the same algorithm
    once more against ITensorMPS.jl's dmrg()/Sweeps/add()) -- see any of
    those docstrings for the derivation. mpscpp2 (itensor_version=2)
    doesn't have this session method yet. There is also no ED
    implementation of this method at all (unlike vev()/gs_energy()/...,
    which all honor self.mode="ED" for cross-validation) -- self.mode="ED"
    is rejected explicitly below rather than silently ignored.

    Non-Hermitian self.hamiltonian is dispatched to a separate solver,
    nhdmrg.py's nhdmrg_generalized() (the non-Hermitian, complex-lambda,
    biorthogonal-quotient generalization of NH-DMRG, mirroring how this
    function itself generalizes plain gs_energy()) -- implemented on the
    same "python"/3/"julia_live" set as the Hermitian path above (no
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
    if self.itensor_version not in (3,"python","julia_live"):
        raise NotImplementedError(
            "gs_energy_generalized is only implemented for "
            "itensor_version=3, 'python' or 'julia_live' so far -- call "
            "chain.setup_cpp(version=3), chain.setup_python() or "
            "chain.setup_julia() first")
    if self._session is None and self.itensor_version!="julia_live":
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
    if self.itensor_version=="python" and self.ns<2:
        # pyitensor's sweeps are two-site as well, so a one-site chain has
        # no update to make and the state never moves off its random start
        # -- but unlike plain gs_energy(), which mode.py routes to ED for
        # this size, there is no ED fallback here, and the outer
        # self-consistent iteration still returns a lambda (the Rayleigh
        # quotient of that untouched state), i.e. a silently wrong number
        # rather than an obvious failure. Confirmed directly: a 1-site
        # chain returned -0.3049 for an exact -0.5. Placed *before* the
        # non-Hermitian dispatch below, unlike the itensor_version==3
        # guard further down: NH-DMRG escapes v3's abort because it never
        # calls ITensor's dmrg(), but it does not escape this one, its own
        # sweep being two-site too (checked, same wrong-number symptom).
        raise RuntimeError(
            "gs_energy_generalized: pyitensor's two-site DMRG can't handle "
            "a chain this short (n=%d < 2 sites) -- use mode=\"ED\" for the "
            "plain ground state, or a longer chain"%self.ns)
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
    if self.itensor_version=="julia_live":
        # This backend has no self._session at all (its state lives in the
        # live Julia session instead), so the whole Lagrange-multiplier
        # iteration runs in one Julia call -- see
        # mpsjulialive/generalized.jl. Everything after it is the same
        # bookkeeping as the session backends below.
        from .mpsjulialive.generalized import gs_energy_generalized as gsg_jl
        lam,wf0 = gsg_jl(self,A,lam0=lam0)
        self.e0 = lam
        self.wf0 = wf0.copy() # the solve's own result: assigned, not injected
        self.computed_gs = True
        return lam
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
    # The solve's own result, which is the session's state already, so it
    # is assigned rather than injected (mark_injected() is for the public
    # setters only). It used to go through set_initial_wf(), which reset
    # computed_gs=False, and an earlier version of this function set
    # computed_gs=True *before* that call, so the next plain gs_energy()
    # re-ran an ordinary ground-state solve over the generalized state.
    self.wf0 = mps.MPS(MBO=self,cpp_handle=self._session.gs_wavefunction()).copy()
    self._gs_injected = None
    self._gs_supplied = False
    self.computed_gs = True
    self._gs_solver_key = solver_key(self) # see gs_is_current
    return lam




def get_gs_manifold(MBO,n=2,tol=1e-3,**kwargs):
    """Return the ground state manifold, i.e. all the states with the
    lowest energy"""
    (es,wfs) = MBO.get_excited_states(n=n,**kwargs)
    es = np.array(es) # defensive: the boolean mask below needs an ndarray
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
    """Set `wf` as the ground state of the chain.

    On a DMRG backend with a session the state is marked as injected
    (mark_injected()): the next ground-state read hands a copy of it to
    the session and gives it its own energy <wf|H|wf>, unswept, so that
    every consumer -- vev(), the dynamical correlators, whichever of the
    Python-side wf0 or the session's own state they read -- measures the
    state that was set. It used to reach the Python-side wf0 only, and the
    next correlator call put the session's solved state back over it
    (2026-09-24b hole hunt, finding 11)."""
    mode = wf.mode # get the mode
    if mode=="DMRG": # DMRG mode
        if not mark_injected(MBO,wf,reconverge=False):
            MBO.computed_gs = True # no session: the state is current as set
    elif mode=="ED": # ED mode
        MBO.get_ED_obj() # generate the ED object
        MBO.ED_obj.computed_gs = True # comptued GS
        MBO.ED_obj.wf0 = wf.v.copy() # copy the array
        # every ED submode measures this state from its own energy, and
        # submode="ED" reads it rather than the dex manifold (2026-09-24c
        # audit, findings 1 and 2); the next ED solve retires the mark
        MBO.ED_obj._injected_state = True

