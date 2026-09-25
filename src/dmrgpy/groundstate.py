from . import mps
import numpy as np


def best_gs(sc,n=1,**kwargs):
    """Compute many ground states, and retain only the best one.

    Every other keyword goes to each of the n gs_energy() calls, which reads
    it (maxde=, reconverge=, ...) or raises TypeError on it, as
    get_gs(best=False) does: get_gs(best=True, **kwargs) used to hand them
    to this function, which took none, so any keyword next to best=True
    raised "best_gs() got an unexpected keyword argument" (2026-09-25
    record, "Left open"). wf0= is refused by name: the n solves each start
    from the backend's own state, which is what makes the lowest of them
    worth keeping, and n solves from one given state are one solve."""
    if kwargs.get("wf0") is not None:
        raise TypeError("get_gs(best=True) keeps the lowest of n solves, "
                "each from the backend's own starting state, so it takes "
                "no wf0=; call get_gs(wf0=...) without best=True to start "
                "from a given state")
    kwargs.pop("wf0",None)
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        sc.computed_gs = False # initialize
        # force a fresh session solve: with the Hamiltonian unchanged the
        # send-cache in gs_energy_single would otherwise return the
        # previous iteration's cached energy instead of re-running DMRG
        sc._session_ham_cache = None
        e0 = sc.gs_energy(**kwargs) # ground state energy
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

    A state with no recorded key and no injection mark was put there with
    nothing to hand it to, by a backend with no session (a restored
    julia_live snapshot, a set_gs() on a chain that fell back to ED), and
    is returned as-is."""
    if not self.computed_gs: return False
    if pending_injection(self) is not None: return False # see above
    key = getattr(self,"_gs_solver_key",None)
    if key is None: return True # no session to hand it to, see above
    return key==solver_key(self)


# What every reader of a chain with no Hamiltonian raises, in the words of
# Many_Body_Chain.get_hamiltonian(), which uses it too
NO_HAMILTONIAN = ("this chain has no Hamiltonian yet; call set_hamiltonian() "
                  "first")


def require_hamiltonian(self):
    """Raise ValueError naming set_hamiltonian() on a chain that has none.

    Readers used to fail wherever the missing Hamiltonian was first
    touched: "'NoneType' object has no attribute 'get_dagger'" in the
    Hermiticity probe on every DMRG backend, and on mode="ED" "No active
    exception to reraise" (Spin_Chain), "'NoneType' object has no
    attribute 'op'" (the fermion chains) or 'shape'/'T' (boson and
    parafermion chains) inside the ED builders (2026-09-25b hole hunt,
    finding 6). Many_Body_Chain.is_hermitian() calls it for every DMRG
    reader, which all probe the Hamiltonian first, and
    Many_Body_Chain._ed_reader() for the ED readers."""
    if getattr(self,"hamiltonian",None) is None:
        raise ValueError(NO_HAMILTONIAN)


def stored_answer_holds(self,kwargs):
    """True when gs_energy()/get_gs() may hand back the stored ground state
    (and its energy) for a call with these keywords, without reaching the
    solver. The one condition both entry points read, so that they cannot
    disagree about a call again.

    It used to be "current, and no wf0=", which returned the stored answer
    before any other keyword was read: a misspelled wf=x or reconverg=False
    was swallowed on a solved chain and raised TypeError on a fresh one,
    and maxde= came back unrefined, -3.3468165405 against -3.4061631313 on
    a fresh chain (2026-09-25b hole hunt, finding 5, and the 2026-09-25
    record's maxde= lead). The stored answer now comes back only for a call
    it already answers, and every other call goes to the solver, which
    honours the keyword or raises on it under its own signature, exactly as
    on a chain that is not current. Per route:

    - the session backends, Hermitian (gs_energy_single()) or not
      (gs_energy_nhdmrg()): besides wf0=None, only reconverge=False or
      None, maxde=None and maxdepth=, which does nothing without maxde=.
      maxde= goes to the solver, which refines; reconverge=True asks for a
      sweep from the state; any other keyword goes to the solver, which
      reads it (gs_energy_nhdmrg()'s H=, tol=, ...), refuses it, or, on the
      non-Hermitian route, ignores it after solving again, as it does on a
      chain that is not current (so get_gs_degeneracy(delta=...) on a
      current non-Hermitian chain re-solves; answering it from the stored
      state also answered gs_energy(H=H2) with the stored energy of H).
      reconverge=False and maxdepth= stay here rather than going on: the
      session has no energy of its own for a state taken as it was set,
      so its gs_energy(skip_dmrg=True) would sweep it, and NH-DMRG would
      solve again over it;
    - julia_live (_gs_energy_julia()): wf0=None only, since that solver
      takes neither maxde= nor maxdepth= and raises on reconverge= without
      a state, so any other keyword raises there as on a fresh chain."""
    if not gs_is_current(self): return False
    if kwargs.get("wf0") is not None: return False
    extra = set(kwargs)-{"wf0"}
    if not extra: return True # nothing else asked
    if self.itensor_version=="julia_live": return False # see above
    if extra-{"reconverge","maxde","maxdepth"}: return False # for the solver
    if kwargs.get("reconverge") is True: return False # a sweep, asked for
    return kwargs.get("maxde") is None # maxde=: a refinement, asked for


def ed_ground_state(self,wf0=None,reconverge=None,maxde=None,maxdepth=5):
    """The ED route of gs_energy() and get_gs(): read the keywords
    gs_energy_single() takes the way an exact solve honours them, and
    return the ED object, holding the state they ask for.

    Both entry points used to pass none of them on, on every route that
    resolves to ED, v3's own fallback below three sites included: so
    gs_energy(wf0=x, reconverge=False) left the ED ground state on the
    chain (|<get_gs()|x>|^2 = 0.0732, and vev(Sz0) -0.2287 against
    <x|Sz0|x> = 0.0390, on a 4-site chain), a misspelled keyword was
    swallowed, and get_gs(wf0=x) raised TypeError from EDchain.get_gs()
    (2026-09-25b hole hunt, finding 2).

    - wf0=x with reconverge=False makes x the chain's state through
      set_gs(), so that every later reader measures x, as on DMRG. x must
      be an ED state (random_state() on this route gives one); an MPS
      raises TypeError, since set_gs() would take its DMRG branch.
    - wf0=x as a start, reconverge None or True: a sweep from x ends on
      the ground state, and the exact solve is that state, so x itself is
      not read; a state set by hand before is dropped, as the sweep
      replaces it on DMRG.
    - maxde= and maxdepth= are met by the exact state, whose energy
      fluctuation is zero, so they are accepted and change nothing: the
      same call on DMRG must not raise because mode.py fell back to ED.
    - any other keyword raises TypeError, from this signature, as
      gs_energy_single()'s does on DMRG.

    The energy the caller then reads is the ED object's gs_energy(), the
    lowest eigenvalue, also after wf0=x with reconverge=False: that is what
    gs_energy(mode="ED") returns after set_gs(x), the 2026-09-24c record's
    open choice, which this route does not decide on its own."""
    require_hamiltonian(self)
    ed = self.get_ED_obj()
    if wf0 is None: return ed
    if getattr(wf0,"mode",None)!="ED":
        raise TypeError("wf0= on a chain that answers by ED must be an ED "
                "state, as random_state() gives on it; got %s. (mode.py "
                "routes this chain to ED: mode=\"ED\", itensor_version=3 "
                "below 3 sites, or no compiled extension.)"
                % type(wf0).__name__)
    if reconverge is False: set_gs(self,wf0) # x, as it is
    elif getattr(ed,"_injected_state",False):
        ed.computed_gs = False # the next read solves, replacing the set state
    return ed


def mark_injected(self,wf,reconverge=False,supplied=True):
    """Store `wf` (a unit-norm copy of it, see unit_copy()) as this chain's
    state, marked as injected by the caller rather than produced by a
    solve.

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
    know it exists. julia_live has no session either, but its solver takes
    a start state, so it is marked the same way and gs_energy()'s
    julia_live branch reads the mark (_gs_energy_julia()); before that,
    set_initial_wf() and set_initial_wf_guess() never reached the Julia
    solve, which started from its own random state (|<x|gs>|^2 anywhere
    from 0.0008 to 0.68 over runs, for a doublet member x of a 3-site
    chain), and set_gs() raised AttributeError on the Julia MPS. On a C++ backend that fell back to ED
    there is nothing to hand the state to, and the old behaviour stands:
    set_gs() stores it as current and set_initial_wf() leaves the next
    solve to that backend. Returns whether the state was marked.

    supplied=False is for the library's own re-marks of a state it
    computed (promote_to_dense): the state is taken unswept like any
    other, but does not count as the caller's, which is what
    submode="SECTOR" asks (see state_supplied()).

    The copy is normalized here, and not only when it is taken, because
    some readers use self.wf0 before any ground-state read does
    (evolve_and_measure() without wf=, for one)."""
    self.wf0 = unit_copy(wf)
    if (getattr(self,"_session",None) is None
            and self.itensor_version!="julia_live"):
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


def mark_lower_edge(self):
    """Record that self.e0 is the Hamiltonian's lower band edge: set by a
    plain Hermitian ground-state solve on julia_live (_gs_energy_julia()),
    and by nothing else.

    The julia_live KPM window takes its lower edge from e0 only when this
    holds, and from a solve of its own otherwise (mpsjulialive/dynamics.
    py). The mark holds the state and the energy it was made for and counts
    only while both are the chain's, the way mark_injected()'s does, so
    every other writer of e0 or wf0 -- a setter, gs_energy_generalized(),
    restart(), a take -- retires it without a line of its own. The window
    used to trust e0 unless the state was marked supplied, so a writer
    that knew nothing of the window was trusted by default: after
    gs_energy_generalized() e0 was lambda, and a metric that put lambda
    above E0 (A = 2*Id, A = 1.5+0.4*Sz0) raised "KPM moments diverging"
    where "python" and v3 returned (2026-09-25b hole hunt, finding 10)."""
    self._gs_lower_edge = (self.wf0,self.e0)


def e0_is_lower_edge(self):
    """True while self.e0 is marked as H's lower band edge, see
    mark_lower_edge()."""
    mark = getattr(self,"_gs_lower_edge",None)
    return (mark is not None and mark[0] is self.wf0 and mark[0] is not None
            and mark[1]==self.e0)


def solve_marks(self):
    """What a solve writes next to the chain's state and a helper that
    solves in place (mpsjulialive/dynamics.py's band-edge solves) must put
    back with it: the solver key and the lower-edge mark. See
    restore_solve_marks()."""
    return (getattr(self,"_gs_solver_key",None),
            getattr(self,"_gs_lower_edge",None))


def restore_solve_marks(self,marks):
    """Put back what solve_marks() read. The band-edge solves run at a
    clamped maxm/nsweeps, so the key they leave would make the restored
    state look stale to gs_is_current() at the chain's own parameters,
    and the next read would solve over it, a set state included."""
    self._gs_solver_key,self._gs_lower_edge = marks


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


def unit_copy(wf):
    """A copy of `wf` divided by its norm: the form in which a state becomes
    the chain's (mark_injected(), the wf0= of gs_energy(), and
    _take_injected_state(), which every unswept take goes through).

    A state set by hand is a ray. gs_energy() and the session vev() divided
    by <x|x>, while KPM, CVM, TD, TDZ, ROOTN, evolve_and_measure() without
    wf=, and every mode="ED" and julia_live reader took the vector as it
    was, so after set_gs(2*s) the KPM sum rule was 4 times the vev() of the
    same chain, gs_energy_fluctuation() of an exact eigenstate was 6|E0| on
    ED and julia_live, and get_gs() handed back a vector of norm 2
    (2026-09-25b hole hunt, finding 1). Normalizing once, where the state
    is taken, gives every reader the same unit vector, and get_gs() then
    returns it.

    Only a norm that double precision cannot divide out is refused: zero,
    or not finite (a NaN from an MPO applied to a state it annihilates).
    There is no floor above that. MPS.normalize()'s absolute 1e-8 is what
    it is not to copy: a state written in small units is still a ray
    (finding 22 of the same hunt), and nothing here knows the scale that
    produced the state, so a relative test has no reference either."""
    n2 = np.real(wf.dot(wf))
    if not (np.isfinite(n2) and n2>0.0):
        raise ValueError("the state given has norm^2 = %r, so it names no "
                "state: a set or start state is taken as the ray x/||x||, "
                "and x = 0 (for instance an operator applied to a state it "
                "annihilates) has no direction"%(n2,))
    return wf*(1.0/np.sqrt(n2))


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


def _on_session(self):
    """True when the chain's states live in a DMRG session (v2, v3,
    "python"). Not the bare test on self._session: setup_julia() leaves the
    previous backend's session on the chain, so a chain switched to
    julia_live still has one, and a Julia MPS handed to it fails."""
    return (self.itensor_version in (2,3,"python")
            and getattr(self,"_session",None) is not None)


def _state_energy(self,wf,left=None):
    """<wf|H|wf>/<wf|wf> for the chain's own Hamiltonian, the energy a state
    that is not a solve's is measured from; real for a Hermitian H. With
    left= it is the biorthogonal <left|H|wf>/<left|wf> of a non-Hermitian
    pair. The chain's own aMb/overlap on the session backends, the Julia
    MPS's own algebra on julia_live."""
    bra = wf if left is None else left
    if _on_session(self):
        e = self.aMb(bra,self.hamiltonian,wf)/self.overlap(bra,wf)
    else: e = bra.aMb(self.hamiltonian,wf)/bra.dot(wf) # julia_live
    if left is None and self.is_hermitian(self.hamiltonian):
        e = float(np.real(e))
    return e


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
    manifold the two coincide.

    julia_live has no session, so there is nothing to push: the state is
    the chain's, with the same energy, read with the Julia MPS's own
    <wf|H|wf>, and no solver key (see gs_is_current()): a set state has no
    solver parameters, and stays the chain's until something replaces it.

    The state is taken as the unit vector x/||x|| (unit_copy()), whichever
    route brought it here, the non-Hermitian gs_energy(wf0=x,
    reconverge=False) included; for a state normalized when it was marked
    this is a no-op, and the division of e0 by <x|x> below with it."""
    wf = unit_copy(wf)
    session = _on_session(self)
    if session:
        _session_parameters(self)
        send_hamiltonian(self) # precondition: H on the session
        self._session.set_wavefunction(detached_copy(wf).cpp_handle)
    e = _state_energy(self,wf)
    self.e0 = e
    self.wf0 = wf
    self._gs_injected = None
    self._gs_supplied = supplied
    self.computed_gs = True
    self.sites_from_file = True
    self.gs_from_file = True
    self.skip_dmrg_gs = True
    self._gs_solver_key = solver_key(self) if session else None
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

    - wf0=, an explicit start: a unit-norm copy (unit_copy()) is handed to
      the session and swept from, or taken unswept with reconverge=False;
    - a state the caller injected (mark_injected()): taken unswept, with
      e0 = <wf|H|wf>, after set_gs()/set_initial_wf(), and swept from after
      set_initial_wf_guess();
    - otherwise the session's own state: its cached energy when it has one
      under the current Hamiltonian and parameters (skip_dmrg_gs, which
      reconverge=True overrides), a sweep from it when it does not.

    maxde= is a tolerance on the energy fluctuation PER SITE,
    ||(H-<H>)|psi>||/ns, the intensive quantity, so gs_energy_fluctuation(),
    which returns the total ||(H-<H>)|psi>||, is to be divided by ns before
    it is compared with it. While the per-site value exceeds maxde, maxm is
    doubled and the state refined with two sweeps, at most maxdepth times;
    the refined energy is returned and the refined state kept.
    """
    supplied = True
    if wf0 is not None:
        mode = "skip" if reconverge is False else "reconverge"
        start = unit_copy(wf0) # the ray, as a set state is (finding 1)
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
          # labelled: the bare "Energy fluctuation =" read as the total
          # gs_energy_fluctuation() returns, ns times this number
          print("Energy fluctuation per site = ",de,maxm)
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


def _gs_energy_julia(self,ishermitian=True,wf0=None,reconverge=None):
    """gs_energy() on julia_live, with gs_energy_single()'s precedence: an
    explicit wf0= (swept from, or taken unswept with reconverge=False),
    then a state the caller injected (taken unswept, with e0 = <wf|H|wf>,
    after set_gs()/set_initial_wf(), swept from after
    set_initial_wf_guess()), then a solve from that backend's own random
    start. Both branches used to call get_gs_dmrg() directly, so the mark
    was never read and every setter was dropped (2026-09-24b hole hunt,
    finding 12, recorded there as open for julia_live). The caller's state
    is not moved: get_gs_dmrg() hands the solver a Julia-side copy, and
    ITensorMPS's dmrg() sweeps a copy of that again (dmrg.jl, psi =
    copy(psi0)). reconverge= means something only next to wf0=, since
    there is no session state to re-sweep here.

    A solve records its solver key, as gs_energy_single() does, so that a
    stored state is current only under the parameters it was solved with:
    none was recorded here, and gs_is_current() reads a missing key as
    "nothing to hand it to", so a convergence ramp over maxm on one chain
    returned the first energy every time (-3.194321 at maxm 2, 4 and 16 on
    an 8-site Heisenberg chain; 2026-09-25 open items, "Left open"). A
    Hermitian solve also marks its energy as H's lower band edge
    (mark_lower_edge()), which the julia_live KPM window reads."""
    if wf0 is None and reconverge is not None:
        raise TypeError("gs_energy(reconverge=...) on julia_live needs a "
                "state to take or sweep from, given as wf0=")
    supplied = True
    if wf0 is not None:
        mode = "skip" if reconverge is False else "reconverge"
        start = unit_copy(wf0) # the ray, as a set state is (finding 1)
    else:
        mode = pending_injection(self)
        start = self.wf0
        if mode is not None and len(self._gs_injected)>2:
            supplied = self._gs_injected[2]
    if mode=="skip": return _take_injected_state(self,start,supplied=supplied)
    from .mpsjulialive.groundstate import get_gs_dmrg
    e0,wf = get_gs_dmrg(self,ishermitian=ishermitian,
            wf0=start if mode=="reconverge" else None)
    # the solve's own result: assigned, not injected
    self.wf0 = wf
    self._gs_injected = None
    self._gs_supplied = False
    self._gs_solver_key = solver_key(self) # see gs_is_current
    if ishermitian: mark_lower_edge(self)
    return e0


def gs_energy(self,**kwargs):
    # a chain with no Hamiltonian is refused, naming set_hamiltonian(),
    # inside is_hermitian() (require_hamiltonian())
    if self.is_hermitian(self.hamiltonian): # put a check for Hermitian
        if self.itensor_version in (2,3,"python"): # C++ or pure-Python version
            return gs_energy_single(self,**kwargs)
        elif self.itensor_version=="julia_live":
            return _gs_energy_julia(self,**kwargs)
        else: raise
    else: # non-Hermitian excited states
        # Julia version has its own function
        if self.itensor_version=="julia_live":
            return _gs_energy_julia(self,ishermitian=False,**kwargs)
        elif self.itensor_version in (2,3,"python"): # real non-Hermitian DMRG
            if kwargs.get("wf0") is not None:
                # NH-DMRG takes no start state, so an explicit wf0= can
                # only be taken as it is; gs_energy_nhdmrg() used to drop
                # it as an unknown keyword and re-solve, so the call
                # returned the NH-DMRG energy with |<x|wf0>|^2 = 0.02 to
                # 0.11 for a random x on a 4-site chain
                if kwargs.get("reconverge") is not False or len(kwargs)>2:
                    raise TypeError("gs_energy(wf0=...) on a non-Hermitian "
                            "Hamiltonian: NH-DMRG takes no start state, so "
                            "a state to sweep from cannot be honoured; pass "
                            "reconverge=False (and nothing else) to take it "
                            "unswept (normalized), as set_gs() does")
                return _take_injected_state(self,kwargs["wf0"].copy())
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

    Returns lambda, which is also kept as self.lam_generalized. The chain's
    state afterwards is the generalized eigenvector wg, and its energy
    self.e0 is wg's own, <wg|H|wg>/<wg|wg> (the biorthogonal
    <psil|H|psir>/<psil|psir> on the non-Hermitian route), the energy a
    state set with set_gs(wg) gets. lambda is not an energy of wg, nor in
    general an energy at all (A = c*1 gives lambda = E0/c), and e0 used to
    be lambda: KPM, CVM, ROOTN and TDZ, which measure from e0, put wg's
    lines lambda - <wg|H|wg> away from where TD and EX, which measure from
    <wg|H|wg>, put them, so for A = 2*Id, where wg is the plain ground
    state, the default KPM had a line at negative frequency (2026-09-25b
    hole hunt, finding 8). A bare gs_energy() afterwards returns that
    <wg|H|wg> (it returned lambda before that fix).

    CAVEAT: self.wf0 afterwards is an eigenvector of the *shifted* problem
    H-lambda*A (or its biorthogonal NH counterpart), not a plain eigenstate
    of self.hamiltonian alone, and every method that reads the chain's
    state -- get_excited_states() (its overlap-penalty anchor), any
    dynamical/KPM correlator, NH-KPM (nonhermitian/kpm.py) -- measures wg,
    exactly as after set_gs(wg). That is what happens on every chain: the
    solve records H as sent, on both the Hermitian and the non-Hermitian
    route, so a correlator reads the generalized state whether or not the
    chain was solved before (it used to re-solve a plain ground state over
    it on a chain whose send-cache was empty, and on a non-Hermitian H on
    every chain, 2026-09-25 fixes). For the plain ground state instead,
    restart() and then gs_energy() before using any other method that
    reads self.wf0."""
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
        self.wf0 = wf0.copy() # the solve's own result: assigned, not injected
        self.e0 = _state_energy(self,self.wf0) # not lam, see above
        self.lam_generalized = lam
        self._gs_injected = None
        self._gs_supplied = False
        self.computed_gs = True
        self._gs_solver_key = solver_key(self) # see gs_is_current
        # no mark_lower_edge(): e0 is the state's energy, not H's lower
        # edge, so the julia_live KPM window solves for that edge itself
        return lam
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    # Through the send-cache, so the chain records that H is on the session.
    # A bare session.set_hamiltonian() here left the cache as it was, so on
    # a chain whose cache was empty the next correlator found H "not on the
    # session", re-solved a plain ground state over the generalized one and
    # measured that (|<wg|wf0>|^2 = 0.60 and <Sz0> -0.48 -> 0 on a 6-site
    # chain), while a chain solved beforehand read the generalized state, as
    # the CAVEAT above says. Skipping an identical re-send loses nothing:
    # both sessions drop their own energy cache after this solve anyway.
    send_hamiltonian(self)
    if self.itensor_version=="python": # pyitensor accepts lam0=None directly
        session_lam0 = lam0
    else: # the compiled v3 binding takes a plain float, NaN meaning "unset"
        session_lam0 = float('nan') if lam0 is None else lam0
    lam = self._session.gs_energy_generalized(A.to_terms(),lam0=session_lam0)
    # The solve's own result, which is the session's state already, so it
    # is assigned rather than injected (mark_injected() is for the public
    # setters only). It used to go through set_initial_wf(), which reset
    # computed_gs=False, and an earlier version of this function set
    # computed_gs=True *before* that call, so the next plain gs_energy()
    # re-ran an ordinary ground-state solve over the generalized state.
    self.wf0 = mps.MPS(MBO=self,cpp_handle=self._session.gs_wavefunction()).copy()
    self.e0 = _state_energy(self,self.wf0) # the state's own energy, see above
    self.lam_generalized = lam
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
    (2026-09-24b hole hunt, finding 11). julia_live is marked the same
    way, and read by _gs_energy_julia().

    The state is set as the ray it names, x/||x||, on every backend, ED
    included, so get_gs() hands back the unit vector (see unit_copy() for
    why, and for the zero state, which raises ValueError)."""
    mode = wf.mode # get the mode
    if mode=="DMRG": # DMRG mode
        if not mark_injected(MBO,wf,reconverge=False):
            MBO.computed_gs = True # no session: the state is current as set
    elif mode=="ED": # ED mode
        v = unit_copy(wf).v # the ray; every ED reader used to read x as given
        MBO.get_ED_obj() # generate the ED object
        MBO.ED_obj.computed_gs = True # comptued GS
        MBO.ED_obj.wf0 = v
        # every ED submode measures this state from its own energy, and
        # submode="ED" reads it rather than the dex manifold (2026-09-24c
        # audit, findings 1 and 2); the next ED solve retires the mark
        MBO.ED_obj._injected_state = True

