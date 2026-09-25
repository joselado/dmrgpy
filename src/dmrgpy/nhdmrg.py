"""
Non-Hermitian DMRG (NH-DMRG) driver, shared by every DMRG backend that
implements the session-level Chain.nhdmrg method: the compiled ITensor
v3 backend (mpscpp3/chain_session.h's Chain::nhdmrg -- the annotated
original), the compiled ITensor v2 backend (mpscpp2's back-port), and
the pure-Python backend (pyitensor/nhdmrg.py). The live Julia backend
(itensor_version="julia_live") has no self._session at all, so it plugs
into the same drivers one level up instead, through
mpsjulialive/nhdmrg.py's per-attempt functions -- everything else here
(the retry loop, the two-sided eigen-residual certificate) is generic
MultiOperator*MPS algebra and is shared unchanged.

The algorithm is a port of ITensorNHDMRG.jl
(https://github.com/tipfom/ITensorNHDMRG.jl) in its default
configuration: "onesided" local Arnoldi solves of A|x> = lambda|x> and
Adag|y> = conj(lambda)|y> on each two-site block, combined with the
"fidelity" truncation of Yamamoto et al., Phys. Rev. B 105, 205125 (both
MPS truncated with the same isometry from the hermitian average
rho = (rho_l + rho_r)/2 of the left/right reduced density matrices).

The optimization targets the eigenvalue with the smallest real part --
the same "ground state" convention used by the pre-existing MPS Arnoldi
route for non-Hermitian Hamiltonians (mpsalgebra's mode="GS"), which is
now only a fallback for backends without a session (julia_live keeps its
own path in groundstate.py).
"""

import math

from . import mps


def _unit_scale_up(cmax):
    """The power of two mpscpp2/3's mo_terms.h unit_scale_up() returns:
    it brings a largest |coefficient| below 1 into [1,2), and is exactly
    1.0 when that coefficient is 1 or more (or not a positive number), so
    a caller that branches on it runs its unscaled code byte for byte
    there. Multiplying by it and dividing back is exact."""
    if not (cmax>0.0) or cmax>=1.0: return 1.0
    _,e = math.frexp(cmax) # cmax in [2^(e-1),2^e), and e <= 0
    return math.ldexp(1.0,min(1-e,1000)) # cmax*up in [1,2)


def _max_abs_coef(terms):
    """Largest |coefficient| of a to_terms() list, read the way mo_terms.h's
    max_abs_coef() reads the list a session receives (raw, unmerged)."""
    return max([abs(c) for c,_ in terms]+[0.0])


def _is_identity_term(ops):
    """Whether a to_terms() factor list is the identity (every factor Id)."""
    return all(name=="Id" for (name,_) in ops)


def _residual_scale(terms):
    """(c, e_id) for the eigen-residual certificate of an operator given
    by its to_terms() list: e_id the summed coefficient of its pure
    identity terms, and c = min(1, largest |coefficient| among the other
    terms), or 1 when there are none. See nhdmrg()'s docstring."""
    e_id = 0.0
    cnon = 0.0
    for coef,ops in terms:
        if _is_identity_term(ops): e_id = e_id + coef
        else: cnon = max(cnon,abs(coef))
    return (min(1.0,cnon) if cnon>0.0 else 1.0),e_id


def _generalized_residual_scale(terms_h,terms_a):
    """(c_H, shift, s_A) for nhdmrg_generalized()'s certificate, which
    divides by c_H + |lambda - shift|*s_A; see its docstring."""
    c_h,e_id = _residual_scale(terms_h)
    s_a = _max_abs_coef(terms_a)
    if not (s_a>0.0): s_a = 1.0 # A=0 has no eigenproblem; any scale will do
    a_id = _residual_scale(terms_a)[1]
    if abs(a_id)>1e-12*s_a: return c_h,e_id/a_id,s_a
    # A with no identity part: lambda cannot absorb H's constant, which
    # is then just one more of H's coefficients
    cnon = max([abs(c) for c,o in terms_h if not _is_identity_term(o)]+[0.0])
    cmax = max(cnon,abs(e_id))
    return (min(1.0,cmax) if cmax>0.0 else 1.0),0.0,s_a


def nhdmrg(self,H=None,krylovdim=20,restarts=2,tol=1e-4,ntries=5):
    """Run non-Hermitian DMRG on a session-backend chain. Returns
    (energy,psil,psir) with energy the (complex) eigenvalue of smallest
    real part and psil/psir the biorthogonal left/right eigenvector MPS,
    normalized so that <psil|psir> = 1 (each tensor pair shares its site
    and link indices, so both behave as ordinary MPS individually).

    - H: operator to diagonalize (defaults to self.hamiltonian)
    - krylovdim/restarts: per-bond local Arnoldi effort; the outer DMRG
      sweeps (self.nsweeps) do the actual converging, so these stay small.
      itensor_version="julia_live" delegates the sweep to the real
      ITensorNHDMRG.jl package rather than to one of this codebase's own
      ports, where these two knobs don't map one-to-one -- krylovdim is
      forwarded, restarts is ignored (see mpsjulialive/nhdmrg.py's
      nhdmrg_attempt for why)
    - tol/ntries: eigen-residual certificate. The non-Hermitian "energy"
      is not a variational bound, so a (rare) stalled sweep can report a
      spurious value below the true spectrum with nothing else looking
      wrong; the only reliable convergence certificate is the pair of
      residuals ||H|psir> - E|psir>|| and ||Hdag|psil> - E*|psil>||. Both
      are checked: the right residual alone would accept a run whose
      anchored adjoint solve locked psil onto a *different* eigenstate,
      since <psil|H|psir>/<psil|psir> equals E identically whenever psir
      alone is an eigenvector. Each run starts from its own random MPS,
      so runs are re-drawn (up to ntries times) until the worse of the
      two relative residuals drops below tol, and the best run is
      returned regardless (converged runs sit at ~1e-14 while stalls sit
      at ~1e-1, so tol's exact value is uncritical). An attempt that
      fails outright (RuntimeError) is redrawn the same way; only when
      *every* attempt fails does this raise.
      "Relative" is to c + |E - e_id|, with e_id the summed coefficient of
      H's pure identity terms and c = min(1, the largest |coefficient| of
      the others): scale-covariant (H -> s*H scales both the residual and
      the denominator by s) and blind to a constant offset (which moves E
      and e_id together and leaves the residual vector alone), and exactly
      the old 1 + |E| for any H with no identity term and a largest
      coefficient of 1 or more, the case the ~1e-14/~1e-1 calibration
      above was measured on. The old 1 + |E| certified every state below
      s = tol/(2||h||) (2e-5 on a 6-site Heisenberg chain) and every
      unconverged run next to an offset of ~1e4 (2026-09-25b audit,
      finding 26). The warning prints this relative number.

    On itensor_version 2, 3 and "python" the session solves 2^k*H, the
    power of two mo_terms.h's unit_scale_up() takes from the largest
    |coefficient| of H's terms (exactly 1, i.e. the unscaled call, at a
    largest coefficient of 1 or more), and the energy is divided back
    exactly; psil/psir are the same eigenvectors either way. The local
    Arnoldi's thresholds (its 1e-13 breakdown, its 1e-10*(1+|lambda|)
    restart test) are absolute, calibrated at unit scale, and so were
    unconverged below s ~ 1e-13 on v3 (E0/s 0.24 to 2.7 off ED on a
    6-site chain) and in small units generally (finding 27); every other
    dmrg() already ran at that scale since the 2026-09-25 record, item 2.
    The session stores nothing of H: the energy handed back, e0 and every
    later reader (send_hamiltonian, NH-KPM) see H in the caller's units.
    The certificate's residuals are formed at the same 2^k on every
    backend (julia_live's solve itself is left in the caller's units).
    """
    if self.itensor_version not in (2,3,"python","julia_live"):
        raise NotImplementedError("nhdmrg requires itensor_version 2, 3, "
                "\"python\" or \"julia_live\" (got "
                +str(self.itensor_version)+"); use the Arnoldi route "
                "(get_excited_states) for the other backends")
    if self._session is None and self.itensor_version!="julia_live":
        # itensor_version 2/3 requested but the matching extension was
        # never compiled (sites.py's initialize() leaves _session None).
        # The two sibling entry points in this file and
        # groundstate.gs_energy_generalized all raise here; without this,
        # nhdmrg() alone died with an AttributeError on the
        # set_sweep_params call below. (julia_live legitimately has no
        # _session -- its state lives in the live Julia session.)
        raise RuntimeError(
            "nhdmrg needs a compiled ITensor extension for "
            "itensor_version="+str(self.itensor_version)+" but none is "
            "available for this chain -- run `python install.py "
            "--itensor-version=3`, or call chain.setup_python() / "
            "chain.setup_julia() to use a backend that needs no compiler")
    if H is None: H = self.hamiltonian
    Hd = H.get_dagger()
    terms = H.to_terms()
    rscale,e_id = _residual_scale(terms) # the certificate's, see docstring
    # H's unit scale, see the docstring (H's adjoint has the same
    # magnitudes, so one factor serves both); exactly 1.0 at a largest
    # coefficient of 1 or more
    up = _unit_scale_up(_max_abs_coef(terms))
    # the certificate's MPS algebra runs at that scale too, since it is
    # not free of units either: in the caller's units "python"'s own
    # H*psir is off in small units (||H psir||/s 2.662 and 3.195 against
    # 2.604 at s=1e-14 and 1e-16 on a 6-site chain, a residual floor of
    # 6e-9*s already at 1e-8), which read a pair whose energy is exact to
    # 5e-15 as a relative residual of 0.16 at s=1e-14 and 0.54 at 1e-16
    Hc,Hcd = (H,Hd) if up==1.0 else (up*H,up*Hd)
    if self.itensor_version=="julia_live":
        from .mpsjulialive.nhdmrg import nhdmrg_attempt
        attempt = lambda: nhdmrg_attempt(self,H,krylovdim=krylovdim,
                restarts=restarts)
    else:
        self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,
                self.noise)
        self._session.set_verbose(self.verbose)
        self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
        terms_dag = Hd.to_terms()
        # the solve at unit scale, see the docstring
        if up!=1.0:
            terms = [(c*up,o) for c,o in terms]
            terms_dag = [(c*up,o) for c,o in terms_dag]
        def attempt():
            energy,hl,hr = self._session.nhdmrg(terms,terms_dag,
                    int(krylovdim),int(restarts))
            if up!=1.0: energy = energy/up # exact, a power of two
            return (energy,mps.MPS(self,cpp_handle=hl).copy(),
                    mps.MPS(self,cpp_handle=hr).copy())
    best = None
    last_error = None
    for i in range(max(1,int(ntries))):
        # A failed attempt is treated as a bad random draw and redrawn,
        # the same way nhdmrg_generalized()'s own loop below does (see its
        # comment): an unlucky start can leave the two vectors
        # (near-)biorthogonally degenerate, which mpsjulialive/nhdmrg.jl's
        # nh_biorthogonal_pair rejects rather than dividing by ~0 --
        # exactly the class of failure ntries>1 exists for. Observed for
        # real on itensor_version="julia_live" with an asymmetric-hopping
        # chain (tests/test_nh_dmrg.py's nh_asymmetric_hopping_chain);
        # before this, one such draw aborted the whole call even though
        # the very next draw would have converged to ~1e-15.
        try:
            energy,psil,psir = attempt()
        except RuntimeError as e:
            last_error = e
            if self.verbose>0:
                print("nhdmrg attempt",i,"raised",repr(e),
                      "-- retrying with a fresh random start")
            continue
        # at unit scale: up*(H psir - E psir), divided by up*(c + |E-e_id|)
        eu = energy*up
        r = Hc*psir - eu*psir
        l = Hcd*psil - eu.conjugate()*psil
        resid = max(abs(r.dot(r))**0.5,abs(l.dot(l))**0.5)/(
                up*(rscale+abs(energy-e_id)))
        if best is None or resid<best[0]:
            best = (resid,energy,psil,psir)
        if resid<tol: break
        if self.verbose>0:
            print("nhdmrg attempt",i,"did not converge, relative residual",
                  resid)
    if best is None:
        # Carry the last attempt's own message and traceback through
        # (`from last_error`). Every backend raises RuntimeError for its
        # own reasons -- an ITensor Error() surfaces as one through
        # pybind11, pyitensor raises one for its guards -- so a
        # deterministic, non-retryable failure also lands here after
        # ntries identical tries. Stating only this driver's guess at the
        # cause would send the user tuning nsweeps/maxm/krylovdim for a
        # problem none of them affect.
        raise RuntimeError(
            "nhdmrg: every attempt (of "+str(ntries)+") failed. If the "
            "error below is a left/right pair that never became "
            "biorthogonalizable, raising nsweeps, maxm or krylovdim may "
            "help; otherwise it is not a convergence problem. Last "
            "attempt's error: "+repr(last_error)) from last_error
    resid,energy,psil,psir = best
    if resid>=tol:
        print("Warning: nhdmrg did not reach the residual tolerance "
              "after",ntries,"tries (best relative residual "+str(resid)+
              "); consider raising nsweeps, maxm or krylovdim")
    return energy,psil,psir


def nhdmrg_generalized(self,A,H=None,krylovdim=20,restarts=2,tol=1e-4,
        ntries=5,lam0=None):
    """Non-Hermitian generalized-eigenvalue NH-DMRG: solves
    H|psi_R>=lambda*A|psi_R> for a possibly non-Hermitian H (defaults to
    self.hamiltonian) and a Hermitian positive-definite metric operator A
    (a MultiOperator, same calling convention as vev()/
    gs_energy_generalized()). Returns (lambda,psil,psir), the complex
    generalized eigenvalue of smallest real part and the biorthogonal
    left/right eigenvector MPS (<psil|psir>=1) -- same return convention
    as nhdmrg() above, generalizing it exactly the way
    gs_energy_generalized() generalizes gs_energy() (see
    pyitensor/nhdmrg.py's nhdmrg_generalized() for the self-consistent
    Lagrange-multiplier algorithm, now with a complex lambda and
    biorthogonal expectation values).

    Implemented for itensor_version="python" (pyitensor/nhdmrg.py's
    nhdmrg_generalized()), itensor_version=3 (mpscpp3/chain_session.h's
    Chain::nhdmrg_generalized, a line-for-line port of the same algorithm
    against this file's own nhdmrg_one_sweep instead of the hand-rolled
    Python one) and itensor_version="julia_live"
    (mpsjulialive/generalized.jl's get_gs_generalized_nhdmrg, the same
    outer loop wrapped around real ITensorNHDMRG.jl sweeps) -- mpscpp2 has
    no analogous session method yet.

    - krylovdim/restarts: per-bond local Arnoldi effort (same as nhdmrg(),
      including julia_live's own caveat about restarts being ignored there)
    - tol/ntries: eigen-residual certificate, same rationale as nhdmrg():
      the non-Hermitian generalized "eigenvalue" is not a variational
      bound, so both residuals ||H|psi_R>-lambda*A|psi_R>|| and
      ||H^dagger|psi_L>-conj(lambda)*A|psi_L>|| are checked (the right
      residual alone would accept a run whose anchored adjoint solve
      locked psi_L onto a *different* eigenstate, since
      <psi_L|H|psi_R>/<psi_L|A|psi_R> equals lambda identically whenever
      psi_R alone is a genuine eigenvector). Each attempt starts from its
      own fresh random MPS (same rationale as nhdmrg()); the best of up
      to ntries attempts is returned regardless of whether tol was met.
      The residuals carry the units of H and lambda those of H over A, so
      they are divided by c_H + |lambda - e_id/a_id|*s_A: c_H and e_id as
      in nhdmrg() (min(1, H's largest non-identity |coefficient|), and
      H's summed identity coefficient), a_id A's summed identity
      coefficient and s_A its largest |coefficient|. That is covariant
      under H -> s*H and invariant under A -> t*A (lambda -> lambda/t)
      and, for A = a_id*Id, under an offset of H; it is exactly nhdmrg()'s
      certificate at A = Id, and exactly the old 1 + |lambda| for an H
      with no identity term and a largest coefficient of 1 or more
      against an A whose largest coefficient is 1 (A = 1 + 0.2*Sz0, say).
      s_A is not capped at 1 the way c_H is: lambda*A has to come out in
      H's units whatever A's own. An A with no identity part cannot
      absorb H's constant, which then just counts towards c_H.
    - lam0: starting lambda estimate passed through unchanged to every
      attempt (defaults to a data-driven guess seeded from each attempt's
      own fresh random state -- see pyitensor/nhdmrg.py's own default).

    On itensor_version 3 and "python" the session solves 2^k*H against
    the same A, the unit scale of nhdmrg() read from H's terms, with lam0
    carried into and lambda back out of those units exactly (see
    nhdmrg()'s docstring).
    """
    if self.itensor_version not in (3,"python","julia_live"):
        raise NotImplementedError(
            "nhdmrg_generalized is only implemented for "
            "itensor_version=3, 'python' or 'julia_live' so far -- call "
            "chain.setup_cpp(version=3), chain.setup_python() or "
            "chain.setup_julia() first")
    if self._session is None and self.itensor_version!="julia_live":
        # same "itensor_version==3 but no compiled extension" gap
        # gs_energy_generalized() guards against -- see its own comment.
        # (julia_live legitimately has no _session at all: its state lives
        # in the live Julia session instead.)
        raise RuntimeError(
            "nhdmrg_generalized needs a compiled ITensor v3 extension "
            "(itensor_version=3) but none is available for this chain -- "
            "run `python install.py --itensor-version=3`, or call "
            "chain.setup_python() to use the pure-Python backend instead")
    if H is None: H = self.hamiltonian
    if H is None:
        raise RuntimeError("nhdmrg_generalized called before set_hamiltonian")
    from . import multioperator
    A = multioperator.obj2MO(A)
    Hd = H.get_dagger()
    terms = H.to_terms()
    terms_a = A.to_terms()
    # the certificate's scale, see the docstring
    rscale,shift,a_scale = _generalized_residual_scale(terms,terms_a)
    # H (not A) at unit scale, see nhdmrg(); lambda scales with H, and the
    # certificate's algebra runs at that scale too
    up = _unit_scale_up(_max_abs_coef(terms))
    Hc,Hcd = (H,Hd) if up==1.0 else (up*H,up*Hd)
    if self.itensor_version=="julia_live":
        from .mpsjulialive.nhdmrg import nhdmrg_generalized_attempt
        attempt = lambda: nhdmrg_generalized_attempt(self,H,A,
                krylovdim=krylovdim,restarts=restarts,lam0=lam0)
    else:
        self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,
                self.noise)
        self._session.set_verbose(self.verbose)
        self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
        terms_dag = Hd.to_terms()
        if up!=1.0:
            terms = [(c*up,o) for c,o in terms]
            terms_dag = [(c*up,o) for c,o in terms_dag]
        if self.itensor_version=="python": # pyitensor accepts lam0=None directly
            # (a NaN stays NaN, i.e. unset, under the scaling)
            session_lam0 = None if lam0 is None else lam0*up
        else: # the compiled v3 binding takes a plain complex, NaN meaning "unset"
            session_lam0 = complex(float('nan'),0.0) if lam0 is None else lam0*up
        def attempt():
            lam,hl,hr = self._session.nhdmrg_generalized(terms,terms_dag,
                    terms_a,int(krylovdim),int(restarts),lam0=session_lam0)
            if up!=1.0: lam = lam/up # exact, a power of two
            return (lam,mps.MPS(self,cpp_handle=hl).copy(),
                    mps.MPS(self,cpp_handle=hr).copy())
    best = None
    last_error = None
    for i in range(max(1,int(ntries))):
        # Each attempt's fresh random start (see docstring) can, on rare
        # unlucky draws, drive the biorthogonal pair into the metric A's
        # near-null-space -- nhdmrg_generalized()'s own guard against that
        # (both backends) raises RuntimeError rather than returning a
        # meaningless lambda. Treated the same as an ordinary
        # resid>=tol failure below (redraw and retry) rather than letting
        # it abort every remaining attempt -- exactly the class of "bad
        # random draw" ntries>1 exists to route around in the first
        # place; found via code review.
        try:
            lam,psil,psir = attempt()
        except RuntimeError as e:
            last_error = e
            if self.verbose>0:
                print("nhdmrg_generalized attempt",i,"raised",repr(e),
                      "-- retrying with a fresh random start")
            continue
        lu = lam*up # the certificate at unit scale, as in nhdmrg()
        r = Hc*psir - lu*(A*psir)
        l = Hcd*psil - lu.conjugate()*(A*psil)
        resid = max(abs(r.dot(r))**0.5,abs(l.dot(l))**0.5)/(
                up*(rscale+abs(lam-shift)*a_scale))
        if best is None or resid<best[0]:
            best = (resid,lam,psil,psir)
        if resid<tol: break
        if self.verbose>0:
            print("nhdmrg_generalized attempt",i,"did not converge, "
                  "relative residual",resid)
    if best is None:
        # Same reasoning as nhdmrg()'s own all-attempts-failed message:
        # report the last attempt's actual error rather than asserting a
        # cause. The near-null-space guard (A not positive definite) is
        # only one of the ways an attempt can fail -- on julia_live the
        # left/right pair failing to biorthogonalize raises here too, and
        # has nothing to do with A -- so naming A unconditionally sent
        # users off to debug a perfectly good metric operator.
        raise RuntimeError(
            "nhdmrg_generalized: every attempt (of "+str(ntries)+") "
            "failed. If the error below is the near-null-space guard, A "
            "is likely not positive definite for this problem; if it is "
            "a left/right pair that never became biorthogonalizable, the "
            "metric is not the issue. Last attempt's error: "
            +repr(last_error)) from last_error
    resid,lam,psil,psir = best
    if resid>=tol:
        print("Warning: nhdmrg_generalized did not reach the residual "
              "tolerance after",ntries,"tries (best relative residual "+
              str(resid)+"); consider raising nsweeps, maxm or krylovdim")
    return lam,psil,psir


def gs_energy_generalized_nhdmrg(self,A,**kwargs):
    """gs_energy_generalized-style entry point for a non-Hermitian
    self.hamiltonian: run nhdmrg_generalized() and store the right
    eigenvector as the chain's ground state wavefunction, mirroring
    gs_energy_nhdmrg()'s own wf0/nh_left_wf handling (including its unit
    normalization of wf0 -- nhdmrg_generalized()'s own psir carries
    <psil|psir>=1 biorthogonal normalization instead).

    Returns lambda, kept as self.lam_generalized; self.e0 is the pair's own
    biorthogonal energy <psil|H|psir>/<psil|psir>, which NH-KPM measures
    from, as on the Hermitian route (groundstate.gs_energy_generalized's
    docstring, and 2026-09-25b hole hunt, finding 8): e0 used to be lambda,
    which is not an energy of the state."""
    lam,psil,psir = nhdmrg_generalized(self,A,**kwargs)
    from .groundstate import _state_energy
    self.e0 = _state_energy(self,psir,left=psil)
    self.lam_generalized = lam
    wf0 = psir.normalize()
    if wf0 is None: wf0 = psir.copy()
    self.nh_left_wf = psil.copy() # left eigenvector, for biorthogonal use
    # the solve's own result: assigned, not injected (groundstate.
    # mark_injected() is for the public setters, and would make the next
    # gs_energy() hand this state to the session as the caller's)
    self.wf0 = wf0.copy()
    self._nh_left_for = self.wf0 # the right state psil pairs with, see gs_energy_nhdmrg
    self._gs_injected = None
    self._gs_supplied = False
    self.computed_gs = True
    from .groundstate import solver_key
    self._gs_solver_key = solver_key(self) # see groundstate.gs_is_current
    _record_hamiltonian_sent(self)
    return lam


def _record_hamiltonian_sent(self):
    """Put H on the session through groundstate's send-cache after an NH
    solve, so that ground_state_on_session() finds the solve's state
    current. NH-DMRG hands the session its terms per call and never set
    them as the session's Hamiltonian, so the next correlator found H
    missing, reset computed_gs and re-solved plain NH-DMRG over the state:
    a generalized one was discarded even on a chain solved first
    (|<wg|wf0>|^2 = 0.6877, e0 from lambda to the plain -1.596396 on a
    4-site chain), and a plain one cost a second solve. No NH reader needs
    more than this: NH-KPM hands the session its own terms and states, and
    CVM_explicit reads self.wf0, so the session holds no NH state to lose.
    julia_live has no session and no send-cache."""
    if getattr(self,"_session",None) is None: return
    from .groundstate import send_hamiltonian
    send_hamiltonian(self)


def gs_energy_nhdmrg(self,**kwargs):
    """gs_energy-style entry point: run NH-DMRG and store the right
    eigenvector as the chain's ground state wavefunction (the state
    observables like vev() act on), mirroring what the Arnoldi
    non-Hermitian branch of groundstate.gs_energy stores -- including its
    unit normalization of wf0 (nhdmrg()'s own psir carries <psil|psir>=1
    biorthogonal normalization instead, so it is renormalized here; the
    biorthogonal pair as such stays available through nhdmrg()).

    Accepts (and ignores) unknown keyword arguments: gs_energy() forwards
    its **kwargs here for non-Hermitian Hamiltonians, and the previous
    Arnoldi route accepted a different set of solver knobs
    (maxit/delta/nkry_min/... -- see algebra/arnolditk.py's mpsarnoldi),
    so a strict signature would turn previously-working calls like
    get_gs_degeneracy(delta=...) into TypeErrors.

    H= is the exception, refused rather than ignored: the result is stored
    as the chain's state, and a pair solved for another operator is not
    the chain's ground state. It used to be accepted and stored with e0,
    the pair and the solver key, so the chain's own Hamiltonian's NH-KPM,
    finding H on the session, read the other operator's pair and energy
    (-1.836506+0.072051j against -1.596396 on a 4-site chain, 0.587 of the
    peak off; 2026-09-25b hole hunt, finding 9), and the Hermitian route
    already raised TypeError on the same keyword. nhdmrg(H=...) returns
    that pair without storing it."""
    if "H" in kwargs:
        raise TypeError("gs_energy(H=...): the ground state of an operator "
                "other than the chain's own Hamiltonian is not the chain's "
                "ground state, so it is not stored as one; nhdmrg(H=...) "
                "returns its (energy, psil, psir) without touching the "
                "chain, or set_hamiltonian(H) makes it the chain's own")
    known = ("krylovdim","restarts","tol","ntries")
    passed = {k:v for k,v in kwargs.items() if k in known}
    ignored = [k for k in kwargs if k not in known]
    if ignored and self.verbose>0:
        print("nhdmrg: ignoring keyword arguments",ignored,
              "(not NH-DMRG parameters)")
    e0,psil,psir = nhdmrg(self,**passed)
    self.computed_gs = True
    from .groundstate import solver_key
    self._gs_solver_key = solver_key(self) # see groundstate.gs_is_current
    self.e0 = e0
    # unit norm, matching the Arnoldi route's convention (MPS.normalize
    # returns a fresh normalized copy, or None for a degenerate state)
    wf0 = psir.normalize()
    self.wf0 = wf0 if wf0 is not None else psir.copy()
    self.nh_left_wf = psil.copy() # left eigenvector, for biorthogonal use
    # ...of this right state, the object itself, as groundstate.
    # mark_injected() does: anything that replaces self.wf0 (set_gs(), a
    # take of gs_energy(wf0=x, reconverge=False), restart()) unpairs the
    # two, and nonhermitian/kpm.py refuses an unpaired right state instead
    # of pairing it with a left state solved for another one
    self._nh_left_for = self.wf0
    self._gs_injected = None
    self._gs_supplied = False # a solve's state, not the caller's
    _record_hamiltonian_sent(self)
    return self.e0
