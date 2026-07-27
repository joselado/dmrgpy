"""
Non-Hermitian DMRG (NH-DMRG) driver, shared by every DMRG backend that
implements the session-level Chain.nhdmrg method: the compiled ITensor
v3 backend (mpscpp3/chain_session.h's Chain::nhdmrg -- the annotated
original), the compiled ITensor v2 backend (mpscpp2's back-port), and
the pure-Python backend (pyitensor/nhdmrg.py).

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

from . import mps


def nhdmrg(self,H=None,krylovdim=20,restarts=2,tol=1e-4,ntries=5):
    """Run non-Hermitian DMRG on a session-backend chain. Returns
    (energy,psil,psir) with energy the (complex) eigenvalue of smallest
    real part and psil/psir the biorthogonal left/right eigenvector MPS,
    normalized so that <psil|psir> = 1 (each tensor pair shares its site
    and link indices, so both behave as ordinary MPS individually).

    - H: operator to diagonalize (defaults to self.hamiltonian)
    - krylovdim/restarts: per-bond local Arnoldi effort; the outer DMRG
      sweeps (self.nsweeps) do the actual converging, so these stay small
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
      at ~1e-1, so tol's exact value is uncritical).
    """
    if self.itensor_version not in (2,3,"python"):
        raise NotImplementedError("nhdmrg requires itensor_version 2, 3 "
                "or \"python\" (got "+str(self.itensor_version)+"); use "
                "the Arnoldi route (get_excited_states) for the other "
                "backends")
    if H is None: H = self.hamiltonian
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,
            self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    Hd = H.get_dagger()
    terms = H.to_terms()
    terms_dag = Hd.to_terms()
    best = None
    for i in range(max(1,int(ntries))):
        energy,hl,hr = self._session.nhdmrg(terms,terms_dag,
                int(krylovdim),int(restarts))
        psil = mps.MPS(self,cpp_handle=hl).copy()
        psir = mps.MPS(self,cpp_handle=hr).copy()
        r = H*psir - energy*psir
        l = Hd*psil - energy.conjugate()*psil
        resid = max(abs(r.dot(r))**0.5,abs(l.dot(l))**0.5)/(1.0+abs(energy))
        if best is None or resid<best[0]:
            best = (resid,energy,psil,psir)
        if resid<tol: break
        if self.verbose>0:
            print("nhdmrg attempt",i,"did not converge, residual",resid)
    resid,energy,psil,psir = best
    if resid>=tol:
        print("Warning: nhdmrg did not reach the residual tolerance "
              "after",ntries,"tries (best residual "+str(resid)+"); "
              "consider raising nsweeps, maxm or krylovdim")
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

    Only implemented for itensor_version="python" (pyitensor) so far --
    mpscpp2/mpscpp3 have no analogous session method yet (CLAUDE.md's
    "Implement in pyitensor first" scoping, same as
    gs_energy_generalized()'s own original scoping before it gained a v3
    port).

    - krylovdim/restarts: per-bond local Arnoldi effort (same as nhdmrg())
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
    - lam0: starting lambda estimate passed through unchanged to every
      attempt (defaults to a data-driven guess seeded from each attempt's
      own fresh random state -- see pyitensor/nhdmrg.py's own default).
    """
    if self.itensor_version!="python":
        raise NotImplementedError(
            "nhdmrg_generalized is only implemented for "
            "itensor_version='python' so far -- call chain.setup_python() first")
    if H is None: H = self.hamiltonian
    if H is None:
        raise RuntimeError("nhdmrg_generalized called before set_hamiltonian")
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,
            self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    from . import multioperator
    A = multioperator.obj2MO(A)
    Hd = H.get_dagger()
    terms = H.to_terms()
    terms_dag = Hd.to_terms()
    terms_a = A.to_terms()
    best = None
    for i in range(max(1,int(ntries))):
        lam,hl,hr = self._session.nhdmrg_generalized(terms,terms_dag,terms_a,
                int(krylovdim),int(restarts),lam0=lam0)
        psil = mps.MPS(self,cpp_handle=hl).copy()
        psir = mps.MPS(self,cpp_handle=hr).copy()
        r = H*psir - lam*(A*psir)
        l = Hd*psil - lam.conjugate()*(A*psil)
        resid = max(abs(r.dot(r))**0.5,abs(l.dot(l))**0.5)/(1.0+abs(lam))
        if best is None or resid<best[0]:
            best = (resid,lam,psil,psir)
        if resid<tol: break
        if self.verbose>0:
            print("nhdmrg_generalized attempt",i,"did not converge, residual",resid)
    resid,lam,psil,psir = best
    if resid>=tol:
        print("Warning: nhdmrg_generalized did not reach the residual "
              "tolerance after",ntries,"tries (best residual "+str(resid)+
              "); consider raising nsweeps, maxm or krylovdim")
    return lam,psil,psir


def gs_energy_generalized_nhdmrg(self,A,**kwargs):
    """gs_energy_generalized-style entry point for a non-Hermitian
    self.hamiltonian: run nhdmrg_generalized() and store the right
    eigenvector as the chain's ground state wavefunction, mirroring
    gs_energy_nhdmrg()'s own wf0/nh_left_wf handling (including its unit
    normalization of wf0 -- nhdmrg_generalized()'s own psir carries
    <psil|psir>=1 biorthogonal normalization instead)."""
    lam,psil,psir = nhdmrg_generalized(self,A,**kwargs)
    self.e0 = lam
    wf0 = psir.normalize()
    self.wf0 = wf0 if wf0 is not None else psir.copy()
    self.nh_left_wf = psil.copy() # left eigenvector, for biorthogonal use
    self.set_initial_wf(self.wf0) # resets computed_gs=False, so...
    self.computed_gs = True        # ...this must come after, not before
                                    # (see gs_energy_generalized's own fix
                                    # for exactly this ordering bug)
    return lam


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
    get_gs_degeneracy(delta=...) into TypeErrors."""
    known = ("H","krylovdim","restarts","tol","ntries")
    passed = {k:v for k,v in kwargs.items() if k in known}
    ignored = [k for k in kwargs if k not in known]
    if ignored and self.verbose>0:
        print("nhdmrg: ignoring keyword arguments",ignored,
              "(not NH-DMRG parameters)")
    e0,psil,psir = nhdmrg(self,**passed)
    self.computed_gs = True
    self.e0 = e0
    # unit norm, matching the Arnoldi route's convention (MPS.normalize
    # returns a fresh normalized copy, or None for a degenerate state)
    wf0 = psir.normalize()
    self.wf0 = wf0 if wf0 is not None else psir.copy()
    self.nh_left_wf = psil.copy() # left eigenvector, for biorthogonal use
    return self.e0
