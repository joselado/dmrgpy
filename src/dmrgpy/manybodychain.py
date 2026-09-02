from __future__ import print_function
import numpy as np
import os
from . import mps
from . import timedependent
from . import groundstate
from . import operatornames
from . import correlator
from . import densitymatrix
from . import dynamics
from . import funtk
from . import vev
from . import mpsalgebra
from . import entanglement
from . import entropy
from . import excited
from . import effectivehamiltonian
from . import multioperator
from .cppext import DEFAULT_ITENSOR_VERSION

dmrgpath = os.path.dirname(os.path.realpath(__file__)) # path to this package

one = np.matrix(np.identity(3))



class Coupling():
  def __init__(self,i,j,g):
    """Store a two-site coupling constant g between sites i and j"""
    self.i = i
    self.j = j
    self.g = g




class Many_Body_Chain():
  @property
  def maxm(self):
      """Maximum MPS bond dimension kept in the DMRG/TDVP truncations."""
      return self._maxm

  @maxm.setter
  def maxm(self,m):
      """Validated on assignment rather than several calls deep.

      maxm<=0 was checked nowhere: itensor_version=2/3 died inside
      ITensor with an uncatchable SIGABRT ("m > maxdim"), taking the
      user's whole process with them, and itensor_version="python"
      silently returned a bond-dimension-1 (product-state) energy for
      maxm=0 and something else again for maxm<0. Raising here names the
      parameter and the value, and does it at the line that set it."""
      try: mi = int(m)
      except (TypeError,ValueError):
          raise ValueError("maxm must be a positive integer, got "+repr(m))
      if mi<1:
          raise ValueError("maxm must be >= 1, got "+repr(m))
      self._maxm = mi

  def __init__(self,sites,itensor_version=DEFAULT_ITENSOR_VERSION,**kwargs):
      """Create a many-body chain over the given list of sites, and
      initialize its DMRG/ED backend (itensor_version selects which one)"""
      self.sites = sites # list of the sites
      self.Id = self.get_operator("Id",0)
#      self.path = id_generator() # random ID in dmrgpy_tmp
      self.ns = len(sites) # number of sites
      self.mode = None # no mode (use the input parameter)
      self.exchange = 0 # zero
      self.fields = 0 # zero
      self.pairing = 0
      self.hubbard = 0
      self.hopping = 0
      self.resorder = False # reorder the indexes
      self.resordered_indexes = None # reordered indexes
      self.fermionic = False
      self.sites_from_file = False
      self.excited_gram_schmidt = False # it does not seem very effective
      self.hamiltonian = None # Hamiltonian, as a multioperator
      self._dcex_excited_cache = None # cache for dcex.py's excited-state
          # search (submode="EX" dynamical correlator), invalidated in
          # restart()
      self.hubbard_matrix = np.zeros((self.ns,self.ns)) # empty matrix
      self.use_ampo_hamiltonian = False # use ampo Hamiltonian
      # additional arguments
      self.kpmmaxm = 50 # bond dimension in KPM
      self.maxm = 30 # bond dimension in wavefunctions
      self.mpomaxm = 5000 # bond dimension for operators
      self.nsweeps = 15 # number of sweeps
      # White density-matrix noise for the DMRG sweeps, applied at full
      # amplitude for the FIRST HALF of the sweeps and exactly zero for the
      # second (mpscpp3/chain_session.h: `sw <= nsweeps_/2 ? noise_ : 0.0`),
      # so it is an escape-from-a-local-minimum aid that anneals away rather
      # than a permanent perturbation.
      #
      # 1e-7 follows ITensor's own convention: its Sweeps constructor
      # defaults noise to 0, and its canonical DMRG samples lead with
      # 1E-7 and decay (ITensor/sample/dmrg.cc and mixedspin.cc use
      # `sweeps.noise() = 1E-7,1E-8,0.0`; the hubbard_2d samples use
      # 1e-6..1e-8).
      #
      # This default used to be 1e-1 -- six orders larger, and three to four
      # above anything in ITensor's own examples. Measured on a spinful
      # two-orbital open chain (16 spin-orbitals) against an exact
      # diagonalization, at matched maxm/nsweeps/penalties with noise as the
      # only variable: 1e-1 cost +1.29 meV of accuracy and 1.4x runtime at
      # maxm=40. The mechanism is the annealing schedule -- at 1e-1 the
      # second, noise-free half of the sweeps does not fully recover from
      # what the first half did.
      #
      # Among SMALL values the same model mildly prefers no noise at all
      # (maxm=60, error/time): 0 -> 0.043 meV/141s, 1e-8 -> 0.052 meV/170s,
      # 1e-6 -> 0.052 meV/174s, 1e-4 -> 0.082 meV/174s. So 1e-7 is not
      # chosen because it beats 0 here; it costs ~0.01 meV and ~20% runtime
      # against 0 on this model. It is kept because noise is an escape-from-
      # a-local-minimum mechanism and a default of exactly 0 removes it for
      # every caller, including the models that need it -- a trap that no
      # amount of extra sweeps or bond dimension escapes is a real failure
      # mode in this codebase (see docs/known_issue_idmrg_product_state_
      # collapse.md for the iDMRG analogue, which needed noise added, not
      # removed). A caller who has checked their model can still set 0, and
      # groundstate.py's own reconvergence retry already does exactly that.
      self.noise = 1e-7 # noise for dmrg (ITensor's own sample value)
      # Bond-dimension ramp for the ground state (mpscppN/chain_session.h's
      # Chain::make_sweeps_ramped(), pyitensor/chain.py's
      # _make_sweeps_ramped()): rather than running all nsweeps at the
      # full self.maxm, spend the first bond_ramp_fraction of the schedule
      # growing the sweep bond dimension geometrically from
      # bond_ramp_start up to self.maxm, and hold it there for the rest.
      # Since two-site DMRG costs ~O(m^3) per bond, those early sweeps --
      # which mostly just locate the right variational subspace and gain
      # little from a large m -- become much cheaper, and the expensive
      # full-maxm sweeps start from an already-good state. This is the
      # standard ITensor sweep-table idiom (see the example table at the
      # top of ITensor's own sweeps.h).
      #
      # The ramp spans a fixed *fraction* of the schedule rather than
      # growing as fast as it can (doubling per sweep, which is what was
      # tried first): growing fast minimizes the number of cheap sweeps
      # and so gives away nearly all of the speedup: on a 30-site
      # inhomogeneous Heisenberg-Hubbard chain at nsweeps=20, doubling
      # leaves only 4 of 20 sweeps below a maxm of 150. Filling half the
      # schedule instead measures 2.15x at maxm=90 and 1.7-2.0x at
      # maxm=60 on that model, for the same energy to ~1e-8 -- rerun it
      # with examples/groundstate/bond_dimension_ramp.
      # At the default 0.5 the second half of the schedule -- and hence
      # the returned energy -- always runs at the full self.maxm.
      self.bond_ramp = True # ramp the bond dimension across gs sweeps
      self.bond_ramp_start = 10 # bond dimension of the first sweep
      self.bond_ramp_fraction = 0.5 # fraction of the sweeps spent ramping
      # The noise term is most useful while the state is still being
      # built; it decays by this factor per ramping sweep and is off
      # entirely once the schedule reaches the full self.maxm, so the
      # final, converged sweeps are noise-free.
      self.bond_ramp_noise_decay = 0.1
      self.verbose = False # print the C++ DMRG per-sweep progress output
      self.kpmcutoff = 1e-12 # cutoff in KPM
      self.cutoff = 1e-12 # cutoff in ground state
      self.tevol_custom_exp = True # custom exponential function for Tevol
      self.tevol_method = "TDVP" # real-time MPS evolution method: "TDVP"
          # (default, itensor_version in (3,"python"), two-site TDVP via
          # mpscpp3/chain_session.h's or pyitensor's quench_tdvp()/
          # evolve_and_measure_tdvp()), "TDVP_GSE" (one-site TDVP with
          # Krylov global subspace expansion, arXiv:2005.06104 -- same
          # itensor_version support as "TDVP"; see tdvp_gse_* below),
          # "TEBD" (2nd-order-Trotter TEBD, pyitensor/tebd.py or
          # mpscpp3/tebd.h -- itensor_version in (3,"python"), and only
          # for a strictly nearest-neighbor Hamiltonian; gates are built
          # once from the bare bond Hamiltonians and reused every step,
          # so it's cheaper per step than TDVP whenever it applies),
          # "AUTO" (tries "TEBD" first and transparently falls back to
          # "TDVP" if self.hamiltonian turns out not to be strictly
          # nearest-neighbor -- same itensor_version support as "TEBD";
          # see timedependent.py's _tebd_or_tdvp() for the actual
          # try/fallback and why this isn't the default), or
          # "MPO" (the legacy 2nd-order Taylor-expanded evoloperator()
          # backup, the only option for itensor_version=2, since none of
          # TDVP/TDVP_GSE/TEBD exist there -- a v2-API port of TDVP was
          # attempted but reverted after a severe, unresolved performance
          # regression at n>~10 sites, see git history around
          # mpscpp2/TDVP/ if picking this up again). See
          # timedependent.py's evolution_dmrg_DC()/evolve_and_measure_dmrg()
          # for the actual dispatch.
      self.tdvp_gse_sweeps = 3 # "TDVP_GSE" only: number of leading TDVP
          # sweeps preceded by a global-subspace-expansion step; after
          # that, bond dimension is left to whatever GSE already grew it
          # to (one-site TDVP alone never grows it further) and only
          # plain one-site TDVP sweeps run. Mirrors mpscpp3/TDVP/sample/
          # run.cc's own "if(n<3)" -- GSE mainly matters early, while the
          # state's bond dimension is still small.
      self.tdvp_gse_krylov_order = 3 # "TDVP_GSE" only: dimension of the
          # Krylov subspace {phi,H*phi,...} built at each GSE step
          # (TDVP/README.md's "KrylovOrd").
      self.tdvp_gse_cutoff = 1e-8 # "TDVP_GSE" only: SVD cutoff used both
          # for each Krylov-vector MPO application and for the final
          # density-matrix truncation in global_subspace_expand().
      self.cvm_tol = 1e-5 # tolerance for CVM
      self.cvm_nit = 1e3 # iterations for CVM
      self.cvm_patience = 50 # CVM CG early stop: iterations without a
          # meaningful (>0.1% relative) best-residual improvement before
          # concluding the truncation-imposed residual floor is reached
          # (see cvm.py::cvm_correction_vector)
      self.cvm_blowup = 100.0 # CVM CG early stop: running residual this
          # many times above the best one means the truncated recurrence
          # is diverging past the floor
      self.cvm_solver = "cg" # CVM linear solver: "cg" (global conjugate
          # gradient over whole-MPS primitives, every backend, the
          # historical behavior) or "variational" (Jeckelmann's dynamical
          # DMRG: minimize <x|M|x>-2Re<b|x> by two-site sweeping, so the
          # truncation is part of the ansatz instead of an error injected
          # after every operation -- itensor_version="python" only, and
          # only for A=B^dagger; see pyitensor/ddmrg.py)
      self.cvm_nsweeps = 6 # sweeps for cvm_solver="variational"
      self.cvm_maxm = self.maxm # bond dimension for the CVM correction
          # vector, independent of the ground state's own maxm (mirrors
          # kpmmaxm): a correction vector solving [(H-omega)^2+eta^2]xc=b
          # can need a different entanglement structure than the GS.
      self.kpm_scale = 0.7 # scaling of the spectra for KPM
      self.kpm_accelerate = True # set to true
      self.kpm_n_scale = 3 # scaling factor for the number of polynomials
      self.gs_from_file = False # start from a random wavefunction
      self.excited_from_file = False # read excited states
      self.e0 = None # no ground state energy
      self.wf0 = None # no initial WF
      self.skip_dmrg_gs = False # skip the DMRG minimization
      self.computed_gs = False # computed the GS already
      self.fit_td = False # use fitting procedure in time evolution
      self.itensor_version = itensor_version # ITensor version
      self.has_ED_obj = False # ED object has been computed
      self.ED_obj = None # no ED object
      self.kpm_extrapolate = False # use extrapolation
      self.kpm_extrapolate_factor = 2.0 # factor for the extrapolation
      self.kpm_extrapolate_mode = "plain" # mode of the extrapolation
      # KPM energy truncation (Holzner et al., PRB 83, 195115 (2011),
      # Sec. III-B) -- itensor_version "python" or 3 only (mpscpp2/v2 has
      # no port), see pyitensor/kpm_energy_truncation.py +
      # mpscpp3/chain_session.h's kpm_dynamical_correlator_truncated()
      # (an independent method from kpm_dynamical_correlator, not a
      # branch inside it) and kpmdmrg.py's own dispatch between the two.
      # Off by default (existing kpm_scale behavior is unaffected);
      # enables kpm_scale to be pushed below its current safe floor
      # (~0.5 in this codebase's bandwidth-midpoint-centered rescaling
      # convention) for higher KPM spectral resolution, at the cost of
      # the extra truncation-sweep work per Chebyshev vector.
      self.kpm_energy_truncate = False # enable energy truncation
      self.kpm_truncate_dK = 10 # per-site Krylov subspace dimension
      self.kpm_truncate_nsweeps = 3 # number of truncation sweeps
      self.kpm_truncate_threshold = 1.0 # energy threshold (rescaled units)
      # in-process session (mpscpp2/mpscpp3's pybind11 bindings.cc for
      # itensor_version 2/3, or pyitensor.chain.Chain -- a plain Python
      # object, no pybind11 involved -- for itensor_version="python"),
      # built in sites.py::initialize() for itensor_version in
      # (2,3,"python"). Stays None if a C++ extension isn't compiled, in
      # which case mode.py falls back to ED -- there is no file-based DMRG
      # backend left to fall back to. ("python" has no such precondition:
      # it's always available, see cppext.py.)
      self._session = None
      # Conserved-sector (quantum-number) targeting, off by default -- see
      # set_conserved_sector(). None means "no sector": DMRG searches the
      # full Hilbert space, which is what every backend has always done.
      self.conserved_sector = None
      # whether self._session actually carries that sector. False both
      # when no sector is set and when only the ED backend can represent
      # it (a chain whose DMRG sites lack the quantum number, e.g. an Sz
      # sector on the interleaved spinful chain, which mode="ED" targets
      # and DMRG cannot) -- mode.py refuses DMRG in the latter case.
      self._sector_on_session = False
      self.initialize(**kwargs)
      # and initialize the sites
  def initialize(self,**kwargs):
      """Initialize the sites"""
      if self.mode=="ED": return # do nothing
      if self.itensor_version in [2,3,"python","julia"]:
          from .sites import initialize
          initialize(self)
      elif self.itensor_version=="julia_live":
          from .mpsjulialive.sites import initialize
          initialize(self)
      else: raise
  def __deepcopy__(self,memo):
      """Deepcopy that never tries to copy the in-process extension session
      (an opaque Chain object -- pybind11-backed for itensor_version 2/3,
      see mpscpp2/mpscpp3's bindings.cc, or a plain pyitensor.chain.Chain
      for itensor_version="python"): the C++ handle has no pickle/deepcopy
      support, and cloning a live session doesn't have a well-defined
      meaning anyway for either backend. clone() (used by
      bandwidth()/lowest_eigenvalue()) needs a working session of its own
      though, so a fresh Chain is built for the clone instead of copying
      the original one (a Chain only needs the site list to construct; any
      wavefunction/Hamiltonian state is set up again by whatever the clone
      is used for)."""
      from copy import deepcopy
      cls = self.__class__
      out = cls.__new__(cls)
      memo[id(self)] = out
      for k,v in self.__dict__.items():
          if k=="_session": out.__dict__[k] = None
          elif k=="_session_ham_cache": out.__dict__[k] = None
              # groundstate.py's Hamiltonian-send cache holds a reference
              # to the session itself: deepcopying it would try to copy
              # the pybind11 Chain (unsupported), and the clone gets a
              # fresh, empty session below anyway, so it must start with
              # no send-cache
          else: out.__dict__[k] = deepcopy(v,memo)
      if self._session is not None:
          from . import cppext
          out._session = cppext.get_backend(self.itensor_version).Chain(out.sites)
          # a fresh session starts sector-less; re-apply the clone's own
          # sector so it searches the same Hilbert space the original did
          if out.conserved_sector: out._apply_conserved_sector()
      return out
  def _reset_dmrg_state(self):
      """Invalidate any cached ground state computed under the previous
      backend/session. Without this, switching backend on an existing
      chain (setup_cpp/setup_python/setup_julia) would leave
      self.computed_gs=True and self.wf0 pointing at the OLD backend's
      wavefunction: gs_energy() would then silently return the stale
      energy without recomputing on the new backend, and vev()/get_gs()
      would hand the old backend's MPS handle to the new backend's
      session (a hard crash for the C++ backends, since an mpscpp2 MPS
      cpp_handle is not valid input to an mpscpp3 Chain, or vice versa).
      Same fields restart() resets, minus has_ED_obj (the ED backend is
      independent of itensor_version, so it doesn't need invalidating
      just because the DMRG backend changed)."""
      self.computed_gs = False
      self.gs_from_file = False
      self.skip_dmrg_gs = False
      self.wf0 = None
  def set_conserved_sector(self,**qns):
      """Confine every calculation on this chain to one quantum-number sector.

      Called with no arguments, sector targeting is switched off again and
      the chain goes back to searching the full Hilbert space (the default).
      Otherwise each keyword is a conserved quantity and its target value:

          fc.set_conserved_sector(Nf=6)        # exactly 6 particles
          sc.set_conserved_sector(Sz=0)        # total Sz=0
          hc.set_conserved_sector(Nf=8,Sz=0)   # Hubbard: both at once

      Which quantities a chain offers follows from its site types: `Nf`
      (particle number) for spinless and spinful fermions, `Sz` (total spin
      projection, in ITensor's integer 2*Sz units -- Sz=1 is one spin-1/2's
      worth) for every spin chain and for spinful fermions, `Nb` (boson
      number) for boson chains. Parafermion chains have no sector support.
      Asking for only some of what a site type offers is allowed and means
      exactly what it says: a Hubbard chain with Nf alone conserves particle
      number while leaving Sz free, so Sz-breaking terms remain legal.

      Two things follow from turning this on, both enforced rather than left
      to be discovered. Every operator built on the chain -- the Hamiltonian,
      but also anything passed to vev()/correlators -- must itself conserve
      the requested quantities, or a ValueError names the offending operator
      (`Sx` on an Sz-conserving chain, a bare `C` on an Nf-conserving one);
      a Hamiltonian written as Sx.Sx+Sy.Sy+Sz.Sz is fine, since that sum does
      conserve Sz even though its individual terms do not. And the sector
      invalidates the Hamiltonian/ground state already sent to the session,
      so the next gs_energy() re-sends and re-solves.

      Three solvers implement this: DMRG with itensor_version=3 (genuinely
      block-sparse QN tensors), DMRG with itensor_version="python" (a
      charge grading on dense storage plus a charge penalty on the
      variational solves, see pyitensor/sector.py), and mode="ED" (the
      charge is diagonal in the ED product basis, so a sector is a set of
      basis states and every assembled operator is restricted to the
      corresponding submatrix, see edtk/edchain.py). The API, the units
      (Sz in 2*Sz units) and the answers are the same on all three;
      itensor_version=2 and "julia_live" have no quantum numbers at all
      and still refuse rather than answering with the global ground state.

      Two things about the ED implementation differ in kind rather than in
      result. It builds the *full* Hilbert space and then restricts, so a
      sector buys a smaller eigenproblem and not a smaller construction.
      And an unreachable target (Nf=99 on an 8-site chain) surfaces at the
      first calculation rather than here, since the ED object is built
      lazily; the quantum-number *names* are checked immediately either way.
      """
      previous = self.conserved_sector # restored below if the request is rejected
      previous_on_session = self._sector_on_session
      if not qns: # switch sector targeting off
          self.conserved_sector = None
      else:
          for k,v in qns.items():
              if not isinstance(v,(int,np.integer)):
                  raise ValueError("set_conserved_sector: %s must be an integer "
                                   "(Sz is in 2*Sz units), got %s"%(k,repr(v)))
          self.conserved_sector = {k:int(v) for k,v in qns.items()}
      if qns and self.itensor_version not in (3,"python"):
          # DMRG on this backend cannot do it, but ED can -- as long as
          # this chain type defines the requested quantum numbers at all
          try: charges = self.get_sector_charge_operators()
          except NotImplementedError: charges = dict()
          unknown = sorted(k for k in qns if k not in charges)
          if len(charges)==0 or len(unknown)>0:
              self.conserved_sector = previous
              raise NotImplementedError(
                  "set_conserved_sector%s is implemented for "
                  "itensor_version=3, itensor_version=\"python\" and "
                  "mode=\"ED\" (this chain uses %s%s). Switch with "
                  "setup_cpp(version=3) or setup_python(), or run it "
                  "with mode=\"ED\"."
                  %("" if len(unknown)==0 else " with %s"%unknown,
                    repr(self.itensor_version),
                    "" if len(charges)>0 else
                    ", whose ED backend defines no conserved quantities"))
      try:
          self._apply_conserved_sector()
      except Exception: # unreachable target, wrong quantum number for these
          self.conserved_sector = previous # sites, ...: leave the chain as it was
          self._sector_on_session = previous_on_session
          # ...and put the session back on the sector it was already in,
          # so the rejected request costs nothing at all
          try: self._apply_conserved_sector()
          except Exception: self._sector_on_session = previous_on_session
          raise
      self.restart() # the Hamiltonian/ground state were built on the old sites
  def get_sector_charge_operators(self):
      """The conserved quantities this chain offers to
      set_conserved_sector(), as {name: MultiOperator} measuring each one
      over the whole chain -- `Nf` counts particles, `Sz` is in ITensor's
      integer 2*Sz units (so one spin-1/2's worth is 1), `Nb` counts
      bosons.

      This is what the ED backend targets a sector with: each operator is
      diagonal in the ED product basis, so its eigenvalue labels every
      basis state and the sector is the set of states carrying the
      requested labels (see edtk/edchain.py::set_conserved_sector). The
      DMRG backends do not use it -- they get their quantum numbers from
      the site set instead -- but the names and units are deliberately the
      same on both sides.

      Chain types that conserve nothing (parafermions) leave this raising,
      which is what makes set_conserved_sector refuse them.
      """
      raise NotImplementedError(
          "%s defines no conserved quantities for set_conserved_sector"
          %type(self).__name__)
  def _ed_only_quantum_numbers(self):
      """Quantum numbers this chain's ED backend can target but its DMRG
      site set cannot, so that set_conserved_sector knows not to try
      pushing them onto the session (and not to report the session's
      complaint as the answer).

      Empty for almost every chain -- the two sides deliberately agree on
      names and units. The exception is the interleaved spinful fermionic
      chain, whose DMRG representation is 2n spinless fermionic sites that
      know about Nf and nothing about spin.
      """
      return set()
  def _apply_sector_to_ed(self,edobj):
      """Push this chain's conserved sector down to a freshly built ED
      object. Called from every get_ED_obj(): restart() (which
      set_conserved_sector calls) drops the cached ED object, so the
      sector is applied exactly once, when the object is rebuilt."""
      sector = getattr(self,"conserved_sector",None)
      if not sector: return edobj
      edobj.set_conserved_sector(sector,self.get_sector_charge_operators())
      return edobj
  def _apply_conserved_sector(self):
      """Push self.conserved_sector down to the in-process session, which is
      where the site set actually gets rebuilt (chain_session.h's
      Chain::set_conserved_sector). hasattr-guarded so an extension compiled
      before this existed still works as long as no sector is requested."""
      self._sector_on_session = False
      ed_only_ok = self.mode=="ED" # ED targets the sector by itself
      if self.itensor_version not in (3,"python"):
          # this DMRG backend has no quantum numbers at all
          if self.conserved_sector and not ed_only_ok:
              raise NotImplementedError(
                  "this chain targets the conserved sector %s, which DMRG "
                  "only implements for itensor_version=3 and "
                  "itensor_version=\"python\" -- clear it with "
                  "set_conserved_sector(), or run it with mode=\"ED\""
                  %(self.conserved_sector,))
          return
      if self._session is None: return
      if not hasattr(self._session,"set_conserved_sector"):
          if self.conserved_sector and not ed_only_ok:
              raise NotImplementedError(
                  "this compiled mpscpp3 extension predates conserved-sector "
                  "support -- rebuild it with 'python install.py "
                  "--itensor-version=3'")
          return
      if self.conserved_sector and ed_only_ok \
              and len(set(self.conserved_sector) & self._ed_only_quantum_numbers())>0:
          # A chain run through ED may legitimately ask for a quantum
          # number its DMRG site set does not carry (Sz on the interleaved
          # spinful chain, whose DMRG sites are spinless fermions). Leave
          # the session sector-free in that case: the ED backend targets
          # the sector by itself, and mode.py refuses DMRG on this chain
          # rather than letting it answer with the global ground state.
          # Only the quantum numbers named above skip the session -- every
          # other request still goes through it, so an unreachable target
          # is still reported here rather than at the first calculation.
          return
      qns = sorted((self.conserved_sector or {}).items())
      self._session.set_conserved_sector(qns)
      self._sector_on_session = bool(self.conserved_sector)
  def set_sector_penalty(self,lam=None):
      """Strength of the charge penalty the pure-Python backend adds to its
      variational solves while a conserved sector is set
      (itensor_version="python" only -- itensor_version=3 confines the
      calculation with genuinely block-sparse tensors and has no such
      knob). `None` restores the default, derived from the Hamiltonian's
      own coefficient scale.

      The penalty is lambda*sum_k (Q_k-q_k)^2, identically zero on the
      target sector, so this changes no converged number -- only how
      strongly an excursion out of the sector is suppressed on the way
      there. It exists because dense storage lets a truncating SVD leak
      ~1e-16 of amplitude across charge sectors, which a variational sweep
      then amplifies; see pyitensor/sector.py. Setting it to 0 disables
      confinement entirely and is only useful for demonstrating that.
      """
      if self.itensor_version!="python":
          raise NotImplementedError(
              "set_sector_penalty applies to itensor_version=\"python\" only "
              "(this chain uses %s): the other backends with sector support "
              "confine the calculation structurally instead."
              %repr(self.itensor_version))
      self._session.set_sector_penalty(lam)
  def promote_to_dense(self):
      """Leave conserved-sector mode while *keeping* the state computed in it.

      `set_conserved_sector()` with no arguments also switches the sector
      off, but it throws the ground state away along with everything else
      built on the sector's site set. This keeps it, converted exactly to
      its dense equivalent, so that the operations a sector forbids become
      available on a state that was nevertheless obtained inside one:

          fc.set_conserved_sector(Nf=6)   # solve at exactly 6 particles
          fc.gs_energy()
          fc.promote_to_dense()           # ... then leave the sector
          wf = fc.applyoperator(fc.C[3],fc.wf0)   # legal only now

      While a sector is set, every operator built on the chain must
      conserve it -- a bare `C` changes Nf, `Sx` changes Sz -- and that is
      enforced rather than merely discouraged, because ITensor's AutoMPO
      aborts the whole process over a flux-violating term instead of
      raising something catchable. After promoting, the same wavefunction
      lives on ordinary dense indices and all of it is legal again:
      applying `C`/`Cdag`, a pairing or transverse-field quench, a
      dynamical correlator of a charge-changing operator.

      The conversion is exact and involves no re-solve: a QN-conserving MPS
      is the same state as its dense counterpart, only stored block-sparsely
      (ITensor's own `removeQNs`, see `Chain::promote_to_dense` in
      mpscpp3/chain_session.h). What it costs is the block-sparsity speedup,
      which is gone from here on -- so promote after the expensive
      sector-confined part of the calculation, not before it. On
      itensor_version="python" the state was dense all along (see
      pyitensor/sector.py) and the conversion only relabels its site
      indices, so nothing is lost there; what promotion buys is the same on
      both backends -- the charge-changing operations the sector forbids.

      The Hamiltonian and the band-edge/iDMRG/VUMPS caches are dropped (they
      were built on the QN indices) and re-sent on the next call that needs
      them; the ground-state energy and wavefunction are kept, so a bare
      `gs_energy()` afterwards returns the sector's energy rather than
      re-solving unconstrained. Call `restart()` if an unconstrained
      re-solve is what you want.

      Nothing happens if no sector is set. itensor_version=3 and
      itensor_version="python" only, like `set_conserved_sector` itself.
      """
      if not self.conserved_sector: return # nothing to promote
      if self.itensor_version not in (3,"python"):
          raise NotImplementedError(
              "promote_to_dense is implemented for itensor_version=3 and "
              "itensor_version=\"python\" only (this chain uses %s)."
              %repr(self.itensor_version))
      if self._session is None:
          # No live session to promote: the sector was only ever recorded
          # on the Python object (_apply_conserved_sector is a no-op
          # without one), so clearing it here is the whole job.
          self.conserved_sector = None
          return
      if not hasattr(self._session,"promote_to_dense"):
          raise NotImplementedError(
              "this compiled mpscpp3 extension predates conserved-sector "
              "promotion -- rebuild it with 'python install.py "
              "--itensor-version=3'")
      self._session.promote_to_dense()
      self.conserved_sector = None
      # The session dropped its Hamiltonian MPO along with the QN sites, so
      # the send-cache must not claim it is still there. (sector_key() would
      # already force a re-send on its own, since the sector changed -- this
      # keeps it true independently of that.)
      self._session_ham_cache = None
      # self.wf0 is a *separate* C++ MPS from the session's own (see
      # groundstate.gs_energy_single, which copies gs_wavefunction() out),
      # so promoting the session does not reach it.
      self.wf0 = self.promote_mps(self.wf0)
  def promote_mps(self,wf):
      """Convert one wavefunction from a conserved sector's QN-carrying site
      indices to this chain's dense ones, exactly -- the per-wavefunction
      counterpart of `promote_to_dense`, which only reaches the chain's own
      ground state. Any `MPS` obtained while the sector was set (an excited
      state, an `applyoperator` result, a state saved before promoting) needs
      this before it can be contracted against an operator built afterwards.
      Returns a new `MPS`; a no-op on one that is already dense."""
      if wf is None: return None
      if self.itensor_version not in (3,"python") \
              or getattr(wf,"cpp_handle",None) is None:
          return wf
      if not hasattr(self._session,"promote_mps"):
          # Returning wf unconverted would hand back a handle still on the
          # sector's QN indices, whose next contraction aborts the process
          # from inside ITensor rather than raising -- say so instead.
          raise NotImplementedError(
              "this compiled mpscpp3 extension predates conserved-sector "
              "promotion -- rebuild it with 'python install.py "
              "--itensor-version=3'")
      out = wf.copy()
      out.cpp_handle = self._session.promote_mps(wf.cpp_handle)
      return out
  def setup_julia(self):
      """Setup the Julia mode"""
      self.itensor_version = "julia_live"
      self._reset_dmrg_state()
      self.initialize()
  def setup_cpp(self,version=DEFAULT_ITENSOR_VERSION):
      """Setup the C++ mode (version 2 = ITensor v2, 3 = ITensor v3)"""
      self.itensor_version = version
      self._reset_dmrg_state()
      self.initialize()
  def setup_python(self):
      """Setup the pure-Python DMRG backend (pyitensor.chain.Chain, no
      compiler/pybind11 needed) -- see cppext.py."""
      self.itensor_version = "python"
      self._reset_dmrg_state()
      self.initialize()
  def get_mode(self,**kwargs):
      """Resolve the effective calculation mode ("DMRG" or "ED") for this
      chain, see mode.py"""
      from .mode import get_mode
      return get_mode(self,**kwargs)
  def copy(self):
      """Return a copy of the chain (alias for clone())"""
      return self.clone() # clone and create a new one
      from copy import deepcopy
      return deepcopy(self)
  def clone(self):
      """
      Clone the object
      """
      from copy import deepcopy
      name = "dmrgpy_clone_"+str(np.random.randint(10000))
      out = deepcopy(self) # full copy of the object
      out.path = "/tmp/"+name # new path (label only, see sites.py::initialize)
      return out # return new object
  def set_hamiltonian(self,MO,restart=True): 
      """Set the Hamiltonian.

      `MO` is normally a MultiOperator (a symbolic sum of products of named
      site operators). It may instead be an already-built MPO -- a
      StaticOperator, as returned by `toMPO()` and closed under the MPO
      algebra that class exposes (`*` for MPO products, `+` for compressed
      MPO sums, scalar scaling). That route exists for operators that are
      cheap as an MPO and expensive as a term list: a total-spin penalty
      g*S^2_total is O(n^2) symbolic terms but only three squares of
      extensive one-body operators as MPOs, so

          Sx = sum(self.Sx) ; Sy = sum(self.Sy) ; Sz = sum(self.Sz)
          S2 = self.toMPO(Sx)*self.toMPO(Sx) + self.toMPO(Sy)*self.toMPO(Sy) \
               + self.toMPO(Sz)*self.toMPO(Sz)
          self.set_hamiltonian(self.toMPO(H) + g*S2)

      stays compact at any system size. Only itensor_version=3 implements
      it (Chain::set_hamiltonian_mpo); other backends raise. Note it is a
      structural convenience rather than automatically a speedup --
      ITensor's toMPO already compresses a term list well (measured: 6.7x
      the terms cost 1.4x the sweep time on a spinful chain), so reach for
      it when the term list is awkward or quadratic, not on the assumption
      that term count drives the solve."""
      if restart: self.restart() # restar the calculation
      self.hamiltonian = MO
      self.use_ampo_hamiltonian = True # use ampo Hamiltonian
  def get_heff(self,**kwargs):
      """Return effective Hamiltonian"""
      return effectivehamiltonian.get_effective_hamiltonian_coefficients(self,
              **kwargs)
  def bandwidth(self,h,**kwargs):
      """Compute the bandwidth of an Hermitian operator"""
      mbc = self.clone() # clone the object
      mbc.set_hamiltonian(h) ; e0 = mbc.gs_energy(**kwargs)
      mbc.set_hamiltonian(-h) ; e1 = mbc.gs_energy(**kwargs)
      return -e0 -e1
  def lowest_eigenvalue(self,X,**kwargs):
      """Given an operator X, return its smallest eigenvalue"""
      mbc = self.clone() # clone the object
      mbc.set_hamiltonian(X) 
      return mbc.gs_energy(**kwargs)
  def restart(self):
      """Restart the calculation"""
      self.computed_gs = False
      self.gs_from_file = False
      self.has_ED_obj = False # restart ED obj
      self.skip_dmrg_gs = False
      self.wf0 = None # initial file for GS
      self._dcex_excited_cache = None # invalidate cached excited states
          # (dcex.py), tied to the ground state being replaced above
      self._sector_states_cache = None # same, for sectordc.py's cached
          # per-sector solves (which hold a whole clone chain alive, so
          # dropping them here also releases that session)
      self._session_ham_cache = None # invalidate groundstate.py's
          # Hamiltonian-send cache so the next gs_energy() re-sends the
          # Hamiltonian (clearing the session's energy/band-edge caches)
          # even when the terms are unchanged -- restart() promises a
          # genuinely cold recalculation
  def is_hermitian(self,H):
      """Check if an operator is Hermitian.

      Cached by H's own identity (not by value -- MultiOperator terms are
      never mutated in place, see multioperator.py's copy() docstring, so
      identity is a safe proxy for "unchanged"): gs_energy() calls this
      unconditionally on every invocation, so a parameter sweep or
      best_gs()'s repeated gs_energy() calls against the same
      self.hamiltonian would otherwise re-pay the (still nontrivial even
      after mpsalgebra.is_hermitian's own bond-dimension fix) witness
      check on every single call. Only the most recently checked H is
      cached -- callers here always check either self.hamiltonian
      (unchanged across such repeated calls) or a fresh clone's H
      (bandwidth(), lowest_eigenvalue(), each a different object with no
      repeat to benefit from anyway), so a one-entry cache captures the
      case that matters without the bookkeeping of an unbounded one."""
      cache = getattr(self,'_is_hermitian_cache',None)
      if cache is not None and cache[0] is H:
          return cache[1]
      from .mpsalgebra import is_hermitian
      out = is_hermitian(self,H)
      self._is_hermitian_cache = (H,out)
      return out
  def vev_MB(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.vev(self,MO,**kwargs)
  def get_gs_degeneracy(self,**kwargs):
      """Return the degeneracy of the ground state"""
      from . import degeneracy
      return degeneracy.gs_degeneracy(self,**kwargs)
  def vev(self,MO,mode="DMRG",**kwargs):
      """Compute the vacuum expectation value <GS|MO|GS> of an operator,
      dispatching between DMRG and ED"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          if self.itensor_version in (2,3,"python"): # C++ or pure-Python version
              return vev.vev(self,MO,**kwargs)
          elif self.itensor_version=="julia_live": # julia live version
              from .mpsjulialive.vev import vev as vevjl
              return vevjl(self,MO,**kwargs)
          else: raise
      elif mode=="ED":
          MOf = self.toMPO(MO,mode="ED") # fast operator
          return self.get_ED_obj().vev(MOf,**kwargs) # ED object
      else: raise
  def metts_vev(self,MO,T,**kwargs):
      """Finite-temperature <MO> via METTS sampling (White & Stoudenmire,
      arXiv:1002.1305) -- see vevtk/mettsvev.py. Implemented for
      itensor_version='python' and 3 (not 2, which has no TDVP). MO may
      also be a list/tuple of MultiOperators, measured together on one
      shared METTS sample chain (returns a list of (mean,stderr) pairs
      instead of a single pair) -- much cheaper than one metts_vev call
      per operator."""
      from .vevtk.mettsvev import metts_vev as _metts_vev
      return _metts_vev(self,MO,T,**kwargs)
  def metts_dynamical_correlator(self,name,T,**kwargs):
      """Real-time finite-temperature dynamical correlator
      C_AB(t) = <A(t)B>_T via dynamical METTS sampling (Wang, McClarty,
      Dankova, Honecker & Wietek, arXiv:2405.18484, Sec. II) -- see
      vevtk/mettsdynamicalcorrelator.py. Implemented for
      itensor_version='python' and 3 (not 2, which has no TDVP), same
      restriction as metts_vev. Validate against
      self.get_dynamical_correlator(mode="ED", submode="ED", T=T,
      name=name) (edtk/dynamics.py's dynamical_correlator_finite_T), the
      exact Lehmann-sum reference for small systems."""
      from .vevtk.mettsdynamicalcorrelator import metts_dynamical_correlator as _mdc
      return _mdc(self,name,T,**kwargs)
  def test_ED(self):
      """Test the ED object"""
      self.get_ED_obj().test()
  def toMPO(self,H,**kwargs):
      """Transport an operator into a matrix-product operator"""
      return mpsalgebra.toMPO(self,H,**kwargs)
  def excited_vev_MB(self,MO,**kwargs):
      """
      Compute a vacuum expectation value
      """
      return vev.excited_vev(self,MO,**kwargs)
  def excited_vev(self,MO,**kwargs):
      """Compute a vacuum expectation value using excited states"""
      return self.excited_vev_MB(MO,**kwargs)
  def generate_bilinear(self,fun,A,B):
      """Generic bilinear term"""
      fun = funtk.obj2fun(fun) # set function
      # Index by len(A)/len(B), not self.ns: for chain classes whose
      # per-site operator lists are indexed by logical location rather
      # than physical site (e.g. Spinful_Fermionic_Chain,
      # Mixed_Spin_Fermion_Chain), len(A)/len(B) < self.ns, and indexing
      # up to self.ns would raise IndexError. For every chain whose
      # lists are physical-site-indexed already, len(A)==len(B)==self.ns,
      # so this is unchanged there.
      h = multioperator.msum(fun(i,j)*A[i]*B[j]
              for i in range(len(A)) for j in range(len(B)))
      return 0.5*(h + h.get_dagger()) # Hermitian
  def update_hamiltonian(self):
      """Rebuild the Hamiltonian from the stored hopping/hubbard/pairing/
      exchange terms"""
      h = self.hopping + self.hubbard + self.pairing
      h = h + self.exchange
      self.set_hamiltonian(h)
  def get_dagger(self,m):
      """Return the Hermitian conjugate of an operator"""
      return m.get_dagger() # dummy method
  def overlap(self,wf1,wf2,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.overlap(self,wf1,wf2,**kwargs)
  def aMb(self,wf1,M,wf2,**kwargs):
      """Compute the overlap <a|M|b>"""
      return mpsalgebra.overlap_aMb(self,wf1,M,wf2,**kwargs)
  def operator_norm(self,op,**kwargs):
      """Estimate the norm of an operator"""
      return mpsalgebra.operator_norm(self,op,**kwargs)
  def is_zero_operator(self,op,**kwargs):
      """Check if this is the zero operator"""
      out = self.operator_norm(op,**kwargs)
      return out<1e-4
  def exponential(self,h,wf,**kwargs):
      """Compute the overlap"""
      return mpsalgebra.exponential(self,h,wf,**kwargs)
  def applyoperator(self,A,wf,**kwargs):
      """Apply an operator"""
      return mpsalgebra.applyoperator(self,A,wf,**kwargs)
  def applyinverse(self,A,wf,**kwargs):
      """Apply an operator"""
      return mpsalgebra.applyinverse(self,A,wf,**kwargs)
  def summps(self,wf1,wf2,**kwargs):
      """Apply an operator"""
      return mpsalgebra.summps(self,wf1,wf2,**kwargs)
  def scale_mps(self,x,wf):
      """Multiply an MPS by a number (see mpsalgebra.scale_mps)"""
      return mpsalgebra.scale_mps(self,wf,x)
  def trace(self,A,**kwargs):
      """Compute the trace of an operator"""
      return mpsalgebra.trace(self,A,**kwargs)
  def inverse_trace(self,A,**kwargs):
      """Compute the trace of the inverse of an operator"""
      return mpsalgebra.inverse_trace(self,A,**kwargs)
  def set_pairings_MB(self,fun):
      """Generic pairing term"""
      h = self.generate_bilinear(fun,self.C,self.C)
      self.pairing = h
      self.update_hamiltonian()
  def set_hubbard_MB(self,fun):
      """Hubbard term"""
      h = self.generate_bilinear(fun,self.N,self.N)
      self.hubbard = h # store
      self.update_hamiltonian()
  def set_hoppings_MB(self,fun):
      """Hubbard term"""
      h = self.generate_bilinear(fun,self.Cdag,self.C)
      self.hopping = h # store
      self.update_hamiltonian()
  def run(self,**kwargs):
      """Run the calculation, dispatching through mode.py"""
      from .mode import run
      return run(self,**kwargs)
  def get_bond_entropy(self,wf,i,j):
      """Return the entanglement entropy of two sites"""
      return entropy.bond_entropy(self,wf,i,j)
  def get_site_entropy(self,wf,b):
      """Return the entanglement entropy of a site"""
      return entropy.site_entropy(self,wf,b)
  def get_mutual_information(self,wf,i,j):
      """Return the mutual information"""
      return entropy.mutual_information(self,wf,i,j)
  def get_pair_entropy(self,wf,i,j):
      """Return the entanglement entropy of two sites with the
      rest of the system"""
      return entropy.pair_entropy(self,wf,i,j)
  def get_correlation_matrix(self,**kwargs):
      """Return the correlation matrix"""
      return entanglement.get_correlation_matrix(self,**kwargs)
  def get_correlation_eigenvalues(self,**kwargs):
      """Return the eigenvalues of the correlation matrix"""
      return entanglement.get_correlation_eigenvalues(self,**kwargs)
  def get_correlation_entropy(self,**kwargs):
      """Return the entanglement entropy computed from the correlation
      matrix"""
      return entanglement.get_correlation_entropy(self,**kwargs)
  def get_correlated_orbitals(self,**kwargs):
      """Return the correlated orbitals"""
      return entanglement.get_correlated_orbitals(self,**kwargs)
  def get_correlated_density(self,**kwargs):
      """Return the correlated density"""
      return entanglement.get_correlated_density(self,**kwargs)
  def get_dynamical_correlator_MB(self,**kwargs):
      """Return a dynamical correlator, computed with DMRG"""
      return dynamics.get_dynamical_correlator(self,**kwargs)
  def get_dynamical_correlator(self,mode="DMRG",name=None,i=0,j=0,**kwargs):
      """Return a dynamical correlator, dispatching between DMRG and ED"""
      mode = self.get_mode(mode=mode) # overwrite mode
      # Resolve name= here, once, for every submode and both solvers. The
      # documented string form ("ZZ", "cdc", ... together with i=/j=) used
      # to be understood only by the TD/CVM/TDZ submodes under mode="DMRG"
      # -- the ones that happened to call str2MO themselves -- while the
      # default submode="KPM", submode="EX" and every mode="ED" path
      # crashed on it several frames deep. This chain object is also the
      # only thing that *can* resolve it: an EDchain has no Sx/Sz/C
      # attributes of its own.
      if name is not None:
          name = operatornames.str2MO(self,name,i=i,j=j)
          kwargs["name"] = name
      if mode=="DMRG":
          return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED":
          return self.get_ED_obj().get_dynamical_correlator(**kwargs)
  def get_spectral_function(self,*args,**kwargs):
      """Single-particle spectral function A_ij(w) of a fermionic chain,
      assembled from the eigenstates of the N+1 and N-1 particle-number
      sectors (see sectordc.py). Frequencies are measured from the
      chemical potential by default. itensor_version=3 and
      itensor_version="python" only, like every sector-based method."""
      from . import sectordc
      return sectordc.spectral_function(self,*args,**kwargs)
  def get_spin_spectral_function(self,*args,**kwargs):
      """Dynamical spin structure factor S_ij(w), resolved into its Sz
      channels: S^{+-} and S^{-+} from the Sz+/-2 sectors, S^{zz} from the
      ground state's own sector (see sectordc.py). itensor_version=3 and
      itensor_version="python" only."""
      from . import sectordc
      return sectordc.spin_spectral_function(self,*args,**kwargs)
  def get_dynamical_correlator_moments(self,**kwargs):
      """Return the raw KPM Chebyshev moments of a dynamical correlator,
      together with the rescaling data needed to reconstruct a spectrum
      from them: (mus,emin,emax,scale,n,delta).

      This is the expensive half of submode="KPM"'s
      get_dynamical_correlator(); feeding the result to
      kpmdmrg.dynamical_correlator_from_moments() reconstructs a spectrum
      with any kernel ("jackson", "lorentz", "plain", "hodc") without
      repeating the DMRG work, which is what makes a kernel-vs-kernel
      comparison meaningful (same moments, no re-convergence in between).
      DMRG/KPM only -- there is no ED counterpart, since the ED path
      builds spectra by explicit summation rather than from moments."""
      from . import kpmdmrg
      return kpmdmrg.dynamical_correlator_moments(self,**kwargs)
  def get_distribution(self,mode="DMRG",**kwargs):
      """Return the distribution of an operator's spectrum"""
      # get_mode(), like every neighbouring entry point: without it the
      # DMRG branch ran even when mode.py had routed everything else to ED
      # (mode="ED", or itensor_version=3 on a <3-site chain), handing
      # session-only code an ED State and failing with the opaque
      # "'State' object has no attribute 'cpp_handle'"
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          from . import distribution
          return distribution.get_distribution(self,**kwargs)
     #     raise # not implemented
       #   return dynamics.get_dynamical_correlator(self,**kwargs)
      elif mode=="ED": 
          return self.get_ED_obj().get_distribution(**kwargs)
      else: raise
  def get_distribution_moments(self,mode="DMRG",**kwargs):
      """Return the moments of the distribution of an operator's spectrum"""
      mode = self.get_mode(mode=mode) # overwrite mode, see get_distribution
      if mode=="DMRG":
          from . import distribution
          return distribution.get_distribution_moments(self,**kwargs)
      elif mode=="ED":
          raise NotImplementedError(
              "get_distribution_moments has no ED implementation (the ED "
              "path builds spectra by explicit summation, not from KPM "
              "moments). Note mode.py routes to ED on its own when the "
              "requested C++ extension is unavailable, or for "
              "itensor_version=3 on a chain with fewer than 3 sites.")
      else: raise ValueError("Unrecognized mode "+repr(mode))
  def get_excited(self,mode="DMRG",**kwargs):
      """Return excitation energies.

      Both modes return one entry per eigenSTATE, so a degenerate level
      appears once per member and the two agree index by index (up to
      DMRG convergence). `get_gap()` reads es[1]-es[0], which is therefore
      ~0 -- correctly -- whenever the ground level is degenerate; ask for
      enough levels that the multiplet AND the state you want both fit
      (with a three-fold degenerate ground state, the first genuine
      excitation is n=4).

      `mode="ED"` above `algebra.maxsize` used to drop members of
      degenerate levels, which made ED look like it returned distinct
      levels and made DMRG look like it duplicated states; that was an
      ARPACK defect, fixed in algebra._deflated_lowest_hermitian. See
      docs/ed_vs_dmrg_degenerate_multiplets.md.

      To tell apart states whose energies are too close to separate,
      identify them by a quantum number rather than by position: on a
      spin-rotation invariant model evaluate <S^2> per state (2 for a
      triplet member, 0 for a singlet)."""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          return excited.get_excited(self,**kwargs) # return excitation energies
      elif mode=="ED": 
          return self.get_ED_obj().get_excited(**kwargs) # ED
  def get_full_matrix(self,name):
      """Return the full matrix of a named operator, via ED"""
      return self.get_ED_obj().get_operator(name) # get the full operator
  def get_excited_states(self,mode="DMRG",**kwargs):
      """Return the excited STATES themselves (not just their energies).

      Like `get_excited`, both modes return one entry per eigenstate, so a
      degenerate level appears once per member -- see its docstring and
      docs/ed_vs_dmrg_degenerate_multiplets.md."""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG":
          return excited.get_excited_states(self,**kwargs) # return es and waves
      elif mode=="ED": 
          return self.get_ED_obj().get_excited_states(**kwargs) # ED
  def get_gap(self,**kwargs):
    """Return the gap"""
    es = self.get_excited(n=2,**kwargs)
    return es[1] -es[0]
  def get_hamiltonian():
      """Return the Hamiltonian as a multioperator"""
      return self.hamiltonian
  def nhdmrg(self,**kwargs):
      """Non-Hermitian DMRG (itensor_version 2, 3 or "python"): return
      (energy,psil,psir), the eigenvalue with smallest real part of the
      Hamiltonian (or of an operator passed as H=...) together with the
      biorthogonal left/right eigenvector MPS. See nhdmrg.py"""
      from .nhdmrg import nhdmrg
      return nhdmrg(self,**kwargs)
  def gs_energy_fluctuation(self,**kwargs):
      """Compute the energy fluctuations"""
      h = self.get_hamiltonian()
      e = self.vev(h)
      e2 = self.vev(h,npow=2)
      return np.sqrt(np.abs(e2-e**2))
  def set_initial_wf_guess(self,wf):
      """Set the initial guess, and perform the DMRG GS calculation"""
      self.set_initial_wf(wf,reconverge=True)
  def set_initial_wf(self,wf,reconverge=False):
      """Use a certain wavefunction as initial guess"""
      self.computed_gs = False
      if wf is None:
        self.gs_from_file = False # use a wavefunction from a file
      else:
        self.gs_from_file = True # use a wavefunction from a file
        self.wf0 = wf.copy() # name of the wavefunction
        if reconverge: self.skip_dmrg_gs = False # reconverge the calculation
        else: self.skip_dmrg_gs = True # reconverge the calculation
  def set_gs(self,wf):
      """Set the ground state"""
      from .groundstate import set_gs 
      groundstate.set_gs(self,wf) # set this as ground state
  def get_gs(self,best=False,n=1,mode="DMRG",**kwargs):
      """Return the ground state"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": # DMRG mode
        if groundstate.gs_is_current(self): # stored, and still valid
#            self.wf0.write(name=self.wf0.name,path=self.path)
            return self.wf0
        if best: groundstate.best_gs(self,n=n,**kwargs) # best ground state
        else: self.gs_energy(**kwargs) # perform a ground state calculation
        return self.wf0 # return wavefunction
      elif mode=="ED": return self.get_ED_obj().get_gs(**kwargs)
  def get_gs_manifold(self,**kwargs):
      """Return the ground-state manifold"""
      return groundstate.get_gs_manifold(self,**kwargs)
  def gs_energy(self,mode="DMRG",**kwargs):
      """Return the ground state energy"""
      mode = self.get_mode(mode=mode) # overwrite mode
      if mode=="DMRG": 
          # not just computed_gs: a stored energy is only an answer to
          # this call if it was computed under the solver parameters in
          # force now (groundstate.gs_is_current)
          if groundstate.gs_is_current(self): return self.e0
          return groundstate.gs_energy(self,**kwargs)
      elif mode=="ED": return self.get_ED_obj().gs_energy() # ED object
      else: raise
  def gs_energy_generalized(self,A,**kwargs):
      """Smallest generalized eigenvalue lambda solving
      H|psi>=lambda*A|psi> (H the chain's own Hamiltonian, A a Hermitian
      positive-definite metric MultiOperator). See
      groundstate.gs_energy_generalized's docstring -- itensor_version
      "python" or 3 only, for now, and its own CAVEAT paragraph: afterward
      self.wf0 is an eigenstate of the shifted problem, not a plain
      eigenstate of self.hamiltonian, which other methods that read
      self.wf0 (get_excited_states, dynamical/KPM correlators, ...) have
      no way to detect."""
      return groundstate.gs_energy_generalized(self,A,**kwargs)
  def get_correlator(self,pairs=[],**kwargs):
      """Return a correlator, default one"""
      print("Method get_correlator is deprecated, use vev instead")
      from . import spinchain
      from . import fermionchain
      if type(self)==spinchain.Spin_Chain: 
          ops = [self.Sz[i]*self.Sz[j] for (i,j) in pairs]
      elif type(self)==fermionchain.Fermionic_Chain: 
          ops = [self.Cdag[i]*self.C[j] for (i,j) in pairs]
      else: raise
      return np.array([self.vev(o,**kwargs) for o in ops])
  # Many_Body_Chain.evolution() used to live here. It called
  # timedependent.evolution(), which does not exist and (per git history)
  # never did in this form, so every call raised AttributeError on every
  # backend and in every mode -- dead API, removed rather than
  # resurrected. The real entry points are evolution_ABA()/
  # evolve_and_measure() (timedependent.py) and
  # get_dynamical_correlator(submode="TD").
  def get_rdm(self,**kwargs):
      """
      Compute the reduced density matrix
      """
      return densitymatrix.reduced_dm(self,**kwargs) # return DM
  def get_kpm_scale(self):
      """
      Return an estimate of the bandwidth
      """
      return 3*self.ns # estimated bandwidth
  def random_state(self,mode="DMRG",orthogonal=None):
      """Generate a random MPS"""
      if self.mode is not None: mode = self.mode # redefine
      if mode in ["DMRG","MPS"]:
         if self.itensor_version in (2,3,"python"): # C++ or pure-Python version
             from . import mps
             if orthogonal is None: return mps.random_mps(self)
             else: return mps.orthogonal_random_mps(self,orthogonal)
         elif self.itensor_version=="julia_live": # Julia version
             from .mpsjulialive import mps
             return mps.random_mps(self)
      elif mode=="ED":
          return self.get_ED_obj().random_state()
      else: raise
  def random_mps(self,**kwargs):
      """Generate a random MPS (alias for random_state)"""
      return self.random_state(**kwargs)
  def get_operator(self,name,i=None):
      """Return a certain multioperator"""
      from . import multioperator
      return multioperator.obj2MO([[name,i]]) # return operator



#from fermionchain import Fermionic_Hamiltonian
#from spinchain import Spin_Hamiltonian






def setup_sweep(self,mode="default"):
  """Setup the sweep parameters"""
  sweep = dict() # dictionary
  sweep["cutoff"] = 1e-06
  if mode=="default": # default mode
    sweep["n"] = "20"
    sweep["maxm"] = "100" 
  elif mode=="fast": # default mode
    sweep["n"] = "3"
    sweep["maxm"] = "20" 
  elif mode=="accurate": # default mode
    sweep["n"] = "10"
    sweep["maxm"] = "50" 
  else: raise
  sweep["n"] = self.nsweeps
  sweep["maxm"] = self.maxm
  sweep["cutoff"] = self.cutoff
  self.sweep = sweep # initialize
#  write_sweeps(self) # write the sweeps




Many_Body_Hamiltonian = Many_Body_Chain # temporal workaround


import string
import random

def id_generator(size=20, chars=string.ascii_uppercase + string.digits):
   out = ''.join(random.choice(chars) for _ in range(size))
   return "/tmp/dmrgpy_tmp/"+out # temporal folder




