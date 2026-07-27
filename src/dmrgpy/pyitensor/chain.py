"""Chain: the session/handle facade, mirroring mpscpp3/chain_session.h's
Chain class method-by-method (and mpscpp3/bindings.cc's pybind11 surface,
minus the pybind11 plumbing itself -- there is no C++ boundary here, so
every method just returns native Python/numpy objects directly).

Args("Cutoff",...,"MaxDim",...) bags become explicit cutoff=/maxdim=
keyword arguments throughout, since Python already has real keyword
arguments and reimplementing a stringly-typed Args class would add
nothing. Every other design choice, quirk, and cross-reference below
points back to the specific chain_session.h comment it mirrors -- this
file is a close transcription, not a redesign, precisely so it inherits
that file's already-debugged behavior rather than re-deriving it.
"""

import itertools
from functools import lru_cache

import numpy as np

from .autompo import AutoMPO
from .dmrg import dmrg, dmrg_excited, dmrg_generalized
from .mpobuilder import to_mpo
from .mpsalgebra import applyMPO, inner, nmultMPO
from .mpsalgebra import sum as mps_sum
from .mpsalgebra import randomMPS, traceC
from .mpsalgebra import _fresh_link_copy
from .metts import metts_thermal_average as _metts_thermal_average
from .metts import metts_dynamical_correlator as _metts_dynamical_correlator
from .sites import SiteX
from .svd import svd
from .sweeps import Sweeps
from .tdvp import tdvp_step as _tdvp_step_fn
from .gse import global_subspace_expand as _global_subspace_expand_fn
from .tensor import commonIndex, contract_many, dag, delta, noPrime, prime, swapPrime

_BUILD_CUTOFF = 1e-14  # mo_terms.h's build_mpo() never exposes a cutoff knob at all


@lru_cache(maxsize=1)
def _four_pt_perm_table():
    """Sign and index-role data for all 4! ways a strictly increasing site
    quadruple (a,b,c,d) can be assigned to (i,j,k,l) of <Cdag_i C_j Cdag_k
    C_l> -- a direct port of the same table built in
    Chain::four_correlation_tensor_sweep (mpscpp3/chain_session.h; see its
    docstring for the full derivation, including the empirically-found
    extra sign for non-alternating Cdag/C type patterns). Returns a tuple
    of (mask, irank, jrank, krank, lrank, sign) rows."""
    table = []
    for p in itertools.permutations(range(4)):
        inv = [0, 0, 0, 0]
        for r in range(4):
            inv[p[r]] = r
        inversions = sum(1 for x in range(4) for y in range(x + 1, 4) if p[x] > p[y])
        mask = 0
        for r in range(4):
            if p[r] % 2 == 0:
                mask |= (1 << r)
        sign = 1 if inversions % 2 == 0 else -1
        if mask not in (5, 10):
            sign = -sign
        table.append((mask, inv[0], inv[1], inv[2], inv[3], sign))
    return tuple(table)


class Chain:
    def __init__(self, site_types):
        self.sites = SiteX(site_types)
        self.H = None
        self.have_H = False
        self.wf0 = None
        self._wf0_energy = None
        self._bandwidth_min = None
        self._bandwidth_max = None

        self.maxm = 30
        self.nsweeps = 15
        self.cutoff = 1e-12
        self.noise = 1e-1
        self.mpomaxm = 5000
        self.verbose = False

    def num_sites(self):
        return self.sites.length()

    def site_dim(self, site):
        return self.sites.dim(site)

    def random_mps(self):
        return self._default_mps()

    def set_sweep_params(self, maxm, nsweeps, cutoff, noise):
        self.maxm = maxm
        self.nsweeps = nsweeps
        self.cutoff = cutoff
        self.noise = noise

    def set_mpomaxm(self, mpomaxm):
        self.mpomaxm = mpomaxm

    def set_verbose(self, verbose):
        self.verbose = verbose

    def set_hamiltonian(self, terms):
        ampo = AutoMPO.from_terms(self.sites, terms)
        self.H = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        self.have_H = True
        self._wf0_energy = None  # any cached energy is now stale
        self._bandwidth_min = None  # ...and so is any cached bandwidth
        self._bandwidth_max = None

    def gs_energy(self, skip_dmrg=False):
        if not self.have_H:
            raise RuntimeError("Chain.gs_energy called before set_hamiltonian")
        if skip_dmrg and self.wf0 is not None and self._wf0_energy is not None:
            return self._wf0_energy
        if self.wf0 is None:
            self.wf0 = self._default_mps()
        sweeps = self._make_sweeps()
        energy = dmrg(self.wf0, self.H, sweeps, quiet=not self.verbose)
        self._wf0_energy = energy
        return energy

    def gs_wavefunction(self):
        if self.wf0 is None:
            raise RuntimeError("Chain.gs_wavefunction called before gs_energy")
        return self.wf0

    def set_wavefunction(self, wf):
        self.wf0 = wf
        self._wf0_energy = None  # energy no longer matches wf0

    def gs_energy_generalized(self, terms_a, lam0=None):
        """Smallest generalized eigenvalue lambda solving
        H|psi>=lambda*A|psi> for this chain's own Hamiltonian H (already
        built via set_hamiltonian()) and a metric operator A given by
        terms_a. A must be Hermitian positive definite for the
        self-consistent Lagrange-multiplier iteration in
        dmrg_generalized() (dmrg.py) to converge to the smallest
        generalized eigenvalue -- not checked here, same convention as
        set_hamiltonian() not checking H is Hermitian either. Mutates
        and stores self.wf0 like gs_energy(), but doesn't share its
        plain-<H> energy cache (skip_dmrg=True on a later gs_energy()
        call would otherwise return this generalized run's lambda,
        which is not <psi|H|psi>) -- invalidated below instead."""
        if not self.have_H:
            raise RuntimeError("Chain.gs_energy_generalized called before set_hamiltonian")
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_a), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        if self.wf0 is None:
            self.wf0 = self._default_mps()
        sweeps = self._make_sweeps()
        lam = dmrg_generalized(self.wf0, self.H, A, sweeps, quiet=not self.verbose, lam0=lam0)
        self._wf0_energy = None  # stale: no longer <wf0|H|wf0>, see docstring
        return lam

    def excited_states(self, n, scale_lagrange=1.0, do_gram_schmidt=False):
        if not self.have_H:
            raise RuntimeError("Chain.excited_states called before set_hamiltonian")
        if self.wf0 is None:
            self.gs_energy()  # ensure a ground state is available
        sweeps = self._make_sweeps()
        psi0 = self.wf0.copy()
        psi0.normalize()
        wfs = [psi0]
        weight = self._bandwidth() * scale_lagrange
        for _ in range(1, n):
            psi1 = self._default_mps()
            dmrg_excited(psi1, self.H, wfs, weight, sweeps, quiet=not self.verbose)
            psi1.normalize()
            wfs.append(psi1)
        if do_gram_schmidt:
            wfs = self._gram_schmidt(wfs)
        energies, fluctuations, wavefunctions = [], [], []
        for wf in wfs:
            fluctuations.append(self._energy_fluctuation(wf, self.H))
            energies.append(inner(wf, self.H, wf).real)
            wavefunctions.append(wf)
        return energies, fluctuations, wavefunctions

    def vev(self, terms, wf, npow=1):
        A = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf.copy()
        nrm = inner(psi, psi).real ** 0.5
        psi = psi * (1.0 / nrm)
        if npow == 1:
            return inner(psi, A, psi)
        psi1 = psi
        for _ in range(npow - 1):
            psi1 = self._apply_mpo(A, psi1)
        return inner(psi, A, psi1)

    def apply_operator(self, terms, wf):
        A = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        return self._apply_mpo(A, wf)

    def overlap(self, wf1, wf2):
        return inner(wf1, wf2)

    def overlap_aMb(self, wf1, terms, wf2):
        A = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        return inner(wf1, A, wf2)

    def sum_mps(self, wf1, wf2):
        return mps_sum(wf1, wf2, cutoff=self.cutoff, maxdim=self.maxm)

    def conjugate(self, wf):
        out = wf.copy()
        for i in range(1, out.length() + 1):
            out.set_A(i, dag(out.A(i)))
        return out

    def reduced_dm(self, wf, site):
        # psi /= innerC(psi,psi).real() -- divides by the norm *squared*,
        # not its square root; preserved exactly as chain_session.h has it
        # (its own comment traces this back to a note in mpscpp2's version).
        # commonIndex(psi.A(site),psi.A(site+1)) is called unconditionally
        # in the original too, i.e. it was never defended against site
        # being the last site either -- not special-cased here either.
        psi = wf.copy()
        nrm2 = inner(psi, psi).real
        psi = psi * (1.0 / nrm2)
        psi.position(site)
        ir = commonIndex(psi.A(site), psi.A(site + 1))
        s = self.sites.si(site)
        rho = psi.A(site) * dag(prime(psi.A(site), s, ir))
        for k in range(site + 1, psi.length() + 1):
            rho = rho * psi.A(k)
            rho = rho * dag(prime(psi.A(k), "Link"))
        out = []
        rho.visit(lambda z: out.append(z))
        dim = self.sites.dim(site)
        return np.array(out, dtype=complex).reshape(dim, dim)

    def exponential_apply(self, terms, wf, tau, nsteps):
        H = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        taui = tau / float(nsteps)
        expH = self._custom_exp(H, taui)
        psi1 = wf
        for _ in range(nsteps):
            psi1 = self._apply_mpo(expH, psi1)
        return psi1

    def build_operator(self, terms):
        return to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)

    def nhdmrg(self, terms_h, terms_hadj, krylovdim=20, restarts=2):
        """Non-Hermitian DMRG: optimize a biorthogonal left/right
        eigenpair of the non-Hermitian operator given by terms_h,
        targeting the eigenvalue with smallest real part; terms_hadj must
        be the adjoint operator's terms (MultiOperator.get_dagger() on
        the Python side). Port of mpscpp3/chain_session.h's Chain::nhdmrg
        (the annotated original) -- see nhdmrg.py in this package.
        Returns (energy, psil, psir)."""
        from .nhdmrg import nhdmrg as _nhdmrg
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h),
                   cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        HA = to_mpo(AutoMPO.from_terms(self.sites, terms_hadj),
                    cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        # fresh random start every run (never wf0): stalled runs are
        # detected by the caller's eigen-residual check and re-drawn
        # (see dmrgpy's nhdmrg.py retry loop)
        psi0 = self._default_mps()
        sweeps = self._make_sweeps()
        return _nhdmrg(H, HA, psi0, sweeps, krylovdim=krylovdim,
                       restarts=restarts, quiet=not self.verbose)

    def nhdmrg_generalized(self, terms_h, terms_hadj, terms_a,
            krylovdim=20, restarts=2, lam0=None):
        """Non-Hermitian generalized-eigenvalue NH-DMRG: solves
        H|psi_R>=lambda*A|psi_R> for a possibly non-Hermitian operator
        (terms_h, with terms_hadj its adjoint -- same
        MultiOperator.get_dagger() convention as nhdmrg() above) and a
        Hermitian positive-definite metric operator A (terms_a). See
        nhdmrg.py's nhdmrg_generalized() for the self-consistent
        Lagrange-multiplier algorithm (the non-Hermitian counterpart of
        gs_energy_generalized()'s own, generalized to a complex lambda
        and biorthogonal expectation values). Returns
        (lambda, psil, psir) with <psil|psir> = 1, same convention as
        nhdmrg()."""
        from .nhdmrg import nhdmrg_generalized as _nhdmrg_generalized
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h),
                   cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        HA = to_mpo(AutoMPO.from_terms(self.sites, terms_hadj),
                    cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_a),
                   cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        # fresh random start every run (never wf0), same rationale as
        # nhdmrg() above: the non-Hermitian "energy" is not a variational
        # bound, so a stalled run can only be detected/cured by the
        # caller's eigen-residual check and a redraw.
        psi0 = self._default_mps()
        sweeps = self._make_sweeps()
        return _nhdmrg_generalized(H, HA, A, psi0, sweeps,
                krylovdim=krylovdim, restarts=restarts,
                quiet=not self.verbose, lam0=lam0)

    def apply_pure_operator(self, A, wf):
        return self._apply_mpo(A, wf)

    def tdvp_step(self, H, wf, dt, num_center=2):
        """One TDVP step of size dt, given an already-built MPO H (e.g.
        from build_operator()) and a raw MPS handle wf -- lets a caller
        (tdz.py's TDZ submode) drive the evolution one variable-sized step
        at a time, unlike quench_tdvp/evolve_and_measure_tdvp above,
        which loop internally over a fixed number of equal, real dt steps
        and thus can't be reused for a per-step-varying complex time
        increment. dt may be any complex number: tdvp.py's tdvp_step()
        forward/backward coefficients are already generic to complex dt
        (the backward half-sweep's coefficient is just the negative of
        the forward one, not its conjugate -- the real-time-only case
        just happens to make those coincide), so real time (dt purely
        imaginary by this module's own -i*dt convention, matching
        mpscpp3's chain_session.h) and complex time (TDZ) share this same
        code path unchanged. num_center=2 (default) is two-site TDVP;
        num_center=1 is one-site TDVP, which doesn't grow bond dimension
        on its own -- pair with global_subspace_expand() below. NOTE: this
        method's own (H, wf, dt) argument order (matching mpscpp3's
        Chain::tdvp_step(H, psi, dt, ...)) is the REVERSE of the
        module-level tdvp.py::tdvp_step(psi, H, dt, ...) it calls below --
        any new call site added here should double check which of the two
        orderings it means to match."""
        return _tdvp_step_fn(wf.copy(), H, dt, cutoff=self.cutoff,
                maxdim=self.maxm, niter=50, num_center=num_center)

    def global_subspace_expand(self, H, phi, krylov_order, cutoff, maxdim=0):
        """Krylov-subspace global subspace expansion (arXiv:2005.06104),
        growing phi's bond dimension using H so one-site TDVP can keep up
        with two-site TDVP's own SVD-driven growth. maxdim=0 (matching
        the mpscpp3 binding's own sentinel) means "uncapped" for the
        per-Krylov-vector applyMPO steps; the *enlarged* bond dimension
        itself is always hard-capped at self.maxm (matching
        Chain::global_subspace_expand()'s own "MaxDim",maxm_ on the
        v3/mpscpp3 side -- see gse.py's own comment for why this is
        needed regardless of maxdim)."""
        return _global_subspace_expand_fn(H, phi, krylov_order, cutoff,
                maxdim=(maxdim if maxdim > 0 else None),
                bond_maxdim=self.maxm)

    def evolve_taylor_step(self, H, wf, z):
        """Applies one Taylor-expanded exp(z*H) step (_evoloperator()
        above, z may be complex) to an already-built MPO H and MPS wf --
        the MPO-Taylor (non-TDVP) analogue of tdvp_step() above, used by
        tdz.py's "TDZ" submode as a cross-check / non-TDVP alternative
        here, matching mpscpp2/mpscpp3's own evolve_taylor_step()
        bindings (mpscpp2 has no TDVP, so this is its only route to
        TDZ)."""
        expH = self._evoloperator(H, z)
        return self._apply_mpo(expH, wf)

    def multiply_operators(self, A, B):
        return self._mult_mpo(A, B)

    def sum_operators(self, A, B):
        return self._sum_mpo(A, B)

    def scale_operator(self, A, z):
        return A * z

    def trace_operator(self, A):
        return traceC(A)

    def hermitian_operator(self, A):
        out = A.copy()
        for j in range(1, out.length() + 1):
            out.set_A(j, dag(swapPrime(out.A(j), 0, 1, "Site")))
        return out

    def overlap_aMb_operator(self, wf1, A, wf2):
        return inner(wf1, A, wf2)

    def bond_entropy(self, wf, b):
        psi = wf.copy()
        psi.position(b)
        twosite = psi.A(b) * psi.A(b + 1)
        left_link = commonIndex(psi.A(b), psi.A(b - 1)) if b > 1 else None
        s_b = self.sites.si(b)
        left_inds = ([left_link] if left_link else []) + [s_b]
        _, _, _, spectrum = svd(twosite, left_inds, cutoff=0.0, maxdim=None)
        SvN = 0.0
        for p in spectrum.eigs():
            if p > 1e-12:
                SvN += -p * np.log(p)
        return SvN

    def quench(self, terms_h, terms_i, terms_j, nt, dt, fit_td=True):
        if self.wf0 is None:
            self.gs_energy()
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        H = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        EGS = inner(self.wf0, H, self.wf0).real / inner(self.wf0, self.wf0).real
        ampo_h.add(-EGS, "Id", 1)
        expH = self._evoloperator(to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm), dt)
        A1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi1 = self._apply_mpo(A1, self.wf0)
        psi2 = self._apply_mpo(A2, self.wf0)
        norm0 = np.sqrt(inner(psi1, psi1))
        correlator = []
        for _ in range(nt):
            psi1 = self._apply_mpo(expH, psi1, x0=psi1) if fit_td else self._apply_mpo(expH, psi1)
            psi1.normalize()
            psi1 = psi1 * norm0
            correlator.append(inner(psi2, psi1))
        return correlator, psi1

    def evolve_and_measure(self, terms_h, terms_op, wf, nt, dt, fit_td=True):
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        expH = self._evoloperator(to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm), dt)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf
        correlator = []
        for _ in range(nt):
            psi = self._apply_mpo(expH, psi, x0=psi) if fit_td else self._apply_mpo(expH, psi)
            correlator.append(inner(psi, A, psi))
        return correlator, psi

    def quench_tdvp(self, terms_h, terms_i, terms_j, nt, dt):
        if self.wf0 is None:
            self.gs_energy()
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        H = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        EGS = inner(self.wf0, H, self.wf0).real / inner(self.wf0, self.wf0).real
        ampo_h.add(-EGS, "Id", 1)
        Hshift = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi1 = self._apply_mpo(A1, self.wf0)
        psi2 = self._apply_mpo(A2, self.wf0)
        norm0 = np.sqrt(inner(psi1, psi1))
        correlator = []
        for _ in range(nt):
            psi1 = _tdvp_step_fn(psi1, Hshift, dt, cutoff=self.cutoff, maxdim=self.maxm, niter=50)
            psi1.normalize()
            psi1 = psi1 * norm0
            correlator.append(inner(psi2, psi1))
        return correlator, psi1

    def evolve_and_measure_tdvp(self, terms_h, terms_op, wf, nt, dt):
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf
        correlator = []
        for _ in range(nt):
            psi = _tdvp_step_fn(psi, H, dt, cutoff=self.cutoff, maxdim=self.maxm, niter=50)
            correlator.append(inner(psi, A, psi))
        return correlator, psi

    def quench_tdvp_gse(self, terms_h, terms_i, terms_j, nt, dt, gse_sweeps,
            krylov_order, gse_cutoff):
        """GSE counterpart of quench_tdvp() above: identical setup/
        measurement, but each per-step evolution is one-site TDVP
        (num_center=1) preceded by a global_subspace_expand() call for
        the first gse_sweeps steps -- mirrors
        Chain::quench_tdvp_gse()/mpscpp3/chain_session.h."""
        if self.wf0 is None:
            self.gs_energy()
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        H = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        EGS = inner(self.wf0, H, self.wf0).real / inner(self.wf0, self.wf0).real
        ampo_h.add(-EGS, "Id", 1)
        Hshift = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi1 = self._apply_mpo(A1, self.wf0)
        psi2 = self._apply_mpo(A2, self.wf0)
        norm0 = np.sqrt(inner(psi1, psi1))
        correlator = []
        for it in range(nt):
            if it < gse_sweeps:
                psi1 = self.global_subspace_expand(Hshift, psi1, krylov_order, gse_cutoff)
            psi1 = _tdvp_step_fn(psi1, Hshift, dt, cutoff=self.cutoff, maxdim=self.maxm,
                    niter=50, num_center=1)
            psi1.normalize()
            psi1 = psi1 * norm0
            correlator.append(inner(psi2, psi1))
        return correlator, psi1

    def evolve_and_measure_tdvp_gse(self, terms_h, terms_op, wf, nt, dt, gse_sweeps,
            krylov_order, gse_cutoff):
        """GSE counterpart of evolve_and_measure_tdvp() above -- see
        quench_tdvp_gse()'s docstring."""
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf
        correlator = []
        for it in range(nt):
            if it < gse_sweeps:
                psi = self.global_subspace_expand(H, psi, krylov_order, gse_cutoff)
            psi = _tdvp_step_fn(psi, H, dt, cutoff=self.cutoff, maxdim=self.maxm,
                    niter=50, num_center=1)
            correlator.append(inner(psi, A, psi))
        return correlator, psi

    def cvm_dynamical_correlator(self, terms_i, terms_j, omega, eta, energy, tol, max_it):
        if self.wf0 is None:
            raise RuntimeError("Chain.cvm_dynamical_correlator called before gs_energy")
        S1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        S2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        z = complex(omega + energy, eta)
        ampo = AutoMPO(self.sites)
        ampo.add(z, "Id", 1)
        zId = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = mps_sum(zId, self.H * (-1.0), cutoff=self.cutoff, maxdim=self.maxm)
        b = self._apply_mpo(S2, self.wf0)
        x = self._bicstab(A, b, tol, max_it)
        G = inner(self.wf0, S1, x)
        return -G.imag / np.pi

    def apply_inverse(self, terms, wf, tol, max_it):
        A = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        return self._bicstab(A, wf, tol, max_it)

    def correlation_matrix(self, wf):
        n = self.sites.length()
        out = np.zeros((n, n), dtype=complex)
        for i in range(n):
            for j in range(i, n):
                ampo = AutoMPO(self.sites)
                ampo.add(1.0, "Cdag", i + 1, "C", j + 1)
                op = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
                c = inner(wf, op, wf)
                out[i, j] = c
                out[j, i] = np.conj(c)
        return out

    def four_correlation_tensor(self, wf, accelerate=True):
        n = self.sites.length()
        out = np.zeros((n, n, n, n), dtype=complex)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        current, conj_idx = (i, j, k, l), (l, k, j, i)
                        if accelerate and current > conj_idx:
                            continue
                        ampo = AutoMPO(self.sites)
                        ampo.add(1.0, "Cdag", i + 1, "C", j + 1, "Cdag", k + 1, "C", l + 1)
                        op = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
                        c = inner(wf, op, wf)
                        out[i, j, k, l] = c
                        if current != conj_idx or not accelerate:
                            out[l, k, j, i] = np.conj(c)
        return out

    def four_correlation_tensor_sweep(self, wf, accelerate=True):
        """Single-sweep, environment-reuse four-point correlator tensor
        <Cdag_i C_j Cdag_k C_l> -- a close port of
        Chain::four_correlation_tensor_sweep (mpscpp3/chain_session.h; see
        its docstring for the full algorithm and the derivation of the
        sign table in _four_pt_perm_table() above), so itensor_version=
        "python" gets the same speedup over four_correlation_tensor()'s
        independent per-tuple AutoMPO builds that itensor_version=3 does.

        Only the pairwise-distinct-index entries go through the fast
        sweep; the remaining (subdominant) repeated-index entries fall
        back to four_correlation_tensor()'s own per-tuple AutoMPO method,
        same as the C++ version.

        The "prime the bra's Link tags" trick chain_session.h uses to keep
        a self-overlap's bra and ket legs from colliding doesn't apply
        here (this engine never primes Link indices for that purpose, see
        mpsalgebra.py's inner()) -- instead, following inner()'s own
        convention, the bra side is built once per outer site 'a' via
        _fresh_link_copy(psi) (freshly relabeled Link indices, physical
        indices untouched), and the running environment carries one leg
        from the *ket*'s own (unrelabeled) links and one from this fresh
        bra copy's links.

        `accelerate` only gates the (subdominant) repeated-index fallback
        loop below, unlike four_correlation_tensor()'s own accelerate,
        which skips ~half of the *dominant* per-tuple AutoMPO builds via
        conjugate-pair symmetry. The pairwise-distinct fast-sweep body's
        dominant cost (the shared environment sweep, and the six
        per-pattern scalar() evaluations) is paid once regardless of how
        many output entries a given leaf value gets scattered into, so
        there is no equivalent saving to skip there -- don't expect
        accelerate=True to give the same speedup here it gives
        four_correlation_tensor().
        """
        n = self.sites.length()
        out = np.zeros((n, n, n, n), dtype=complex)

        # Deliberately NOT renormalized: four_correlation_tensor() computes
        # the raw inner(wf,op,wf) with no enforced normalization, and the
        # repeated-index fallback loop below calls inner(wf,...) on this
        # same, unmodified wf -- both must agree on that convention, or the
        # two halves of one output tensor would carry different,
        # inconsistent overall scales for any non-unit-norm wf (an earlier
        # version of this method normalized psi here, mirroring
        # chain_session.h's reduced_dm() convention, which is harmless for
        # an already-unit-norm ground state but silently mis-scaled the
        # pairwise-distinct entries against the repeated-index fallback for
        # any wf that wasn't already unit-normalized, e.g. wf*c). psi is
        # still copied (not aliased) because position() mutates its gauge
        # in place.
        psi = wf.copy()

        def apply_local(T, site, name):
            if name is None:
                return T
            return noPrime(T * self.sites.op(name, site), "Site")

        def fold(E, bra_tensors, site, apply_f, opname):
            T = psi.A(site)
            if apply_f:
                T = apply_local(T, site, "F")
            T = apply_local(T, site, opname)
            piece = dag(bra_tensors[site - 1])
            pieces = [p for p in (E, T, piece) if p is not None]
            return contract_many(pieces)

        for ai in range(0, n - 3):
            a = ai + 1
            psi.position(a)
            bra_tensors = _fresh_link_copy(psi)
            lbase = None
            if a > 1:
                ket_link = commonIndex(psi.A(a - 1), psi.A(a))
                bra_link = commonIndex(bra_tensors[a - 2], bra_tensors[a - 1])
                lbase = delta(ket_link, bra_link)
            for opa_i in range(2):
                opa = "Cdag" if opa_i == 0 else "C"
                la = fold(lbase, bra_tensors, a, False, opa)
                lab_running = la
                for bi in range(ai + 1, n - 2):
                    b = bi + 1
                    if bi > ai + 1:  # gap (a,b): 1 op so far -> F
                        lab_running = fold(lab_running, bra_tensors, b - 1, True, None)
                    for opb_i in range(2):
                        opb = "Cdag" if opb_i == 0 else "C"
                        lab = fold(lab_running, bra_tensors, b, False, opb)
                        labc_running = lab
                        for ci in range(bi + 1, n - 1):
                            c = ci + 1
                            if ci > bi + 1:  # gap (b,c): 2 ops -> no F
                                labc_running = fold(labc_running, bra_tensors, c - 1, False, None)
                            for opc_i in range(2):
                                opc = "Cdag" if opc_i == 0 else "C"
                                labc = fold(labc_running, bra_tensors, c, False, opc)
                                labcd_running = labc
                                for di in range(ci + 1, n):
                                    d = di + 1
                                    if di > ci + 1:  # gap (c,d): 3 ops -> F
                                        labcd_running = fold(labcd_running, bra_tensors, d - 1, True, None)
                                    for opd_i in range(2):
                                        mask = ((1 if opa_i == 0 else 0) | (2 if opb_i == 0 else 0) |
                                                (4 if opc_i == 0 else 0) | (8 if opd_i == 0 else 0))
                                        if bin(mask).count("1") != 2:
                                            continue  # not a Cdag,C,Cdag,C pattern
                                        opd = "Cdag" if opd_i == 0 else "C"
                                        labcd = fold(labcd_running, bra_tensors, d, False, opd)
                                        if d < n:
                                            ket_link = commonIndex(psi.A(d), psi.A(d + 1))
                                            bra_link = commonIndex(bra_tensors[d - 1], bra_tensors[d])
                                            labcd = labcd * delta(ket_link, bra_link)
                                        val = labcd.scalar()
                                        pos = (a - 1, b - 1, c - 1, d - 1)
                                        for (pmask, ir, jr, kr, lr, sign) in _four_pt_perm_table():
                                            if pmask != mask:
                                                continue
                                            i, j, k, l = pos[ir], pos[jr], pos[kr], pos[lr]
                                            out[i, j, k, l] = sign * val

        # Repeated-index entries (subdominant, not covered above): same
        # per-tuple AutoMPO method as four_correlation_tensor().
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        if len({i, j, k, l}) == 4:
                            continue
                        current, conj_idx = (i, j, k, l), (l, k, j, i)
                        if accelerate and current > conj_idx:
                            continue
                        ampo = AutoMPO(self.sites)
                        ampo.add(1.0, "Cdag", i + 1, "C", j + 1, "Cdag", k + 1, "C", l + 1)
                        op_ = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
                        corr_val = inner(wf, op_, wf)
                        out[i, j, k, l] = corr_val
                        if current != conj_idx or not accelerate:
                            out[l, k, j, i] = np.conj(corr_val)
        return out

    def kpm_dynamical_correlator(self, terms_i, terms_j, kpmmaxm, kpm_scale, kpm_accelerate,
                                  kpm_n_scale, delta, kpm_cutoff):
        if not self.have_H:
            raise RuntimeError("Chain.kpm_dynamical_correlator called before set_hamiltonian")
        if self.wf0 is None:
            self.gs_energy()
        scaled_H, emin, emax, scale = self._scaled_hamiltonian(kpm_scale)
        n = int(round((emax - emin) / delta)) * kpm_n_scale
        m1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        m2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi1 = self._apply_mpo_with(m1, self.wf0, kpm_cutoff, kpmmaxm)
        psi2 = self._apply_mpo_with(m2, self.wf0, kpm_cutoff, kpmmaxm)
        moments = self._kpm_moments(scaled_H, psi1, psi2, n, kpmmaxm, kpm_cutoff, kpm_accelerate)
        return moments, emin, emax, scale, n

    def general_kpm(self, terms_x, wfa, wfb, kpmmaxm, kpm_accelerate, num_polynomials, kpm_cutoff):
        m = to_mpo(AutoMPO.from_terms(self.sites, terms_x), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        return self._kpm_moments(m, wfa, wfb, num_polynomials, kpmmaxm, kpm_cutoff, kpm_accelerate)

    def metts_vev(self, terms_ops, T, nsamples=200, nwarmup=20,
                  dbeta_half_step=0.05, basis_ops=("Sz", "Sx"), seed=None, niter=30,
                  njobs=1):
        """Finite-temperature <A> for one or more operators A via METTS
        sampling (White & Stoudenmire, arXiv:1002.1305 -- see metts.py for
        the full algorithm). Reuses self.H (already built by
        set_hamiltonian) and this chain's own cutoff/maxm as the TDVP
        truncation controls, same convention as every other method here.
        njobs is a keyword-only *addition*, python-only -- vevtk/
        mettsvev.py never forwards it to the compiled v3 session, see that
        module's own njobs guard.

        terms_ops: a list of MultiOperator.to_terms() outputs, one per
        observable. All operators are measured on the *same* sampled
        Markov chain of METTS states -- the (nwarmup+nsamples) imaginary-
        time evolutions are the expensive part of this method, and
        metts_thermal_average() already shares them across every op in
        `ops`, so this just exposes that instead of forcing one full
        sampling run per operator (see vevtk/mettsvev.py, the
        Many_Body_Chain-facing wrapper, for the single-operator
        convenience path built on top of this).

        Returns (means, stderrs), each a list matching terms_ops -- see
        metts_thermal_average()'s own docstring for what these mean and
        their (Markov-correlated, so optimistic) error bars.
        """
        if not self.have_H:
            raise RuntimeError("Chain.metts_vev called before set_hamiltonian")
        ops = [to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF,
                      maxdim=self.mpomaxm) for terms_op in terms_ops]
        beta = 1.0 / T
        return _metts_thermal_average(
            self.H, self.sites, ops, beta, nsamples, nwarmup=nwarmup,
            dbeta_half_step=dbeta_half_step, cutoff=self.cutoff, maxdim=self.maxm,
            basis_ops=basis_ops, seed=seed, niter=niter, njobs=njobs)

    def metts_dynamical_correlator(self, terms_a, terms_b, T, nt, dt, nsamples,
                                    nwarmup=20, dbeta_half_step=0.05,
                                    basis_ops=("Sz", "Sx"), seed=None, niter=30,
                                    tdvp_cutoff=None, tdvp_maxdim=None,
                                    tdvp_niter=50, njobs=1):
        """Real-time finite-temperature correlator C_AB(t)=<A(t)B>_T via
        dynamical METTS sampling (arXiv:2405.18484, Sec. II) -- see
        metts.py's metts_dynamical_correlator() for the algorithm itself.
        Reuses self.H (already built by set_hamiltonian) and this chain's
        own cutoff/maxm for the imaginary-time METTS sampling; tdvp_cutoff/
        tdvp_maxdim (default None -> this chain's own cutoff/maxm, same as
        metts.py's own default) separately control the real-time evolution
        of |v_i(t)>/|w_i(t)> -- see metts.py's own tdvp_cutoff/tdvp_maxdim
        docstring for why that evolution generally wants a looser cutoff/
        larger bond dimension than the imaginary-time sampling step.

        terms_a, terms_b: MultiOperator.to_terms() outputs for A and B --
        used exactly as given (no dagger), same convention as
        vev()/apply_operator() above.

        Returns (means, stderrs), each a length-nt array -- see
        metts_dynamical_correlator()'s own docstring for their meaning.
        """
        if not self.have_H:
            raise RuntimeError("Chain.metts_dynamical_correlator called before set_hamiltonian")
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_a), cutoff=_BUILD_CUTOFF,
                   maxdim=self.mpomaxm)
        B = to_mpo(AutoMPO.from_terms(self.sites, terms_b), cutoff=_BUILD_CUTOFF,
                   maxdim=self.mpomaxm)
        beta = 1.0 / T
        return _metts_dynamical_correlator(
            self.H, self.sites, A, B, beta, nt, dt, nsamples, nwarmup=nwarmup,
            dbeta_half_step=dbeta_half_step, cutoff=self.cutoff, maxdim=self.maxm,
            tdvp_cutoff=tdvp_cutoff, tdvp_maxdim=tdvp_maxdim,
            basis_ops=basis_ops, seed=seed, niter=niter, tdvp_niter=tdvp_niter,
            njobs=njobs)

    def nhkpm_moments(self, terms_hs, terms_hs_dag, wfa, wfb, n, kpmmaxm, kpmcutoff):
        """Non-Hermitian KPM (port of NHKPM.jl, see
        mpscpp3/chain_session.h's Chain::nhkpm_moments for the full
        algorithm comment and nonhermitian/kpm.py for the calling
        convention): biorthogonal Chebyshev moments from a coupled
        forward/adjoint recursion built from hs=(z*Id-H)/E_max and its
        dagger (terms_hs/terms_hs_dag, built at the MultiOperator level in
        Python -- one call per requested frequency z). Returns mu_n,
        length n, with mu_n[k] = <wfb|vn_{2k-1}> in the 1-based notation
        of the reference."""
        if n < 1:
            raise ValueError("Chain.nhkpm_moments: n must be >= 1")
        hs = to_mpo(AutoMPO.from_terms(self.sites, terms_hs),
                    cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        hsd = to_mpo(AutoMPO.from_terms(self.sites, terms_hs_dag),
                     cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)

        v = wfa * 1.0
        alpha_prev2 = self._apply_mpo_with(hsd, v, kpmcutoff, kpmmaxm)  # alpha[1]
        alpha_prev1 = mps_sum(self._apply_mpo_with(hs, alpha_prev2, kpmcutoff, kpmmaxm) * 2.0,
                               v * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)  # alpha[2]

        vn_prev2 = v * 1.0  # vn[1]
        vn_prev1 = self._apply_mpo_with(hs, vn_prev2, kpmcutoff, kpmmaxm) * 2.0  # vn[2]

        mu = [inner(wfb, vn_prev2)]
        for _ in range(1, n):
            alpha_x = mps_sum(self._apply_mpo_with(hsd, alpha_prev1, kpmcutoff, kpmmaxm) * 2.0,
                               alpha_prev2 * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)
            alpha_y = mps_sum(self._apply_mpo_with(hs, alpha_x, kpmcutoff, kpmmaxm) * 2.0,
                               alpha_prev1 * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)
            t = mps_sum(self._apply_mpo_with(hsd, vn_prev1, kpmcutoff, kpmmaxm) * 2.0,
                        vn_prev2 * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)
            vn_x = mps_sum(alpha_prev1 * 2.0, t, cutoff=kpmcutoff, maxdim=kpmmaxm)
            vn_y = mps_sum(self._apply_mpo_with(hs, vn_x, kpmcutoff, kpmmaxm) * 2.0,
                           vn_prev1 * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)

            mu.append(inner(wfb, vn_x))

            alpha_prev2, alpha_prev1 = alpha_x, alpha_y
            vn_prev2, vn_prev1 = vn_x, vn_y
        return mu

    # -- private helpers, mirroring chain_session.h's own private section --

    def _default_mps(self):
        return randomMPS(self.sites, self.maxm)

    def _apply_mpo(self, K, x, x0=None):
        out = applyMPO(K, x, x0=x0, cutoff=self.cutoff, maxdim=self.maxm)
        out.noPrime("Site")
        return out

    def _apply_mpo_with(self, K, x, cutoff, maxdim, x0=None):
        out = applyMPO(K, x, x0=x0, cutoff=cutoff, maxdim=maxdim)
        out.noPrime("Site")
        return out

    def _make_sweeps(self, ns=None, maxdim=None):
        if ns is None:
            ns = self.nsweeps
        if maxdim is None:
            maxdim = self.maxm
        sweeps = Sweeps(ns)
        sweeps.maxdim = maxdim
        sweeps.cutoff = self.cutoff
        sweeps.noise = self.noise
        for i in range(ns // 2, ns):
            sweeps.setnoise(i, 0.0)
        return sweeps

    def _minimum_energy(self):
        if self._bandwidth_min is None:
            self._bandwidth_min = self.gs_energy(skip_dmrg=True)
        return self._bandwidth_min

    def _maximum_energy(self):
        if self._bandwidth_max is None:
            # Reduced-effort DMRG on -H, mirroring the compiled
            # backends' maximum_energy(): only a spectral *bound*
            # (Chebyshev window / excited-state penalty weight), never a
            # physical result. A variational underestimate is tolerated
            # by kpm_scale's margin (~bandwidth/6 of headroom) and only
            # shrinks the KPM moment count; a too-tight bound is caught
            # loudly by _check_kpm_moment.
            psi = self._default_mps()
            sweeps = self._make_sweeps(ns=min(self.nsweeps, 5),
                                       maxdim=min(self.maxm, 20))
            neg_H = self.H * (-1.0)
            self._bandwidth_max = -dmrg(psi, neg_H, sweeps, quiet=not self.verbose)
        return self._bandwidth_max

    def _bandwidth(self):
        return self._maximum_energy() - self._minimum_energy()

    def _energy_fluctuation(self, psi1, H):
        psi1 = psi1.copy()
        psi1.normalize()
        psi2 = self._apply_mpo(H, psi1)
        de = inner(psi1, psi2).real
        de = inner(psi2, psi2).real - de * de
        return de

    def _gram_schmidt(self, wfs):
        wfs = list(wfs)
        for i in range(1, len(wfs)):
            for j in range(i):
                proj = wfs[j] * inner(wfs[j], wfs[i])
                wf = mps_sum(wfs[i], proj * (-1.0), cutoff=self.cutoff, maxdim=self.maxm)
                wf.normalize()
                wfs[i] = wf
        return wfs

    def _scaled_hamiltonian(self, kpm_scale):
        emin = self._minimum_energy()
        emax = self._maximum_energy()
        shift = -(emin + emax) / 2.0
        ampo = AutoMPO(self.sites)
        ampo.add(shift, "Id", 1)
        shift_mpo = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        m = mps_sum(self.H, shift_mpo, cutoff=self.cutoff, maxdim=self.mpomaxm)
        scale = 1.0 / ((emax - emin) * kpm_scale)
        m = m * scale
        return m, emin, emax, scale

    def _same_mps(self, vi, vj, maxm, cutoff):
        d = mps_sum(vi * 1.0, vj * (-1.0), cutoff=cutoff, maxdim=maxm)
        dd = np.sqrt(inner(d, d).real)
        return dd < 1e-10

    def _sum_mpo(self, A1, A2):
        return mps_sum(A1, A2, cutoff=self.cutoff, maxdim=self.maxm)

    def _mult_mpo(self, A1, A2):
        out = nmultMPO(A1, A2.copy().prime(), cutoff=self.cutoff, maxdim=self.mpomaxm)
        out.mapPrime(2, 1)
        return out

    def _bicstab(self, A, b, tol, max_it):
        x = b.copy()
        r_old = mps_sum(b, self._apply_mpo(A, x) * (-1.0), cutoff=self.cutoff, maxdim=self.maxm)
        r_ = r_old.copy()
        p = r_old.copy()
        k = 0
        while k < max_it:
            Ap = self._apply_mpo(A, p)
            alpha = inner(self.conjugate(r_old), r_) / inner(self.conjugate(Ap), r_)
            s = mps_sum(r_old, Ap * (-alpha), cutoff=self.cutoff, maxdim=self.maxm)
            As = self._apply_mpo(A, s)
            w = inner(self.conjugate(As), s) / inner(self.conjugate(As), As)
            x = mps_sum(x, mps_sum(p * alpha, s * w, cutoff=self.cutoff, maxdim=self.maxm),
                        cutoff=self.cutoff, maxdim=self.maxm)
            r_new = mps_sum(s, As * (-w), cutoff=self.cutoff, maxdim=self.maxm)
            res = np.sqrt(abs(inner(self.conjugate(r_new), r_new).real))
            if res <= tol:
                break
            beta = (alpha / w) * inner(self.conjugate(r_new), r_) / inner(self.conjugate(r_old), r_)
            p = mps_sum(r_new, mps_sum(p, Ap * (-w), cutoff=self.cutoff, maxdim=self.maxm) * beta,
                        cutoff=self.cutoff, maxdim=self.maxm)
            r_old = r_new
            k += 1
        return x

    def _custom_exp(self, H, z):
        Iden = self._identity_mpo()
        out = self._sum_mpo(Iden, H * z)
        H2 = self._mult_mpo(H, H)
        out = self._sum_mpo(out, H2 * (0.5 * z * z))
        return out

    def _evoloperator(self, H, dt):
        # exp(-i*dt*H) Taylor-expanded to (nominally) 3rd order -- verbatim
        # port of mpscpp3/chain_session.h's evoloperator(), including its
        # deliberately-reproduced quirk: H3 is computed but the z^3/6 term
        # multiplies H2 again, not H3. See that file's own comment for why
        # this is preserved rather than fixed. dt may be complex (see
        # evolve_taylor_step() below, used by tdz.py's "TDZ" submode):
        # z=-1j*dt generalizes complex(0.0,-dt) correctly for that case
        # (the two coincide whenever dt is real).
        Iden = self._identity_mpo()
        z = -1j * dt
        out = self._sum_mpo(Iden, H * z)
        H2 = self._mult_mpo(H, H)
        H3 = self._mult_mpo(H, H2)  # computed to match the original; unused below, see note
        del H3
        out = self._sum_mpo(out, H2 * (0.5 * z * z))
        out = self._sum_mpo(out, H2 * (z * z * z / 6.0))  # NOTE: original uses H2 here, not H3
        return out

    def _identity_mpo(self):
        ampo = AutoMPO(self.sites)
        ampo.add(1.0, "Id", 1)
        return to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)

    def _kpm_moments_full(self, m, vi, vj, n, kpmmaxm, kpmcutoff):
        out = []
        v = vi * 1.0
        am = vi * 1.0
        a = self._apply_mpo_with(m, v, kpmcutoff, kpmmaxm)
        # legitimate moments <vj|T_k|vi> are bounded by ||vi||*||vj||
        # (NOT by the zeroth moment <vj|vi>, which can be ~0 for a
        # near-orthogonal cross-correlator pair)
        bound = np.sqrt(abs(inner(vi, vi).real * inner(vj, vj).real))
        out.append(inner(vj, v))
        out.append(inner(vj, a))
        for _ in range(n):
            ap = self._apply_mpo_with(m, a, kpmcutoff, kpmmaxm)
            ap = mps_sum(ap * 2.0, am * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)
            out.append(inner(vj, ap))
            self._check_kpm_moment(out, bound)
            am = a * 1.0
            a = ap * 1.0
        return out

    def _kpm_moments_accelerated(self, m, vi, n, kpmmaxm, kpmcutoff):
        out = []
        a = self._apply_mpo_with(m, vi, kpmcutoff, kpmmaxm)
        am = vi * 1.0
        mu0 = inner(vi, vi)
        mu1 = inner(vi, a)
        # here vi==vj, so the moment bound ||vi||*||vj|| is just mu0
        bound = abs(mu0)
        out.append(mu0)
        out.append(mu1)
        for _ in range(n // 2):
            ap = self._apply_mpo_with(m, a, kpmcutoff, kpmmaxm)
            ap = mps_sum(ap * 2.0, am * (-1.0), cutoff=kpmcutoff, maxdim=kpmmaxm)
            bk = 2.0 * inner(a, a) - mu0
            bk1 = 2.0 * inner(a, ap) - mu1
            out.append(bk)
            out.append(bk1)
            self._check_kpm_moment(out, bound)
            am = a * 1.0
            a = ap * 1.0
        return out

    @staticmethod
    def _check_kpm_moment(out, bound):
        # Chebyshev moments of a correctly scaled Hamiltonian (spectrum
        # inside [-1,1]) satisfy |<vj|T_k|vi>| <= ||vi||*||vj|| = bound;
        # exponential growth beyond it means the scaled spectrum leaked
        # outside [-1,1] (band-edge estimate too tight for the chosen
        # kpm_scale) and every subsequent moment is garbage. Mirrors the
        # compiled backends' check_kpm_moment; the +1.0 keeps the
        # threshold meaningful when both norms are tiny.
        if abs(out[-1]) > 1e3 * (bound + 1.0):
            raise RuntimeError(
                "KPM moments diverging: scaled spectrum outside [-1,1] "
                "(band-edge estimate too tight; increasing kpm_scale "
                "widens the safety margin)")

    def _kpm_moments(self, m, vi, vj, n, kpmmaxm, kpmcutoff, accelerate):
        if accelerate and self._same_mps(vi, vj, self.maxm, self.cutoff):
            return self._kpm_moments_accelerated(m, vi, n, kpmmaxm, kpmcutoff)
        return self._kpm_moments_full(m, vi, vj, n, kpmmaxm, kpmcutoff)
