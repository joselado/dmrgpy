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
from .mpsalgebra import applyMPO, inner, nmultMPO, BatchedBras
from .mpsalgebra import sum as mps_sum
from .mpsalgebra import randomMPS, traceC
from .mpsalgebra import _fresh_link_copy
from .metts import metts_thermal_average as _metts_thermal_average
from .metts import metts_dynamical_correlator as _metts_dynamical_correlator
from .sites import SiteX
from .svd import svd
from .sweeps import Sweeps
from .mpscontainer import _link_at
from .tdvp import tdvp_step as _tdvp_step_fn
from .tebd import TEBDEvolver as _TEBDEvolver
from .gse import global_subspace_expand as _global_subspace_expand_fn
from .kpm_energy_truncation import energy_truncate as _kpm_energy_truncate
from .tensor import ITensor, commonIndex, contract_many, dag, delta, noPrime, prime, swapPrime

_BUILD_CUTOFF = 1e-14  # mo_terms.h's build_mpo() never exposes a cutoff knob at all


def _mps_arrays_lpr(psi):
    """Every site tensor of `psi` as a plain (chi_left, d, chi_right) NumPy
    array, with length-1 axes inserted at the two open ends.

    Index order on an MPS tensor is not guaranteed, so each is transposed
    explicitly against its own (left link, Site, right link) rather than
    assumed -- `transpose_to` raises rather than silently mis-ordering if the
    tensor is not what this expects."""
    n = psi.length()
    out = []
    for i in range(1, n + 1):
        T = psi.A(i)
        phys = next(ind for ind in T.inds if ind.hastags("Site"))
        left = commonIndex(psi.A(i - 1), T) if i > 1 else None
        right = commonIndex(T, psi.A(i + 1)) if i < n else None
        order = [ind for ind in (left, phys, right) if ind is not None]
        arr = T.transpose_to(order)
        if left is None:
            arr = arr.reshape((1,) + arr.shape)
        if right is None:
            arr = arr.reshape(arr.shape + (1,))
        out.append(np.ascontiguousarray(arr))
    return out


def _four_pt_site_sort_sign(site_list):
    """(-1)**(number of inversions) of `site_list` -- the fermionic sign
    picked up by reordering a product of operators at those sites into
    site-sorted order.

    `autompo.HTerm.resolve()` does not supply this, because by the time it
    runs the factors are already site-sorted: `HTerm.add()` inserts each
    factor in site order and applies the anticommutation sign incrementally
    as it goes. Callers that bypass `HTerm` -- like the local folds here,
    which resolve per-site matrices directly -- have to supply it
    themselves, and this is that sign.

    It is deliberately the *simple* form of `HTerm.add()`'s rule, valid
    because every factor in a four-point <Cdag C Cdag C> string is
    fermionic: then the correction is just the parity of the sorting
    permutation, counting pairs whose sites differ (same-site pairs keep
    their relative order and contribute nothing). `HTerm.add()` is the
    general version, which also has to check *which* of the operators it
    passes are fermionic. `test_four_pt_sign_agrees_with_hterm` pins the two
    against each other so this cannot drift.

    Verified against the AutoMPO+to_mpo+inner pipeline these folds replace
    on 67 tuples (7 hand-picked plus 60 random) -- without it, tuples such
    as Cdag_1 C_1 Cdag_3 C_2 come out with a flipped sign."""
    inv = 0
    for x in range(len(site_list)):
        for y in range(x + 1, len(site_list)):
            if site_list[x] > site_list[y]:
                inv += 1
    return -1 if inv % 2 else 1


def _four_pt_repeated_tuples(n):
    """Every (i,j,k,l) in [0,n)^4 whose four indices are NOT pairwise
    distinct, enumerated directly rather than by scanning all n^4 and
    discarding the distinct ones.

    Scanning is O(n^4) pure-Python iterations just to *find* an O(n^3) set;
    at n=40 that is 2.5M wasted loop steps. Enumerating instead by (choice
    of the m<=3 distinct sites) x (surjective assignment of the four
    positions onto them) touches exactly the tuples that are yielded:
    n + 14*C(n,2) + 36*C(n,3)."""
    for m in (1, 2, 3):
        for combo in itertools.combinations(range(n), m):
            for assign in itertools.product(range(m), repeat=4):
                if len(set(assign)) != m:
                    continue  # must actually use all m sites
                yield tuple(combo[a] for a in assign)


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



def _max_link_dim(psi):
    """Largest bond dimension of an MPS, i.e. ITensor's own
    maxLinkDim(psi). Used by Chain.gs_energy() to keep the bond-dimension
    ramp from truncating a warm starting wavefunction."""
    out = 0
    for i in range(1, psi.length()):
        link = _link_at(psi, i, i + 1)
        if link is not None and link.dim > out:
            out = link.dim
    return out


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

        # Ground-state bond-dimension ramp, see set_bond_ramp()/
        # _make_sweeps_ramped(). Defaults mirror both the C++ backends'
        # and Many_Body_Chain's own (manybodychain.py), which is what
        # actually drives them in practice.
        self.bond_ramp = True
        self.bond_ramp_start = 10
        self.bond_ramp_fraction = 0.5
        self.bond_ramp_noise_decay = 0.1

        # KPM energy truncation (Holzner et al. PRB 83, 195115 (2011),
        # Sec. III-B), see kpm_energy_truncation.py -- off by default, so
        # kpm_scale/kpm_dynamical_correlator()'s existing behavior is
        # unaffected unless a caller opts in via set_kpm_energy_truncation.
        self.kpm_energy_truncate = False
        self.kpm_truncate_dK = 10
        self.kpm_truncate_nsweeps = 3
        self.kpm_truncate_threshold = 1.0

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

    def set_bond_ramp(self, enabled, start, fraction, noise_decay):
        """Ground-state bond-dimension ramp -- the counterpart of
        mpscppN/chain_session.h's Chain::set_bond_ramp(). See
        _make_sweeps_ramped() for what the schedule looks like;
        enabled=False restores the original flat-maxm schedule exactly."""
        self.bond_ramp = enabled
        self.bond_ramp_start = start
        self.bond_ramp_fraction = fraction
        self.bond_ramp_noise_decay = noise_decay

    def set_verbose(self, verbose):
        self.verbose = verbose

    def set_kpm_energy_truncation(self, enabled, dK, n_sweeps, threshold):
        """Python-only KPM control knob (no C++ v2/v3 counterpart -- see
        kpmdmrg.py, which only calls this when itensor_version=="python").
        See kpm_energy_truncation.energy_truncate() for what dK/n_sweeps/
        threshold mean; enabled=False (the default) reproduces this
        class's original kpm_dynamical_correlator()/general_kpm()
        behavior exactly."""
        self.kpm_energy_truncate = enabled
        self.kpm_truncate_dK = dK
        self.kpm_truncate_nsweeps = n_sweeps
        self.kpm_truncate_threshold = threshold

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
        # A warm start (set_wavefunction, or simply the previous
        # gs_energy() call's own solution) must not be truncated by the
        # ramp's early low-maxdim sweeps -- pass its bond dimension as the
        # ramp floor. A fresh start is built directly at the ramp's first
        # maxdim rather than at the full maxm only to be truncated.
        floor_dim = 0 if self.wf0 is None else _max_link_dim(self.wf0)
        sweeps = self._make_sweeps_ramped(self.nsweeps, self.maxm, floor_dim)
        if self.wf0 is None:
            self.wf0 = self._default_mps(sweeps.at(0)[0])
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

    def batched_bras(self, wfs):
        """Prepare a fixed family of bra MPS for repeated overlap against
        many different kets, returning an opaque handle whose .overlaps(ket)
        gives every <bra_b|ket> from a single batched sweep.

        Exists for tdz.py's complex-time correlator, which overlaps one
        evolving state against the same n_max+1 bras at every time step;
        see mpsalgebra.BatchedBras for what it does and why the batch is
        what makes those overlaps worth putting on a device. Backends
        without this method fall back to a loop over overlap(), which is
        what tdz.py does when hasattr() says so."""
        return BatchedBras(wfs)

    def overlap_aMb(self, wf1, terms, wf2):
        A = to_mpo(AutoMPO.from_terms(self.sites, terms), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        return inner(wf1, A, wf2)

    def sum_mps(self, wf1, wf2):
        return mps_sum(wf1, wf2, cutoff=self.cutoff, maxdim=self.maxm)

    def scale_mps(self, wf, z):
        """Rescale an MPS by a complex number, mirroring
        mpscppN/chain_session.h's scale_mps(). mpscontainer.MPS.__mul__
        already scales a single site tensor, so this is O(chi^2 d) with no
        contraction and no truncation -- see the C++ comment for why
        mps.py needs it (the alternative is a full identity-MPO sweep per
        scalar multiplication)."""
        return wf * complex(z)

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
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi2, psi1))
            psi1 = self._apply_mpo(expH, psi1, x0=psi1) if fit_td else self._apply_mpo(expH, psi1)
            psi1.normalize()
            psi1 = psi1 * norm0
        return correlator, psi1

    def evolve_and_measure(self, terms_h, terms_op, wf, nt, dt, fit_td=True):
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        expH = self._evoloperator(to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm), dt)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf
        correlator = []
        for _ in range(nt):
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi, A, psi))
            psi = self._apply_mpo(expH, psi, x0=psi) if fit_td else self._apply_mpo(expH, psi)
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
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi2, psi1))
            psi1 = _tdvp_step_fn(psi1, Hshift, dt, cutoff=self.cutoff, maxdim=self.maxm, niter=50)
            psi1.normalize()
            psi1 = psi1 * norm0
        return correlator, psi1

    def evolve_and_measure_tdvp(self, terms_h, terms_op, wf, nt, dt):
        """`wf` is copied before evolving: _tdvp_step_fn mutates its input
        MPS's tensor list in place (via set_A()), so evolving it directly
        would silently corrupt the caller's own wavefunction object --
        including self.wf0 itself on the common wf=None default path in
        timedependent.py's evolve_and_measure_dmrg(). Confirmed directly:
        without this copy, a call with the default wf changes what
        self.get_gs() already computed, so a second, unrelated
        measurement on "the same" ground state silently sees a partially
        time-evolved state instead."""
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf.copy()
        correlator = []
        for _ in range(nt):
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi, A, psi))
            psi = _tdvp_step_fn(psi, H, dt, cutoff=self.cutoff, maxdim=self.maxm, niter=50)
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
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi2, psi1))
            if it < gse_sweeps:
                psi1 = self.global_subspace_expand(Hshift, psi1, krylov_order, gse_cutoff)
            psi1 = _tdvp_step_fn(psi1, Hshift, dt, cutoff=self.cutoff, maxdim=self.maxm,
                    niter=50, num_center=1)
            psi1.normalize()
            psi1 = psi1 * norm0
        return correlator, psi1

    def evolve_and_measure_tdvp_gse(self, terms_h, terms_op, wf, nt, dt, gse_sweeps,
            krylov_order, gse_cutoff):
        """GSE counterpart of evolve_and_measure_tdvp() above -- see
        quench_tdvp_gse()'s docstring and evolve_and_measure_tdvp()'s own
        docstring for why `wf` is copied here too."""
        H = to_mpo(AutoMPO.from_terms(self.sites, terms_h), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf.copy()
        correlator = []
        for it in range(nt):
            # Measure before evolving, so correlator[k] is C(k*dt),
            # matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            # pairs it with -- see evolution_dmrg_DC()'s own comment there.
            correlator.append(inner(psi, A, psi))
            if it < gse_sweeps:
                psi = self.global_subspace_expand(H, psi, krylov_order, gse_cutoff)
            psi = _tdvp_step_fn(psi, H, dt, cutoff=self.cutoff, maxdim=self.maxm,
                    niter=50, num_center=1)
        return correlator, psi

    def quench_tebd(self, terms_h, terms_i, terms_j, nt, dt):
        """TEBD counterpart of quench_tdvp() above: identical setup/
        measurement, but each per-step evolution is a 2nd-order Trotter
        TEBD step (tebd.py's TEBDEvolver, gates built once up front) in
        place of TDVP's per-step Krylov exponentiation -- only valid for
        a strictly nearest-neighbor terms_h (TEBDEvolver/
        bond_hamiltonians() raise NotImplementedError otherwise). Folds
        the same ground-energy shift into the Hamiltonian passed to
        TEBDEvolver as quench_tdvp() applies to Hshift, for the same
        reason (removing the state's own e^{-i*EGS*t} global phase)."""
        if self.wf0 is None:
            self.gs_energy()
        ampo_h = AutoMPO.from_terms(self.sites, terms_h)
        H = to_mpo(ampo_h, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        EGS = inner(self.wf0, H, self.wf0).real / inner(self.wf0, self.wf0).real
        terms_shifted = list(terms_h) + [(-EGS, [("Id", 1)])]
        evolver = _TEBDEvolver(self.sites, terms_shifted, dt, cutoff=self.cutoff, maxdim=self.maxm)
        A1 = to_mpo(AutoMPO.from_terms(self.sites, terms_i), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        A2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi1 = self._apply_mpo(A1, self.wf0)
        psi2 = self._apply_mpo(A2, self.wf0)
        norm0 = np.sqrt(inner(psi1, psi1))
        correlator = []
        for _ in range(nt):
            # Measure before evolving -- see quench()'s own comment.
            correlator.append(inner(psi2, psi1))
            psi1 = evolver.step(psi1)
            psi1.normalize()
            psi1 = psi1 * norm0
        return correlator, psi1

    def evolve_and_measure_tebd(self, terms_h, terms_op, wf, nt, dt):
        """TEBD counterpart of evolve_and_measure_tdvp() above -- see
        quench_tebd()'s docstring and evolve_and_measure_tdvp()'s own
        docstring for why `wf` is copied here too (TEBDEvolver.step()
        mutates its input in place, same as _tdvp_step_fn)."""
        evolver = _TEBDEvolver(self.sites, terms_h, dt, cutoff=self.cutoff, maxdim=self.maxm)
        A = to_mpo(AutoMPO.from_terms(self.sites, terms_op), cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        psi = wf.copy()
        correlator = []
        for _ in range(nt):
            # Measure before evolving -- see quench()'s own comment.
            correlator.append(inner(psi, A, psi))
            psi = evolver.step(psi)
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

    def ddmrg_correction_vector(self, terms_i, terms_j, omega, eta, energy,
                                maxdim, nsweeps, x0=None):
        """Correction-vector dynamical correlator by variational sweeping
        (ddmrg.py) rather than by cvm.py's global conjugate gradient.

        Returns (value, x), the spectral weight -W_min/(pi*eta) and the
        converged correction vector -- the latter so a caller sweeping a
        frequency window can seed the next point with it (see
        ddmrg.correction_vector's note on why warm-starting is safe here
        and is not for the global CG).

        terms_i/terms_j are the two operators of <GS|A (z-H)^-1 B|GS>;
        this solver requires A = B^dagger and cvm.py checks that before
        calling, but both are taken so the signature matches
        cvm_dynamical_correlator's and so B is built from terms_j exactly
        as that path does."""
        from . import ddmrg
        if self.wf0 is None:
            raise RuntimeError("Chain.ddmrg_correction_vector called before gs_energy")
        S2 = to_mpo(AutoMPO.from_terms(self.sites, terms_j),
                    cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        # (H - E0 - omega) as an MPO, then M = that squared plus eta^2.
        ampo = AutoMPO(self.sites)
        ampo.add(-(energy + omega), "Id", 1)
        shift = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        Hsh = mps_sum(self.H, shift, cutoff=self.cutoff, maxdim=self.mpomaxm)
        ampo2 = AutoMPO(self.sites)
        ampo2.add(eta * eta, "Id", 1)
        eta2 = to_mpo(ampo2, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        M = mps_sum(self._mult_mpo(Hsh, Hsh), eta2,
                    cutoff=self.cutoff, maxdim=self.mpomaxm)
        b = self._apply_mpo(S2, self.wf0) * (-eta)
        x, W = ddmrg.correction_vector(M, b, maxdim, nsweeps, self.cutoff,
                                        x0=x0)
        return -W / (np.pi * eta), x

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
        sweep. The remaining repeated-index entries -- there are O(n^3) of
        them, but each is a *local* operator on at most 3 distinct sites --
        are folded directly by `_four_pt_fill_repeated`. They used to fall
        back to four_correlation_tensor()'s per-tuple AutoMPO method, as
        the C++ version still does; that was measured at 96% of this
        method's total runtime at n=12, so "subdominant" (as this docstring
        used to call it) was exactly backwards.

        The "prime the bra's Link tags" trick chain_session.h uses to keep
        a self-overlap's bra and ket legs from colliding doesn't apply
        here (this engine never primes Link indices for that purpose, see
        mpsalgebra.py's inner()) -- instead, following inner()'s own
        convention, the bra side is built once per outer site 'a' via
        _fresh_link_copy(psi) (freshly relabeled Link indices, physical
        indices untouched), and the running environment carries one leg
        from the *ket*'s own (unrelabeled) links and one from this fresh
        bra copy's links.

        `accelerate` only gates the repeated-index entries, unlike
        four_correlation_tensor()'s own accelerate,
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

        # Repeated-index entries (not covered by the sweep above). These
        # used to go through the same per-tuple AutoMPO+to_mpo+inner method
        # as four_correlation_tensor(), described as "subdominant" because
        # there are O(n^3) of them against the sweep's O(n^4). Measured,
        # that reasoning was exactly backwards: fewer tuples times a far
        # more expensive per-tuple cost dominated everything. At n=12 the
        # fallback was 96% of the total runtime (24.2s of 25.3s: 4500
        # to_mpo builds at 15.9s plus 4500 full-chain inner products at
        # 8.4s), against 1.1s for the actual sweep.
        #
        # Each of these operators is a product of four Cdag/C factors on at
        # most 3 distinct sites, i.e. a *local* operator -- there is no need
        # to compile an MPO over the whole chain and sweep it. Resolve it to
        # per-site matrices once (Jordan-Wigner threading included) and fold
        # it directly, over [mn,mx] only, with the same mixed-canonical
        # shortcut the sweep body uses at its two ends.
        self._four_pt_fill_repeated(out, psi, n, accelerate)
        return out

    def _four_pt_fill_repeated(self, out, psi, n, accelerate):
        """Fill `out`'s repeated-index entries by local operator folds.

        Two subtleties, both load-bearing:

        * `HTerm.resolve()` supplies the per-site matrices *and* the JW F
          strings, but drops the fermionic sign for reordering operators
          across different sites -- see `_four_pt_site_sort_sign`, which
          supplies it. Getting this wrong flips the sign of entries like
          (i,j,k,l)=(0,0,2,1) with no other symptom.
        * The fold must stay the *raw* <wf|Op|wf>, matching the sweep half
          of this same tensor (see the "deliberately NOT renormalized"
          comment above). `position()` is a gauge transformation and
          preserves it, so folding against the already-positioned `psi` is
          consistent -- verified on a deliberately non-unit-norm wf, the
          case the earlier normalization bug was found with.

        Work is grouped by the tuple's minimum site so one `position()` call
        serves every tuple starting there.

        The fold runs on plain NumPy arrays rather than ITensor objects.
        That is not premature: with the operator content resolved, every
        shape and axis order here is already known, so routing each step
        through `ITensor.__mul__` only pays for Index matching that answers a
        question we have the answer to. Profiled at n=12, the ITensor form
        spent ~80% of its time in `mul_plan`/`_pairwise_result_size`/
        `Index.__hash__` rather than in arithmetic. It also drops the
        `_fresh_link_copy` bra entirely -- a self-overlap's bra is just
        `conj()` of the same array, and the Link-identity collision that copy
        exists to avoid cannot arise when the contraction is positional."""
        from .sites.base import is_fermionic
        by_min = {}
        for (i, j, k, l) in _four_pt_repeated_tuples(n):
            if accelerate and (i, j, k, l) > (l, k, j, i):
                continue
            by_min.setdefault(min(i, j, k, l), []).append((i, j, k, l))

        # per-site operator matrices are fixed for the whole run
        mats = {}

        def site_matrix(site, name):
            key = (site, name)
            if key not in mats:
                mats[key] = self.sites.site_type(site).matrix(name)
            return mats[key]

        for mn0 in sorted(by_min):
            mn = mn0 + 1  # 1-based
            psi.position(mn)
            arrays = _mps_arrays_lpr(psi)   # (chi_l, d, chi_r) per site
            for (i, j, k, l) in by_min[mn0]:
                mx = max(i, j, k, l) + 1
                by_site = {}
                for nm, site in (("Cdag", i + 1), ("C", j + 1),
                                  ("Cdag", k + 1), ("C", l + 1)):
                    by_site.setdefault(site, []).append(nm)
                E = np.eye(arrays[mn - 1].shape[0], dtype=complex)
                carry = False
                for site in range(mn, mx + 1):
                    names = by_site.get(site)
                    if names is None:
                        mat = site_matrix(site, "F") if carry else None
                    else:
                        odd = sum(1 for nm in names if is_fermionic(nm)) % 2 == 1
                        mat = site_matrix(site, names[0])
                        for nm in names[1:]:
                            mat = site_matrix(site, nm) @ mat
                        if carry != odd:
                            mat = site_matrix(site, "F") @ mat
                        carry = carry != odd
                    A = arrays[site - 1]
                    # tensordot, not einsum: at these sizes np.einsum's
                    # `optimize=True` re-derives a contraction path on every
                    # call and costs more than the contraction itself (it
                    # made this loop measurably slower than the ITensor form
                    # it replaced), while `optimize=False` falls back to a
                    # slow C loop for 3-operand contractions.
                    if mat is not None:
                        A_op = np.tensordot(A, mat, axes=([1], [0]))   # (l,r,o)
                        A_op = np.ascontiguousarray(A_op.transpose(0, 2, 1))
                    else:
                        A_op = A
                    X = np.tensordot(E, A_op, axes=([0], [0]))         # (L,p,r)
                    E = np.tensordot(X, np.conj(A), axes=([0, 1], [0, 1]))  # (r,R)
                val = _four_pt_site_sort_sign((i, j, k, l)) * complex(np.trace(E))
                out[i, j, k, l] = val
                if (i, j, k, l) != (l, k, j, i) or not accelerate:
                    out[l, k, j, i] = np.conj(val)

    def four_correlation_tensor_fold(self, wf, cdag_ops, c_ops, accelerate=True):
        """<Cdag_i C_j Cdag_k C_l> for flat fermionic *modes*, by local
        operator folds -- no MPO is built for any tuple.

        `cdag_ops`/`c_ops` are one `(operator_name, site_1based)` pair per
        mode, so this covers both a spinless chain (mode == site, names
        "Cdag"/"C") and a native spinful one (two modes per site, names
        "Cdagup"/"Cup"/"Cdagdn"/"Cdn" on an ElectronSite). Nothing here is
        flavor-aware beyond those names: Jordan-Wigner threading and
        same-site composition come from the site type's own operator
        matrices and "F", exactly as `autompo.HTerm.resolve()` derives them,
        and the cross-site reordering sign from
        `_four_pt_site_sort_sign` (which `resolve()` does not supply --
        see its docstring).

        This exists because the only tensor available for native spinful
        sites was `ctmode="explicit"`, which builds an MPO and sweeps the
        whole chain per tuple; profiled at 4 sites, 55% of its time was in
        `to_mpo` and 38% in `inner`, i.e. 93% doing work a local operator
        does not need. Unlike `four_correlation_tensor_sweep` this does not
        reuse environments *across* tuples, so it is O(n_modes^4 * n_sites)
        rather than O(n_modes^4); it is a large constant-factor win over
        building an MPO per tuple, not an asymptotic one."""
        from .sites.base import is_fermionic
        nm = len(c_ops)
        if len(cdag_ops) != nm:
            raise ValueError("four_correlation_tensor_fold: cdag_ops and "
                              "c_ops must describe the same modes")
        out = np.zeros((nm, nm, nm, nm), dtype=complex)
        n = self.sites.length()
        psi = wf.copy()
        mats = {}

        def site_matrix(site, name):
            key = (site, name)
            if key not in mats:
                mats[key] = self.sites.site_type(site).matrix(name)
            return mats[key]

        by_min = {}
        for i in range(nm):
            for j in range(nm):
                for k in range(nm):
                    for l in range(nm):
                        if accelerate and (i, j, k, l) > (l, k, j, i):
                            continue
                        sites4 = (cdag_ops[i][1], c_ops[j][1],
                                  cdag_ops[k][1], c_ops[l][1])
                        by_min.setdefault(min(sites4), []).append((i, j, k, l))

        for mn in sorted(by_min):
            psi.position(mn)
            arrays = _mps_arrays_lpr(psi)
            for (i, j, k, l) in by_min[mn]:
                factors = [cdag_ops[i], c_ops[j], cdag_ops[k], c_ops[l]]
                site_list = [f[1] for f in factors]
                mx = max(site_list)
                by_site = {}
                for name, site in factors:
                    by_site.setdefault(site, []).append(name)
                E = np.eye(arrays[mn - 1].shape[0], dtype=complex)
                carry = False
                for site in range(mn, mx + 1):
                    names = by_site.get(site)
                    if names is None:
                        mat = site_matrix(site, "F") if carry else None
                    else:
                        odd = sum(1 for x in names if is_fermionic(x)) % 2 == 1
                        mat = site_matrix(site, names[0])
                        for x in names[1:]:
                            mat = site_matrix(site, x) @ mat
                        if carry != odd:
                            mat = site_matrix(site, "F") @ mat
                        carry = carry != odd
                    A = arrays[site - 1]
                    if mat is not None:
                        A_op = np.tensordot(A, mat, axes=([1], [0]))
                        A_op = np.ascontiguousarray(A_op.transpose(0, 2, 1))
                    else:
                        A_op = A
                    X = np.tensordot(E, A_op, axes=([0], [0]))
                    E = np.tensordot(X, np.conj(A), axes=([0, 1], [0, 1]))
                val = _four_pt_site_sort_sign(tuple(site_list)) * complex(np.trace(E))
                out[i, j, k, l] = val
                if (i, j, k, l) != (l, k, j, i) or not accelerate:
                    out[l, k, j, i] = np.conj(val)
        return out

    def four_correlation_tensor_batched(self, wf, cdag_ops, c_ops,
                                        accelerate=True):
        """<Cdag_i C_j Cdag_k C_l> with every transfer contraction batched
        over the site combinations that share it -- see
        `pyitensor/fourpoint.py` for the trie this rests on.

        Same signature and same answer as `four_correlation_tensor_fold`,
        whose conventions it inherits wholesale. `accelerate` is accepted
        and ignored: it exists on the other two implementations to skip the
        conjugate half of the tuples, and here the cost is in the shared
        environments rather than per tuple, so there is nothing to skip."""
        from .fourpoint import four_correlation_tensor_batched
        return four_correlation_tensor_batched(self, wf, cdag_ops, c_ops)

    def kpm_dynamical_correlator(self, terms_i, terms_j, kpmmaxm, kpm_scale, kpm_accelerate,
                                  kpm_n_scale, delta, kpm_cutoff):
        if not self.have_H:
            raise RuntimeError("Chain.kpm_dynamical_correlator called before set_hamiltonian")
        if self.wf0 is None:
            self.gs_energy()
        if self.kpm_energy_truncate:
            scaled_H, emin, emax, scale = self._scaled_hamiltonian_gs_anchored(kpm_scale)
        else:
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

    def _default_mps(self, maxm=None):
        """Random starting MPS. maxm defaults to self.maxm; gs_energy()
        passes the ramp's first sweep dimension instead (see
        _make_sweeps_ramped())."""
        if maxm is None:
            maxm = self.maxm
        return randomMPS(self.sites, max(1, int(maxm)))

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

    def _make_sweeps_ramped(self, ns, maxdim, floor_dim):
        """Sweep schedule with a bond-dimension ramp: spend the first
        bond_ramp_fraction of the schedule growing the bond dimension
        geometrically from bond_ramp_start up to the target maxdim, then
        hold it there, instead of running every sweep at the full target.
        The early sweeps are much cheaper (two-site DMRG is ~O(m^3)) while
        still doing most of the work of finding the right variational
        subspace, so the expensive full-maxdim sweeps start from an
        already-good state -- the standard ITensor sweep-table idiom.

        Growing as fast as possible instead (doubling every sweep until
        maxdim is reached) minimizes the number of *cheap* sweeps and so
        gives away almost all of the speedup -- see
        mpscpp3/chain_session.h's make_sweeps_ramped(), of which this is
        the exact counterpart, for the measurement.

        floor_dim is the bond dimension of an already-present starting
        wavefunction, or 0 when starting fresh: ramping down from a
        converged warm state would truncate exactly the information the
        caller wanted kept, so the ramp never starts below it.

        Noise decays by bond_ramp_noise_decay per ramping sweep and is off
        entirely once the schedule reaches full maxdim, so the final,
        converged sweeps are noise-free. (Note the pure-Python DMRG in
        dmrg.py has no noise term at all -- see its docstring -- so the
        noise part of the schedule is carried for parity with the compiled
        backends rather than acted on here.)
        """
        if not self.bond_ramp or ns < 2 or maxdim <= 1:
            return self._make_sweeps(ns, maxdim)
        start = self.bond_ramp_start if self.bond_ramp_start > 0 else 10
        start = max(start, floor_dim)
        if start >= maxdim:
            return self._make_sweeps(ns, maxdim)
        nr = int(ns * self.bond_ramp_fraction)
        nr = min(max(nr, 1), ns - 1)  # always ends at full maxdim
        sweeps = Sweeps(ns)
        sweeps.cutoff = self.cutoff
        for i in range(ns):
            if i < nr:
                # floor(x+0.5), matching C++'s std::llround exactly
                # (Python's round() is banker's rounding, which differs
                # by one on an exact .5)
                m = int(np.floor(start * (maxdim / start) ** (i / nr) + 0.5))
                m = min(max(m, 1), maxdim)
                nz = self.noise * self.bond_ramp_noise_decay ** i
            else:
                m = maxdim
                nz = 0.0
            sweeps.setmaxdim(i, m)
            sweeps.setnoise(i, nz)
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

    def _scaled_hamiltonian_gs_anchored(self, kpm_scale, safety=0.025):
        """Ground-state-anchored rescaling (Holzner et al., PRB 83,
        195115 (2011), Eqs. (20)-(21b)): places the ground state energy
        E0 at -(1-safety/2) and E0+Ws at +(1-safety/2), where
        Ws=(emax-emin)*kpm_scale is the *same* window-width formula
        _scaled_hamiltonian() uses (so kpm_scale keeps the same meaning/
        units in both) -- but anchored at E0 rather than centered on the
        full many-body bandwidth's midpoint.

        Used by kpm_dynamical_correlator() instead of
        _scaled_hamiltonian() whenever kpm_energy_truncate is enabled:
        shrinking kpm_scale only delivers Sec. III-B's promised
        resolution gain if the narrowed window also sits where a typical
        correlator's spectral weight actually lives -- just above the
        ground state, since its Chebyshev vectors are built by acting
        operators on |0>, not around the geometric middle of the entire
        many-body spectrum. Confirmed directly: with
        _scaled_hamiltonian()'s own bandwidth-midpoint centering, any
        kpm_scale below the safe ~0.5 floor clips the ground state
        itself out of the window first (before ever reaching a genuinely
        useful, narrower-but-still-correlator-covering regime), since
        the ground state already sits at the spectrum's own lower edge,
        not near its middle.

        Returns (m, emin, emax, scale) with the *same* meaning
        downstream (kpm_dynamical_correlator()/kpmdmrg.py's moment-count
        and energy-axis reconstruction) as _scaled_hamiltonian()'s own
        return value, even though emin/emax are no longer literally the
        true many-body band edges: emin=E0 (still correct -- the energy
        axis kpmdmrg.py reconstructs is "excitation energy above the
        ground state", i.e. E_phys-E0, unchanged), and
        emax=E0+Ws (the actual retained window's top, which is what
        should set the polynomial order N via kpmdmrg.py's
        n=(emax-emin)/delta, not the full, unused many-body bandwidth).
        """
        e0 = self._minimum_energy()
        w = self._bandwidth()
        ws = w * kpm_scale
        wp = 1.0 - safety / 2.0
        a = ws / (2.0 * wp)
        scale = 1.0 / a
        shift = -e0 - wp * a
        ampo = AutoMPO(self.sites)
        ampo.add(shift, "Id", 1)
        shift_mpo = to_mpo(ampo, cutoff=_BUILD_CUTOFF, maxdim=self.mpomaxm)
        m = mps_sum(self.H, shift_mpo, cutoff=self.cutoff, maxdim=self.mpomaxm)
        m = m * scale
        return m, e0, e0 + ws, scale

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

    def _maybe_energy_truncate(self, psi, m):
        """Opt-in post-processing of a freshly-formed Chebyshev vector
        `psi` (mutated in place), per set_kpm_energy_truncation() /
        kpm_energy_truncation.py -- a no-op unless explicitly enabled, so
        existing kpm_scale callers see no behavior change. `m` must be
        the same rescaled/shifted Hamiltonian MPO used to build `psi`'s
        own recursion."""
        if not self.kpm_energy_truncate:
            return
        _kpm_energy_truncate(psi, m, dK=self.kpm_truncate_dK,
                              n_sweeps=self.kpm_truncate_nsweeps,
                              threshold=self.kpm_truncate_threshold)

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
            self._maybe_energy_truncate(ap, m)
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
            self._maybe_energy_truncate(ap, m)
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
