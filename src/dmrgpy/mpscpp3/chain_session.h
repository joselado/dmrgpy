#include <tuple>
#include <map> // Chain::MettsEigCache (Chain::metts_vev)
#include <algorithm> // std::next_permutation (four_correlation_tensor_sweep)
#include <random> // std::mt19937_64/std::discrete_distribution (Chain::metts_vev)
#include <cmath> // std::isnan (Chain::gs_energy_generalized's lam0 sentinel)
#include <limits> // std::numeric_limits<double>::quiet_NaN() (ditto)
#include <functional> // std::function (Chain::vx_build_linear_map's own action callback)

// TDVP/tdvp.h and TDVP/basisextension.h are quoted includes, so they
// resolve relative to this file's own directory (mpscpp3/) with no
// Makefile/include-path change needed. mpscpp3-only: TDVP/ has no v2
// counterpart (mpscpp2's Chain has no TDVP methods at all), so
// tdvp_step()/quench_tdvp()/evolve_and_measure_tdvp()/
// global_subspace_expand() below only exist in this file. tdvp_step()'s
// num_center parameter selects one-site (1) or two-site (2, the default)
// TDVP; two-site TDVP grows bond dimension via SVD the same way two-site
// DMRG does, so global_subspace_expand() (TDVP/basisextension.h's
// addBasis(), implementing the Krylov-subspace basis-enrichment scheme of
// arXiv:2005.06104/Phys. Rev. B 102, 094315) is what lets one-site TDVP
// grow bond dimension instead -- see TDVP/README.md for the algorithm.
#include "TDVP/tdvp.h"
#include "TDVP/basisextension.h"

// Port of mpscpp2/chain_session.h to the ITensor v3 API. The Chain class's
// public methods, semantics, and preserved-not-fixed quirks (see
// evoloperator()'s note below) are unchanged from the v2 session/handle
// model -- only the ITensor-library calls themselves are updated:
//
//  - IQTensor/IQIndex no longer exist in v3 (merged into ITensor/Index, QN
//    blocks now live directly on the Index); this actually *simplifies* a
//    couple of spots that relied on IQTensor<->ITensor conversions in v2
//    (e.g. bond_entropy's "ITensor twosite = psi.A(b)*psi.A(b+1)" is a
//    trivial same-type multiplication now, not a cross-type one).
//  - exactApplyMPO(K,x,args)/fitApplyMPO(psi,K,res,args) were unified into
//    a single applyMPO(): applyMPO(K,x,args) (Method="DensityMatrix" by
//    default) replaces exactApplyMPO, and applyMPO(K,x,x0,args) (Method
//    defaults to "Fit" when a guess x0 is given) replaces fitApplyMPO's
//    in-place res-as-initial-guess pattern used here (every fitApplyMPO
//    call in this file is of the form fitApplyMPO(psi,K,psi,args), i.e.
//    the guess is the same state being updated, so it maps directly to
//    applyMPO(K,psi,psi,args)). v2's "older interface" exactApplyMPO(x,K,args)
//    (reversed argument order) has no v3 equivalent, so every such call
//    site below has its arguments swapped back to the canonical (K,x) order.
//  - overlap()/overlapC() still exist in v3 but only as deprecated wrappers
//    around inner()/innerC() (same signatures); this file uses inner/innerC
//    directly -- except plain inner() (unlike v2's overlap(), and unlike
//    innerC()) throws outright ("Cannot use inner(...) with complex
//    MPS/MPO, use innerC(...) instead") the moment either operand is
//    complex-valued. Since MOTerm::coef is always Cplx, any MPO/MPS built
//    from Python-supplied terms is complex-typed even when every
//    coefficient's imaginary part is exactly 0 (confirmed directly: the
//    KPM dynamical correlator path aborts here in same_mps() once wf0_
//    picks up a complex-but-real-valued MPO application upstream). Every
//    "Re[<x|y>]" use below goes through innerC(...).real() instead of
//    plain inner(...) so it works regardless of realness.
//  - Args "Maxm" still works (auto-forwarded to "MaxDim" with a deprecation
//    warning) but this file uses "MaxDim" directly.
//  - IndexType Site/Link selectors for prime()/swapPrime() don't exist in
//    v3 (tags replaced IndexType); this file selects "Site" by the concrete
//    physical Index (sites_.si(site), always tagged "Site" by every site
//    class used here) and "Link" by TagSet("Link").
//  - trace_operator() no longer needs to build an explicit identity MPO to
//    compute <Id|A> = Tr[A]; v3 provides traceC(A) directly.
//
// See the KPMResult/ExcitedResult/TimeEvolutionResult comments in
// mpscpp2/chain_session.h for what each struct mirrors on the old
// file-based backend -- unchanged here, not repeated.
struct KPMResult
    {
    std::vector<std::complex<double>> moments;
    double emin, emax, scale;
    int num_polynomials;
    };

// Diagnostics from Chain::kpm_energy_truncate() -- see that method's own
// comment for what each field means (Eqs. 40-41 of the paper cited
// there).
struct KPMTruncStats
    {
    double avg_truncated_weight;
    double state_change_norm;
    };

struct ExcitedResult
    {
    std::vector<double> energies;
    std::vector<double> fluctuations;
    std::vector<MPS> wavefunctions;
    };

struct TimeEvolutionResult
    {
    std::vector<std::complex<double>> correlator;
    MPS final_wf;
    };

struct NHDMRGResult
    {
    std::complex<double> energy;
    MPS psil, psir;
    };

// C++ port of pyitensor/idmrg.py's IDMRGResult -- see
// Chain::idmrg_ground_state's own doc comment for the algorithm and the
// scope this port covers (ground-state energy density only for now,
// static correlators not yet ported).
struct IdmrgResult
    {
    double density = 0.0; // converged (or best-so-far) energy per site
    bool converged = false;
    int niter_done = 0; // macro-iterations actually run
    };

// C++ port of pyitensor/vumps.py's VUMPSResult -- the converged (or
// best-effort) mixed-gauge ground state at a fixed bond dimension D, in
// the thermodynamic limit. See Chain::vumps_ground_state's own doc
// comment (and this file's own "-- VUMPS / excitation ansatz private
// helpers --" section, near idmrg_make_W_end) for why this whole
// feature -- unlike idmrg_ground_state above -- is implemented as plain
// dense row-major arrays closed over LAPACK (itensor::zgeev_wrapper/
// zheev_wrapper/zgesv_wrapper/zgesvd_wrapper) rather than ITensor
// tensor-network objects: d_g (the grouped supersite's own physical
// dimension) and D are always small enough here (n_uc<=2, reach<=1
// bonds -- see vumps_check_reach_one) that this is both simpler and
// far less risky than re-deriving VUMPS's already extremely subtle
// fixed-point/gauge bookkeeping (see pyitensor/vumps.py's and
// pyitensor/idmrg_excitations.py's own module docstrings -- eight
// independent investigation passes were needed to get the tangent-space
// excitation ansatz right against ITensor-style Index-based tensors the
// first time) against ITensor v3's Index/IndexSet machinery a second
// time. Only the ground-state energy density is exposed here -- no
// vev/correlator (VUMPSResult has no per-sublattice U_list, same
// documented limitation pyitensor.vumps.VUMPSResult already carries).
struct VumpsResult
    {
    double e0 = 0.0; // converged (or best-so-far) energy per PHYSICAL site
    bool converged = false;
    int niter_done = 0; // VUMPS outer iterations actually run, at the FINAL bond dimension
    double gauge_mismatch = 0.0; // final ||AC-AL@C||+||AC-C@AR||, normalized -- see vumps_gauge_mismatch
    };

// A converged/summed VUMPS-mixed-gauge uniform iMPS with no ground-state-
// specific bookkeeping -- C++ analogue of pyitensor/vumps.py's own
// UniformMPS, returned by Chain::vumps_apply_mpo/vumps_imps_sum. Load into
// a Chain's own snapshot via Chain::vumps_load_uniform_state to make
// vumps_onsite_expectation/vumps_two_point_correlator (or a further
// apply_mpo/imps_sum call) see it -- see those methods' own comments.
struct VumpsUniformResult
    {
    int D = 0, d_g = 0;
    std::vector<Cplx> AL, AR, C, AC; // (D,d_g,D)/(D,d_g,D)/(D,D)/(D,d_g,D), row-major
    Cplx eta = Cplx(0,0); // norm diagnostic, see vx_canonicalize_n1's own comment
    };

// -- Internal (never exposed to Python) VUMPS/excitation-ansatz plumbing
// structs -- see Chain's own "-- VUMPS / excitation ansatz private
// helpers --" section (near idmrg_make_W_end) for how each is used; kept
// as plain top-level structs, exactly like IdmrgChan/IdmrgBond/IdmrgOnsite/
// IdmrgAutomatonRow above, rather than nested classes.

// One reach-1 bond's own automaton content after grouping: mat_a (the
// S->pending transition, applied at the bond's "earlier" site) and mat_b
// (pending->F, "later" site), both (d_g,d_g) row-major, M[i,o] convention
// -- C++ analogue of pyitensor idmrg_excitations._pending_channels' own
// per-channel tuple.
struct PendingChan
    {
    std::vector<Cplx> mat_a, mat_b;
    };

// A PendingChan plus its own precomputed one-more-site bond content
// (Lvec_a/Rvec_b, each (D,D)) -- C++ analogue of pyitensor/vumps.py's own
// _precompute_bond_environments per-channel tuple.
struct PendingBondEnv
    {
    std::vector<Cplx> mat_a, mat_b, Lvec_a, Rvec_b;
    };

// A VUMPS starting point (AL,AR: (D,d_g,D); C: (D,D)) -- see
// Chain::vumps_random_init/vumps_grow_init.
struct VumpsInit
    {
    std::vector<Cplx> AL, AR, C;
    };

// One full vumps_single_run() attempt's own result.
struct VumpsRunResult
    {
    std::vector<Cplx> AL, AR, C, GL, GR;
    double e_cell = 0.0;
    bool converged = false;
    int niter = 0;
    double mismatch = 0.0;
    };

// One outer-iteration's own environment build (Chain::vumps_build_environments).
struct VumpsEnv
    {
    std::vector<Cplx> GL, GR;
    double e_cell = 0.0;
    std::vector<PendingBondEnv> bond_envs;
    };

// A channel-resolved (D,D)-per-channel map {S,F,pending[0],pending[1],...}
// -- C++ analogue of pyitensor idmrg_excitations._build_GBL/_build_GBR's
// own dict return value, kept as an explicit struct here (no dict
// container needed: only ever S/F/pending, exactly the channel layout
// vumps_pending_channels already enumerates).
struct VumpsChannelMap
    {
    std::vector<Cplx> S, F;
    std::vector<std::vector<Cplx>> pending;
    };

// Chain::td_dynamical_correlator_window's own return value: S(x,t) =
// <psi0|A_x exp(-iHt)B_0|psi0> for every t=ts[it], x=xs[ix] -- S is flat
// row-major, S[it*xs.size()+ix]. The k,omega Fourier transform is done on
// the Python side (infinitechain.py), reusing timedependent.py's own
// _fourier_transform_correlator exactly as the pyitensor backend already
// does for this same submode -- see td_dynamical_correlator_window's own
// comment for why.
struct TdWindowResult
    {
    std::vector<double> ts;
    std::vector<int> xs;
    std::vector<std::complex<double>> S;
    };

// A finite-state-automaton "channel" at one sublattice position, in the
// exact sense idmrg.py's own _build_automaton docstring documents: "S"
// (start, kind=0) self-loops via Id and is the only channel a new 2-site
// term may start from; "F" (accumulator, kind=1) self-loops via Id
// (picking up onsite terms) and receives completed pending channels;
// kind=2 is a "pending" channel for a still-open 2-site term, identified
// by (bond_index, r) exactly as in the Python port.
struct IdmrgChan
    {
    int kind; // 0=S, 1=F, 2=pending
    int bond_index; // valid iff kind==2
    int r; // valid iff kind==2
    };

// A classified 2-site term (idmrg.py's "bonds" dict) -- mat_a/mat_b are
// dense, row-major (d*d)-sized operator matrices in the same (in,out)
// convention as ITensor's own SiteSet::op() (see idmrg_op_dense).
// Fermionic terms (carry_ferm) are not supported by this C++ port yet
// (always false here) -- see Chain::idmrg_ground_state's own comment.
struct IdmrgBond
    {
    int rel_a, rel_b;
    std::vector<Cplx> mat_a, mat_b;
    bool carry_ferm;
    Cplx coef;
    };

// A classified 1-site (onsite) term (idmrg.py's "onsite" list entries).
struct IdmrgOnsite
    {
    int rel;
    Cplx coef;
    std::vector<Cplx> mat;
    };

// The precomputed, dense transition-matrix content for one sublattice
// position p's automaton tensor -- idmrg.py's own W_bulk[p], but kept as
// plain data (never a persistent ITensor, see idmrg_ground_state's own
// comment for why) so a fresh ITensor with fresh Index objects can be
// built from it on demand. flat[(li*right_n+ri)*d*d + s_in*d + s_out] is
// the (s_in,s_out) matrix element of the li->ri channel transition (all
// zero, i.e. absent from the automaton, unless idmrg_build_row's own
// transition rules set it).
struct IdmrgAutomatonRow
    {
    int p, d, left_n, right_n;
    std::vector<Cplx> flat;
    };

class Chain
    {
    public:
    explicit Chain(std::vector<int> const& site_types)
        : sites_(SpinX(site_types))
        { }

    int
    num_sites() const { return sites_.length(); }

    // Local Hilbert-space dimension at a (1-based) site, e.g. for reshaping
    // reduced_dm()'s flat output without any float-precision guessing.
    int
    site_dim(int site) const { return dim(sites_.si(site)); }

    MPS
    random_mps() const { return default_mps(); }

    void
    set_sweep_params(int maxm, int nsweeps, double cutoff, double noise)
        {
        maxm_ = maxm;
        nsweeps_ = nsweeps;
        cutoff_ = cutoff;
        noise_ = noise;
        }

    void
    set_mpomaxm(int mpomaxm) { mpomaxm_ = mpomaxm; }

    void
    set_verbose(bool verbose) { verbose_ = verbose; }

    void
    set_hamiltonian(std::vector<MOTerm> const& terms)
        {
        H_ = build_mpo(sites_,terms,mpomaxm_);
        have_H_ = true;
        have_wf0_energy_ = false; // any cached energy is now stale
        have_bandwidth_min_ = false; // ...and so is any cached bandwidth
        have_bandwidth_max_ = false;
        }

    double
    gs_energy(bool skip_dmrg=false)
        {
        if (!have_H_) Error("Chain::gs_energy called before set_hamiltonian");
        if (skip_dmrg && have_wf0_ && have_wf0_energy_) return wf0_energy_;
        if (!have_wf0_)
            {
            wf0_ = default_mps();
            have_wf0_ = true;
            }
        auto sweeps = make_sweeps();
        double energy = dmrg(wf0_,H_,sweeps,dmrg_args());
        wf0_energy_ = energy;
        have_wf0_energy_ = true;
        return energy;
        }

    MPS const&
    gs_wavefunction() const
        {
        if (!have_wf0_) Error("Chain::gs_wavefunction called before gs_energy");
        return wf0_;
        }

    void
    set_wavefunction(MPS const& wf)
        {
        wf0_ = wf;
        have_wf0_ = true;
        have_wf0_energy_ = false; // energy no longer matches wf0_
        }

    // Generalized-eigenvalue DMRG: solves H|psi>=lambda*A|psi> for a
    // Hermitian positive-definite metric operator A (A=identity reduces
    // this exactly to plain gs_energy()) via a self-consistent Lagrange-
    // multiplier iteration -- a line-for-line port of pyitensor/dmrg.py's
    // dmrg_generalized() against ITensor v3's own dmrg()/Sweeps/sum()
    // instead of the hand-rolled two-site sweep pyitensor uses (see that
    // function's own docstring for the full derivation). Minimizing
    // <psi|H|psi> subject to the metric normalization <psi|A|psi>=1 has
    // stationarity condition (H-lambda*A)|psi>=mu|psi> for Lagrange
    // multiplier lambda -- the *ordinary* (plain-normalized) eigenproblem
    // of the shifted operator H-lambda*A, exactly what a standard DMRG
    // sweep already finds; at mu=0 this is precisely H|psi>=lambda*A|psi>.
    // Each outer iteration therefore (i) builds the MPO H-lambda*A from
    // the current lambda estimate, (ii) runs *one* ordinary DMRG sweep
    // against it (a length-1 Sweeps object, reusing wf0_ as the warm
    // start every iteration -- ITensor's own dmrg() always starts from
    // whatever MPS it's handed, exactly like gs_energy()'s own re-use of
    // wf0_ across calls), then (iii) updates lambda to the freshly-swept
    // state's generalized Rayleigh quotient <psi|H|psi>/<psi|A|psi>. One
    // outer iteration per nsweeps_, with noise halved off partway through
    // exactly like nhdmrg()'s own per-sweep loop above (this can't just
    // hand the whole make_sweeps() schedule to one dmrg() call the way
    // gs_energy() does, since lambda must be updated between sweeps).
    // H-lambda*A is rebuilt from the *original* H_, A each outer
    // iteration (never accumulated), so its own bond dimension never
    // grows across outer iterations -- always at most
    // bondDim(H_)+bondDim(A), summed with this file's usual
    // "MaxDim",mpomaxm_,"Cutoff",cutoff_ convention for MPO+MPO sums
    // (scaled_hamiltonian()'s shift_mpo sum, sum_mpo(), ...). Measured
    // directly: forcing an exact (Cutoff=0) sum here instead -- closer
    // to pyitensor/dmrg.py's own mps_sum(H,A*(-lam)), which does default
    // to a lossless sum -- costs a *lot* more than it buys: on an 8-site
    // benchmark case this made every outer sweep run against a
    // meaningfully higher-bond-dimension Heff (AutoMPO's own term-built
    // H_ is not bond-dimension-minimal, so cutoff_'s ordinary truncation
    // is doing real compression here, not just absorbing floating-point
    // noise) for a ~30x slowdown (0.13s -> 4.2s) and *no* accuracy
    // improvement (the exact-sum run's error against the ED ground truth
    // was not smaller). cutoff_ defaults tiny (1e-12) and mpomaxm_ large,
    // so this convention already agrees with pyitensor to near machine
    // precision in practice -- see
    // examples/groundstate/dmrg_generalized_benchmark. lam0 (optional,
    // NaN meaning "unset"): starting lambda estimate; NaN (the default)
    // instead seeds it from wf0_'s own current generalized Rayleigh
    // quotient.
    //
    // H_ is assumed Hermitian (checked by groundstate.py's Python-level
    // wrapper before this is ever called, not here); A is assumed
    // Hermitian positive definite (not checked at all, same
    // documented-but-unenforced precondition arpacktk.py's
    // mpsiram_generalized has for its own M).
    // ITensor v3's dmrg() aborts the whole process (SIGABRT, "LocalOp is
    // default constructed") rather than raising a catchable exception
    // for chains with fewer than 3 sites (see mode.py's own get_mode(),
    // which exists specifically to route gs_energy()/every other DMRG
    // entry point around this by falling back to ED) -- there is no ED
    // fallback for this method, so that case is rejected explicitly
    // below instead of silently crashing the interpreter.
    double
    gs_energy_generalized(std::vector<MOTerm> const& terms_a,
                           double lam0=std::numeric_limits<double>::quiet_NaN())
        {
        if (!have_H_) Error("Chain::gs_energy_generalized called before set_hamiltonian");
        if (sites_.length()<3)
            Error("Chain::gs_energy_generalized: ITensor v3's two-site dmrg() "
                  "aborts the whole process for chains shorter than 3 sites "
                  "(see mode.py's own itensor_version==3 guard) -- use "
                  "itensor_version=\"python\" for short chains instead");
        auto A = build_mpo(sites_,terms_a,mpomaxm_);
        if (!have_wf0_)
            {
            wf0_ = default_mps();
            have_wf0_ = true;
            }
        double lam = lam0;
        if (std::isnan(lam))
            {
            double a0 = innerC(wf0_,A,wf0_).real();
            double h0 = innerC(wf0_,H_,wf0_).real();
            lam = (std::abs(a0)>1e-14) ? h0/a0 : 0.0;
            }
        for (int sw=1;sw<=nsweeps_;++sw)
            {
            auto Heff = sum(H_,(-lam)*A,{"MaxDim",mpomaxm_,"Cutoff",cutoff_});
            auto sweeps1 = Sweeps(1);
            sweeps1.maxdim() = maxm_;
            sweeps1.cutoff() = cutoff_;
            sweeps1.noise() = (sw<=nsweeps_/2) ? noise_ : 0.0;
            dmrg(wf0_,Heff,sweeps1,dmrg_args());
            double a_psi = innerC(wf0_,A,wf0_).real();
            if (!(std::abs(a_psi)>1e-14))
                // "!(>tol)", not "<tol": a NaN a_psi fails *both*
                // comparisons, so the more obvious "<tol" form would let
                // a NaN silently slip past this guard -- found via code
                // review, same fix applied to pyitensor/dmrg.py's
                // identical guard.
                Error("Chain::gs_energy_generalized: <psi|A|psi> collapsed to "
                      "~0 (or NaN) mid-iteration (A may not be positive "
                      "definite, or the sweep drove psi toward A's "
                      "near-null-space) -- cannot form a meaningful "
                      "generalized Rayleigh quotient");
            double h_psi = innerC(wf0_,H_,wf0_).real();
            lam = h_psi/a_psi;
            }
        have_wf0_energy_ = false; // stale: wf0_ no longer optimizes plain <H_>
        return lam;
        }

    ExcitedResult
    excited_states(int n, double scale_lagrange=1.0, bool do_gram_schmidt=false)
        {
        if (!have_H_) Error("Chain::excited_states called before set_hamiltonian");
        if (!have_wf0_) gs_energy(); // ensure a ground state is available
        auto sweeps = make_sweeps();
        std::vector<MPS> wfs;
        auto psi0 = wf0_;
        psi0.normalize();
        wfs.push_back(psi0);
        double weight = bandwidth()*scale_lagrange;
        for (int i=1;i<n;i++)
            {
            MPS psi1 = default_mps();
            auto args = dmrg_args();
            args.add("Weight",weight);
            dmrg(psi1,H_,wfs,sweeps,args);
            psi1.normalize();
            wfs.push_back(psi1);
            }
        if (do_gram_schmidt) wfs = gram_schmidt(wfs);
        ExcitedResult out;
        for (auto& wf : wfs)
            {
            out.fluctuations.push_back(energy_fluctuation(wf,H_));
            out.energies.push_back(innerC(wf,H_,wf).real());
            out.wavefunctions.push_back(wf);
            }
        return out;
        }

    // Non-Hermitian DMRG (NH-DMRG), mpscpp3-only: variationally optimizes a
    // biorthogonal left/right eigenpair (<psil|, |psir>) of a non-Hermitian
    // H, targeting the eigenvalue with the smallest real part (dmrgpy's
    // "ground state" convention for non-Hermitian Hamiltonians, matching
    // mpsalgebra's Arnoldi mode="GS"). Port of ITensorNHDMRG.jl
    // (https://github.com/tipfom/ITensorNHDMRG.jl) in its default
    // configuration: the "onesided" local solver (independent Arnoldi
    // solves of A|x>=lambda|x> on the right block and
    // Adag|y>=conj(lambda)|y> on the left block -- ordering by real part
    // is identical for both spectra, so the two solves target the same
    // eigenpair) combined with the "fidelity" truncation of Yamamoto et
    // al., Phys. Rev. B 105, 205125: both MPS are truncated with the
    // *same* isometry, obtained from the hermitian average
    // rho=(rho_l+rho_r)/2 of the left and right two-site reduced density
    // matrices. That shared isometry is what keeps psil and psir on
    // identical site *and* link Index objects throughout the sweep, which
    // in turn makes the projected two-site eigenproblem's input and
    // output spaces literally the same index space (the alternative
    // "biorthoblock" truncation of arXiv:2401.15000 needs a non-unitary
    // Schur/Sylvester transform per bond and is not ported). The adjoint
    // Hamiltonian arrives as its own term list, built on the Python side
    // via MultiOperator.get_dagger(), rather than reconstructing
    // dag(swapPrime(H)) index-by-index here -- it reuses build_mpo
    // unchanged and sidesteps v3's prime-level conventions entirely.
    // Starting both MPS from the *same* normalized state makes the pair
    // trivially biorthonormal (<psil|psir>=1 exactly), replacing the
    // reference's initial biorthogonalize! pass (a no-op for identical
    // inputs); the two states then diverge from the first bond update on.
    // One full forward-then-backward NH-DMRG sweep against a fixed
    // (H, HA=H^dagger) pair, mutating psil/psir in place. Builds its own
    // two-sided environments from scratch every call (self-contained,
    // like dmrg.py's own _dmrg_one_sweep/pyitensor's
    // nhdmrg.py::_nhdmrg_one_sweep) rather than assuming any carried over
    // from a previous call -- needed by nhdmrg_generalized() below, whose
    // (H,HA) change every outer iteration, so environments from a
    // previous pair would be wrong. For nhdmrg()'s own multi-sweep loop
    // (same H,HA throughout), rebuilding every sweep instead of reusing
    // the previous sweep's own end-of-sweep environments is numerically
    // equivalent (both are the same deterministic function of the
    // current psil/psir, which don't change between sweep calls) at the
    // cost of one redundant O(N) environment build per sweep -- cheap
    // next to the O(N) per-bond Arnoldi solves each sweep already pays.
    // Returns the final bond's local eigenvalue (nhdmrg()'s own
    // within-sweep tie-break anchor; see Sel::SRTieBreak), starting the
    // anchor from (energy_hint, have_energy_hint) instead of always
    // (0, false) so a caller can continue anchoring across sweeps that
    // share the same (H, HA) -- nhdmrg() itself does,
    // nhdmrg_generalized() below deliberately does not (see its own
    // comment).
    Cplx
    nhdmrg_one_sweep(MPS& psil, MPS& psir, MPO const& H, MPO const& HA,
                      int krylovdim, int restarts, double noise,
                      Cplx energy_hint, bool have_energy_hint)
        {
        int N = sites_.length();
        Cplx energy = energy_hint;
        bool have_energy = have_energy_hint;
        // Two-sided environments (bra psil', ket psir for H; roles swapped
        // for the adjoint problem), one tensor per site; a
        // default-constructed ITensor stands for the scalar 1 past the
        // chain ends. L[i] covers sites 1..i, R[i] covers i..N.
        std::vector<ITensor> Lh(N+2), Rh(N+2), La(N+2), Ra(N+2);
        auto make_env = [&](std::vector<ITensor>& E, MPO const& W,
                            MPS const& ket, MPS const& bra, int i, int prev)
            {
            auto T = E[prev] ? E[prev]*ket.A(i) : ket.A(i);
            T *= W.A(i);
            T *= dag(prime(bra.A(i)));
            E[i] = T;
            };
        for (int i=N;i>=3;--i)
            {
            make_env(Rh,H,psir,psil,i,i+1);
            make_env(Ra,HA,psil,psir,i,i+1);
            }
        // Projected two-site operator at bond b: env legs and MPO output
        // legs all come out primed, so a single noPrime() maps the result
        // back onto the input's own index space (valid because psil and
        // psir share link indices -- see the fidelity note above).
        auto apply_proj = [&](MPO const& W, std::vector<ITensor> const& L,
                              std::vector<ITensor> const& R, int b, ITensor v)
            {
            if (L[b-1]) v *= L[b-1];
            v *= W.A(b);
            v *= W.A(b+1);
            if (R[b+2]) v *= R[b+2];
            v.noPrime();
            return v;
            };
        for (int ha=1;ha<=2;++ha)
        for (int bi=0;bi<N-1;++bi)
            {
            int b = (ha==1) ? 1+bi : N-1-bi;
            auto thl = psil.A(b)*psil.A(b+1);
            auto thr = psir.A(b)*psir.A(b+1);
            // smallest real part, with Re-degenerate candidates
            // tie-broken toward the previous bond's eigenvalue (see
            // arnoldi_smallest_real's Sel comment)
            auto er_thr = arnoldi_smallest_real(
                [&](ITensor const& v) { return apply_proj(H,Lh,Rh,b,v); },
                thr,krylovdim,restarts,
                have_energy ? Sel::SRTieBreak : Sel::SR,energy);
            // anchor the adjoint solve to the right solve's eigenvalue
            auto el_thl = arnoldi_smallest_real(
                [&](ITensor const& v) { return apply_proj(HA,La,Ra,b,v); },
                thl,krylovdim,restarts,
                Sel::Closest,std::conj(er_thr.first));
            energy = er_thr.first;
            have_energy = true;
            thr = er_thr.second;
            thl = el_thl.second;
            thl /= norm(thl);
            thr /= norm(thr);
            // rescale so <thl|thr> = 1 (split between the two states,
            // with the reference's separate real branch: for a real
            // negative overlap the complex branch's sqrt(conj(z))
            // lands on the wrong side of std::sqrt's branch cut and
            // would leave the overlap at -1 instead of +1)
            auto ov = eltC(dag(thl)*thr);
            if (std::abs(ov)>1e-12)
                {
                if (std::abs(ov.imag())<1e-14*std::abs(ov))
                    {
                    double sq = std::sqrt(std::abs(ov.real()));
                    thl /= sq;
                    thr /= (ov.real()<0 ? -sq : sq);
                    }
                else
                    {
                    thl /= std::sqrt(std::conj(ov));
                    thr /= std::sqrt(ov);
                    }
                }
            // fidelity truncation: hermitian average of the left and
            // right reduced density matrices over the indices kept on
            // the orthogonality-moving side of the bond
            auto keep = (ha==1) ? commonInds(psil.A(b),thl)
                                : commonInds(psil.A(b+1),thl);
            auto rl = thl;
            auto rr = thr;
            for (auto const& I : keep) { rl = prime(rl,I); rr = prime(rr,I); }
            rl *= dag(thl);
            rr *= dag(thr);
            auto rho = 0.5*(rl+rr);
            if (noise>0)
                {
                // reference's noiseterm(): cross term between the left
                // and right blocks, hermitized here since rho feeds a
                // hermitian eigensolver (the Julia fidelity path adds
                // it unhermitized and relies on eigen(ishermitian=true)
                // only reading one triangle)
                ITensor X = (ha==1) ? H.A(b) : H.A(b+1);
                if (ha==1 && Lh[b-1]) X = Lh[b-1]*X;
                if (ha==2 && Rh[b+2]) X = X*Rh[b+2];
                auto nt = (X*thl)*dag(noPrime(X*thr));
                rho += (noise/2.0)*(nt+dag(swapPrime(nt,0,1)));
                }
            ITensor U,D;
            diagPosSemiDef(rho,U,D,{"MaxDim",maxm_,"Cutoff",cutoff_,
                                    "Tags","Link,NH"});
            if (ha==1)
                {
                psil.ref(b) = U;
                psil.ref(b+1) = dag(U)*thl;
                psir.ref(b) = U;
                psir.ref(b+1) = dag(U)*thr;
                if (b<N-1)
                    {
                    make_env(Lh,H,psir,psil,b,b-1);
                    make_env(La,HA,psil,psir,b,b-1);
                    }
                }
            else
                {
                psil.ref(b+1) = U;
                psil.ref(b) = thl*dag(U);
                psir.ref(b+1) = U;
                psir.ref(b) = thr*dag(U);
                if (b>1)
                    {
                    make_env(Rh,H,psir,psil,b+1,b+2);
                    make_env(Ra,HA,psil,psir,b+1,b+2);
                    }
                }
            }
        return energy;
        }

    // Non-Hermitian DMRG (NH-DMRG): the algorithm itself -- ITensorNHDMRG.jl
    // port, "onesided" local solver, "fidelity" truncation, the
    // biorthogonal-start rationale -- is documented on nhdmrg_one_sweep()
    // just above (this method is now a thin per-sweep driver around it, see
    // its own comment for why nhdmrg_one_sweep rebuilds environments fresh
    // every call rather than reusing them across sweeps here).
    NHDMRGResult
    nhdmrg(std::vector<MOTerm> const& terms_h,
           std::vector<MOTerm> const& terms_hadj,
           int krylovdim=20, int restarts=2)
        {
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto HA = build_mpo(sites_,terms_hadj,mpomaxm_);
        // always a fresh random start (never wf0_): the non-Hermitian
        // energy is not a variational bound, so a rare stalled run can
        // only be detected by the caller's eigen-residual check and cured
        // by re-running -- which requires every run to draw its own
        // random initial state (see nhdmrg.py's retry loop)
        MPS psir = default_mps();
        psir.position(1);
        psir.normalize();
        MPS psil = psir;
        Cplx energy = 0;
        bool have_energy = false;
        for (int sw=1;sw<=nsweeps_;++sw)
            {
            // mirrors make_sweeps(): noise only in the first half-schedule
            double noise = (sw<=nsweeps_/2) ? noise_ : 0.0;
            energy = nhdmrg_one_sweep(psil,psir,H,HA,krylovdim,restarts,
                    noise,energy,have_energy);
            have_energy = true;
            if (verbose_)
                printfln("NH-DMRG sweep %d energy = %.12f + %.12f i",
                         sw,energy.real(),energy.imag());
            }
        // both sweeps end at bond 1 with ortho "right": center at site 1
        psil.leftLim(0); psil.rightLim(2);
        psir.leftLim(0); psir.rightLim(2);
        NHDMRGResult out;
        // definitive energy: the biorthogonal Rayleigh quotient of the
        // final pair (the last local eigenvalue tracks it but lags half a
        // bond update behind)
        out.energy = innerC(psil,H,psir)/innerC(psil,psir);
        out.psil = psil;
        out.psir = psir;
        return out;
        }

    // Non-Hermitian generalized-eigenvalue NH-DMRG: finds the smallest-
    // real-part lambda solving H|psi_R>=lambda*A|psi_R> (equivalently
    // H^dagger|psi_L>=conj(lambda)*A|psi_L>, since A is Hermitian) for a
    // possibly non-Hermitian H and a Hermitian positive-definite metric A
    // -- the non-Hermitian counterpart of gs_energy_generalized() above,
    // and a line-for-line port of pyitensor/nhdmrg.py's
    // nhdmrg_generalized() against this file's own nhdmrg_one_sweep
    // instead of the hand-rolled Python one. Same self-consistent
    // Lagrange-multiplier trick as gs_energy_generalized(), generalized
    // to a complex lambda and biorthogonal (rather than plain)
    // expectation values: minimizing (in the NH-DMRG sense -- smallest
    // real part, not a variational bound) subject to the metric
    // biorthogonal normalization <psi_L|A|psi_R>=1 gives stationarity
    // condition (H-lambda*A)|psi_R>=mu|psi_R>,
    // (H-lambda*A)^dagger|psi_L>=conj(mu)|psi_L> -- the *ordinary*
    // (<psi_L|psi_R>=1-normalized) NH-DMRG eigenproblem of the shifted
    // pair (H-lambda*A, HA-conj(lambda)*A) (note the adjoint shift uses
    // conj(lambda), not lambda: (lambda*A)^dagger=conj(lambda)*A^dagger=
    // conj(lambda)*A since A is Hermitian) -- exactly what one ordinary
    // NH-DMRG sweep already finds. At mu=0 this is precisely
    // H|psi_R>=lambda*A|psi_R>, so each outer iteration here (i) builds
    // Heff=H-lambda*A and HAeff=HA-conj(lambda)*A from the current lambda
    // estimate, (ii) runs one ordinary NH-DMRG sweep (nhdmrg_one_sweep,
    // unmodified) against them, then (iii) updates lambda to the
    // freshly-swept pair's generalized biorthogonal Rayleigh quotient
    // <psi_L|H|psi_R>/<psi_L|A|psi_R>. A=identity reduces this exactly to
    // plain nhdmrg().
    //
    // Each outer sweep's own within-sweep SRTieBreak anchor
    // (nhdmrg_one_sweep's energy_hint/have_energy_hint) restarts fresh
    // (target 0, no anchor) rather than carrying over from the previous
    // outer iteration -- same rationale as the Python port: (H,HA)
    // changes every outer iteration (unlike nhdmrg()'s own sweep loop,
    // which reuses the same fixed H/HA throughout and so *does* carry the
    // anchor across sweeps), so an anchor built for a different lambda's
    // shifted spectrum has no particular reason to still be meaningful
    // for the next one.
    //
    // No ITensor v3 dmrg() call is made anywhere in nhdmrg_one_sweep (the
    // two-site sweep is hand-rolled directly against arnoldi_smallest_real
    // and manual ITensor contractions, unlike gs_energy_generalized()'s
    // dmrg() call), so the short-chain SIGABRT crash gs_energy_generalized
    // guards against does not apply here -- confirmed directly, a 2-site
    // chain runs this method without aborting.
    //
    // lam0 (optional, NaN real part meaning "unset"): starting lambda
    // estimate; NaN (the default) instead seeds it from the fresh random
    // start's own generalized biorthogonal Rayleigh quotient.
    NHDMRGResult
    nhdmrg_generalized(std::vector<MOTerm> const& terms_h,
                        std::vector<MOTerm> const& terms_hadj,
                        std::vector<MOTerm> const& terms_a,
                        int krylovdim=20, int restarts=2,
                        Cplx lam0=Cplx(std::numeric_limits<double>::quiet_NaN(),0.0))
        {
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto HA = build_mpo(sites_,terms_hadj,mpomaxm_);
        auto A = build_mpo(sites_,terms_a,mpomaxm_);
        MPS psir = default_mps();
        psir.position(1);
        psir.normalize();
        MPS psil = psir;

        Cplx lam = lam0;
        if (std::isnan(lam.real()))
            {
            Cplx a0 = innerC(psil,A,psir);
            Cplx h0 = innerC(psil,H,psir);
            lam = (std::abs(a0)>1e-14) ? h0/a0 : Cplx(0,0);
            }

        for (int sw=1;sw<=nsweeps_;++sw)
            {
            double noise = (sw<=nsweeps_/2) ? noise_ : 0.0;
            auto Heff = sum(H,(-lam)*A,{"MaxDim",mpomaxm_,"Cutoff",cutoff_});
            auto HAeff = sum(HA,(-std::conj(lam))*A,
                              {"MaxDim",mpomaxm_,"Cutoff",cutoff_});
            nhdmrg_one_sweep(psil,psir,Heff,HAeff,krylovdim,restarts,noise,
                              Cplx(0,0),false);
            Cplx a_psi = innerC(psil,A,psir);
            if (!(std::abs(a_psi)>1e-14))
                // "!(>tol)", not "<tol": see gs_energy_generalized's
                // identical guard for why the negated form is needed to
                // also catch NaN.
                Error("Chain::nhdmrg_generalized: <psi_L|A|psi_R> collapsed "
                      "to ~0 (or NaN) mid-iteration (A may not be positive "
                      "definite, or the sweep drove the biorthogonal pair "
                      "toward A's near-null-space) -- cannot form a "
                      "meaningful generalized Rayleigh quotient");
            Cplx h_psi = innerC(psil,H,psir);
            lam = h_psi/a_psi;
            if (verbose_)
                printfln("NH-DMRG-generalized sweep %d lambda = %.12f + %.12f i",
                         sw,lam.real(),lam.imag());
            }
        psil.leftLim(0); psil.rightLim(2);
        psir.leftLim(0); psir.rightLim(2);
        NHDMRGResult out;
        out.energy = lam;
        out.psil = psil;
        out.psir = psir;
        return out;
        }

    std::complex<double>
    vev(std::vector<MOTerm> const& terms, MPS const& wf, int npow=1)
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto psi = wf;
        psi /= sqrt(innerC(psi,psi).real()); // normalize
        Cplx c = 0;
        if (npow==1) { c = innerC(psi,A,psi); }
        if (npow>1)
            {
            auto psi1 = psi;
            for (int i=0;i<npow-1;i++)
                psi1 = apply_mpo(A,psi1,{"MaxDim",maxm_,"Cutoff",cutoff_});
            c = innerC(psi,A,psi1);
            }
        return c;
        }

    MPS
    apply_operator(std::vector<MOTerm> const& terms, MPS const& wf)
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        return apply_mpo(A,wf,args);
        }

    std::complex<double>
    overlap_mps(MPS const& wf1, MPS const& wf2) const
        {
        return innerC(wf1,wf2);
        }

    std::complex<double>
    overlap_aMb(MPS const& wf1, std::vector<MOTerm> const& terms, MPS const& wf2) const
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        return innerC(wf1,A,wf2);
        }

    MPS
    sum_mps(MPS const& wf1, MPS const& wf2) const
        {
        return sum(wf1,wf2,{"MaxDim",maxm_,"Cutoff",cutoff_});
        }

    MPS
    conjugate(MPS const& wf) const
        {
        auto out = wf;
        for (int i=0;i<out.length();++i) out.Aref(i+1).conj();
        return out;
        }

    std::vector<std::complex<double>>
    reduced_dm(MPS const& wf, int site) const
        {
        auto psi = wf;
        psi /= innerC(psi,psi).real(); // normalize (see note in mpscpp2/chain_session.h)
        psi.position(site);
        auto ir = commonIndex(psi.A(site),psi.A(site+1));
        auto rho = psi.A(site)*dag(prime(psi.A(site),sites_.si(site),ir));
        for (int k=site+1;k<=psi.length();++k)
            {
            rho *= psi.A(k);
            rho *= dag(prime(psi.A(k),TagSet("Link")));
            }
        std::vector<std::complex<double>> out;
        auto collect = [&out](Cplx z) { out.push_back(z); };
        rho.visit(collect);
        return out;
        }

    MPS
    exponential_apply(std::vector<MOTerm> const& terms, MPS const& wf,
                       std::complex<double> tau, int nsteps) const
        {
        auto H = build_mpo(sites_,terms,mpomaxm_);
        auto taui = tau/double(nsteps);
        auto expH = custom_exp(H,taui);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        auto psi1 = wf;
        for (int it=1;it<=nsteps;it++) psi1 = apply_mpo(expH,psi1,args);
        return psi1;
        }

    MPO
    build_operator(std::vector<MOTerm> const& terms) const
        {
        return build_mpo(sites_,terms,mpomaxm_);
        }

    MPS
    apply_pure_operator(MPO const& A, MPS const& wf) const
        {
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        return apply_mpo(A,wf,args);
        }

    MPO
    multiply_operators(MPO const& A, MPO const& B) const
        {
        return mult_mpo(A,B);
        }

    // Public wrapper around the private sum_mpo() helper below (mirrors
    // multiply_operators()/mult_mpo() just above), exposed so StaticOperator
    // on the Python side can add two already-built MPOs directly. sum_mpo()
    // itself is just ITensor v3's own sum(MPO,MPO,args) -- a compressed
    // direct sum, algorithmically the same construction as ITensorMPS.jl's
    // `+(::MPO, ::MPO)` (abstractmps.jl's default "densitymatrix" algorithm).
    MPO
    sum_operators(MPO const& A, MPO const& B) const
        {
        return sum_mpo(A,B);
        }

    // Scalar multiple of an already-built MPO -- a plain tensor rescale
    // (multiplies one site's tensor by z), not a contraction, so unlike
    // multiply_operators()/sum_mpo() this doesn't touch bond dimension.
    // Exposed so StaticOperator can implement negation/subtraction on top
    // of sum_operators(), mirroring how Julia's `-(A,B) = +(A,-B)` for
    // MPS/MPO in ITensorMPS.jl's abstractmps.jl reduces to `+` plus a
    // scalar multiple.
    MPO
    scale_operator(MPO const& A, Cplx z) const
        {
        return z*A;
        }

    // Mirrors multmpo.h's trace_mpo_operator() task / operators.h's
    // trace_mpo() (Tr[A] = <Id|A>); v3 provides this directly via traceC()
    // instead of needing to build an explicit identity MPO first.
    std::complex<double>
    trace_operator(MPO const& A) const
        {
        return traceC(A);
        }

    MPO
    hermitian_operator(MPO const& A) const
        {
        auto out = A;
        for (auto j : range1(out.length())) out.Aref(j) = dag(swapPrime(out.A(j),0,1,TagSet("Site")));
        return out;
        }

    std::complex<double>
    overlap_aMb_operator(MPS const& wf1, MPO const& A, MPS const& wf2) const
        {
        return innerC(wf1,A,wf2);
        }

    double
    bond_entropy(MPS const& wf, int b) const
        {
        auto psi = wf;
        psi.position(b);
        ITensor twosite = psi.A(b)*psi.A(b+1);
        auto U = psi.A(b);
        ITensor S,V;
        auto spectrum = svd(twosite,U,S,V);
        double SvN = 0.0;
        for (auto p : spectrum.eigs()) if (p>1E-12) SvN += -p*std::log(p);
        return SvN;
        }

    TimeEvolutionResult
    quench(std::vector<MOTerm> const& terms_h,
           std::vector<MOTerm> const& terms_i,
           std::vector<MOTerm> const& terms_j,
           int nt, double dt, bool fit_td=true)
        {
        if (!have_wf0_) gs_energy();
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = build_ampo(sites_,terms_h);
        ampo += -EGS,"Id",1;
        auto expH = evoloperator(toMPO(ampo),dt);
        auto A1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto A2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (fit_td) psi1 = apply_mpo(expH,psi1,psi1,args);
            else psi1 = apply_mpo(expH,psi1,args);
            psi1.normalize();
            psi1 *= norm0;
            out.correlator.push_back(innerC(psi2,psi1));
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure(std::vector<MOTerm> const& terms_h,
                        std::vector<MOTerm> const& terms_op,
                        MPS const& wf, int nt, double dt, bool fit_td=true)
        {
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        auto ampo = build_ampo(sites_,terms_h);
        auto expH = evoloperator(toMPO(ampo),dt);
        auto A = build_mpo(sites_,terms_op,mpomaxm_);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (fit_td) psi = apply_mpo(expH,psi,psi,args);
            else psi = apply_mpo(expH,psi,args);
            out.correlator.push_back(innerC(psi,A,psi));
            }
        out.final_wf = psi;
        return out;
        }

    // Applies one Taylor-expanded exp(z*H) step (evoloperator() below,
    // see its own comment for the pre-existing z^3/6-uses-H2-not-H3
    // quirk this preserves unchanged) to wf, given an already-built MPO
    // H and a possibly-complex z -- exposed for callers (tdz.py's "TDZ"
    // complex-time-evolution dynamical correlator) that need a
    // per-step-varying complex increment, unlike quench()/
    // evolve_and_measure() above, which use one fixed real dt for their
    // whole internal loop. mpscpp2 (which has no TDVP) uses this same
    // method as its only route to TDZ; here it's a cross-check /
    // non-TDVP alternative to the public tdvp_step() below.
    MPS
    evolve_taylor_step(MPO const& H, MPS const& wf, Cplx z) const
        {
        auto expH = evoloperator(H,z);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        return apply_mpo(expH,wf,args);
        }

    // TDVP counterparts of quench()/evolve_and_measure() above: same
    // physics (real-time evolution under H_, same EGS-shift convention in
    // quench_tdvp() so results are directly comparable to the MPO-Taylor
    // backup), but each per-step evolution is done with two-site TDVP
    // (tdvp_step(), see its comment below) instead of applying the
    // Taylor-expanded evoloperator() MPO. No fit_td flag: TDVP has no
    // MPO-fit variant.
    TimeEvolutionResult
    quench_tdvp(std::vector<MOTerm> const& terms_h,
                std::vector<MOTerm> const& terms_i,
                std::vector<MOTerm> const& terms_j,
                int nt, double dt)
        {
        if (!have_wf0_) gs_energy();
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = build_ampo(sites_,terms_h);
        ampo += -EGS,"Id",1;
        auto Hshift = toMPO(ampo);
        auto A1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto A2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            psi1 = tdvp_step(Hshift,psi1,dt);
            psi1.normalize();
            psi1 *= norm0;
            out.correlator.push_back(innerC(psi2,psi1));
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure_tdvp(std::vector<MOTerm> const& terms_h,
                             std::vector<MOTerm> const& terms_op,
                             MPS const& wf, int nt, double dt)
        {
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto A = build_mpo(sites_,terms_op,mpomaxm_);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            psi = tdvp_step(H,psi,dt);
            out.correlator.push_back(innerC(psi,A,psi));
            }
        out.final_wf = psi;
        return out;
        }

    // GSE counterparts of quench_tdvp()/evolve_and_measure_tdvp() just
    // above: identical setup/measurement, but each per-step evolution is
    // one-site TDVP (tdvp_step(...,1)) preceded by a global_subspace_expand()
    // call for the first gse_sweeps steps -- the Yang-White scheme
    // (arXiv:2005.06104) that lets one-site TDVP's bond dimension keep up
    // with two-site TDVP's own SVD-driven growth. The driver loop stays in
    // C++ (like quench_tdvp/evolve_and_measure_tdvp) rather than in Python,
    // both for consistency with those and to avoid nt Python<->C++ round
    // trips for the (typically large) default nt.
    TimeEvolutionResult
    quench_tdvp_gse(std::vector<MOTerm> const& terms_h,
                     std::vector<MOTerm> const& terms_i,
                     std::vector<MOTerm> const& terms_j,
                     int nt, double dt, int gse_sweeps, int krylov_order,
                     double gse_cutoff)
        {
        if (!have_wf0_) gs_energy();
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = build_ampo(sites_,terms_h);
        ampo += -EGS,"Id",1;
        auto Hshift = toMPO(ampo);
        auto A1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto A2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (it<gse_sweeps)
                psi1 = global_subspace_expand(Hshift,psi1,krylov_order,gse_cutoff,0);
            psi1 = tdvp_step(Hshift,psi1,dt,1);
            psi1.normalize();
            psi1 *= norm0;
            out.correlator.push_back(innerC(psi2,psi1));
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure_tdvp_gse(std::vector<MOTerm> const& terms_h,
                                 std::vector<MOTerm> const& terms_op,
                                 MPS const& wf, int nt, double dt,
                                 int gse_sweeps, int krylov_order,
                                 double gse_cutoff)
        {
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto A = build_mpo(sites_,terms_op,mpomaxm_);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (it<gse_sweeps)
                psi = global_subspace_expand(H,psi,krylov_order,gse_cutoff,0);
            psi = tdvp_step(H,psi,dt,1);
            out.correlator.push_back(innerC(psi,A,psi));
            }
        out.final_wf = psi;
        return out;
        }

    // TEBD counterparts of quench()/quench_tdvp() above: same physics and
    // EGS-shift convention as quench_tdvp() (so results are directly
    // comparable to it and to the MPO-Taylor backup), but each per-step
    // evolution is a 2nd-order-Trotter TEBD step (tebd.h's
    // build_tebd_gates()/tebd_step(), gates built once up front from the
    // bare bond Hamiltonians and reused unchanged every step) instead of
    // TDVP's per-step Krylov exponentiation of an effective Hamiltonian --
    // only valid for a strictly nearest-neighbor terms_h
    // (bond_hamiltonians() throws ITError otherwise; see tebd.h's own
    // comment). Mirrors pyitensor/chain.py's Chain.quench_tebd() (the
    // pure-Python backend's own TEBD, which this file's terms_h/dt
    // conventions match exactly) -- kept as a near-duplicate of
    // quench_tdvp() rather than sharing a driver loop with it, for the
    // same reason quench_tdvp()/quench() themselves aren't unified: each
    // per-step evolution primitive has its own setup (tdvp_step() takes
    // an MPO, tebd_step() takes precomputed BondGates).
    TimeEvolutionResult
    quench_tebd(std::vector<MOTerm> const& terms_h,
                std::vector<MOTerm> const& terms_i,
                std::vector<MOTerm> const& terms_j,
                int nt, double dt)
        {
        if (!have_wf0_) gs_energy();
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto terms_shifted = terms_h;
        MOTerm shift_term;
        shift_term.coef = Cplx(-EGS,0.0);
        shift_term.factors = {{"Id",1}};
        terms_shifted.push_back(shift_term);
        auto h_bonds = bond_hamiltonians(sites_,terms_shifted);
        auto gates = build_tebd_gates(sites_,h_bonds,dt);
        auto A1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto A2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            tebd_step(psi1,gates,cutoff_,maxm_);
            psi1.normalize();
            psi1 *= norm0;
            out.correlator.push_back(innerC(psi2,psi1));
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure_tebd(std::vector<MOTerm> const& terms_h,
                             std::vector<MOTerm> const& terms_op,
                             MPS const& wf, int nt, double dt)
        {
        auto h_bonds = bond_hamiltonians(sites_,terms_h);
        auto gates = build_tebd_gates(sites_,h_bonds,dt);
        auto A = build_mpo(sites_,terms_op,mpomaxm_);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            tebd_step(psi,gates,cutoff_,maxm_);
            out.correlator.push_back(innerC(psi,A,psi));
            }
        out.final_wf = psi;
        return out;
        }

    // One step exp(dt*H) of two-site TDVP (TDVP/tdvp.h), given an
    // already-built MPO H (e.g. from build_operator()) and MPS psi --
    // exposed publicly (moved here from a private helper of the same
    // name used only by quench_tdvp()/evolve_and_measure_tdvp() above)
    // so a caller can drive the evolution one variable-sized step at a
    // time, unlike those two methods, which loop internally over a fixed
    // number of equal, real dt steps and thus can't be reused for a
    // per-step-varying complex time increment (see tdz.py's "TDZ"
    // dynamical-correlator submode, which needs exactly that: the
    // per-step contour increment dz_k = exp(-i*alpha0*f(t_k))*dt varies
    // with t_k). dt may be any complex number: TDVP/README.md documents
    // its own "t" argument (here dt, since tdvp()'s own convention is
    // U=exp(t*H)) as "real, imaginary, or complex" -- real time
    // evolution (dt purely imaginary here, by this method's own
    // -i*dt convention) and complex time evolution (TDZ) share this same
    // code path unchanged, nothing backend-specific needed for either.
    // One Sweeps(1) TDVP "sweep" is a left-to-right pass with a dt/2 step
    // followed by a right-to-left pass with another dt/2 step, i.e.
    // exactly one full step of size dt -- matching how
    // quench()/evolve_and_measure() apply their evoloperator() MPO once
    // per iteration. niter=50 bounds the Lanczos iterations used to
    // solve each local effective TDVP equation (tdvp.h's "MaxIter"); not
    // exposed as a knob since neither the MPO-Taylor path this replaces
    // has an equivalent one. num_center selects one-site (1) or two-site
    // (2, the default -- matches every existing caller, which all predate
    // this parameter) TDVP; "Truncate" follows TDVP/README.md's own
    // stated per-NumCenter default (off for one-site, on for two-site)
    // since one-site TDVP conserves bond dimension exactly and truncating
    // it would only ever shrink it, never grow it -- growth for one-site
    // comes from global_subspace_expand() below instead.
    // niter=50 matches every pre-existing caller (none passed anything
    // else); metts_vev() below is the first caller that needs a
    // different value, hence the added (default-preserving) parameter
    // rather than a second near-duplicate method.
    MPS
    tdvp_step(MPO const& H, MPS psi, Cplx dt, int num_center=2, int niter=50) const
        {
        Cplx t = Cplx(0.0,-1.0)*dt;
        auto sweeps = Sweeps(1);
        sweeps.maxdim() = maxm_;
        sweeps.cutoff() = cutoff_;
        sweeps.niter() = niter;
        tdvp(psi,H,t,sweeps,{"Quiet",!verbose_,"Silent",!verbose_,
                              "NumCenter",num_center,
                              "Truncate",num_center==2,
                              "DoNormalize",true});
        return psi;
        }

    // METTS (Minimally Entangled Typical Thermal States; E.M. Stoudenmire
    // and S.R. White, arXiv:1002.1305) -- a direct port of
    // pyitensor/metts.py's algorithm onto real ITensor v3 (see that
    // module's own docstring for the full derivation this only summarizes
    // here). Samples a Markov chain of classical product states (CPS):
    // each is imaginary-time evolved by beta/2 via tdvp_step() above
    // (dt=Cplx(0,-1)*step -- that method's own -i*dt convention already
    // makes a purely real dt give genuine decay rather than rotation, so
    // no new evolution primitive is needed here), then collapsed back
    // down to a new CPS by collapse_to_cps() (private, below) using
    // diagHermitian() for the per-site basis rotation, alternating the
    // collapse basis through basis_ops from one sample to the next for
    // ergodicity (paper, Sec. 3.3). A plain sample average of
    // <phi|A|phi> over the resulting chain converges to the true thermal
    // average, because this specific Markov chain's stationary
    // distribution already carries the correct Boltzmann weight -- no
    // importance reweighting needed (paper, Eq. (3)-(5)).
    //
    // Returns (means, stderrs), one entry per op in terms_ops: stderrs[k]
    // is a naive i.i.d. estimate over the nsamples retained samples (after
    // nwarmup discarded equilibration steps), so -- since consecutive
    // METTS samples are Markov-correlated -- it's likely optimistic unless
    // dbeta_half_step/nwarmup are generous enough that samples actually
    // decorrelate. All ops in terms_ops are measured on the same sampled
    // Markov chain of METTS states -- the (nwarmup+nsamples) imaginary-
    // time evolutions below are the expensive part, so sharing them across
    // every requested op (rather than resampling once per op) is the
    // whole point of taking a list here instead of a single op (mirrors
    // pyitensor/metts.py's metts_thermal_average, which already worked
    // this way).
    std::pair<std::vector<Cplx>,std::vector<double>>
    metts_vev(std::vector<std::vector<MOTerm>> const& terms_ops, double T, int nsamples,
              int nwarmup, double dbeta_half_step,
              std::vector<std::string> const& basis_ops,
              unsigned long seed, int niter=30) const
        {
        if (!have_H_) Error("Chain::metts_vev called before set_hamiltonian");
        if (T<=0) Error("Chain::metts_vev: T must be > 0");
        if (basis_ops.empty()) Error("Chain::metts_vev: basis_ops must be non-empty");
        if (nsamples<1) Error("Chain::metts_vev: nsamples must be >= 1");
        if (terms_ops.empty()) Error("Chain::metts_vev: terms_ops must be non-empty");
        std::vector<MPO> ops;
        ops.reserve(terms_ops.size());
        for (auto const& terms_op : terms_ops) ops.push_back(build_mpo(sites_,terms_op,mpomaxm_));
        int nops = (int)ops.size();
        double beta_half = 1.0/(2.0*T);
        int nsteps = std::max(1,(int)std::ceil(beta_half/dbeta_half_step));
        Cplx dt = Cplx(0.0,-1.0)*(beta_half/double(nsteps));

        std::mt19937_64 rng(seed);
        int nbasis = (int)basis_ops.size();
        auto eigcache = metts_build_eigcache(basis_ops);
        MPS cps = random_cps(basis_ops[0],rng,eigcache);

        std::vector<std::vector<Cplx>> samples(nops);
        for (auto& s : samples) s.reserve(nsamples);
        int total_iters = nwarmup+nsamples;
        for (int it=0; it<total_iters; ++it)
            {
            MPS phi = cps;
            for (int s=0; s<nsteps; ++s)
                {
                phi = tdvp_step(H_,phi,dt,2,niter);
                phi /= sqrt(innerC(phi,phi).real());
                }
            if (it>=nwarmup)
                {
                for (int k=0; k<nops; ++k) samples[k].push_back(innerC(phi,ops[k],phi));
                }
            cps = collapse_to_cps(phi,basis_ops[it % nbasis],rng,eigcache);
            }

        // nwarmup>=total_iters (e.g. nwarmup alone consuming every
        // iteration) would otherwise leave samples empty and divide by
        // zero below -- not reachable through the Python-facing default
        // kwargs, but not guarded against a caller passing an oversized
        // nwarmup either, so check explicitly rather than silently
        // returning NaN.
        if (samples[0].empty())
            Error("Chain::metts_vev: no samples were retained (nwarmup >= nwarmup+nsamples?)");
        std::vector<Cplx> means(nops,Cplx(0.0,0.0));
        std::vector<double> stderrs(nops,0.0);
        for (int k=0; k<nops; ++k)
            {
            Cplx mean = 0;
            for (auto& v : samples[k]) mean += v;
            mean /= double(samples[k].size());
            double var = 0.0;
            for (auto& v : samples[k])
                {
                double d = std::abs(v-mean);
                var += d*d;
                }
            double stderr_ = samples[k].size()>1
                ? std::sqrt(var/double(samples[k].size()-1)/double(samples[k].size())) : 0.0;
            means[k] = mean;
            stderrs[k] = stderr_;
            }
        return {means,stderrs};
        }

    // Dynamical METTS (real-time finite-T correlators), arXiv:2405.18484
    // (Wang, McClarty, Dankova, Honecker & Wietek, "Spectroscopy and
    // complex-time correlations using minimally entangled typical thermal
    // states"), Sec. II / Algorithm 1 -- a direct port of
    // pyitensor/metts.py's metts_dynamical_correlator() onto real ITensor
    // v3, following the exact same METTS Markov chain metts_vev() above
    // does (imaginary-time tdvp_step + collapse_to_cps), generalized from
    // the static <A>_beta average to the two-time correlator
    // C_AB(t) = <A(t)B>_beta = <e^{iHt} A e^{-iHt} B>_beta: for every
    // retained sample |phi>, set |v(0)>=B|phi>, |w(0)>=|phi>, measure
    // C(t_k)=<w(t_k)|A|v(t_k)>, then real-time-evolve both |v> and |w>
    // independently under H_ via tdvp_step (dt real here -- that
    // method's own -i*dt convention already makes a purely real dt give
    // genuine rotation, i.e. real time, exactly as quench_tdvp() above
    // already relies on). A plain sample average over every retained
    // METTS sample converges to C_AB(t), for the same reason
    // metts_vev()'s own plain average converges to <A>_beta -- see that
    // method's comment above.
    //
    // Returns (means,stderrs), each a length-nt array over
    // t=0,dt,...,(nt-1)*dt -- same (means,stderrs) convention as
    // metts_vev() above, just array-valued per time step instead of
    // per-operator; the ts array itself is left for the Python-facing
    // wrapper to construct from (nt,dt), same as
    // pyitensor.chain.Chain.metts_dynamical_correlator's own split.
    std::pair<std::vector<Cplx>,std::vector<double>>
    metts_dynamical_correlator(std::vector<MOTerm> const& terms_a,
              std::vector<MOTerm> const& terms_b,
              double T, int nt, double dt, int nsamples,
              int nwarmup, double dbeta_half_step,
              std::vector<std::string> const& basis_ops,
              unsigned long seed, int niter=30, int tdvp_niter=50) const
        {
        if (!have_H_) Error("Chain::metts_dynamical_correlator called before set_hamiltonian");
        if (T<=0) Error("Chain::metts_dynamical_correlator: T must be > 0");
        if (basis_ops.empty()) Error("Chain::metts_dynamical_correlator: basis_ops must be non-empty");
        if (nsamples<1) Error("Chain::metts_dynamical_correlator: nsamples must be >= 1");
        if (nt<1) Error("Chain::metts_dynamical_correlator: nt must be >= 1");

        auto A = build_mpo(sites_,terms_a,mpomaxm_);
        auto B = build_mpo(sites_,terms_b,mpomaxm_);
        double beta_half = 1.0/(2.0*T);
        int nsteps = std::max(1,(int)std::ceil(beta_half/dbeta_half_step));
        Cplx dtau = Cplx(0.0,-1.0)*(beta_half/double(nsteps));
        Cplx dtr = Cplx(dt,0.0); // real time step, tdvp_step's own -i*dt convention

        std::mt19937_64 rng(seed);
        int nbasis = (int)basis_ops.size();
        auto eigcache = metts_build_eigcache(basis_ops);
        MPS cps = random_cps(basis_ops[0],rng,eigcache);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);

        std::vector<std::vector<Cplx>> samples(nt); // samples[k] over retained METTS samples
        for (auto& s : samples) s.reserve(nsamples);

        int total_iters = nwarmup+nsamples;
        for (int it=0; it<total_iters; ++it)
            {
            MPS phi = cps;
            for (int s=0; s<nsteps; ++s)
                {
                phi = tdvp_step(H_,phi,dtau,2,niter);
                phi /= sqrt(innerC(phi,phi).real());
                }
            if (it>=nwarmup)
                {
                MPS v = apply_mpo(B,phi,args);
                // |v(t)> is deliberately *not* unit-norm (paper's own
                // Eq. right after (eq:auxstatewi): "while |w_i(t)> is
                // normalized, we generically have <v_i(t)|v_i(t)> =
                // constant != 1") -- but tdvp_step()'s underlying
                // ITensorTDVP tdvp() call always passes "DoNormalize",
                // true, which force-renormalizes its *input* MPS to unit
                // norm on every single call, silently corrupting v's own
                // (should-stay-constant, generically !=1) norm at the
                // very first step and pinning it at 1 forever after --
                // confirmed directly (a factor-of-2 discrepancy against
                // the exact ED reference for a v(0) with norm 0.5). Same
                // fix quench_tdvp()/quench_tdvp_gse() above already apply
                // to their own similarly-unnormalized psi1: save the true
                // norm once, then restore it after every tdvp_step() call
                // (undoing that forced renormalization to 1) -- w needs
                // no such correction since it starts at unit norm (a
                // normalized METTS sample) and real-time evolution
                // preserves that exactly, making DoNormalize's forced
                // renormalization to 1 a no-op for it.
                Cplx normv0 = std::sqrt(innerC(v,v));
                MPS w = phi;
                for (int k=0; k<nt; ++k)
                    {
                    samples[k].push_back(innerC(w,A,v));
                    if (k<nt-1)
                        {
                        v = tdvp_step(H_,v,dtr,2,tdvp_niter);
                        v.normalize();
                        v *= normv0;
                        w = tdvp_step(H_,w,dtr,2,tdvp_niter);
                        }
                    }
                }
            cps = collapse_to_cps(phi,basis_ops[it % nbasis],rng,eigcache);
            }

        // Same reachable-only-with-an-oversized-nwarmup guard as
        // metts_vev() above.
        if (samples[0].empty())
            Error("Chain::metts_dynamical_correlator: no samples were retained (nwarmup >= nwarmup+nsamples?)");
        std::vector<Cplx> means(nt,Cplx(0.0,0.0));
        std::vector<double> stderrs(nt,0.0);
        for (int k=0; k<nt; ++k)
            {
            Cplx mean = 0;
            for (auto& v : samples[k]) mean += v;
            mean /= double(samples[k].size());
            double var = 0.0;
            for (auto& v : samples[k])
                {
                double d = std::abs(v-mean);
                var += d*d;
                }
            double stderr_ = samples[k].size()>1
                ? std::sqrt(var/double(samples[k].size()-1)/double(samples[k].size())) : 0.0;
            means[k] = mean;
            stderrs[k] = stderr_;
            }
        return {means,stderrs};
        }

    // Global subspace expansion (TDVP/basisextension.h's addBasis()):
    // enriches phi's local bases with a Krylov subspace {phi, H*phi,
    // H^2*phi, ...} of dimension krylov_order, then discards the least
    // significant directions via density-matrix truncation -- the
    // Yang-White scheme (arXiv:2005.06104) that lets one-site TDVP grow
    // bond dimension the way two-site TDVP does via SVD. Mirrors
    // TDVP/sample/run.cc's own call: a flat per-order cutoff (truncK, one
    // entry per one of the krylov_order-1 MPO applications) is simpler to
    // expose as a single Python-level knob than a per-order vector, and
    // is what the sample itself uses (epsilonK = {1E-12, 1E-12}). When
    // maxdim>0 the maxdimK overload is used instead (README's "typical
    // strategy" for when phi's bond dimension is no longer small: cap
    // each Krylov application at a fixed bond dimension rather than a
    // fixed truncation error) -- exposed for that future use case, but
    // every current caller (quench_tdvp_gse/evolve_and_measure_tdvp_gse
    // below, and pyitensor/gse.py's mirrored maxdim parameter) always
    // passes maxdim=0, so this branch is unreached in practice today.
    MPS
    global_subspace_expand(MPO const& H, MPS phi, int krylov_order,
                            double cutoff, int maxdim) const
        {
        // krylov_order<=1 means zero Krylov companions (krylov_order-1
        // MPO applications) -- a no-op, mirrored by pyitensor/gse.py's own
        // identical guard. Without this, krylov_order<=0 would make
        // krylov_order-1 negative, which as the std::vector<...> size
        // argument just below underflows to a huge unsigned value and
        // throws std::length_error instead of doing nothing.
        if (krylov_order<=1) return phi;
        auto args = Args("Cutoff",cutoff,"Method","DensityMatrix",
                          "KrylovOrd",krylov_order,"DoNormalize",true,
                          "Quiet",!verbose_);
        if (maxdim>0)
            {
            auto maxdimK = std::vector<int>(krylov_order-1,maxdim);
            addBasis(phi,H,maxdimK,args);
            }
        else
            {
            auto truncK = std::vector<Real>(krylov_order-1,cutoff);
            addBasis(phi,H,truncK,args);
            }
        // addBasis()/denmatSumDecomp() (basisextension.h, vendored
        // unmodified from upstream) only bounds each bond's *new*
        // directions by cutoff -- there's no per-call cap on the
        // resulting *total* bond dimension, so for a large-enough system
        // with genuinely entangled Krylov vectors, every one of the
        // requested gse_sweeps calls can keep adding more, with no
        // ceiling (confirmed directly: n=20, default settings, bond
        // dimension ran to 370 -- 9x this chain's own maxm=40 -- before
        // being killed). Passing "MaxDim",maxm_ into addBasis()'s own
        // Args doesn't fix this either: it caps the *new* directions
        // alone (denmatSumDecomp's diag_hermitian call) at maxm_, not
        // "new+existing" -- so growth is merely slower (+40 per call
        // instead of unbounded) rather than actually capped.
        //
        // A separate, explicit truncating sweep afterward fixes the
        // *count*, but critically must use "Cutoff",0.0 here, NOT
        // `cutoff` -- GSE's whole point is that the directions it just
        // added carry ~zero weight in phi's own state (that's what makes
        // it state-preserving), *before* any subsequent one-site TDVP
        // evolution has a chance to rotate/populate them. A Cutoff-based
        // truncation right after GSE therefore always discards exactly
        // those directions again, which was confirmed directly to
        // silently turn every later one-site TDVP step into a pure
        // global-phase no-op (bond dimension right back to phi's
        // original rank, and a rank-1 local tensor under one-site TDVP
        // literally cannot do anything but rotate its own phase). A
        // MaxDim-only (Cutoff=0) truncation only trims anything at all
        // once the *count* actually exceeds maxm_, leaving GSE's added
        // directions alone whenever there's room -- exactly the desired
        // "hard cap, but don't undo GSE" behavior.
        phi.position(length(phi),{"Cutoff",0.0,"MaxDim",maxm_});
        phi.position(1,{"Cutoff",0.0,"MaxDim",maxm_});
        return phi;
        }

    double
    cvm_dynamical_correlator(std::vector<MOTerm> const& terms_i,
                             std::vector<MOTerm> const& terms_j,
                             double omega, double eta, double energy,
                             double tol, int max_it) const
        {
        if (!have_wf0_) Error("Chain::cvm_dynamical_correlator called before gs_energy");
        auto S1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto S2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        const Cplx z(omega+energy,eta);
        auto ampo = AutoMPO(sites_);
        ampo += z,"Id",1;
        auto zId = toMPO(ampo);
        auto A = sum(zId,(-1.0)*H_,args);
        auto b = apply_mpo(S2,wf0_,args);
        auto x = bicstab(A,b,tol,max_it,args);
        Cplx G = innerC(wf0_,S1,x);
        return -G.imag()/M_PI;
        }

    MPS
    apply_inverse(std::vector<MOTerm> const& terms, MPS const& wf,
                  double tol, int max_it) const
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        return bicstab(A,wf,tol,max_it,args);
        }

    std::vector<std::complex<double>>
    correlation_matrix(MPS const& wf) const
        {
        int N = sites_.length();
        std::vector<std::complex<double>> out(N*N);
        for (int i=0;i<N;i++)
            {
            for (int j=i;j<N;j++)
                {
                auto ampo = AutoMPO(sites_);
                ampo += 1.0,"Cdag",i+1,"C",j+1;
                auto op = toMPO(ampo);
                auto c = innerC(wf,op,wf);
                out[i*N+j] = c;
                out[j*N+i] = std::conj(c);
                }
            }
        return out;
        }

    std::vector<std::complex<double>>
    four_correlation_tensor(MPS const& wf, bool accelerate=true) const
        {
        int N = sites_.length();
        std::vector<std::complex<double>> out(N*N*N*N,0.0);
        auto idx = [N](int i,int j,int k,int l) { return ((i*N+j)*N+k)*N+l; };
        if (accelerate)
            {
            for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
            for (int l=0;l<N;l++)
                {
                std::tuple<int,int,int,int> current{i,j,k,l};
                std::tuple<int,int,int,int> conjugate{l,k,j,i};
                if (current<=conjugate)
                    {
                    auto ampo = AutoMPO(sites_);
                    ampo += 1.0,"Cdag",i+1,"C",j+1,"Cdag",k+1,"C",l+1;
                    auto op = toMPO(ampo);
                    auto c = innerC(wf,op,wf);
                    out[idx(i,j,k,l)] = c;
                    if (current!=conjugate) out[idx(l,k,j,i)] = std::conj(c);
                    }
                }
            }
        else
            {
            for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
            for (int l=0;l<N;l++)
                {
                auto ampo = AutoMPO(sites_);
                ampo += 1.0,"Cdag",i+1,"C",j+1,"Cdag",k+1,"C",l+1;
                auto op = toMPO(ampo);
                auto c = innerC(wf,op,wf);
                out[idx(i,j,k,l)] = c;
                out[idx(l,k,j,i)] = std::conj(c);
                }
            }
        return out;
        }

    // Native-spinful-site (Electron/Hubbard, site-type code 1) version of
    // four_correlation_tensor() above: returns the same
    // <Cdag_i C_j Cdag_k C_l> tensor, but indexed over the 2*N flat
    // fermionic modes (mode 2*s=up, 2*s+1=down at physical site s) that
    // this chain's N native sites represent -- matching exactly the flat
    // indexing fermionchain.Spinful_Fermionic_Chain_Native's own
    // self.C/self.Cdag use (see jordanwigner_spinful.py's module
    // docstring for the convention). Unlike every other Jordan-Wigner
    // string in this codebase (threaded explicitly at the Python level,
    // see mo_terms.h/build_ampo()), this one instead hands ITensor's own
    // AutoMPO the literal flavor-resolved "Cdagup"/"Cup"/"Cdagdn"/"Cdn"
    // names directly and lets its built-in isFermionic()/fermionicTerm()
    // machinery (autompo.cc) do the threading -- those names all start
    // with 'C' (ITensor's own trigger for automatic fermionic handling)
    // and are exactly the ones ElectronSite defines, unlike the bare
    // "Cdag"/"C" the non-spinful four_correlation_tensor() above relies
    // on, which ElectronSite has no operator for at all. This is safe to
    // do only for this one, self-contained, always-freshly-built-AutoMPO
    // calculation -- it is not a change to the general Hamiltonian/MPO
    // pipeline (mo_terms.h), which still always receives pre-dressed
    // Python-level Jordan-Wigner terms, for consistency across backends.
    std::vector<std::complex<double>>
    four_correlation_tensor_spinful(MPS const& wf, bool accelerate=true) const
        {
        int Ns = sites_.length(); // number of physical (Electron) sites
        int N = 2*Ns; // number of flat fermionic modes
        std::vector<std::complex<double>> out(N*N*N*N,0.0);
        auto idx = [N](int i,int j,int k,int l) { return ((i*N+j)*N+k)*N+l; };
        auto cdagname = [](int flat) { return (flat%2==0) ? "Cdagup" : "Cdagdn"; };
        auto cname    = [](int flat) { return (flat%2==0) ? "Cup" : "Cdn"; };
        auto site1based = [](int flat) { return flat/2+1; };
        auto build_mpo = [&](int i,int j,int k,int l)
            {
            auto ampo = AutoMPO(sites_);
            ampo += 1.0,cdagname(i),site1based(i),cname(j),site1based(j),
                        cdagname(k),site1based(k),cname(l),site1based(l);
            return toMPO(ampo);
            };
        if (accelerate)
            {
            for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
            for (int l=0;l<N;l++)
                {
                std::tuple<int,int,int,int> current{i,j,k,l};
                std::tuple<int,int,int,int> conjugate{l,k,j,i};
                if (current<=conjugate)
                    {
                    auto op = build_mpo(i,j,k,l);
                    auto c = innerC(wf,op,wf);
                    out[idx(i,j,k,l)] = c;
                    if (current!=conjugate) out[idx(l,k,j,i)] = std::conj(c);
                    }
                }
            }
        else
            {
            for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
            for (int l=0;l<N;l++)
                {
                auto op = build_mpo(i,j,k,l);
                auto c = innerC(wf,op,wf);
                out[idx(i,j,k,l)] = c;
                out[idx(l,k,j,i)] = std::conj(c);
                }
            }
        return out;
        }

    // Fast single-sweep four-point correlator tensor <Cdag_i C_j Cdag_k C_l>
    // for the plain (non-spinful) case, following the algorithmic idea of
    // ITensorCorrelators.jl (https://github.com/ITensor/ITensorCorrelators.jl):
    // instead of building and applying a fresh AutoMPO for every (i,j,k,l)
    // tuple independently (four_correlation_tensor() above -- O(N^4)
    // AutoMPO/MPO builds, each an O(N)-ish term compilation + MPO
    // compression + full inner-product sweep), reuse a single left-to-right
    // sweep's worth of partial contractions ("environments") across the
    // whole tensor.
    //
    // Scope: only the N(N-1)(N-2)(N-3) entries with i,j,k,l *pairwise
    // distinct* go through the fast sweep below; the remaining (subdominant,
    // O(N^3) of the O(N^4) total) entries with a repeated index fall back to
    // the same per-tuple AutoMPO method as four_correlation_tensor() (see
    // the tail of this function). Repeats need same-site multi-operator
    // products (e.g. Cdag_i C_i = N_i, or Cdag_i Cdag_i = 0) that the sweep
    // below -- keyed on four *strictly increasing* site positions -- isn't
    // built to produce, and since they're a lower order of the total work a
    // slower fallback there doesn't cost the overall asymptotic win.
    //
    // Algorithm, for four strictly increasing sites a<b<c<d: any assignment
    // of these to (i,j,k,l) is some permutation of which of
    // {Cdag_i,C_j,Cdag_k,C_l} sits at which of {a,b,c,d}. Reordering
    // fermionic operators from the abstract (i,j,k,l) order into the
    // physical site order a<b<c<d picks up a sign equal to the parity of
    // that permutation -- standard fermion anticommutation: swapping any two
    // operators at *different* sites flips the sign of the expectation
    // value, regardless of which specific operator types they are. The
    // value in site order is then obtained by threading a Jordan-Wigner "F"
    // (fermion parity) string through every *gap* site strictly between two
    // consecutively-applied operators, alternating: F in the first gap
    // (a,b) (one operator applied so far = odd), none in the second gap
    // (b,c) (two applied = even), F in the third gap (c,d) (three applied =
    // odd). This "odd/even by operator count" rule depends only on *how
    // many* operators have been swept past, not on which types they are, so
    // it is identical for all six possible Cdag/C type patterns at a given
    // (a,b,c,d) -- which is exactly what lets one (a,b,c,d) sweep serve
    // every (i,j,k,l) permutation that maps onto it. This mirrors
    // ITensorCorrelators.jl's add_operator_fermi (site-sorted operator
    // application plus a "par = 1-2*parity(sortperm(...))" sign correction),
    // specialized to a fixed 4-operator Cdag/C/Cdag/C string rather than
    // porting its fully generic n-operator/repeated-site recursion engine.
    //
    // Environment reuse: the nested loops over a<b<c<d each carry a
    // "running" environment extended by exactly one new site per loop step
    // (never rebuilt from scratch), so the O(N) cost of threading the
    // F-string through a gap is paid once per gap endpoint pair, not once
    // per (i,j,k,l) tuple -- turning the O(N^4) independent O(N)-cost MPO
    // builds of four_correlation_tensor() into O(N^4) total *cheap* local
    // tensor contractions (four O(N) nested loops, O(1) amortized work per
    // step, times a small constant from trying both Cdag/C at each of the
    // four operator sites). The two ends outside [a,d] use the standard
    // mixed-canonical-form shortcut instead of being folded in site by
    // site: psi.position(a) makes everything left of a collapse to an
    // identity on the bond just before a, and -- since d>a -- everything
    // right of d is *also* already right-orthonormal in that same gauge, so
    // it collapses too (delta on the bond just after d). The three
    // *interior* gaps do need an explicit per-site fold, since only one
    // side of the chain is canonical relative to the fixed orthogonality
    // center 'a' at a time.
    //
    // `accelerate` here only gates the (subdominant, O(N^3)) repeated-
    // index fallback loop below, unlike four_correlation_tensor()/
    // four_correlation_tensor_spinful() where it skips ~half of the
    // *dominant* per-tuple AutoMPO builds via the (i,j,k,l)<->(l,k,j,i)
    // conjugate-pair symmetry. There is no equivalent saving available in
    // the pairwise-distinct fast-sweep body above: its dominant cost is
    // the shared environment sweep across (a,b,c,d) and the six per-
    // pattern eltC() evaluations, both already paid once regardless of
    // how many of the (up to 24) output entries a given leaf value gets
    // scattered into, so skipping half of those *output writes* would not
    // skip any of the actual work -- confirmed directly (identical wall-
    // clock and identical results for accelerate=true vs false on the
    // dominant part). Don't read `accelerate=true` here as implying the
    // same roughly-2x win it gives four_correlation_tensor().
    std::vector<std::complex<double>>
    four_correlation_tensor_sweep(MPS const& wf, bool accelerate=true) const
        {
        int N = sites_.length();
        std::vector<std::complex<double>> out(N*N*N*N,0.0);
        auto idx = [N](int i,int j,int k,int l) { return ((i*N+j)*N+k)*N+l; };

        // Deliberately NOT renormalized (unlike reduced_dm()'s psi above):
        // four_correlation_tensor() computes the raw innerC(wf,op,wf) with
        // no enforced normalization, and the repeated-index fallback loop
        // below calls innerC(wf,...) on this same, unmodified wf -- both
        // must agree on that convention, or the two halves of one output
        // tensor would carry different, inconsistent overall scales for
        // any non-unit-norm wf (confirmed directly: an earlier version of
        // this method normalized psi here, which is harmless for an
        // already-unit-norm ground state but silently returned a tensor
        // scaled by 1/norm^4 on the pairwise-distinct entries against a
        // scaled-by-norm^2 fallback on the repeated-index entries for any
        // wf that wasn't already unit-normalized -- e.g. wf*c for any
        // scalar c). psi is still copied (not aliased) because position()
        // mutates its gauge in place.
        auto psi = wf;

        // Apply a single-site operator (by name; "" = no-op) to a ket
        // tensor and restore the standard unprimed physical-index
        // convention -- same idiom as apply_mpo()'s noPrime(TagSet("Site")).
        auto apply_local = [&](ITensor T, int site1based, std::string const& name) -> ITensor
            {
            if (name.empty()) return T;
            T = T*op(sites_,name,site1based);
            T.noPrime(TagSet("Site"));
            return T;
            };
        // Fold one more site into a running left environment: apply an
        // optional operator (or "F" for a plain JW-string gap site) to the
        // ket, then contract against the bra (Link-tagged indices primed,
        // so the physical/site index directly traces out for a bare fold
        // and threads through the operator otherwise) -- same pattern as
        // reduced_dm()'s trace sweep above.
        auto fold = [&](ITensor E, int site1based, bool applyF, std::string const& opname) -> ITensor
            {
            ITensor T = psi.A(site1based);
            if (applyF) T = apply_local(T,site1based,"F");
            T = apply_local(T,site1based,opname);
            T = (E ? E*T : T);
            return T*dag(prime(psi.A(site1based),TagSet("Link")));
            };

        // Every strictly increasing quadruple of sites (a,b,c,d) can be
        // assigned to (i,j,k,l) in 4! = 24 ways; precompute, once, the sign
        // and the "which rank supplies which of i,j,k,l" data for all 24 --
        // see the algorithm comment above for the sign rule. mask's bit r
        // (r=0..3, standing for a,b,c,d respectively) is set iff the
        // operator landing at that rank is Cdag-typed (i.e. came from slot
        // 0=i or slot 2=k of the original Cdag_i C_j Cdag_k C_l string).
        struct PermEntry { int mask; int irank,jrank,krank,lrank; int sign; };
        static const std::vector<PermEntry> perms = [] {
            std::vector<PermEntry> t;
            int p[4] = {0,1,2,3}; // p[rank] = which original slot (i,j,k,l) lands at that rank
            do {
                int inv[4]; for (int r=0;r<4;++r) inv[p[r]] = r; // inv[slot] = rank
                int inversions = 0;
                for (int x=0;x<4;++x) for (int y=x+1;y<4;++y) if (p[x]>p[y]) ++inversions;
                int mask = 0; for (int r=0;r<4;++r) if (p[r]%2==0) mask |= (1<<r);
                int sign = (inversions%2==0) ? 1 : -1;
                // Extra correction beyond plain reordering parity: of the
                // six weight-2 masks, only 5 (0b0101, Cdag,C,Cdag,C) and 10
                // (0b1010, C,Cdag,C,Cdag) strictly alternate type at every
                // rank; the other four (3,6,9,12, which all have at least
                // one pair of *adjacent same-type* operators) need one more
                // sign flip on top of the permutation parity above. This
                // isn't a plain fermion-anticommutation fact (that's fully
                // captured by `inversions` already, for genuinely distinct
                // sites): it traces back to the same "-1" that
                // jordanwigner.py's CC()/CdagCdag() carry but CdagC()/
                // CCdag() don't -- reordering two *same-type* operators
                // (Cdag,Cdag or C,C) into JW-string-canonical form takes one
                // extra local F/Adag or F/A anticommutation beyond what a
                // mixed Cdag,C pair needs. Confirmed empirically (not just
                // derived): every mask-3/6/9/12 entry disagreed with
                // four_correlation_tensor() by exactly this sign, for every
                // (a,b,c,d) and N tried, before this correction was added.
                if (mask!=5 && mask!=10) sign = -sign;
                t.push_back({mask,inv[0],inv[1],inv[2],inv[3],sign});
                } while (std::next_permutation(p,p+4));
            return t;
            }();

        for (int ai=0; ai<=N-4; ++ai)
            {
            int a = ai+1; // 1-indexed ITensor site
            psi.position(a);
            ITensor Lbase; // default-constructed = trivial "1" when a==1
            if (a>1)
                {
                auto ln = commonIndex(psi.A(a-1),psi.A(a));
                Lbase = delta(dag(prime(ln)),ln);
                }
            for (int opAi=0; opAi<2; ++opAi)
                {
                std::string opA = (opAi==0) ? "Cdag" : "C";
                ITensor La = fold(Lbase,a,false,opA);
                ITensor LabRunning = La; // extended incrementally as bi grows
                for (int bi=ai+1; bi<=N-3; ++bi)
                    {
                    int b = bi+1;
                    if (bi>ai+1) LabRunning = fold(LabRunning,b-1,true,""); // gap (a,b): 1 op so far -> F
                    for (int opBi=0; opBi<2; ++opBi)
                        {
                        std::string opB = (opBi==0) ? "Cdag" : "C";
                        ITensor Lab = fold(LabRunning,b,false,opB);
                        ITensor LabcRunning = Lab;
                        for (int ci=bi+1; ci<=N-2; ++ci)
                            {
                            int c = ci+1;
                            if (ci>bi+1) LabcRunning = fold(LabcRunning,c-1,false,""); // gap (b,c): 2 ops -> no F
                            for (int opCi=0; opCi<2; ++opCi)
                                {
                                std::string opC = (opCi==0) ? "Cdag" : "C";
                                ITensor Labc = fold(LabcRunning,c,false,opC);
                                ITensor LabcdRunning = Labc;
                                for (int di=ci+1; di<N; ++di)
                                    {
                                    int d = di+1;
                                    if (di>ci+1) LabcdRunning = fold(LabcdRunning,d-1,true,""); // gap (c,d): 3 ops -> F
                                    for (int opDi=0; opDi<2; ++opDi)
                                        {
                                        int mask = (opAi==0?1:0)|(opBi==0?2:0)|
                                                   (opCi==0?4:0)|(opDi==0?8:0);
                                        if (__builtin_popcount((unsigned)mask)!=2) continue; // not a Cdag,C,Cdag,C pattern
                                        std::string opD = (opDi==0) ? "Cdag" : "C";
                                        ITensor Labcd = fold(LabcdRunning,d,false,opD);
                                        if (d<N)
                                            {
                                            auto rn = commonIndex(psi.A(d),psi.A(d+1));
                                            Labcd = Labcd*delta(dag(prime(rn)),rn);
                                            }
                                        Cplx val = eltC(Labcd);
                                        int posArr[4] = {a-1,b-1,c-1,d-1}; // back to 0-indexed sites
                                        for (auto const& e : perms)
                                            {
                                            if (e.mask!=mask) continue;
                                            int i = posArr[e.irank], j = posArr[e.jrank],
                                                k = posArr[e.krank], l = posArr[e.lrank];
                                            out[idx(i,j,k,l)] = double(e.sign)*val;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        // Repeated-index entries (subdominant, not covered above): same
        // per-tuple AutoMPO method as four_correlation_tensor(), applied to
        // the *original* (unnormalized-assumption-only, un-mutated) wf.
        auto build_mpo = [&](int i,int j,int k,int l)
            {
            auto ampo = AutoMPO(sites_);
            ampo += 1.0,"Cdag",i+1,"C",j+1,"Cdag",k+1,"C",l+1;
            return toMPO(ampo);
            };
        for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
        for (int k=0;k<N;k++)
        for (int l=0;l<N;l++)
            {
            bool distinct = (i!=j && i!=k && i!=l && j!=k && j!=l && k!=l);
            if (distinct) continue;
            std::tuple<int,int,int,int> current{i,j,k,l};
            std::tuple<int,int,int,int> conjugate{l,k,j,i};
            if (accelerate && current>conjugate) continue;
            auto op_ = build_mpo(i,j,k,l);
            auto c = innerC(wf,op_,wf);
            out[idx(i,j,k,l)] = c;
            if (!accelerate || current!=conjugate) out[idx(l,k,j,i)] = std::conj(c);
            }
        return out;
        }

    KPMResult
    kpm_dynamical_correlator(std::vector<MOTerm> const& terms_i,
                             std::vector<MOTerm> const& terms_j,
                             int kpmmaxm, double kpm_scale, bool kpm_accelerate,
                             int kpm_n_scale, double delta, double kpm_cutoff)
        {
        if (!have_H_) Error("Chain::kpm_dynamical_correlator called before set_hamiltonian");
        if (!have_wf0_) gs_energy(); // ensure a ground state is available
        auto hs = scaled_hamiltonian(kpm_scale);
        int n = int(std::round((hs.emax-hs.emin)/delta))*kpm_n_scale;
        auto m1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto m2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(m1,wf0_,{"MaxDim",kpmmaxm,"Cutoff",kpm_cutoff});
        auto psi2 = apply_mpo(m2,wf0_,{"MaxDim",kpmmaxm,"Cutoff",kpm_cutoff});
        KPMResult out;
        out.moments = kpm_moments(hs.scaled_H,psi1,psi2,n,kpmmaxm,kpm_cutoff,kpm_accelerate);
        out.emin = hs.emin; out.emax = hs.emax; out.scale = hs.scale;
        out.num_polynomials = n;
        return out;
        }

    std::vector<std::complex<double>>
    general_kpm(std::vector<MOTerm> const& terms_x, MPS const& wfa, MPS const& wfb,
                int kpmmaxm, bool kpm_accelerate, int num_polynomials, double kpm_cutoff)
        {
        auto m = build_mpo(sites_,terms_x,mpomaxm_);
        return kpm_moments(m,wfa,wfb,num_polynomials,kpmmaxm,kpm_cutoff,kpm_accelerate);
        }

    // Non-Hermitian KPM (port of NHKPM.jl, https://github.com/GUANGZECHEN/
    // NHKPM.jl, itself implementing the method of Phys. Rev. Lett. 130,
    // 100401): biorthogonal Chebyshev moments from a *coupled*
    // forward/adjoint recursion. A genuinely non-Hermitian scaled operator
    // hs=(z*Id-H)/E_max has no real spectrum, so the plain single-operator
    // Chebyshev recursion used by kpm_moments_full (valid only once H is
    // rescaled onto the real interval [-1,1]) doesn't apply; the reference
    // instead builds a second sequence (alpha) from hs_dag alongside the
    // usual Chebyshev-like vectors (vn) of hs, which makes vn a valid
    // expansion for the complex-shifted operator. terms_hs/terms_hs_dag
    // are hs and its conjugate transpose, built at the Python/
    // MultiOperator level (see nonhermitian/kpm.py) -- unlike ordinary
    // KPM, the expansion operator itself depends on the frequency z, so
    // this is called once per requested z rather than once for the whole
    // spectrum. wfa/wfb are the ket/bra-side states, already dressed by
    // whichever operators define the correlator. Returns mu_n, length n,
    // with mu_n[k] = <wfb|vn_{2k-1}> in the 1-based notation of the
    // reference (only the odd-order vn ever get dotted with wfb; see
    // nonhermitian/kpm.py's spec_from_moments_nh for the final
    // reconstruction, which is what actually needs those odd orders).
    std::vector<std::complex<double>>
    nhkpm_moments(std::vector<MOTerm> const& terms_hs,
                  std::vector<MOTerm> const& terms_hs_dag,
                  MPS const& wfa, MPS const& wfb,
                  int n, int kpmmaxm, double kpmcutoff) const
        {
        if (n<1) Error("Chain::nhkpm_moments: n must be >= 1");
        auto hs = build_mpo(sites_,terms_hs,mpomaxm_);
        auto hsd = build_mpo(sites_,terms_hs_dag,mpomaxm_);

        auto v = 1.0*wfa;
        auto alpha_prev2 = apply_mpo(hsd,v,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}); // alpha[1]
        auto alpha_prev1 = sum(2.0*apply_mpo(hs,alpha_prev2,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}),
                                -1.0*v,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}); // alpha[2]

        auto vn_prev2 = 1.0*v; // vn[1]
        auto vn_prev1 = 2.0*apply_mpo(hs,vn_prev2,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}); // vn[2]

        std::vector<std::complex<double>> mu;
        mu.reserve(n);
        mu.push_back(innerC(wfb,vn_prev2));
        for (int k=1;k<n;k++)
            {
            auto alpha_x = sum(2.0*apply_mpo(hsd,alpha_prev1,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}),
                                -1.0*alpha_prev2,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            auto alpha_y = sum(2.0*apply_mpo(hs,alpha_x,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}),
                                -1.0*alpha_prev1,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            auto t = sum(2.0*apply_mpo(hsd,vn_prev1,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}),
                         -1.0*vn_prev2,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            auto vn_x = sum(2.0*alpha_prev1,t,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            auto vn_y = sum(2.0*apply_mpo(hs,vn_x,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff}),
                            -1.0*vn_prev1,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});

            mu.push_back(innerC(wfb,vn_x));

            alpha_prev2 = alpha_x; alpha_prev1 = alpha_y;
            vn_prev2 = vn_x; vn_prev1 = vn_y;
            }
        return mu;
        }

    // KPM energy truncation (Holzner, Weichselbaum, McCulloch & von Delft,
    // "Chebyshev matrix product state approach for spectral functions",
    // PRB 83, 195115 (2011), Sec. III-B) -- native ITensor v3 port of
    // pyitensor/kpm_energy_truncation.py's energy_truncate(), wired in as
    // a wholly independent method from kpm_dynamical_correlator()/
    // kpm_moments_full()/kpm_moments_accelerated()/scaled_hamiltonian()
    // above: none of those are read or modified by any of what follows,
    // so the existing (always window~full-bandwidth-safe) KPM path is
    // completely unaffected.
    //
    // Why this exists: the KPM window Ws can be chosen much smaller than
    // the true many-body bandwidth W (Ws only needs to cover the
    // correlator's own, usually much narrower, spectral width), which is
    // what buys spectral resolution. But shrinking Ws makes it easy for
    // the rescaled Hamiltonian H' to have eigenstates with
    // |eigenvalue| > 1 -- outside the domain where Chebyshev polynomials
    // are bounded. Since the recursion is carried out via MPS
    // compression, not in H's eigenbasis, numerical noise alone is
    // enough to leak a little high-energy weight into a Chebyshev
    // vector, and since Chebyshev polynomials grow without bound outside
    // [-1,1], that weight blows up exponentially over subsequent
    // recursion steps -- check_kpm_moment() above only *detects* this
    // (aborts once moments diverge); the machinery below is the actual
    // fix, applied once per Chebyshev vector before it is used for a
    // moment or fed into the next recursion step.
    KPMResult
    kpm_dynamical_correlator_truncated(std::vector<MOTerm> const& terms_i,
                             std::vector<MOTerm> const& terms_j,
                             int kpmmaxm, double kpm_scale, bool kpm_accelerate,
                             int kpm_n_scale, double delta, double kpm_cutoff,
                             int trunc_dK, int trunc_nsweeps, double trunc_threshold)
        {
        if (!have_H_) Error("Chain::kpm_dynamical_correlator_truncated called before set_hamiltonian");
        if (!have_wf0_) gs_energy(); // ensure a ground state is available
        auto hs = scaled_hamiltonian_gs_anchored(kpm_scale);
        int n = int(std::round((hs.emax-hs.emin)/delta))*kpm_n_scale;
        auto m1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto m2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = apply_mpo(m1,wf0_,{"MaxDim",kpmmaxm,"Cutoff",kpm_cutoff});
        auto psi2 = apply_mpo(m2,wf0_,{"MaxDim",kpmmaxm,"Cutoff",kpm_cutoff});
        KPMResult out;
        out.moments = kpm_moments_truncated(hs.scaled_H,psi1,psi2,n,kpmmaxm,kpm_cutoff,
                                             kpm_accelerate,trunc_dK,trunc_nsweeps,trunc_threshold);
        out.emin = hs.emin; out.emax = hs.emax; out.scale = hs.scale;
        out.num_polynomials = n;
        return out;
        }

    // Infinite DMRG (iDMRG): C++ port of pyitensor/idmrg.py's growing/
    // infinite-size algorithm (White 1992, generalized to an n_uc-site
    // unit cell per McCulloch 2008) -- see that module's own docstring
    // for the full derivation and bug history; this is a line-for-line
    // translation of its logic (automaton construction, growth loop,
    // local 2-site solve, HL/HR extension), not an independent
    // re-derivation, specifically to avoid reintroducing bugs that were
    // already found and fixed there over an extensive debugging session.
    //
    // This Chain instance must be constructed with site_types = the
    // *n_uc-site unit cell* (not a full chain) -- infinitechain.py's
    // itensor_version=3 path does exactly this
    // (cppext.get_backend(3).Chain(self.site_types)), so sites_ here
    // directly serves as the unit cell's own SiteSet (sites_.si(p+1),
    // sites_.op(name,p+1) for p=0..n_uc-1), matching pyitensor's own
    // sites_uc = SiteX(site_types). n_uc = sites_.length().
    //
    // Scope, matching infinitechain.py's v1: Hermitian only (no
    // Hermiticity check here, same documented-but-unenforced precondition
    // gs_energy_generalized has for its own A), n_uc<=2 only (rejected
    // explicitly below -- see idmrg.py's own module docstring for why
    // n_uc>=3 needs a genuine redesign of the growth loop's sublattice
    // pairing), and fermionic (Jordan-Wigner-threaded) terms not
    // supported yet (idmrg_classify_terms always sets carry_ferm=false;
    // fine for spin/boson-type models, silently wrong for a fermionic
    // one -- no detection/guard for that case yet). Ground-state energy
    // density only -- static correlators (onsite_expectation/
    // two_point_correlator) are not ported to C++ yet; run gs_energy()
    // via this method for the energy, then reuse pyitensor's own
    // idmrg.py correlator functions against a separately-run
    // itensor_version="python" IDMRGResult if correlators are needed.
    //
    // terms_intra/terms_inter: 1-based-site MOTerm lists in exactly
    // infinitechain.py's own h_intra/h_inter split (see
    // _canonicalize_hamiltonian), but shifted from that module's 0-based
    // MultiOperator.op site convention to this codebase's usual 1-based
    // AutoMPO/HTerm convention -- the Python-side caller does this shift
    // (a trivial +1 per site, no Jordan-Wigner transform, since idmrg.py
    // never uses MultiOperator's own global JW machinery -- it threads
    // fermionic strings itself per bond, not implemented here yet, see
    // above), mirroring MultiOperator.to_terms()'s own site-shift
    // convention (used by every *other* method in this file) without its
    // JW transform.
    //
    // Real ITensor v3 is *stricter* than pyitensor's own array-backed
    // ITensor reimplementation about Index identity: IndexSet::checkIndexSet()
    // throws ITError("Duplicate indices in index set") the instant a
    // tensor would carry the same Index object on two different axes --
    // exactly the construction idmrg.py's own W_bulk[p] persistent-tensor
    // design relies on for n_uc=1 (a single sublattice's left and right
    // automaton legs are literally the same Index there), and exactly
    // the root cause of that module's own worst bug (HL's and HR's own
    // mpo axes colliding across macro-iterations, see idmrg.py's
    // idmrg_ground_state comment). This port sidesteps the whole bug
    // class from the start rather than retrofitting a fix: it never
    // builds a persistent, reusable W_bulk[p] ITensor at all. Instead
    // idmrg_build_row precomputes each sublattice's channel-transition
    // content as a plain dense array (IdmrgAutomatonRow), and a fresh
    // ITensor is minted from that data -- with fresh Index objects on
    // every leg that needs independence -- every single time it's used,
    // both for the local-solve W_pL/W_pR and for the HL/HR extension
    // step's own re-minted mpo axis (idmrg_make_W/idmrg_make_W_start/
    // idmrg_make_W_end below), i.e. this port is built from the start
    // with the exact fix idmrg.py's own idmrg_ground_state() had to
    // retrofit after finding that collision, not a naive first attempt
    // at the naive (persistent-tensor) design.
    IdmrgResult
    idmrg_ground_state(std::vector<MOTerm> const& terms_intra,
                        std::vector<MOTerm> const& terms_inter,
                        int maxm, double cutoff, int maxiter, double etol,
                        int krylovdim, int restarts)
        {
        int n_uc = sites_.length();
        if (n_uc>2)
            Error("Chain::idmrg_ground_state: n_uc>2 is not supported yet "
                  "(see this method's own comment)");
        if (maxiter<2)
            Error("Chain::idmrg_ground_state: maxiter must be >= 2 (energy "
                  "density is a finite difference between two "
                  "macro-iterations)");

        std::vector<IdmrgOnsite> onsite;
        std::vector<IdmrgBond> bonds;
        idmrg_classify_terms(terms_intra,terms_inter,n_uc,onsite,bonds);

        std::vector<std::vector<IdmrgChan>> chans(n_uc);
        for (int p=0;p<n_uc;++p) chans[p] = idmrg_channels_at(bonds,n_uc,p);

        std::vector<IdmrgAutomatonRow> rows;
        rows.reserve(n_uc);
        for (int p=0;p<n_uc;++p)
            rows.push_back(idmrg_build_row(p,n_uc,bonds,onsite,
                                            chans[p],chans[(p+1)%n_uc]));

        ITensor HL, HR; // default-constructed = "None" (see this file's
                        // own porting-notes comment on ITensor's explicit
                        // operator bool())
        Index HL_bra, HL_ket, HR_bra, HR_ket, HL_mpo, HR_mpo;
        bool have_HL_ket=false, have_HR_ket=false;

        double energy = 0.0;
        bool have_prev_energy=false, have_prev_density=false;
        double prev_energy=0.0, prev_density=0.0;
        bool converged=false;
        int macro_iter=0;

        // Per-sublattice converged ket tensor (idmrg_U_[p_L]=U, set every
        // mstep -- overwritten every macro-iteration, so it holds the
        // *last executed* macro-iteration's own tensors once the loop
        // exits) plus the natural (left,right) Link Index each one carries
        // as solved -- td_dynamical_correlator_window()'s own window
        // builder needs these explicitly (see this method's own snapshot
        // comment below) rather than re-deriving them from tags, since two
        // Link-tagged indices (the "old"/already-absorbed one and the
        // freshly minted SVD bond) would otherwise be ambiguous to tell
        // apart after the fact.
        std::vector<ITensor> snap_U(n_uc);
        std::vector<Index> snap_U_left(n_uc), snap_U_right(n_uc);
        // Snapshot of (HL,HL_bra,HL_ket,HL_mpo,HR,HR_bra,HR_ket,HR_mpo) as
        // they stand *before* each macro-iteration's own mstep loop
        // mutates them -- overwritten every macro-iteration, so whatever
        // it holds once the loop exits is exactly the environment entering
        // the *last executed* macro-iteration, consistent by construction
        // with that same macro-iteration's own snap_U (U_list[0]'s own
        // left leg is exactly this snapshot's HL_ket, since mstep=0 always
        // solves against whatever HL/HL_ket stood at this snapshot point)
        // -- direct C++ port of pyitensor/idmrg.py's own
        // env_window_boundary, see idmrg_ground_state's own module-level
        // comment there for the full derivation. maxiter>=2 (checked
        // above) guarantees the loop always fully executes macro_iter=0
        // before density (hence any possible break) is even computable, so
        // by the time any snapshot below is taken from the *last executed*
        // macro-iteration, HL/HR are never still "None" -- unlike
        // pyitensor's own env_window_boundary (which starts at `(None,)*8`
        // and can in principle be returned as such if maxiter's own
        // minimum were ever bypassed), this port's snap_have_HL_/
        // snap_have_HR_ are therefore always true by the time they are
        // read after the loop, kept as explicit bools anyway for
        // documentation/defense-in-depth rather than assumed silently.
        ITensor snap_HL, snap_HR;
        Index snap_HL_bra, snap_HL_ket, snap_HL_mpo;
        Index snap_HR_bra, snap_HR_ket, snap_HR_mpo;
        bool snap_have_HL=false, snap_have_HR=false;

        // prev_local[mstep]: the flattened local ground vector produced at
        // this unit-cell position by the *previous* macro-iteration,
        // threaded back in as the next macro-iteration's own Arnoldi warm
        // start -- direct C++ port of pyitensor/idmrg.py's own prev_local
        // (see idmrg_local_solve's own comment for why this is not
        // cosmetic). Persists across macro-iterations (declared outside
        // the loop below), starts empty (no warm start available yet).
        std::vector<std::vector<Cplx>> prev_local(n_uc);

        for (macro_iter=0; macro_iter<maxiter; ++macro_iter)
            {
            snap_HL=HL; snap_HL_bra=HL_bra; snap_HL_ket=HL_ket; snap_HL_mpo=HL_mpo;
            snap_HR=HR; snap_HR_bra=HR_bra; snap_HR_ket=HR_ket; snap_HR_mpo=HR_mpo;
            snap_have_HL=bool(HL); snap_have_HR=bool(HR);
            for (int mstep=0; mstep<n_uc; ++mstep)
                {
                int p_L = mstep;
                int p_R = n_uc-1-mstep;
                bool first_step = (macro_iter==0 && mstep==0);
                bool same_phys = (p_L==p_R);

                Index phys_L_in = sites_.si(p_L+1);
                Index phys_L_out = prime(phys_L_in);
                Index phys_R_in = same_phys ? sim(sites_.si(p_R+1)) : sites_.si(p_R+1);
                Index phys_R_out = prime(phys_R_in);

                // Crossing bond between W_pL and W_pR within *this*
                // micro-step's own local solve -- guaranteed the same
                // dimension on both sides since n_uc<=2 keeps
                // (p_L+1)%n_uc == p_R (see this method's own comment on
                // why n_uc>=3 breaks this pairing). link_tag is a
                // function-local static (see arnoldi_smallest_real's own
                // comment on why): avoids re-parsing "Link" into a
                // TagSet on every micro-step.
                static const TagSet link_tag("Link");
                Index cross_idx = Index(rows[p_L].right_n,link_tag);

                ITensor W_pL, W_pR;
                if (first_step)
                    {
                    W_pL = idmrg_make_W_start(rows[p_L],cross_idx,phys_L_in,phys_L_out);
                    W_pR = idmrg_make_W_end(rows[p_R],cross_idx,phys_R_in,phys_R_out);
                    }
                else
                    {
                    W_pL = idmrg_make_W(rows[p_L],HL_mpo,cross_idx,phys_L_in,phys_L_out);
                    W_pR = idmrg_make_W(rows[p_R],cross_idx,HR_mpo,phys_R_in,phys_R_out);
                    }

                auto [energy_here,U,V,new_bond_u,new_bond_v,theta_flat] = idmrg_local_solve(
                    HL,W_pL,phys_L_in,W_pR,phys_R_in,HR,
                    HL_bra,HL_ket,have_HL_ket,HR_bra,HR_ket,have_HR_ket,
                    cutoff,maxm,krylovdim,restarts,prev_local[mstep]);
                energy = energy_here;
                if (!theta_flat.empty()) prev_local[mstep] = std::move(theta_flat);

                Index left_ket_old = HL_ket; bool have_left_ket_old = have_HL_ket;
                Index right_ket_old = HR_ket; bool have_right_ket_old = have_HR_ket;

                // Fresh, independently-minted mpo-axis identities for the
                // *extend* step only (not used in the solve above, so
                // this doesn't disturb the W_pL/W_pR crossing bond the
                // solve just relied on) -- guarantees the new HL's and
                // HR's own mpo axes can never collide with each other or
                // with a future reuse of rows[p], regardless of how
                // small n_uc is (see this method's own top comment).
                Index new_HL_mpo = sim(cross_idx);
                ITensor W_pL_ext = replaceInds(W_pL,{cross_idx},{new_HL_mpo});
                Index new_HR_mpo = sim(cross_idx);
                ITensor W_pR_ext = replaceInds(W_pR,{cross_idx},{new_HR_mpo});

                bool have_HL = bool(HL), have_HR = bool(HR);
                auto [newHL,newHL_bra] = idmrg_extend_HL(
                    HL,HL_bra,have_HL,W_pL_ext,U,
                    left_ket_old,have_left_ket_old,new_bond_u);
                HL = newHL; HL_bra = newHL_bra;
                HL_ket = new_bond_u; have_HL_ket = true;

                auto [newHR,newHR_bra] = idmrg_extend_HR(
                    HR,HR_bra,have_HR,W_pR_ext,V,
                    right_ket_old,have_right_ket_old,new_bond_v);
                HR = newHR; HR_bra = newHR_bra;
                HR_ket = new_bond_v; have_HR_ket = true;

                HL_mpo = new_HL_mpo;
                HR_mpo = new_HR_mpo;

                snap_U[p_L] = U; snap_U_left[p_L] = left_ket_old; snap_U_right[p_L] = new_bond_u;
                }

            bool have_density = have_prev_energy;
            double density = have_density ? (energy-prev_energy)/(2.0*n_uc) : 0.0;
            if (verbose_)
                println("idmrg macro-iter ",macro_iter,": E=",energy,
                        " density=",(have_density?density:0.0));
            if (have_density && have_prev_density &&
                std::abs(density-prev_density)<etol)
                converged = true;
            prev_energy = energy; have_prev_energy = true;
            if (have_density) { prev_density = density; have_prev_density = true; }
            if (converged) break;
            }

        IdmrgResult out;
        out.density = prev_density;
        out.converged = converged;
        // C++'s own for-loop increments macro_iter to maxiter itself (not
        // maxiter-1) before the loop condition fails and exits, unlike
        // Python's "for macro_iter in range(maxiter)" whose loop variable
        // retains its *last-assigned* value (maxiter-1) when not broken
        // out of early -- min() here reconciles the two conventions so
        // niter_done matches Python's own idmrg_ground_state exactly in
        // both the converged (break, macro_iter<maxiter) and
        // ran-to-completion (macro_iter==maxiter) cases.
        out.niter_done = std::min(macro_iter+1,maxiter);

        idmrg_n_uc_ = n_uc;
        idmrg_rows_ = rows;
        idmrg_U_ = snap_U;
        idmrg_U_left_ = snap_U_left;
        idmrg_U_right_ = snap_U_right;
        idmrg_HL_ = snap_HL; idmrg_HL_bra_ = snap_HL_bra;
        idmrg_HL_ket_ = snap_HL_ket; idmrg_HL_mpo_ = snap_HL_mpo;
        idmrg_HR_ = snap_HR; idmrg_HR_bra_ = snap_HR_bra;
        idmrg_HR_ket_ = snap_HR_ket; idmrg_HR_mpo_ = snap_HR_mpo;
        idmrg_have_HL_ = snap_have_HL; idmrg_have_HR_ = snap_have_HR;
        have_idmrg_snapshot_ = true;

        return out;
        }

    // VUMPS ground state (Zauner-Stauber et al., arXiv:1701.07035;
    // Vanderstraeten/Haegeman/Verstraete, arXiv:1810.07006, Algorithm 4)
    // -- C++ port of pyitensor/vumps.py's vumps_ground_state, itself an
    // alternative to idmrg_ground_state's growing/infinite-size algorithm:
    // solves directly, at a FIXED target bond dimension D and in the
    // thermodynamic limit, for the mixed-gauge {AL,AR,C} fixed point of a
    // single (grouped) supersite's transfer/Hamiltonian environment. See
    // pyitensor/vumps.py's own module docstring for the full algorithm
    // derivation (regularized GL/GR environment solves, H_AC/H_C
    // eigenproblems, the orthogonal-Procrustes AL/AR update) and for the
    // documented D>1 convergence-robustness caveats this port shares
    // (single-attempt VUMPS from a random start is not reliable for D>1;
    // mitigated the same way, a D-ramp with multiple restarts per step --
    // see vumps_ground_state's own private helper below for the one
    // simplification this port takes relative to pyitensor's own driver:
    // no "beat a known-good smaller-D energy" safety-net budget, see that
    // helper's own comment).
    //
    // Same scope as idmrg_ground_state: n_uc<=2 (this Chain's own
    // site_types), Hermitian only, and -- specific to VUMPS/the
    // excitation ansatz -- every 2-site term must have "reach" (grouped-
    // supersite separation) exactly 1 (checked below, throws ITError
    // otherwise; TFIM/Heisenberg/XX-style nearest-neighbor models satisfy
    // this trivially once grouped into a <=2-site supersite).
    //
    // Unlike idmrg_ground_state (which works entirely with genuine
    // ITensor tensor-network objects), this method and its private
    // helpers (see "-- VUMPS / excitation ansatz private helpers --"
    // below, near idmrg_make_W_end) build a dense (Dw,Dw,d_g,d_g)
    // finite-state automaton array directly from idmrg_classify_terms/
    // idmrg_channels_at/idmrg_build_row's own already-validated dense
    // per-sublattice output (reused verbatim, not re-derived), then group
    // adjacent sublattices' automaton rows into one supersite exactly
    // mirroring pyitensor/vumps.py's own _group_automaton -- see
    // vumps_group_automaton's own comment.
    //
    // On return, this Chain's own private snapshot (vumps_AL_/vumps_AR_/
    // vumps_C_/vumps_GL_/vumps_GR_/vumps_W_) is set so
    // vumps_excitation_energies() can be called afterward on the same
    // Chain instance -- mirrors idmrg_ground_state's own
    // have_idmrg_snapshot_ pattern.
    VumpsResult
    vumps_ground_state(std::vector<MOTerm> const& terms_intra,
                        std::vector<MOTerm> const& terms_inter,
                        int D, double tol, int maxiter, int nrestarts)
        {
        int n_uc = sites_.length();
        // User-facing, recoverable input validation -- throw ITError (a
        // catchable C++ exception pybind11 turns into a Python exception),
        // NOT the Error(...) macro (which aborts the whole process, no
        // exception to catch -- confirmed directly: an earlier version of
        // this guard used Error(...) here and a plain invalid D=0 call
        // took down the entire Python interpreter instead of raising).
        // n_uc>2 is left as idmrg_ground_state's own Error(...) precedent
        // (that branch is never reachable through the public API in
        // practice -- infinitechain.py's constructor already rejects
        // n_uc>2 before any Chain is ever built).
        if (n_uc>2)
            Error("Chain::vumps_ground_state: n_uc>2 is not supported yet "
                  "(see idmrg_ground_state's own comment for the same "
                  "restriction and its rationale)");
        if (D<1)
            throw ITError("Chain::vumps_ground_state: D must be >= 1");
        if (nrestarts<1)
            throw ITError("Chain::vumps_ground_state: nrestarts must be >= 1");

        std::vector<IdmrgOnsite> onsite;
        std::vector<IdmrgBond> bonds;
        idmrg_classify_terms(terms_intra,terms_inter,n_uc,onsite,bonds);
        std::vector<std::vector<IdmrgChan>> chans(n_uc);
        for (int p=0;p<n_uc;++p) chans[p] = idmrg_channels_at(bonds,n_uc,p);
        std::vector<IdmrgAutomatonRow> rows;
        rows.reserve(n_uc);
        for (int p=0;p<n_uc;++p)
            rows.push_back(idmrg_build_row(p,n_uc,bonds,onsite,
                                            chans[p],chans[(p+1)%n_uc]));

        int Dw=0, d_g=0;
        std::vector<Cplx> W;
        vumps_group_automaton(rows,n_uc,Dw,d_g,W);
        vumps_check_reach_one(W,Dw,d_g);

        auto h1 = vumps_onsite_matrix(W,Dw,d_g);
        auto pending = vumps_pending_channels(W,Dw,d_g);

        std::mt19937_64 rng(std::random_device{}());

        auto better = [](VumpsRunResult const& a, VumpsRunResult const& b)
            {
            if (a.converged != b.converged) return a.converged;
            return a.e_cell < b.e_cell;
            };

        VumpsRunResult best; bool have_best=false;
        std::vector<Cplx> prev_AL, prev_AR; int prev_D=0; bool have_prev=false;
        for (int D_cur=1; D_cur<=D; ++D_cur)
            {
            int n_here = (D_cur==D) ? nrestarts : std::min(nrestarts,3);
            VumpsRunResult local_best; bool have_local=false;
            for (int attempt=0; attempt<n_here; ++attempt)
                {
                VumpsInit init; bool has_init=false;
                if (attempt==0 && have_prev)
                    {
                    init = vumps_grow_init(D_cur,d_g,prev_D,prev_AL,prev_AR,rng);
                    has_init=true;
                    }
                try
                    {
                    auto r = vumps_single_run(D_cur,d_g,h1,pending,tol,maxiter,
                                               has_init?&init:nullptr,rng);
                    if (verbose_)
                        println("vumps D=",D_cur," attempt=",attempt,": e0=",
                                r.e_cell/n_uc," converged=",r.converged);
                    if (!have_local || better(r,local_best)) { local_best=r; have_local=true; }
                    }
                catch (ITError const& e)
                    {
                    // A degenerate/near-degenerate transfer-matrix spectrum
                    // (vx_dominant_*_fixed_point's own guard, mirroring
                    // pyitensor's _check_dominant_eigenvalue_nondegenerate)
                    // or a singular regularized environment solve -- both
                    // recoverable by simply retrying with a different
                    // random start, exactly like pyitensor's own
                    // vumps_ground_state driver (its own "one_attempt"
                    // catches RuntimeError the same way).
                    if (verbose_)
                        println("vumps D=",D_cur," attempt=",attempt,": failed (",e.what(),")");
                    }
                }
            if (!have_local)
                throw ITError("Chain::vumps_ground_state: every attempt at D="+
                              std::to_string(D_cur)+" failed (degenerate transfer-matrix "
                              "spectrum, or a singular regularized environment solve) -- "
                              "try increasing nrestarts");
            prev_AL = local_best.AL; prev_AR = local_best.AR; prev_D = D_cur; have_prev = true;
            best = local_best; have_best = true;
            }
        (void)have_best; // always true here: the D_cur loop runs at least once (D>=1 checked above)

        VumpsResult out;
        out.e0 = best.e_cell / n_uc;
        out.converged = best.converged;
        out.niter_done = best.niter;
        out.gauge_mismatch = best.mismatch;

        vumps_D_ = D; vumps_dg_ = d_g; vumps_Dw_ = Dw;
        vumps_AL_ = best.AL; vumps_AR_ = best.AR; vumps_C_ = best.C;
        vumps_GL_ = best.GL; vumps_GR_ = best.GR; vumps_W_ = W;
        have_vumps_snapshot_ = true;
        have_vumps_exc_env_ = false; // stale -- rebuilt lazily by vumps_excitation_energies
        return out;
        }

    // Tangent-space/quasiparticle excitation ansatz (Haegeman et al.;
    // Vanderstraeten/Haegeman/Verstraete, arXiv:1810.07006, Sec.6) on top
    // of a converged vumps_ground_state() -- C++ port of
    // pyitensor/idmrg_excitations.py's excitation_energies, mirroring
    // MPSKit.jl's own channel-resolved GBL(k)/GBR(k) construction and
    // 3-term H_eff(k) exactly as that module's own module docstring
    // documents (its "History" section records the eight independent
    // investigation passes needed to get this right on the Python side
    // first -- this port follows that already-validated algorithm
    // line-for-line, in dense arrays, rather than re-deriving anything).
    //
    // Requires vumps_ground_state() to have been called first on this
    // same Chain. The excitation environment (V_L, the channel-resolved
    // background GL_full/GR_full, the mixed AL/AR transfer fixed points,
    // and H_AC's own Rayleigh quotient lam_AC -- see
    // vumps_build_excitation_environment's own comment) is built once,
    // lazily, on the first call, and cached until the next
    // vumps_ground_state() call -- exactly mirroring
    // pyitensor.infinitechain.Infinite_Many_Body_Chain's own
    // _excitation_env caching.
    //
    // Returns the lowest `n` excitation energies (above the ground state)
    // at momentum `k` (radians per unit cell) -- an ordinary (not
    // generalized) Hermitian eigenproblem of size Dx*D (Dx=D*(d_g-1)),
    // solved directly via itensor::zheev_wrapper.
    std::vector<double>
    vumps_excitation_energies(double k, int n)
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_excitation_energies: called before "
                           "vumps_ground_state (no converged VUMPS snapshot)");
        int D = vumps_D_, d_g = vumps_dg_;
        int Dx = D*(d_g-1);
        if (Dx<=0)
            throw ITError("Chain::vumps_excitation_energies: d_g<=1 (trivial "
                           "physical dimension after grouping) -- no nontrivial "
                           "tangent-space excitation exists");
        if (n<1)
            throw ITError("Chain::vumps_excitation_energies: n must be >= 1");
        if (!have_vumps_exc_env_) vumps_build_excitation_environment();

        auto Hmat = vumps_build_h_eff_dense(k);
        int nH = Dx*D;
        // Hermitize (H_eff(k) is Hermitian by construction -- this only
        // cleans up numerical noise, same convention idmrg_excitations.py's
        // own excitation_energies uses).
        for (int i=0;i<nH;++i)
        for (int j=0;j<nH;++j)
            {
            Cplx v = (Hmat[i*nH+j] + std::conj(Hmat[j*nH+i]))/2.0;
            Hmat[i*nH+j] = v;
            }
        auto evals = vx_hermitian_eigvals(Hmat,nH); // ascending
        int take = std::min(n,(int)evals.size());
        std::vector<double> out(take);
        for (int i=0;i<take;++i) out[i] = evals[i] - vumps_lam_AC_;
        return out;
        }

    // <opname> at sub-site p (0..n_uc-1) of a converged vumps_ground_state()
    // -- C++ port of pyitensor/vumps.py's onsite_expectation, computed
    // directly from AC (Vanderstraeten, Haegeman, Verstraete,
    // arXiv:1810.07006, Eq.(34)): AC already IS the correctly normalized
    // single-(super)site reduced state by construction of the mixed
    // canonical gauge, so no eigenproblem is needed here at all (unlike
    // idmrg_ground_state's own would-be dominant-right-fixed-point
    // machinery, which this backend does not implement -- see this file's
    // module-level comment at IdmrgResult).
    Cplx
    vumps_onsite_expectation(std::string const& opname, int p)
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_onsite_expectation: called before "
                           "vumps_ground_state (no converged VUMPS snapshot)");
        int n_uc = sites_.length();
        if (p<0 || p>=n_uc)
            throw ITError("Chain::vumps_onsite_expectation: p must be in 0.."+
                           std::to_string(n_uc-1)+" (n_uc-1), got "+std::to_string(p));
        int D = vumps_D_, d_g = vumps_dg_;
        auto AC = vumps_compose_AL_C(vumps_AL_,vumps_C_,D,d_g);
        auto M = vumps_embed_group_operator({{p, idmrg_op_dense(p,opname)}});
        auto AC_op = vx_apply_op_ket(M,AC,D,d_g);
        Cplx val = vx_dot_conj(AC,AC_op);
        double norm = vx_dot_conj(AC,AC).real();
        return val/norm;
        }

    // <opname_i(site p_i) opname_j(site p_i + r)> of a converged
    // vumps_ground_state()'s infinite chain, r measured in physical sites
    // (r>=0) -- C++ port of pyitensor/vumps.py's two_point_correlator; same
    // r=0 same-site convention (Mj@Mi, via idmrg_matmul) and n_uc-
    // periodicity, built from the mixed-gauge {AC,AR} exactly as that
    // function's own docstring derives (AC placed at the unit cell
    // containing p_i, AR -- exactly right-orthonormal by construction --
    // at every cell strictly to its right, so the right closure needs no
    // eigenproblem, only a direct trace after zero or more plain-AR
    // transfer-tensor applications).
    Cplx
    vumps_two_point_correlator(std::string const& opname_i, int p_i,
                                std::string const& opname_j, int r)
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_two_point_correlator: called before "
                           "vumps_ground_state (no converged VUMPS snapshot)");
        if (r<0)
            throw ITError("Chain::vumps_two_point_correlator: r must be >= 0");
        int n_uc = sites_.length();
        if (p_i<0 || p_i>=n_uc)
            throw ITError("Chain::vumps_two_point_correlator: p_i must be in 0.."+
                           std::to_string(n_uc-1)+" (n_uc-1), got "+std::to_string(p_i));
        int D = vumps_D_, d_g = vumps_dg_;
        auto AC = vumps_compose_AL_C(vumps_AL_,vumps_C_,D,d_g);
        double norm = vx_dot_conj(AC,AC).real();

        int cell_offset = (p_i + r)/n_uc;
        int p_j = (p_i + r)%n_uc;

        if (cell_offset==0)
            {
            std::vector<Cplx> M;
            if (p_j==p_i)
                {
                int d = dim(sites_.si(p_i+1));
                auto Mi = idmrg_op_dense(p_i,opname_i);
                auto Mj = idmrg_op_dense(p_i,opname_j);
                M = vumps_embed_group_operator({{p_i, idmrg_matmul(Mj,Mi,d)}});
                }
            else
                {
                M = vumps_embed_group_operator({{p_i, idmrg_op_dense(p_i,opname_i)},
                                                 {p_j, idmrg_op_dense(p_j,opname_j)}});
                }
            auto AC_op = vx_apply_op_ket(M,AC,D,d_g);
            return vx_dot_conj(AC,AC_op)/norm;
            }

        auto Mi_embed = vumps_embed_group_operator({{p_i, idmrg_op_dense(p_i,opname_i)}});
        auto AC_op = vx_apply_op_ket(Mi_embed,AC,D,d_g);
        // Open right-bond object: bra/ket both AC, operator on the ket side
        // only -- this already IS the full left closure (AC's own left leg
        // is summed away here), leaving just the (ket-bond, bra-bond) legs
        // open.
        auto X = vx_left_close(AC_op,AC,D,d_g);
        for (auto& x : X) x /= norm;

        if (cell_offset>1)
            {
            auto E_AR = vx_op_transfer_matrix(vumps_AR_,D,d_g,vumps_AR_,false,{});
            for (int i=0;i<cell_offset-1;++i)
                X = vx_apply_transfer_from_left(E_AR,D,X);
            }

        auto Mj_embed = vumps_embed_group_operator({{p_j, idmrg_op_dense(p_j,opname_j)}});
        auto AR_op = vx_apply_op_ket(Mj_embed,vumps_AR_,D,d_g);
        auto E_AR_op = vx_op_transfer_matrix(AR_op,D,d_g,vumps_AR_,false,{});
        X = vx_apply_transfer_from_left(E_AR_op,D,X);

        Cplx tr(0,0);
        for (int i=0;i<D;++i) tr += X[i*D+i];
        return tr;
        }

    // Apply a periodic (bounded) MPO to `this` Chain's own converged VUMPS
    // snapshot (vumps_AL_ etc, set by a prior vumps_ground_state() or
    // vumps_load_uniform_state() call), returning a NEW mixed-gauge
    // uniform iMPS -- C++ port of pyitensor/vumps.py's own apply_mpo (see
    // this file's own "-- apply_mpo / imps_sum private helpers --" section
    // for the algorithm: vx_grow_by_mpo_n1, then vx_canonicalize_n1, then
    // vumps_complete_mixed_gauge). `this` Chain's own snapshot is left
    // untouched (matches Python's own pure-function semantics) -- call
    // vumps_load_uniform_state() on the result (on this Chain or any other
    // sharing the same site_types) to make vumps_onsite_expectation/
    // vumps_two_point_correlator/a further apply_mpo/imps_sum call see it.
    //
    // W_bulk_flat[p] (p=0..n_uc-1): a dense, row-major (Left,in,out,Right)
    // rank-4 operator tensor at unit-cell sublattice p, size
    // Dw_left[p]*d_p*d_p*Dw_right[p] where d_p=site_dim(p+1) -- the SAME
    // convention pyitensor/vumps.py's own apply_mpo takes (a list of n_uc
    // rank-4 ITensors), just as plain flat arrays instead of ITensor
    // objects. W_bulk_flat must represent a genuinely BOUNDED (non-
    // extensive) periodic operator -- the Hamiltonian's own accumulator
    // automaton is explicitly out of scope, see pyitensor/idmrg.py's own
    // "Applying a (bounded) MPO to the converged iMPS" section docstring
    // for the specific failure mode of feeding it in anyway. maxdim<=0
    // means "no cap" (Python's maxdim=None).
    VumpsUniformResult
    vumps_apply_mpo(std::vector<std::vector<Cplx>> const& W_bulk_flat,
                     std::vector<int> const& Dw_left, std::vector<int> const& Dw_right,
                     double cutoff, int maxdim)
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_apply_mpo: called before vumps_ground_state/"
                           "vumps_load_uniform_state (no converged VUMPS snapshot)");
        int n_uc = sites_.length();
        if ((int)W_bulk_flat.size()!=n_uc || (int)Dw_left.size()!=n_uc || (int)Dw_right.size()!=n_uc)
            throw ITError("Chain::vumps_apply_mpo: expected "+std::to_string(n_uc)+
                           " unit-cell sites in W_bulk_flat/Dw_left/Dw_right");
        std::vector<IdmrgAutomatonRow> rows(n_uc);
        for (int p=0;p<n_uc;++p)
            {
            int d = dim(sites_.si(p+1));
            int left_n=Dw_left[p], right_n=Dw_right[p];
            if ((size_t)left_n*d*d*right_n != W_bulk_flat[p].size())
                throw ITError("Chain::vumps_apply_mpo: W_bulk_flat["+std::to_string(p)+
                               "] size does not match Dw_left*d*d*Dw_right");
            std::vector<Cplx> flat((size_t)left_n*right_n*d*d,Cplx(0,0));
            for (int l=0;l<left_n;++l)
            for (int si=0;si<d;++si)
            for (int so=0;so<d;++so)
            for (int r=0;r<right_n;++r)
                flat[((size_t)l*right_n+r)*d*d+si*d+so] =
                    W_bulk_flat[p][((size_t)l*d+si)*d*right_n + (size_t)so*right_n + r];
            rows[p] = IdmrgAutomatonRow{p,d,left_n,right_n,std::move(flat)};
            }
        int Dw, d_g; std::vector<Cplx> W;
        vumps_group_automaton(rows,n_uc,Dw,d_g,W);
        if (d_g != vumps_dg_)
            throw ITError("Chain::vumps_apply_mpo: W_bulk's own grouped physical "
                           "dimension ("+std::to_string(d_g)+") does not match this "
                           "Chain's own vumps snapshot d_g ("+std::to_string(vumps_dg_)+")");

        int D = vumps_D_;
        auto B = vx_grow_by_mpo_n1(W,Dw,d_g,vumps_AL_,D);
        auto canon = vx_canonicalize_n1(B,Dw*D,d_g,cutoff,maxdim);
        auto mg = vumps_complete_mixed_gauge(canon.AL,canon.D,d_g);

        VumpsUniformResult out;
        out.D = canon.D; out.d_g = d_g; out.AL = canon.AL;
        out.AR = mg.AR; out.C = mg.C; out.AC = mg.AC; out.eta = canon.eta;
        return out;
        }

    // Direct sum of `this` Chain's own converged VUMPS snapshot (state
    // "a") with a second, externally-supplied converged state "b" (AL_b,
    // bond dimension D_b, sharing this Chain's own d_g/n_uc -- state "b"
    // need not have come from a vumps_ground_state() run on this SAME
    // Chain instance; the main intended use is two independently-converged
    // Chains, e.g. two symmetry-related ground states) -- C++ port of
    // pyitensor/vumps.py's own imps_sum. Reliably throws ITError for two
    // *ordinary* VUMPSResults (both individually normalized to eta=1 on
    // both the left AND right transfer eigenvalues by mixed-gauge
    // construction, hence a genuinely degenerate combined dominant
    // eigenvalue -- caught by vx_canonicalize_n1's own vx_dominant_*_
    // fixed_point calls) -- see pyitensor/vumps.py's own "Summing two
    // converged VUMPS iMPS" section docstring for the physical derivation;
    // only two states with a genuine per-site norm mismatch (e.g. one
    // deliberately rescaled) have a well-posed sum. `this` Chain's own
    // snapshot is left untouched.
    VumpsUniformResult
    vumps_imps_sum(int D_b, std::vector<Cplx> const& AL_b, double cutoff, int maxdim)
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_imps_sum: called before vumps_ground_state/"
                           "vumps_load_uniform_state (no converged VUMPS snapshot)");
        int d_g = vumps_dg_;
        if ((size_t)D_b*d_g*D_b != AL_b.size())
            throw ITError("Chain::vumps_imps_sum: AL_b's size does not match D_b*d_g*D_b "
                           "(this Chain's own vumps snapshot d_g="+std::to_string(d_g)+")");
        int Da = vumps_D_, D = Da+D_b;
        std::vector<Cplx> raw((size_t)D*d_g*D,Cplx(0,0));
        for (int p=0;p<d_g;++p)
            {
            for (int i=0;i<Da;++i) for (int j=0;j<Da;++j)
                raw[(size_t)(i*d_g+p)*D+j] = vumps_AL_[(i*d_g+p)*Da+j];
            for (int i=0;i<D_b;++i) for (int j=0;j<D_b;++j)
                raw[(size_t)((Da+i)*d_g+p)*D+(Da+j)] = AL_b[(i*d_g+p)*D_b+j];
            }
        auto canon = vx_canonicalize_n1(raw,D,d_g,cutoff,maxdim);
        auto mg = vumps_complete_mixed_gauge(canon.AL,canon.D,d_g);

        VumpsUniformResult out;
        out.D = canon.D; out.d_g = d_g; out.AL = canon.AL;
        out.AR = mg.AR; out.C = mg.C; out.AC = mg.AC; out.eta = canon.eta;
        return out;
        }

    // Loads an externally-computed mixed-gauge uniform iMPS (e.g. the
    // return value of vumps_apply_mpo/vumps_imps_sum, from this OR another
    // Chain sharing the same site_types) into `this` Chain's own VUMPS
    // snapshot, so vumps_onsite_expectation/vumps_two_point_correlator
    // (and further vumps_apply_mpo/vumps_imps_sum calls) see it -- C++
    // analogue of pyitensor/vumps.py's own duck typing (VUMPSResult and
    // UniformMPS are interchangeable there since onsite_expectation/
    // two_point_correlator only ever read .sites_uc/.n_uc/.AC/.AR).
    // GL/GR/W/e0/converged/niter_done/gauge_mismatch are NOT restored (no
    // meaning for a summed/MPO-applied state, matching pyitensor's own
    // UniformMPS, which drops them too) -- vumps_excitation_energies is
    // therefore not meaningful after this call without a fresh
    // vumps_ground_state().
    void
    vumps_load_uniform_state(int D, int d_g, std::vector<Cplx> const& AL,
                              std::vector<Cplx> const& AR, std::vector<Cplx> const& C)
        {
        if ((size_t)D*d_g*D != AL.size() || (size_t)D*d_g*D != AR.size() || (size_t)D*D != C.size())
            throw ITError("Chain::vumps_load_uniform_state: AL/AR/C sizes inconsistent with D,d_g");
        vumps_D_ = D; vumps_dg_ = d_g;
        vumps_AL_ = AL; vumps_AR_ = AR; vumps_C_ = C;
        have_vumps_snapshot_ = true;
        have_vumps_exc_env_ = false;
        }

    // Read back `this` Chain's own current VUMPS snapshot AL/AR/C -- e.g.
    // to rescale/feed a state into vumps_imps_sum's own AL_b argument on a
    // DIFFERENT Chain, or to inspect a vumps_load_uniform_state()-loaded
    // state -- requires vumps_ground_state/vumps_load_uniform_state to
    // have been called first on this same Chain.
    void
    vumps_get_snapshot(int& D, int& d_g, std::vector<Cplx>& AL,
                        std::vector<Cplx>& AR, std::vector<Cplx>& C) const
        {
        if (!have_vumps_snapshot_)
            throw ITError("Chain::vumps_get_snapshot: called before vumps_ground_state/"
                           "vumps_load_uniform_state (no converged VUMPS snapshot)");
        D = vumps_D_; d_g = vumps_dg_; AL = vumps_AL_; AR = vumps_AR_; C = vumps_C_;
        }

    // Real-time IBC-window dynamical correlator S(x,t) =
    // <psi0|A_x exp(-iHt)B_0|psi0> of the infinite chain (Milsted/
    // Vanderstraeten, arXiv:1804.09163, Sec. V.1) -- native ITensor v3
    // counterpart of pyitensor/idmrg_window.py's own
    // dynamical_correlator_td (see that module's own extensive docstring
    // for the algorithm and derivation this follows). Requires
    // idmrg_ground_state() to have been called first on this same Chain
    // (have_idmrg_snapshot_ set there).
    //
    // Unlike pyitensor's own window TDVP (which had to hand-roll every
    // window-aware environment/effective-Hamiltonian function, since
    // pyitensor's generic tdvp.py infers a site's Link via a same-Index
    // chain-neighbor lookup that cannot see a window's extra boundary
    // legs -- see idmrg_window.py's own module docstring), the *real*,
    // vendored ITensor v3 TDVP library (TDVP/tdvp.h) already ships an
    // overload taking explicit left/right boundary tensors
    // (tdvp(psi,H,t,LH,RH,sweeps,args), backed by
    // LocalMPO(H,LH,RH,args)) -- so this port only needs to (1) build an
    // ordinary, explicit finite window MPS/MPO by tiling the converged
    // idmrg_U_/idmrg_rows_ n_window times, and (2) hand idmrg_HL_/idmrg_HR_
    // to that existing overload directly, with no bespoke sweep/
    // environment machinery of its own.
    //
    // idmrg_HL_/idmrg_HR_'s own "bra" leg (idmrg_HL_bra_/idmrg_HR_bra_) is
    // minted via sim() in idmrg_extend_HL/HR (an independent, freshly-named
    // Index, unrelated by construction to idmrg_HL_ket_/idmrg_HR_ket_) --
    // but ITensor's own environment-propagation convention, used
    // everywhere else in this codebase (LocalMPO::makeL/makeR, both here
    // and in the vendored library: `dag(prime(psi(i)))`), instead expects
    // a boundary tensor's own bra leg to be *literally* `prime()` of its
    // own ket leg. relabel_bra_to_prime_ket below performs that one-time
    // conversion (a pure Index relabeling, `replaceInds`) so idmrg_HL_/
    // idmrg_HR_ can be handed to the ordinary tdvp()/LocalMPO machinery
    // unchanged -- confirmed necessary and sufficient by hand-derivation of
    // LocalMPO::makeL's own contraction (see this method's own commit
    // history for the derivation): LH is only ever used *directly* (no
    // propagation needed) at the window's own two open ends, so this
    // relabeling is the only adaptation needed, not a deeper redesign.
    //
    // Every precondition/recoverable-failure check in this method and its
    // private helpers below (idmrg_build_window, idmrg_close_array_chain,
    // ...) uses `throw ITError(...)`, never the `Error(...)` macro --
    // Error() calls itensor::error() (itensor/util/error.h), which prints
    // and unconditionally abort()s the whole process without ever
    // actually throwing, so pybind11 never gets a chance to translate it
    // into a catchable Python exception (see check_kpm_moment's own
    // extensive comment above for the fuller derivation and the
    // Python-side symptom this avoids: a hard interpreter-killing SIGABRT
    // instead of a catchable RuntimeError). Every failure mode here (no
    // prior idmrg_ground_state call, x_values reaching outside the
    // window, a bond-dimension mismatch from poor iDMRG convergence) is
    // exactly this kind of recoverable, caller-facing condition.
    TdWindowResult
    td_dynamical_correlator_window(int n_window, std::string const& opname_A,
                                    std::string const& opname_B,
                                    double dt, int nt,
                                    std::vector<int> const& x_values,
                                    int maxdim, double cutoff, int niter,
                                    bool connected, int p_i)
        {
        if (!have_idmrg_snapshot_)
            throw ITError("Chain::td_dynamical_correlator_window called before "
                  "idmrg_ground_state (no converged environment snapshot "
                  "to cap a window with)");
        int n_uc = idmrg_n_uc_;
        if (!(0<=p_i && p_i<n_uc))
            throw ITError("Chain::td_dynamical_correlator_window: p_i out of range");
        if (n_window<1)
            throw ITError("Chain::td_dynamical_correlator_window: n_window must be >= 1");

        auto win = idmrg_build_window(n_window);
        int n = win.n;
        int center = idmrg_window_center(n_window,n_uc,p_i);

        // ITensor's own environment-propagation convention (see this
        // method's own top comment) -- convert once, outside the time
        // loop. Computed *before* perturbing win.psi below, since eshift
        // (just underneath) needs it evaluated on the unperturbed ground
        // window too.
        ITensor LH = idmrg_relabel_bra_to_prime_ket(idmrg_HL_,idmrg_HL_bra_,idmrg_HL_ket_);
        ITensor RH = idmrg_relabel_bra_to_prime_ket(idmrg_HR_,idmrg_HR_bra_,idmrg_HR_ket_);

        // eshift: the window's own total energy (baseline + genuine
        // window physics), measured on the *unperturbed* ground window --
        // mirrors pyitensor/idmrg_window.py's own window_tdvp_step
        // `eshift` fix exactly (see that function's own docstring for the
        // full derivation/justification, and window_total_energy's own
        // docstring point 2). idmrg_HL_/idmrg_HR_ are, just like their
        // pyitensor counterparts, not energy-baseline-subtracted -- they
        // carry a large, macro-iteration-count-dependent additive
        // constant left over from idmrg_ground_state's own growth, so
        // TDVP evolution of the (perturbed) window under the unshifted
        // LH/RH picks up a spurious global phase exp(-i*const*t) that
        // varies run to run for the *same* physical, equally-converged
        // ground state -- confirmed directly on the pyitensor side (see
        // window_tdvp_step's own docstring); this v3 path shares the
        // identical unshifted-LH/RH construction, so it has the same bug.
        // Measured via a throwaway copy of win.psi and a null (t=0) tdvp()
        // step purely to read off its own returned Rayleigh-quotient
        // energy (TDVPWorker's own per-bond expectation value, which
        // already correctly includes the LH/RH boundary caps) -- avoids
        // hand-rolling a fresh LH*mpo*RH sandwich contraction. The
        // *actual* win.psi below is left untouched by this measurement.
        //
        // Deliberately does NOT reuse the caller's own maxdim/cutoff here
        // (unlike the real evolution sweep below): TDVPWorker's per-bond
        // step still runs an SVD split (with the requested truncation)
        // even at t=0 -- exp(0*Heff)=Id makes the *local update* a no-op,
        // but the split itself is not skipped, so reusing a caller-chosen
        // maxdim smaller than the window's own already-converged bond
        // dimension would silently truncate this throwaway copy *before*
        // its energy is read off, biasing eshift (and therefore every
        // reported S(x,t), via the post-hoc exp(+i*eshift*t) correction
        // below) by exactly the discarded weight -- confirmed as a real
        // risk, not hypothetical, for any caller passing a smaller
        // maxdim to td_dynamical_correlator_window than the window's own
        // natural bond dimension (a common pattern: a cheaper maxdim for
        // the time-evolution sweep than the ground-state solve used).
        // cutoff=0 and a generous maxdim make this split lossless (up to
        // floating point), matching pyitensor's own window_total_energy,
        // which contracts exactly with no truncation at all.
        double eshift;
        {
        MPS psi_for_energy = win.psi;
        psi_for_energy.position(1);
        auto sweeps_e = Sweeps(1);
        sweeps_e.maxdim() = 100000;
        sweeps_e.cutoff() = 0.0;
        sweeps_e.niter() = niter;
        Args args_e("Quiet",true,"Silent",true,"NumCenter",2,"Truncate",true,"DoNormalize",true);
        eshift = tdvp(psi_for_energy,win.mpo,Cplx(0,0),LH,RH,sweeps_e,args_e);
        }
        // Documented-not-fixed (code review): this path has no dense-
        // matrix exact cross-check analogous to pyitensor's own
        // test_window_tdvp_step_eshift_matches_exact_dense_evolution
        // (tests/test_idmrg_window_free_fermion.py) -- building one here
        // would need exposing this window's own tensors to Python (or a
        // C++-side dense comparison) that doesn't exist yet. The fix
        // above (untruncated eshift measurement) is verified only by
        // tests/test_idmrg_window_v3.py's own loose reproducibility/
        // consistency checks, not by a strong regression test: manually
        // reintroducing the truncation-coupling bug this block fixes
        // (reusing the caller's own maxdim/cutoff here instead) did NOT
        // reliably fail those checks at the test module's own modest
        // maxm=8 -- the resulting eshift bias was too small there to
        // separate cleanly from ordinary evolution-truncation noise at
        // low maxdim. Isolating it would need a larger maxm and a more
        // extreme maxdim mismatch than currently used, adding real
        // runtime to an already-slow test module; judged disproportionate
        // for now given the fix itself is unambiguously more correct by
        // construction (exact vs. truncated measurement) regardless of
        // whether a dramatic before/after test currently demonstrates it.

        // Sec. V.1 step 3: perturb (B_0|psi>), in place.
        idmrg_window_apply_local_op(win,center,opname_B);
        // Establish a genuine orthogonality center now (via ITensor's own
        // robust QR-sweep-based position(), which -- unlike inner()/
        // innerC() below -- correctly handles win.psi's own two extra
        // dangling boundary legs at sites 1/n): TDVPWorker (TDVP/tdvp.h)
        // would otherwise do this on its own first call, but the
        // norm-save/restore fix just below needs it established *before*
        // that first call too (norm(MPS) requires isOrtho()).
        win.psi.position(1);

        Cplx mean_B = Cplx(0,0);
        std::map<int,Cplx> background;
        if (connected)
            {
            int p_B = (center-1)%n_uc;
            mean_B = idmrg_onsite_expectation(p_B,opname_B);
            for (int x : x_values)
                {
                int p_A = ((center+x-1)%n_uc+n_uc)%n_uc;
                background[x] = mean_B*idmrg_onsite_expectation(p_A,opname_A);
                }
            }

        // The converged unit cell's own dominant right transfer-matrix
        // fixed point (idmrg_window_snapshot_correlator's own right-edge
        // closing weight) depends only on the static, unperturbed ground
        // state -- never on win's own (TDVP-evolving) state -- so it is
        // computed once here, outside the per-time-step loop below,
        // rather than by every one of that loop's own nt calls (each of
        // which would otherwise redo the same power iteration, capped at
        // 2000 steps).
        int p_right = (win.n-1)%n_uc;
        auto rho_after = idmrg_all_right_fixed_points();
        Index rl = idmrg_U_right_[p_right];
        auto rho_flat = idmrg_matrix_to_array(rho_after[p_right],rl,prime(rl));
        int chi_right = dim(rl);

        auto sweeps = Sweeps(1);
        sweeps.maxdim() = maxdim;
        sweeps.cutoff() = cutoff;
        sweeps.niter() = niter;
        Args args("Quiet",!verbose_,"Silent",!verbose_,"NumCenter",2,
                   "Truncate",true,"DoNormalize",true);
        Cplx t = Cplx(0.0,-1.0)*Cplx(dt,0.0);

        TdWindowResult out;
        out.xs = x_values;
        out.ts.resize(nt);
        out.S.assign((size_t)nt*x_values.size(),Cplx(0,0));
        for (int it=0; it<nt; ++it)
            {
            out.ts[it] = it*dt;
            auto snap = idmrg_window_snapshot_correlator(win,opname_A,x_values,center,
                                                           rho_flat,chi_right);
            // Undo the spurious global phase exp(-i*eshift*t) win.psi has
            // picked up from evolving under unshifted LH/RH (see eshift's
            // own comment above) -- applied to the *raw* snapshot value,
            // before background subtraction, since `background` is
            // computed from the static, un-evolved ground state and never
            // carries this phase to begin with.
            Cplx phase = std::exp(Cplx(0.0,1.0)*eshift*out.ts[it]);
            for (size_t ix=0; ix<x_values.size(); ++ix)
                {
                Cplx val = snap[ix]*phase;
                if (connected) val -= background[x_values[ix]];
                out.S[(size_t)it*x_values.size()+ix] = val;
                }
            if (it<nt-1)
                {
                // win.psi is generally not unit-norm (a perturbed,
                // IBC-capped window has no reason to be) -- but
                // ITensorTDVP's own tdvp() always force-renormalizes its
                // input to 1 ("DoNormalize",true, no way to turn off
                // without also disabling the -- needed -- per-bond SVD
                // truncation, see TDVP/tdvp.h's own TDVPWorker). Save/
                // restore the true norm around the call, exactly the fix
                // already established elsewhere in this file for the same
                // issue (quench_tdvp/quench_tdvp_gse/metts_vev's own v(t),
                // see metts_vev's own comment for the full derivation) --
                // but using norm(win.psi) here, *not* those callers' own
                // sqrt(innerC(psi,psi)): innerC()/inner() assume an
                // ordinary open MPS with no dangling legs at sites 1/n,
                // silently mis-contracting (confirmed directly: aborts
                // inside ITensor's own findIndex, "more than one Index
                // found") given win.psi's own two extra boundary legs
                // (idmrg_HL_ket_/idmrg_HR_ket_). norm(MPS), by contrast,
                // just takes the Frobenius norm of the current
                // orthogonality-center tensor (isOrtho() required --
                // hence win.psi.position(1) above, and why this stays
                // correct after every subsequent tdvp() call too, which
                // always leaves win.psi orthogonal, see TDVPWorker's own
                // tail) -- correct regardless of any extra dangling legs,
                // since it never needs to contract them away at all.
                Cplx norm0 = Cplx(norm(win.psi),0.0);
                tdvp(win.psi,win.mpo,t,LH,RH,sweeps,args);
                win.psi.normalize();
                win.psi *= norm0;
                }
            }
        return out;
        }

    private:
    // v2's plain "MPS(sites_)" (no InitState) doesn't carry over as a bare
    // MPS(SiteSet)/InitState-based product state: with sites_ built
    // ConserveQNs=false (see get_sites.h's comment for why), a *product*
    // starting state (every tensor bond dimension 1, e.g. all sites in
    // their first basis state) leaves DMRG's local two-site Davidson step
    // completely unable to find the true, entangled ground state -- the
    // classic "DMRG stuck at an exact product eigenstate" trap, confirmed
    // directly (a Heisenberg chain from this kind of start converges to
    // the trivial product energy no matter how large a noise term is
    // added; the same Hamiltonian from an actual random MPS converges to
    // the correct, entangled ground state immediately). randomMPS(sites_,m)
    // (m = maxm_, i.e. the wavefunction's own target bond dimension) is
    // what actually reproduces v2's effective behavior of an unconstrained
    // search that isn't trapped at a symmetric starting point.
    MPS
    default_mps() const { return randomMPS(sites_,maxm_); }

    // -- METTS helpers (Chain::metts_vev, public, above) --

    // A CPS (classical product state: bond dimension 1 throughout) built
    // from per-site rank-1 ITensors over each site's bare physical index
    // (no Link legs yet). Fresh dim-1 Link indices are attached here via
    // setElt() -- the same "outer product with a one-hot tensor" idiom
    // MPS's own init_tensors() (mps.cc, backing the MPS(InitState)
    // constructor) uses to build a product state from named basis
    // states, generalized here to an arbitrary per-site vector (a
    // METTS collapse outcome is a diagHermitian() eigenvector, not
    // necessarily one of a site type's named states). leftLim/rightLim
    // are set to match MPS(InitState)'s own convention (orthogonality
    // center at site 1) since a bond-dimension-1 chain is trivially
    // canonical everywhere once each site vector is unit-normalized (true
    // here: every vector handed in is either a diagHermitian()
    // eigenvector or built from one).
    MPS
    build_product_state(std::vector<ITensor> const& site_vectors) const
        {
        int N = sites_.length();
        MPS psi(N);
        std::vector<Index> links;
        links.reserve(N-1);
        for (int k=1; k<N; ++k) links.push_back(Index(1,TagSet("Link,l="+std::to_string(k))));
        for (int i=1; i<=N; ++i)
            {
            ITensor T = site_vectors[i-1];
            if (i>1) T *= setElt(links[i-2](1));
            if (i<N) T *= setElt(links[i-1](1));
            psi.set(i,T);
            }
        psi.leftLim(0);
        psi.rightLim(2);
        return psi;
        }

    // Per-(opname,site) diagHermitian() result, cached by
    // metts_build_eigcache() below so random_cps()/collapse_to_cps()
    // don't recompute the same site operator's eigendecomposition from
    // scratch at every one of a METTS run's nwarmup+nsamples iterations
    // -- a given (opname,i)'s eigenbasis never changes across the whole
    // sampling run. D (the eigenvalue tensor diagHermitian also produces)
    // is never read anywhere in this file, so only U and its eig index
    // are kept.
    struct MettsEigBasis { ITensor U; Index eig; };
    using MettsEigCache = std::map<std::pair<std::string,int>,MettsEigBasis>;

    MettsEigCache
    metts_build_eigcache(std::vector<std::string> const& basis_ops) const
        {
        MettsEigCache cache;
        int N = sites_.length();
        for (auto const& opname : basis_ops)
            for (int i=1; i<=N; ++i)
                {
                auto key = std::make_pair(opname,i);
                if (cache.count(key)) continue;
                auto opT = sites_.op(opname,i);
                ITensor U,D;
                diagHermitian(opT,U,D);
                cache[key] = MettsEigBasis{U,uniqueIndex(U,opT)};
                }
        return cache;
        }

    // A random CPS: at each site, uniformly picks one eigenstate of the
    // named single-site Hermitian operator opname (e.g. "Sz") -- a valid
    // seed for the METTS Markov chain regardless of choice (the chain's
    // stationary distribution doesn't depend on the starting CPS, only
    // how many warmup steps are needed to reach it -- White & Stoudenmire,
    // arXiv:1002.1305, Sec. 2).
    MPS
    random_cps(std::string const& opname, std::mt19937_64& rng,
               MettsEigCache const& eigcache) const
        {
        int N = sites_.length();
        std::vector<ITensor> vectors(N);
        for (int i=1; i<=N; ++i)
            {
            auto const& eb = eigcache.at(std::make_pair(opname,i));
            std::uniform_int_distribution<int> dist(1,itensor::dim(eb.eig));
            int k = dist(rng);
            vectors[i-1] = eb.U*setElt(eb.eig(k));
            }
        return build_product_state(vectors);
        }

    // Sequential ("perfect sampling") collapse of MPS wf onto a new CPS
    // built from eigenstates of the single-site Hermitian operator
    // opname -- direct port of pyitensor/metts.py's collapse_to_cps(),
    // see its own docstring for the full derivation; the paper's "quantum
    // measurement" collapse step (Sec. 2.2) generalized from a single
    // qubit measurement to the whole chain via one left-to-right sweep.
    //
    // Canonicalizing to site 1 first (wf.position(1)) makes every site
    // i>1 right-orthonormal; contracting the running "collapsed-so-far"
    // amplitude L into site i's tensor and rotating into the eigenbasis
    // via dag(U) gives `rot`, whose diagonal marginal probabilities
    // (summed over whatever's left of the chain to the right) fall out
    // directly as the diagonal of rot*dag(prime(rot,eig)) -- no explicit
    // partial trace needed, since right-orthonormality already performs
    // it, and (unlike the manual-index-bookkeeping pyitensor port) no
    // special-casing of the last site either: ITensor's own automatic
    // index contraction handles "no right link left" (i==N) the same way
    // as any other site.
    MPS
    collapse_to_cps(MPS wf, std::string const& opname, std::mt19937_64& rng,
                     MettsEigCache const& eigcache) const
        {
        int N = sites_.length();
        wf.position(1);
        std::vector<ITensor> vectors(N);
        ITensor L; // default-constructed/"null" (operator bool()==false) until i>1's first assignment
        for (int i=1; i<=N; ++i)
            {
            auto const& eb = eigcache.at(std::make_pair(opname,i));
            auto const& U = eb.U;
            auto const& eig = eb.eig;
            int dim_i = itensor::dim(eig);

            ITensor T = wf.A(i);
            if (L) T *= L;
            ITensor rot = dag(U)*T;
            ITensor probs = rot*dag(prime(rot,eig));

            std::vector<double> probs_raw(dim_i), p(dim_i);
            double total = 0.0;
            for (int k=1; k<=dim_i; ++k)
                {
                probs_raw[k-1] = probs.eltC(eig(k),prime(eig)(k)).real();
                total += probs_raw[k-1];
                }
            for (int k=0; k<dim_i; ++k) p[k] = probs_raw[k]/total;
            std::discrete_distribution<int> dist(p.begin(),p.end());
            int k = dist(rng)+1; // 1-based

            vectors[i-1] = U*setElt(eig(k));
            if (i<N)
                {
                ITensor Lnew = rot*setElt(eig(k));
                L = Lnew*(1.0/std::sqrt(probs_raw[k-1]));
                }
            }
        return build_product_state(vectors);
        }

    // Thin wrappers around v3's applyMPO() that always restore the result's
    // physical index to the standard, unprimed convention. Plain applyMPO()
    // doesn't do this itself: it contracts K against whichever leg of K
    // actually matches x's own physical index, so if x already carries a
    // primed physical leg (e.g. the result of an *earlier* unwrapped
    // MPO application), the *new* result comes out the other way round
    // (unprimed in this case, but in general whatever K's other leg is) --
    // numerically consistent for a single use, but silently breaks the very
    // next sum()/mult of two such results against each other (sum() and
    // operator+= require literally matching index structure). Confirmed
    // directly: cvm_dynamical_correlator()'s S2|wf0_> came back with a
    // primed physical leg while apply_mpo(A,x=that state)'s came back
    // unprimed, and bicstab()'s very first "sum(b,-1*apply_mpo(A,x))"
    // aborted with "different index structure". Every applyMPO() call in
    // this class goes through here instead so that invariant -- "every MPS
    // this class hands back or feeds into sum()/bicstab() has an unprimed
    // physical index" -- always holds, not just in the call sites that
    // happened to get exercised so far.
    //
    // TagSet("Site") re-parses its string argument every call (same cost
    // already profiled and fixed for TagSet("a") in arnoldi_smallest_real
    // below -- see that function's own comment); apply_mpo() is the
    // single most-invoked helper in this file (every MPO application:
    // KPM's Chebyshev recursion, time evolution, the generalized-
    // eigenproblem solver, bicstab/CVM, custom_exp/evoloperator, ...), so
    // a function-local static avoids re-parsing it on every one of those
    // calls instead of just the iDMRG-specific functions the earlier pass
    // covered.
    MPS
    apply_mpo(MPO const& K, MPS const& x, Args const& args) const
        {
        static const TagSet site_tag("Site");
        auto out = applyMPO(K,x,args);
        out.noPrime(site_tag);
        return out;
        }

    MPS
    apply_mpo(MPO const& K, MPS const& x, MPS const& x0, Args const& args) const
        {
        static const TagSet site_tag("Site");
        auto out = applyMPO(K,x,x0,args);
        out.noPrime(site_tag);
        return out;
        }

    // Restarted Arnoldi on a matrix-free operator over ITensor "vectors"
    // (the two-site blocks of nhdmrg()): builds a krylovdim-step Krylov
    // space with modified Gram-Schmidt (plus one re-orthogonalization
    // pass), diagonalizes the small Hessenberg matrix once per build with
    // ITensor's dense non-hermitian eigen(), keeps the Ritz pair selected
    // by `sel`, and restarts from that Ritz vector. This is the C++
    // stand-in for the reference's KrylovKit eigsolve(...;
    // ishermitian=false, which=:SR) local solver -- like there, it runs
    // at low accuracy per bond (the outer DMRG sweeps do the actual
    // converging). Remaining restarts are skipped when the standard
    // Arnoldi residual bound ||A y - lambda y|| = |h(m+1,m)|*|last
    // Ritz-vector component| certifies the build already converged --
    // that check reuses the one end-of-build diagonalization, so it is
    // free, and it halves the matvec count on already-converged bonds.
    // A per-step stop-on-stable-Ritz-value exit (the idea pyitensor's
    // hermitian _lanczos_ground_state uses) was tried and reverted:
    // unlike the symmetric tridiagonal solve there, the non-hermitian
    // Hessenberg eigen() per step costs as much as the matvec it hopes
    // to save at these block sizes -- measured directly on a 20-site
    // chain, the per-step variant REGRESSED nhdmrg from ~30s to ~104s.
    // (The vendored ITensor/itensor/iterativesolvers.h has its own
    // arnoldi(), but it only supports LargestMagnitude/SmallestReal
    // selection -- neither the Closest anchoring nor the SRTieBreak
    // below is expressible with it, and those selections are
    // load-bearing, see the Sel notes.)
    //
    // Ritz-value selection (`sel`):
    //  - SR: smallest real part (the reference's :SR).
    //  - Closest: closest to `target`. nhdmrg() uses this for the adjoint
    //    (left-eigenvector) solve, with target=conj(lambda_right): the
    //    reference selects :SR independently on both sides, but when
    //    several eigenvalues share the smallest real part (e.g. a complex-
    //    conjugate pair of a PT-symmetric Hamiltonian, where Re-ordering
    //    cannot distinguish a+bi from a-bi) the two independent solves can
    //    lock onto *different* members of that degenerate set -- and the
    //    left/right vectors of different eigenstates are mutually
    //    biorthogonal, so <thl|thr> collapses to ~0 and the sweep never
    //    converges (confirmed directly on a staggered-imaginary-field XX
    //    chain). Anchoring the left solve to the right solve's eigenvalue
    //    pins both onto the same eigenpair.
    //  - SRTieBreak: smallest real part, but among Ritz values whose real
    //    parts sit within a small tolerance of the *global* minimum, the
    //    one closest to `target` (the previous bond's eigenvalue). Used
    //    by the right solve from the second bond on, so a Re-degenerate
    //    +-Im pair cannot make consecutive bond solves flip branch
    //    mid-sweep (the same failure mode as above, one level up: each
    //    flip re-targets the truncation onto a different eigenstate, and
    //    the sweep stalls at a non-eigenstate -- also confirmed on the
    //    same XX chain). The global-min-then-candidates formulation is
    //    deliberately identical to pyitensor/nhdmrg.py's _select_ritz
    //    (an earlier greedy running-best variant here could settle on a
    //    different member of a degenerate cluster than the Python port,
    //    silently diverging the backends on exactly the spectra the
    //    tie-break exists for).
    enum class Sel { SR, Closest, SRTieBreak };

    // Ritz-value selection given a Krylov-subspace-size-m eigenvalue list;
    // factored out of arnoldi_smallest_real so its early-exit convergence
    // check (see that method's own early_tol comment) can reuse *exactly*
    // the same selection logic as the final post-loop selection, rather
    // than risking the two silently drifting apart. Pure post-processing
    // of an already-computed eigenvalue list, no ITensor calls -- the
    // global-min-then-candidates SRTieBreak formulation is deliberately
    // identical to pyitensor/nhdmrg.py's _select_ritz (an earlier greedy
    // running-best variant here could settle on a different member of a
    // degenerate cluster than the Python port, silently diverging the
    // backends on exactly the spectra the tie-break exists for).
    static int
    arnoldi_select_kbest(std::vector<Cplx> const& ev, Sel sel, Cplx target)
        {
        int m = (int)ev.size();
        int kbest = 1;
        if (sel==Sel::Closest)
            {
            for (int k=2;k<=m;++k)
                if (std::abs(ev[k-1]-target)<std::abs(ev[kbest-1]-target))
                    kbest = k;
            }
        else
            {
            double remin = ev[0].real();
            for (int k=2;k<=m;++k) remin = std::min(remin,ev[k-1].real());
            if (sel==Sel::SRTieBreak)
                {
                double degtol = 1e-6*(1.0+std::abs(remin));
                kbest = 0;
                for (int k=1;k<=m;++k)
                    if (ev[k-1].real()<remin+degtol)
                        if (kbest==0 ||
                            std::abs(ev[k-1]-target)<std::abs(ev[kbest-1]-target))
                            kbest = k;
                }
            else
                {
                for (int k=2;k<=m;++k)
                    if (ev[k-1].real()<ev[kbest-1].real()) kbest = k;
                }
            }
        return kbest;
        }

    // early_tol<0 (default): disabled, matching this method's original
    // behavior exactly -- every pre-existing caller (nhdmrg_one_sweep's
    // two calls above) always builds the *full* krylovdim-size Krylov
    // subspace before checking convergence, same as before this parameter
    // existed. early_tol>=0 (idmrg_local_solve's own call, below) checks
    // the *same* Arnoldi residual bound this method's between-restart
    // check already trusts to declare "good enough, skip further
    // restarts" -- but after every Krylov-vector addition instead of only
    // once the full krylovdim has been built, so a local solve that's
    // already converged (the common case away from a near-degenerate
    // spectrum) stops paying for further matvecs it doesn't need. This
    // isn't new solver math: the "happy breakdown" branch just above
    // already proves an early-terminated (m<krylovdim) Krylov subspace is
    // handled correctly by the existing post-loop selection/build code
    // (whatever final m is, that code already runs on it) -- early_tol
    // just generalizes that same trigger from "residual is exactly zero"
    // to "residual is below tolerance", reusing arnoldi_select_kbest so
    // the *decision* of which Ritz pair is "best" can never disagree
    // between the early check and the eventual real answer.
    template <typename Fn>
    std::pair<Cplx,ITensor>
    arnoldi_smallest_real(Fn&& A, ITensor x0, int krylovdim, int restarts,
                          Sel sel=Sel::SR, Cplx target=0,
                          double early_tol=-1.0) const
        {
        // TagSet("a") parses its string argument every call (confirmed
        // via profiling to be a real, measurable cost: TagSet parsing +
        // its own strtol calls showed up at several percent of total
        // iDMRG runtime) -- a function-local static parses it exactly
        // once and every subsequent call just copies the already-parsed
        // TagSet. idmrg_extend_HL/HR below apply the same fix to their
        // own TagSet("Site").
        static const TagSet a_tag("a");
        Cplx lambda = 0;
        for (int r=0;r<restarts;++r)
            {
            double nx = norm(x0);
            if (nx<1e-14) break; // degenerate start, keep whatever we have
            x0 /= nx;
            std::vector<ITensor> V;
            V.push_back(x0);
            std::vector<std::vector<Cplx>> h(krylovdim+1,
                    std::vector<Cplx>(krylovdim,Cplx(0,0)));
            int m = 0;
            for (int j=0;j<krylovdim;++j)
                {
                auto w = A(V.at(j));
                for (int i=0;i<=j;++i)
                    {
                    auto c = eltC(dag(V.at(i))*w);
                    h.at(i).at(j) = c;
                    w -= c*V.at(i);
                    }
                for (int i=0;i<=j;++i)
                    {
                    auto c = eltC(dag(V.at(i))*w);
                    h.at(i).at(j) += c;
                    w -= c*V.at(i);
                    }
                m = j+1;
                double nw = norm(w);
                h.at(j+1).at(j) = nw;
                if (nw<1e-13) break; // happy breakdown: invariant subspace
                // Early-exit convergence check (see early_tol's own
                // comment above). Skipped below m=8 (a residual bound
                // from a tiny subspace is not a reliable signal) and only
                // run every 3rd iteration above that: measured directly
                // (S=1/2 Heisenberg iDMRG) that most local solves in the
                // later, expensive macro-iterations (largest bond
                // dimension) do NOT converge before the full krylovdim,
                // so checking every single iteration mostly just pays for
                // an extra small eigen() + Index/TagSet construction with
                // no payoff in the cases that dominate wall time --
                // checking less often keeps that wasted cost small while
                // still catching genuine early convergence (common in the
                // cheap, small-bond-dimension early macro-iterations)
                // within at most 2 extra Krylov vectors of when it
                // actually happens. Always skipped when early_tol<0
                // (every pre-existing caller), so this adds zero cost to
                // nhdmrg_one_sweep's own two arnoldi_smallest_real calls.
                if (early_tol>=0 && m>=8 && m%3==0)
                    {
                    auto a_chk = Index(m,a_tag);
                    auto Hm_chk = ITensor(prime(a_chk),a_chk);
                    for (int ii=0;ii<m;++ii)
                    for (int jj=0;jj<m;++jj)
                        if (ii<=jj+1) Hm_chk.set(prime(a_chk)(ii+1),a_chk(jj+1),h.at(ii).at(jj));
                    ITensor Wchk,Dchk;
                    eigen(Hm_chk,Wchk,Dchk);
                    auto cchk = commonIndex(Wchk,Dchk);
                    std::vector<Cplx> ev_chk(m);
                    for (int k=1;k<=m;++k) ev_chk[k-1] = eltC(Dchk,cchk(k),prime(cchk)(k));
                    int kbest_chk = arnoldi_select_kbest(ev_chk,sel,target);
                    Cplx ebest_chk = ev_chk[kbest_chk-1];
                    double resid_chk = nw*std::abs(eltC(Wchk,a_chk(m),cchk(kbest_chk)));
                    if (resid_chk<early_tol*(1.0+std::abs(ebest_chk))) break;
                    }
                if (j+1<krylovdim) V.push_back(w/nw);
                }
            auto a = Index(m,a_tag);
            auto Hm = ITensor(prime(a),a);
            for (int i=0;i<m;++i)
            for (int j=0;j<m;++j)
                if (i<=j+1) Hm.set(prime(a)(i+1),a(j+1),h.at(i).at(j));
            ITensor W,D;
            eigen(Hm,W,D);
            auto c = commonIndex(W,D);
            std::vector<Cplx> ev(m);
            for (int k=1;k<=m;++k) ev[k-1] = eltC(D,c(k),prime(c)(k));
            int kbest = arnoldi_select_kbest(ev,sel,target);
            Cplx ebest = ev[kbest-1];
            ITensor xnew;
            for (int i=0;i<m;++i)
                {
                auto coef = eltC(W,a(i+1),c(kbest));
                if (!xnew) xnew = coef*V.at(i);
                else xnew += coef*V.at(i);
                }
            lambda = ebest;
            x0 = xnew;
            // standard Arnoldi residual bound for the selected Ritz pair:
            // ||A y - lambda y|| = |h(m+1,m)| * |last component of the
            // Hessenberg eigenvector|; when the build already converged,
            // further restarts would just rebuild the same subspace
            double resid_est = h.at(m).at(m-1).real()
                              *std::abs(eltC(W,a(m),c(kbest)));
            if (resid_est<1e-10*(1.0+std::abs(ebest))) break;
            }
        return {lambda,x0};
        }

    // -- iDMRG private helpers (Chain::idmrg_ground_state, public, above) --
    // A direct C++ translation of pyitensor/idmrg.py's _term_site_matrices/
    // _classify_terms/_active_channels_at/_onsite_matrix/_build_automaton --
    // see that module's own docstrings for the algorithm; comments here
    // only note where the C++ port genuinely differs (dense-array-only
    // automaton storage, 1-based sites, no fermionic support yet).

    // Dense (row-major, (in,out)-convention) matrix for a named single-site
    // operator at 0-based sublattice position p0based -- the C++ analogue
    // of pyitensor's SiteType.matrix(name), read directly off
    // SiteSet::op() (unprimed leg "In", primed leg "Out", same convention
    // _term_site_matrices's own docstring cross-references).
    std::vector<Cplx>
    idmrg_op_dense(int p0based, std::string const& name) const
        {
        int site1 = p0based+1;
        int d = dim(sites_.si(site1));
        auto Op = sites_.op(name,site1);
        auto s = sites_.si(site1);
        auto sp = prime(s);
        std::vector<Cplx> M(d*d,Cplx(0,0));
        for (int i=1;i<=d;++i)
        for (int j=1;j<=d;++j)
            M[(i-1)*d+(j-1)] = eltC(Op,s(i),sp(j));
        return M;
        }

    static std::vector<Cplx>
    idmrg_matmul(std::vector<Cplx> const& A, std::vector<Cplx> const& B, int d)
        {
        std::vector<Cplx> C(d*d,Cplx(0,0));
        for (int i=0;i<d;++i)
        for (int j=0;j<d;++j)
            {
            Cplx acc(0,0);
            for (int k=0;k<d;++k) acc += A[i*d+k]*B[k*d+j];
            C[i*d+j] = acc;
            }
        return C;
        }

    // Classifies every term (terms_intra ++ terms_inter, MOTerm's own
    // 1-based sites converted to 0-based here to match idmrg.py's own
    // arithmetic exactly) into 1-site (onsite) and 2-site (bonds) pieces,
    // composing multiple same-site operator factors via reversed matrix
    // product exactly like _term_site_matrices does (mats[-1] first,
    // then left-multiplied by each earlier factor in reverse -- i.e. for
    // factors written [A,B,C] at one site, the composed matrix is C@B@A).
    void
    idmrg_classify_terms(std::vector<MOTerm> const& terms_intra,
                          std::vector<MOTerm> const& terms_inter, int n_uc,
                          std::vector<IdmrgOnsite>& onsite,
                          std::vector<IdmrgBond>& bonds) const
        {
        std::vector<MOTerm> all_terms = terms_intra;
        all_terms.insert(all_terms.end(),terms_inter.begin(),terms_inter.end());
        for (auto const& term : all_terms)
            {
            std::map<int,std::vector<std::string>> by_site;
            for (auto const& f : term.factors) by_site[f.site-1].push_back(f.name);
            std::vector<int> rel_sites;
            for (auto const& kv : by_site) rel_sites.push_back(kv.first);
            std::sort(rel_sites.begin(),rel_sites.end());
            if (rel_sites.size()>2)
                Error("Chain::idmrg_ground_state: a term touches more than "
                      "2 distinct sites -- only 1- and 2-site terms are "
                      "supported");
            std::vector<std::vector<Cplx>> mats;
            for (int rel : rel_sites)
                {
                int p = ((rel % n_uc)+n_uc) % n_uc;
                auto const& names = by_site[rel];
                int d = dim(sites_.si(p+1));
                auto combined = idmrg_op_dense(p,names.back());
                for (int k=(int)names.size()-2;k>=0;--k)
                    combined = idmrg_matmul(combined,idmrg_op_dense(p,names[k]),d);
                mats.push_back(std::move(combined));
                }
            if (rel_sites.size()==1)
                onsite.push_back(IdmrgOnsite{rel_sites[0],term.coef,mats[0]});
            else if (rel_sites.size()==2)
                {
                int a = rel_sites[0], b = rel_sites[1];
                if (a>=n_uc)
                    Error("Chain::idmrg_ground_state: internal error -- "
                          "inter-cell term does not touch the central cell");
                bonds.push_back(IdmrgBond{a,b,mats[0],mats[1],false,term.coef});
                }
            }
        }

    // [(bond_index,r), ...] for every 2-site term with a "pending"
    // instance active just before absorbing sublattice p -- exactly
    // idmrg.py's own _active_channels_at, same formula/derivation (see
    // that function's own docstring for why "+1").
    std::vector<std::pair<int,int>>
    idmrg_active_channels_at(std::vector<IdmrgBond> const& bonds, int n_uc, int p) const
        {
        std::vector<std::pair<int,int>> out;
        for (int bi=0;bi<(int)bonds.size();++bi)
            {
            int reach = bonds[bi].rel_b - bonds[bi].rel_a;
            for (int r=1;r<=reach;++r)
                {
                int m = ((bonds[bi].rel_b - r + 1) % n_uc + n_uc) % n_uc;
                if (m==p) out.push_back({bi,r});
                }
            }
        return out;
        }

    std::vector<IdmrgChan>
    idmrg_channels_at(std::vector<IdmrgBond> const& bonds, int n_uc, int p) const
        {
        std::vector<IdmrgChan> chans;
        chans.push_back(IdmrgChan{0,-1,-1}); // S
        chans.push_back(IdmrgChan{1,-1,-1}); // F
        for (auto const& br : idmrg_active_channels_at(bonds,n_uc,p))
            chans.push_back(IdmrgChan{2,br.first,br.second});
        return chans;
        }

    // Sum of coef*mat over every onsite term attributed to sublattice p --
    // an empty vector (not a zero matrix) means "no onsite term here",
    // exactly like idmrg.py's own _onsite_matrix returning None.
    std::vector<Cplx>
    idmrg_onsite_matrix(std::vector<IdmrgOnsite> const& onsite, int p, int d) const
        {
        std::vector<Cplx> mat;
        for (auto const& o : onsite)
            {
            if (o.rel != p) continue;
            if (mat.empty())
                {
                mat.resize(d*d);
                for (int k=0;k<d*d;++k) mat[k] = o.coef*o.mat[k];
                }
            else
                for (int k=0;k<d*d;++k) mat[k] += o.coef*o.mat[k];
            }
        return mat;
        }

    // The dense transition-matrix content for sublattice p, given its own
    // left/right channel lists -- a direct translation of
    // _build_automaton's own big if/elif transition-rule chain (see that
    // function's docstring for the derivation of each case). Fermionic
    // strings are not supported yet (the last branch always uses Id,
    // never sites_.op("F",...); idmrg_classify_terms always sets
    // carry_ferm=false so this is never exercised differently).
    IdmrgAutomatonRow
    idmrg_build_row(int p, int n_uc, std::vector<IdmrgBond> const& bonds,
                     std::vector<IdmrgOnsite> const& onsite,
                     std::vector<IdmrgChan> const& left,
                     std::vector<IdmrgChan> const& right) const
        {
        int d = dim(sites_.si(p+1));
        auto onsite_mat = idmrg_onsite_matrix(onsite,p,d);
        IdmrgAutomatonRow row;
        row.p=p; row.d=d; row.left_n=(int)left.size(); row.right_n=(int)right.size();
        row.flat.assign((size_t)row.left_n*row.right_n*d*d,Cplx(0,0));
        auto setmat = [&](int li,int ri,std::vector<Cplx> const& mat)
            {
            size_t base = ((size_t)li*row.right_n+ri)*d*d;
            for (int k=0;k<d*d;++k) row.flat[base+k] = mat[k];
            };
        auto eye = [&]()
            {
            std::vector<Cplx> I(d*d,Cplx(0,0));
            for (int i=0;i<d;++i) I[i*d+i] = Cplx(1,0);
            return I;
            };
        for (int li=0;li<row.left_n;++li)
        for (int ri=0;ri<row.right_n;++ri)
            {
            auto const& lch = left[li];
            auto const& rch = right[ri];
            bool l_triv = (lch.kind==0 || lch.kind==1);
            bool r_triv = (rch.kind==0 || rch.kind==1);
            if (lch.kind==0 && rch.kind==0)
                {
                setmat(li,ri,eye());
                }
            else if (lch.kind==1 && rch.kind==1)
                {
                setmat(li,ri,eye());
                }
            else if (lch.kind==0 && rch.kind==1)
                {
                // This site's own onsite term: starts and completes in one
                // step, direct from S into the accumulator F -- NOT added
                // onto F,F's own self-loop. This is a real, confirmed drift
                // bug fix (found auditing this port against the current
                // pyitensor/idmrg.py, whose own _build_periodic_mpo
                // docstring documents finding and fixing exactly this
                // mistake: "an earlier, wrong version of this function did
                // exactly that -- mat = Id; if onsite_mat is not None: mat
                // += onsite_mat on the F,F entry", which silently drops
                // every onsite term for a chain with no bond terms at all
                // (W stays block-diagonal between {S} and {F,pending...},
                // so <S|W^N|F> is identically 0 for every N), and produces
                // an exponential-in-iteration blow-up once at least one
                // bond term does activate F (each further absorbed site
                // re-adds the onsite content into an already-summed F
                // channel instead of summing it in exactly once) -- see
                // that docstring's own two confirmed reproducers. Putting
                // it here (S,F) instead is the standard upper-triangular-
                // in-channel-index automaton-MPO construction for summing
                // an onsite term into a running total exactly once per site.
                if (!onsite_mat.empty()) setmat(li,ri,onsite_mat);
                }
            else if (lch.kind==0 && !r_triv)
                {
                auto const& b = bonds[rch.bond_index];
                if (rch.r==b.rel_b-b.rel_a && b.rel_a==p)
                    {
                    std::vector<Cplx> mat(d*d);
                    for (int k=0;k<d*d;++k) mat[k] = b.coef*b.mat_a[k];
                    setmat(li,ri,mat);
                    }
                }
            else if (rch.kind==1 && !l_triv)
                {
                auto const& b = bonds[lch.bond_index];
                if (lch.r==1) setmat(li,ri,b.mat_b);
                }
            else if (!l_triv && !r_triv)
                {
                if (lch.bond_index==rch.bond_index && rch.r==lch.r-1)
                    setmat(li,ri,eye()); // carry_ferm always false in this port
                }
            }
        return row;
        }

    // Fresh ITensor from row's dense data (bulk case: rank 4). left_idx's
    // dim must equal row.left_n, right_idx's row.right_n, phys_in/out's
    // row.d -- the caller is responsible for supplying Index objects of
    // the right dimension (idmrg_ground_state always does, by
    // construction: cross_idx/HL_mpo/HR_mpo are minted or tracked with
    // exactly these dimensions).
    ITensor
    idmrg_make_W(IdmrgAutomatonRow const& row, Index left_idx, Index right_idx,
                 Index phys_in, Index phys_out) const
        {
        ITensor T(left_idx,phys_in,phys_out,right_idx);
        for (int li=0;li<row.left_n;++li)
        for (int ri=0;ri<row.right_n;++ri)
        for (int si=0;si<row.d;++si)
        for (int so=0;so<row.d;++so)
            {
            Cplx v = row.flat[((size_t)li*row.right_n+ri)*row.d*row.d + si*row.d+so];
            if (v != Cplx(0,0))
                T.set({left_idx(li+1),phys_in(si+1),phys_out(so+1),right_idx(ri+1)},v);
            }
        return T;
        }

    // Boundary tensor for the very first-ever micro-step's left side:
    // row's "S" row only (li=0) -- the C++ analogue of idmrg.py's
    // _project_channel(W_bulk[0],"left",0), but built directly from dense
    // data instead of projecting an existing persistent ITensor (this
    // port never builds one, see idmrg_ground_state's own top comment).
    ITensor
    idmrg_make_W_start(IdmrgAutomatonRow const& row, Index right_idx,
                        Index phys_in, Index phys_out) const
        {
        ITensor T(phys_in,phys_out,right_idx);
        int li = 0; // S
        for (int ri=0;ri<row.right_n;++ri)
        for (int si=0;si<row.d;++si)
        for (int so=0;so<row.d;++so)
            {
            Cplx v = row.flat[((size_t)li*row.right_n+ri)*row.d*row.d + si*row.d+so];
            if (v != Cplx(0,0))
                T.set({phys_in(si+1),phys_out(so+1),right_idx(ri+1)},v);
            }
        return T;
        }

    // Boundary tensor for the very first-ever micro-step's right side:
    // row's "F" column only (ri=1) -- analogue of
    // _project_channel(W_bulk[n_uc-1],"right",1).
    ITensor
    idmrg_make_W_end(IdmrgAutomatonRow const& row, Index left_idx,
                      Index phys_in, Index phys_out) const
        {
        ITensor T(left_idx,phys_in,phys_out);
        int ri = 1; // F
        for (int li=0;li<row.left_n;++li)
        for (int si=0;si<row.d;++si)
        for (int so=0;so<row.d;++so)
            {
            Cplx v = row.flat[((size_t)li*row.right_n+ri)*row.d*row.d + si*row.d+so];
            if (v != Cplx(0,0))
                T.set({left_idx(li+1),phys_in(si+1),phys_out(so+1)},v);
            }
        return T;
        }

    // -- VUMPS / excitation ansatz private helpers (Chain::vumps_ground_state
    // and Chain::vumps_excitation_energies, public, above) --
    //
    // Dense, row-major complex linear algebra throughout (never ITensor
    // tensor-network objects) -- see VumpsResult's own comment for why.
    // Shape conventions (matching pyitensor/vumps.py's and
    // pyitensor/idmrg_excitations.py's own numpy array conventions
    // exactly, so this is a line-for-line port, not a fresh derivation):
    //   AL/AR/AC/B: (D,d_g,D), flattened row-major as (l*d_g+p)*D+r.
    //   C, and every "environment"/"fixed point" (GL,GR,r,l,...): (D,D),
    //   flattened row-major as i*D+j.
    //   A transfer tensor E: (D,D,D,D) axes (l,L,r,R) -- l/r the ket legs,
    //   L/R the bra legs -- flattened row-major as ((l*D+L)*D+r)*D+R.
    //   A grouped automaton W: (Dw,Dw,d_g,d_g) axes (chan_l,chan_r,
    //   phys_in,phys_out), flattened row-major as ((l*Dw+r)*d_g+si)*d_g+so
    //   -- channel 0="S" (start), channel 1="F" (accumulator), channels
    //   2.. "pending" (one per reach-1 bond), same convention
    //   idmrg_channels_at/idmrg_build_row already use.
    // An (D,D,D,D) transfer tensor E, reshaped (D*D,D*D) with row=(l,L),
    // col=(r,R), is exactly the row-major flattening above -- i.e. E
    // itself, reinterpreted with n=D*D, already IS that matrix, needing
    // no data movement (used directly by vx_apply_transfer/
    // vx_dominant_right_fixed_point below).

    static std::vector<Cplx>
    vx_eye(int n)
        {
        std::vector<Cplx> I((size_t)n*n,Cplx(0,0));
        for (int i=0;i<n;++i) I[i*n+i] = Cplx(1,0);
        return I;
        }

    static std::vector<Cplx>
    vx_dagger(std::vector<Cplx> const& A, int m, int n) // A:(m,n) row-major -> (n,m) row-major
        {
        std::vector<Cplx> B((size_t)n*m);
        for (int i=0;i<m;++i)
        for (int j=0;j<n;++j)
            B[j*m+i] = std::conj(A[i*n+j]);
        return B;
        }

    static std::vector<Cplx>
    vx_hermitize(std::vector<Cplx> const& A, int n)
        {
        std::vector<Cplx> H((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            H[i*n+j] = (A[i*n+j] + std::conj(A[j*n+i]))/2.0;
        return H;
        }

    static std::vector<Cplx>
    vx_matmul(std::vector<Cplx> const& A, int am, int ak,
              std::vector<Cplx> const& B, int bk, int bn)
        {
        if (ak!=bk) Error("Chain::vumps: vx_matmul inner dimension mismatch");
        std::vector<Cplx> C((size_t)am*bn,Cplx(0,0));
        for (int i=0;i<am;++i)
        for (int j=0;j<bn;++j)
            {
            Cplx acc(0,0);
            for (int k=0;k<ak;++k) acc += A[i*ak+k]*B[k*bn+j];
            C[i*bn+j] = acc;
            }
        return C;
        }

    // T4 (D,D,D,D) reinterpreted as an (n,n) matrix (n=D*D, row=(l,L),
    // col=(r,R)) applied to a flat length-n vector -- see this section's
    // own top comment for why no reshaping/data-movement is needed.
    static std::vector<Cplx>
    vx_matvec_row(std::vector<Cplx> const& M, int n, std::vector<Cplx> const& v)
        {
        std::vector<Cplx> y(n,Cplx(0,0));
        for (int i=0;i<n;++i)
            {
            Cplx acc(0,0);
            for (int j=0;j<n;++j) acc += M[i*n+j]*v[j];
            y[i] = acc;
            }
        return y;
        }

    // v^T @ M (v a flat length-n row vector) -- mirror of vx_matvec_row,
    // contracting M's own FIRST (row) index against v instead.
    static std::vector<Cplx>
    vx_vecmat_row(std::vector<Cplx> const& v, std::vector<Cplx> const& M, int n)
        {
        std::vector<Cplx> y(n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            Cplx acc(0,0);
            for (int i=0;i<n;++i) acc += v[i]*M[i*n+j];
            y[j] = acc;
            }
        return y;
        }

    // pyitensor idmrg._apply_transfer(E4,rho) = einsum('lLrR,rR->lL',E4,rho)
    // -- E's own row-major (l,L,r,R) flattening is exactly the (D*D,D*D)
    // matrix this contraction needs (row=(l,L), col=(r,R)), so this is a
    // plain matrix-vector product.
    static std::vector<Cplx>
    vx_apply_transfer(std::vector<Cplx> const& E, int D, std::vector<Cplx> const& rho)
        { return vx_matvec_row(E,D*D,rho); }

    // pyitensor idmrg._apply_transfer_from_left(E4,rho) = einsum('lL,lLrR->rR',rho,E4)
    // -- the mirror contraction (rho against E's row index instead of column).
    static std::vector<Cplx>
    vx_apply_transfer_from_left(std::vector<Cplx> const& E, int D, std::vector<Cplx> const& rho)
        { return vx_vecmat_row(rho,E,D*D); }

    // T (D,d_g,D) with operator M (M[i,o] convention: M applied to T's
    // physical leg, M's row=input/T's own physical index, M's col=output)
    // -- pyitensor idmrg_excitations._apply_op_ket's einsum('io,aic->aoc',M,T).
    static std::vector<Cplx>
    vx_apply_op_ket(std::vector<Cplx> const& M, std::vector<Cplx> const& T, int D, int d_g)
        {
        std::vector<Cplx> Y((size_t)D*d_g*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int o=0;o<d_g;++o)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int i=0;i<d_g;++i) acc += M[i*d_g+o]*T[(l*d_g+i)*D+r];
            Y[(l*d_g+o)*D+r] = acc;
            }
        return Y;
        }

    // T (D,d_g,D) with its right bond contracted against R (D,D), R's
    // FIRST index matching T's right bond -- pyitensor idmrg_excitations.
    // _cap_right's einsum('aoc,cb->aob',T,R).
    static std::vector<Cplx>
    vx_cap_right(std::vector<Cplx> const& T, int D, int d_g, std::vector<Cplx> const& R)
        {
        std::vector<Cplx> Y((size_t)D*d_g*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int o=0;o<d_g;++o)
        for (int r2=0;r2<D;++r2)
            {
            Cplx acc(0,0);
            for (int r=0;r<D;++r) acc += T[(l*d_g+o)*D+r]*R[r*D+r2];
            Y[(l*d_g+o)*D+r2] = acc;
            }
        return Y;
        }

    // T's left bond contracted against L (D,D), L's FIRST index matching
    // T's left bond, L's SECOND index becoming the new left index (NOT
    // symmetric with vx_cap_right -- see pyitensor idmrg_excitations.
    // _cap_left's own comment) -- einsum('ba,boc->aoc',L,T).
    static std::vector<Cplx>
    vx_cap_left(std::vector<Cplx> const& L, std::vector<Cplx> const& T, int D, int d_g)
        {
        std::vector<Cplx> Y((size_t)D*d_g*D,Cplx(0,0));
        for (int a=0;a<D;++a)
        for (int o=0;o<d_g;++o)
        for (int c=0;c<D;++c)
            {
            Cplx acc(0,0);
            for (int b=0;b<D;++b) acc += L[b*D+a]*T[(b*d_g+o)*D+c];
            Y[(a*d_g+o)*D+c] = acc;
            }
        return Y;
        }

    // E4[l,L,r,R] for an explicit ket/bra pair (D,d_g,D each, not
    // necessarily the same tensor), with an operator M optionally applied
    // to the ket's physical leg first -- pyitensor idmrg_excitations.
    // _op_transfer_matrix: einsum('lpr,LpR->lLrR', ket', conj(bra)).
    static std::vector<Cplx>
    vx_op_transfer_matrix(std::vector<Cplx> const& ket, int D, int d_g,
                           std::vector<Cplx> const& bra,
                           bool hasM, std::vector<Cplx> const& M)
        {
        std::vector<Cplx> ket2 = hasM ? vx_apply_op_ket(M,ket,D,d_g) : ket;
        std::vector<Cplx> E((size_t)D*D*D*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int Lb=0;Lb<D;++Lb)
        for (int r=0;r<D;++r)
        for (int R=0;R<D;++R)
            {
            Cplx acc(0,0);
            for (int p=0;p<d_g;++p)
                acc += ket2[(l*d_g+p)*D+r]*std::conj(bra[(Lb*d_g+p)*D+R]);
            E[((l*D+Lb)*D+r)*D+R] = acc;
            }
        return E;
        }

    // A (D,D,D,D) tensor of all zeros -- transposed(2,3)->(0,1) NOT needed
    // anywhere; only used as a convenience for "no operator" transfer
    // tensors, kept for symmetry with vx_op_transfer_matrix's own hasM=false
    // path (which never allocates a dummy M at all -- see its own callers).

    // Solve A@x=b (A: n x n row-major, b: length n) via itensor::
    // zgesv_wrapper (column-major internally -- converts on entry, LAPACK
    // solves/overwrites b in place with x).
    static std::vector<Cplx>
    vx_solve(std::vector<Cplx> const& A_row, int n, std::vector<Cplx> b)
        {
        std::vector<Cplx> Acol((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*n] = A_row[i*n+j];
        auto info = zgesv_wrapper(n,1,Acol.data(),b.data());
        if (info!=0)
            throw ITError("Chain::vumps: zgesv_wrapper failed (info="+std::to_string(info)+
                           ") -- the regularized environment linear system is singular "
                           "(a near-critical/gapless transfer spectrum?)");
        return b;
        }

    // Dense (n,n) matrix (n=D*D) representing a linear map action: (D,D)
    // flat -> (D,D) flat, built one standard basis matrix at a time --
    // same style pyitensor's own _dense_linear_map uses.
    static std::vector<Cplx>
    vx_build_linear_map(int D, std::function<std::vector<Cplx>(std::vector<Cplx> const&)> const& action)
        {
        int n = D*D;
        std::vector<Cplx> Mat((size_t)n*n,Cplx(0,0));
        std::vector<Cplx> e(n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            e[j] = Cplx(1,0);
            auto col = action(e);
            for (int i=0;i<n;++i) Mat[i*n+j] = col[i];
            e[j] = Cplx(0,0);
            }
        return Mat;
        }

    // trace(conj(A) @ X) = sum_{i,k} conj(A[i,k])*X[k,i] -- NOTE: conj(A)
    // is A's ELEMENTWISE conjugate, not its conjugate transpose (matches
    // numpy's own `A.conj()` used throughout pyitensor/vumps.py and
    // pyitensor/idmrg_excitations.py at exactly these call sites).
    static Cplx
    vx_trace_conjA_X(std::vector<Cplx> const& A, std::vector<Cplx> const& X, int D)
        {
        Cplx s(0,0);
        for (int i=0;i<D;++i)
        for (int k=0;k<D;++k)
            s += std::conj(A[i*D+k])*X[k*D+i];
        return s;
        }

    // sum_i conj(A[i])*B[i] over two equal-length flat arrays, treated as
    // plain vectors regardless of their own logical shape -- C++ analogue
    // of numpy's `np.einsum('...,...->', np.conj(A), B)` reduction used
    // throughout pyitensor/vumps.py's own static-correlator section
    // (`onsite_expectation`/`two_point_correlator`'s norm/val computations).
    static Cplx
    vx_dot_conj(std::vector<Cplx> const& A, std::vector<Cplx> const& B)
        {
        Cplx s(0,0);
        for (size_t i=0;i<A.size();++i) s += std::conj(A[i])*B[i];
        return s;
        }

    // ket/bra (D,d_g,D) closed over their LEFT bond and physical legs
    // together, leaving the (right-ket, right-bra) pair open: result[r,R] =
    // sum_{l,p} ket[l,p,r]*conj(bra[l,p,R]) -- C++ analogue of pyitensor/
    // vumps.py's own `np.einsum('lor,loR->rR', ket, np.conj(bra))` (used
    // there with `bra` already `np.conj`-applied by the caller; this helper
    // conjugates `bra` itself instead, so callers pass the UN-conjugated
    // tensor, matching every other vx_* helper's convention).
    static std::vector<Cplx>
    vx_left_close(std::vector<Cplx> const& ket, std::vector<Cplx> const& bra, int D, int d_g)
        {
        std::vector<Cplx> X((size_t)D*D,Cplx(0,0));
        for (int r=0;r<D;++r)
        for (int R=0;R<D;++R)
            {
            Cplx acc(0,0);
            for (int l=0;l<D;++l)
            for (int p=0;p<d_g;++p)
                acc += ket[(l*d_g+p)*D+r]*std::conj(bra[(l*d_g+p)*D+R]);
            X[r*D+R] = acc;
            }
        return X;
        }

    // The regularized (I - phase*T + [projector]) linear solve shared by
    // pyitensor's own _solve_left_environment/_solve_right_environment
    // (VUMPS's own GL/GR) AND idmrg_excitations._channel_resolvent
    // (excitation ansatz's own GBL(k)/GBR(k)) -- see vumps.py's and
    // idmrg_excitations.py's own docstrings at those functions for the
    // full derivation of why the projector term is needed exactly when
    // phase==1 (T's own dominant eigenvalue sits exactly at 1 there) and
    // must be omitted otherwise (adding it when not needed would corrupt
    // an otherwise-nonsingular system). from_left selects
    // apply_transfer_from_left (true) or apply_transfer (false); T here
    // is always E's own "identity/no-operator" transfer tensor (E_id for
    // VUMPS's own GL/GR, E_RL/E_LR for the excitation ansatz's GBL/GBR).
    std::vector<Cplx>
    vx_regularized_solve(int D, Cplx phase, std::vector<Cplx> const& E, bool from_left,
                          std::vector<Cplx> const& proj_out, std::vector<Cplx> const& proj_in,
                          std::vector<Cplx> const& rhs) const
        {
        bool add_projector = std::abs(phase-Cplx(1.0,0.0)) < 1e-10;
        auto action = [&](std::vector<Cplx> const& x)->std::vector<Cplx>
            {
            auto Tx = from_left ? vx_apply_transfer_from_left(E,D,x) : vx_apply_transfer(E,D,x);
            std::vector<Cplx> out((size_t)D*D);
            for (int k=0;k<D*D;++k) out[k] = x[k] - phase*Tx[k];
            if (add_projector)
                {
                Cplx s = vx_trace_conjA_X(proj_in,x,D);
                for (int k=0;k<D*D;++k) out[k] += proj_out[k]*s;
                }
            return out;
            };
        auto Mat = vx_build_linear_map(D,action);
        return vx_solve(Mat,D*D,rhs);
        }

    // Eigenvalues/right-eigenvectors of a general (non-Hermitian) complex
    // n x n matrix via itensor::zgeev_wrapper, returning the dominant
    // (largest-|eigenvalue|) pair -- C++ analogue of pyitensor idmrg.py's
    // own eig()+argsort dominant-eigenvalue selection, including the same
    // near-degeneracy guard (_DEGENERACY_RTOL=1e-9, see
    // _check_dominant_eigenvalue_nondegenerate's own extensive docstring
    // for why this threshold, not a looser one, was chosen).
    static std::pair<Cplx,std::vector<Cplx>>
    vx_dominant_eig(std::vector<Cplx> const& T_row, int n)
        {
        std::vector<Cplx> Acol((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*n] = T_row[i*n+j];
        std::vector<Cplx> evals(n), vl(1,Cplx(0,0)), vr((size_t)n*n,Cplx(0,0));
        auto info = zgeev_wrapper('N','V',n,Acol.data(),evals.data(),vl.data(),vr.data());
        if (info!=0)
            throw ITError("Chain::vumps: zgeev_wrapper failed (info="+std::to_string(info)+")");
        std::vector<int> order(n);
        for (int i=0;i<n;++i) order[i]=i;
        std::sort(order.begin(),order.end(),
                  [&](int a,int b){ return std::abs(evals[a]) > std::abs(evals[b]); });
        if (n>1 && std::abs(evals[order[1]]) > (1.0-1e-9)*std::abs(evals[order[0]]))
            throw ITError("Chain::vumps: the transfer matrix's dominant eigenvalue is "
                           "(near-)degenerate -- a single dominant fixed point is not "
                           "well-defined here (see pyitensor/idmrg.py's own "
                           "_check_dominant_eigenvalue_nondegenerate for the physical "
                           "reason this can happen, e.g. a gapless/critical chain)");
        int idx = order[0];
        std::vector<Cplx> vec(n);
        for (int i=0;i<n;++i) vec[i] = vr[i + (size_t)idx*n];
        return {evals[idx], vec};
        }

    // Dominant RIGHT fixed point rho of transfer tensor E (apply_transfer(E,rho)=eta*rho),
    // normalized to trace(rho)=1 -- C++ analogue of pyitensor idmrg.py's
    // _dominant_right_fixed_point. E's own row-major (l,L,r,R) flattening
    // already IS the (D*D,D*D) matrix vx_dominant_eig needs (see this
    // section's own top comment) -- no reshaping.
    static std::pair<std::vector<Cplx>,Cplx>
    vx_dominant_right_fixed_point(std::vector<Cplx> const& E, int D)
        {
        auto [eta,vec] = vx_dominant_eig(E,D*D);
        Cplx tr(0,0);
        for (int i=0;i<D;++i) tr += vec[i*D+i];
        if (std::abs(tr) < 1e-13)
            throw ITError("Chain::vumps: dominant right fixed point has ~zero trace -- "
                           "degenerate/ill-defined normalization");
        for (auto & v : vec) v /= tr;
        return {vec,eta};
        }

    // Dominant LEFT fixed point (the right eigenvector of E's own
    // transpose, as an (n,n) matrix) -- C++ analogue of pyitensor idmrg.py's
    // _dominant_left_fixed_point.
    static std::pair<std::vector<Cplx>,Cplx>
    vx_dominant_left_fixed_point(std::vector<Cplx> const& E, int D)
        {
        int n = D*D;
        std::vector<Cplx> Et((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Et[i*n+j] = E[j*n+i];
        auto [eta,vec] = vx_dominant_eig(Et,n);
        Cplx tr(0,0);
        for (int i=0;i<D;++i) tr += vec[i*D+i];
        if (std::abs(tr) < 1e-13)
            throw ITError("Chain::vumps: dominant left fixed point has ~zero trace -- "
                           "degenerate/ill-defined normalization");
        for (auto & v : vec) v /= tr;
        return {vec,eta};
        }

    // Mixed-transfer fixed points (r,l), normalized so trace(l@r)=1 --
    // C++ analogue of pyitensor idmrg_excitations._mixed_fixed_points,
    // including its own |eta|~=1 convergence sanity check.
    static std::pair<std::vector<Cplx>,std::vector<Cplx>>
    vx_mixed_fixed_points(std::vector<Cplx> const& E, int D)
        {
        auto [r,eta_r] = vx_dominant_right_fixed_point(E,D);
        auto [l,eta_l] = vx_dominant_left_fixed_point(E,D);
        for (Cplx eta : {eta_r,eta_l})
            if (std::abs(std::abs(eta)-1.0) > 1e-6)
                throw ITError("Chain::vumps: the mixed transfer tensor's dominant "
                               "eigenvalue magnitude is not close to 1 -- the VUMPS "
                               "ground state this excitation environment was built "
                               "from is not well converged");
        Cplx norm(0,0);
        for (int i=0;i<D;++i)
        for (int k=0;k<D;++k)
            norm += l[i*D+k]*r[k*D+i]; // trace(l@r)
        std::vector<Cplx> l2 = l;
        for (auto & v : l2) v /= norm;
        return {r,l2};
        }

    // Economy ("thin") SVD A=U*diag(s)*Vt (A: m x n row-major, k=min(m,n),
    // U: m x k, Vt: k x n, both row-major on return) via itensor::
    // zgesvd_wrapper('S',...).
    static std::pair<std::vector<Cplx>,std::vector<Cplx>>
    vx_economy_svd(std::vector<Cplx> const& A_row, int m, int n)
        {
        int k = std::min(m,n);
        std::vector<Cplx> Acol((size_t)m*n);
        for (int i=0;i<m;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*m] = A_row[i*n+j];
        std::vector<LAPACK_REAL> s(k);
        std::vector<Cplx> Ucol((size_t)m*k), Vtcol((size_t)k*n);
        char jobz = 'S';
        LAPACK_INT mm=m, nn=n, info=0;
        zgesvd_wrapper(&jobz,&mm,&nn,Acol.data(),s.data(),Ucol.data(),Vtcol.data(),&info);
        if (info!=0)
            throw ITError("Chain::vumps: zgesvd_wrapper failed (info="+std::to_string(info)+")");
        std::vector<Cplx> U_row((size_t)m*k), Vt_row((size_t)k*n);
        for (int i=0;i<m;++i) for (int j=0;j<k;++j) U_row[i*k+j] = Ucol[i+(size_t)j*m];
        for (int i=0;i<k;++i) for (int j=0;j<n;++j) Vt_row[i*n+j] = Vtcol[i+(size_t)j*k];
        return {U_row,Vt_row};
        }

    // Full SVD's U factor only (m x m, row-major) via itensor::
    // zgesvd_wrapper('A',...) -- requires m>=n (only ever called on a
    // (D*d_g,D) matrix with d_g>=1, so this always holds); used only to
    // extract a null-space isometry (vumps_null_space_left below) from
    // U's own columns beyond the first n.
    static std::vector<Cplx>
    vx_full_svd_U(std::vector<Cplx> const& A_row, int m, int n)
        {
        if (m<n) Error("Chain::vumps: vx_full_svd_U requires m>=n (internal precondition)");
        std::vector<Cplx> Acol((size_t)m*n);
        for (int i=0;i<m;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*m] = A_row[i*n+j];
        int l = std::min(m,n);
        std::vector<LAPACK_REAL> s(l);
        std::vector<Cplx> Ucol((size_t)m*m), Vtcol((size_t)l*n);
        char jobz = 'A';
        LAPACK_INT mm=m, nn=n, info=0;
        zgesvd_wrapper(&jobz,&mm,&nn,Acol.data(),s.data(),Ucol.data(),Vtcol.data(),&info);
        if (info!=0)
            throw ITError("Chain::vumps: zgesvd_wrapper (full U) failed (info="+
                           std::to_string(info)+")");
        std::vector<Cplx> U_row((size_t)m*m);
        for (int i=0;i<m;++i) for (int j=0;j<m;++j) U_row[i*m+j] = Ucol[i+(size_t)j*m];
        return U_row;
        }

    // Lowest eigenpair of a Hermitian n x n matrix (row-major) via
    // itensor::zheev_wrapper -- used in place of pyitensor's own Lanczos
    // (dmrg._lanczos_ground_state) for VUMPS's own H_AC/H_C solves: D and
    // d_g are always small enough here (n_uc<=2, reach<=1 bonds) that
    // exact dense diagonalization every outer iteration is both simpler
    // and at least as robust as a truncated Krylov solve (no separate
    // Krylov-dimension knob to tune), at the cost of O(n^3) instead of
    // O(n^2) per outer iteration -- negligible at the bond dimensions
    // this port targets.
    static std::pair<double,std::vector<Cplx>>
    vx_hermitian_ground_state(std::vector<Cplx> const& H_row, int n)
        {
        std::vector<Cplx> Acol((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*n] = H_row[i*n+j];
        std::vector<LAPACK_REAL> d(n);
        auto info = zheev_wrapper(n,Acol.data(),d.data());
        if (info!=0)
            throw ITError("Chain::vumps: zheev_wrapper failed (info="+std::to_string(info)+")");
        std::vector<Cplx> evec(n);
        for (int i=0;i<n;++i) evec[i] = Acol[i]; // column 0 (ascending eigenvalues)
        return {d[0], evec};
        }

    // All eigenvalues (ascending) of a Hermitian n x n matrix.
    static std::vector<double>
    vx_hermitian_eigvals(std::vector<Cplx> const& H_row, int n)
        {
        std::vector<Cplx> Acol((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Acol[i+j*n] = H_row[i*n+j];
        std::vector<LAPACK_REAL> d(n);
        auto info = zheev_wrapper(n,Acol.data(),d.data());
        if (info!=0)
            throw ITError("Chain::vumps: zheev_wrapper failed (info="+std::to_string(info)+")");
        return d;
        }

    // Modified Gram-Schmidt orthonormalization of A's own columns (A: m x
    // n row-major, m>=n, assumed full column rank) -- used only to build
    // an exactly isometric VUMPS initial guess. This is a DELIBERATE
    // simplification relative to pyitensor/vumps.py's own
    // _random_initial_state/_grow_initial_state (which reuse idmrg.py's
    // fixed-point-based _canonicalize_periodic): any exactly isometric
    // starting AL0/AR0 is an equally valid VUMPS starting point -- the
    // VUMPS iteration itself (not the initialization method) is what
    // drives (AL,AR,C) to a self-consistent, Hamiltonian-optimal fixed
    // point (see pyitensor/vumps.py's own module docstring, step 5, and
    // its "crude, mutually-INconsistent starting point" comment on
    // _random_initial_state itself) -- so this port does not need to
    // reproduce _canonicalize_periodic's own considerably more involved
    // two-sided fixed-point/PSD-square-root/SVD-truncation machinery
    // (whose only role, for VUMPS init, would be to ALSO produce an
    // isometric tensor -- ordinary Gram-Schmidt already guarantees that,
    // exactly, without needing any of the transfer-matrix fixed-point
    // machinery _canonicalize_periodic uses for a different purpose --
    // general-n_uc apply_mpo/imps_sum bond truncation, not relevant here
    // since VUMPS builds AL/AR from scratch at the target D directly).
    static std::vector<Cplx>
    vx_gram_schmidt_columns(std::vector<Cplx> A, int m, int n)
        {
        for (int j=0;j<n;++j)
            {
            for (int p=0;p<j;++p)
                {
                Cplx dot(0,0);
                for (int i=0;i<m;++i) dot += std::conj(A[i*n+p])*A[i*n+j];
                for (int i=0;i<m;++i) A[i*n+j] -= dot*A[i*n+p];
                }
            double nrm=0.0;
            for (int i=0;i<m;++i) nrm += std::norm(A[i*n+j]);
            nrm = std::sqrt(nrm);
            if (nrm < 1e-12)
                throw ITError("Chain::vumps: a random initial tensor's column became "
                               "~0 under Gram-Schmidt (extremely unlikely for a genuine "
                               "random start) -- retry");
            for (int i=0;i<m;++i) A[i*n+j] /= nrm;
            }
        return A;
        }

    // (D,d_g,D) with its (left,right) bond legs swapped -- pyitensor's
    // own transpose(A,(2,1,0)) trick (see _random_initial_state's own
    // comment for why this turns a left-canonicalization into a
    // right-canonical tensor once swapped back).
    static std::vector<Cplx>
    vx_transpose_lr(std::vector<Cplx> const& A, int D, int d_g)
        {
        std::vector<Cplx> Y((size_t)D*d_g*D);
        for (int l=0;l<D;++l)
        for (int p=0;p<d_g;++p)
        for (int r=0;r<D;++r)
            Y[(r*d_g+p)*D+l] = A[(l*d_g+p)*D+r];
        return Y;
        }

    // A raw (D,d_g,D) random tensor, optionally seeded by embedding a
    // smaller (D_old,d_g,D_old) already-good tensor into its top-left
    // block plus small noise everywhere (including the embedded block) --
    // pyitensor's own _random_raw_tensor.
    static std::vector<Cplx>
    vx_random_raw_tensor(int D, int d_g, std::mt19937_64& rng,
                          std::vector<Cplx> const* seed=nullptr, int D_old=0,
                          double noise=0.05)
        {
        std::normal_distribution<double> nd(0.0,1.0);
        std::vector<Cplx> A0((size_t)D*d_g*D);
        for (auto & v : A0) v = noise*Cplx(nd(rng),nd(rng));
        if (seed)
            for (int l=0;l<D_old;++l)
            for (int p=0;p<d_g;++p)
            for (int r=0;r<D_old;++r)
                A0[(l*d_g+p)*D+r] += (*seed)[(l*d_g+p)*D_old+r];
        return A0;
        }

    // Fresh (AL0,AR0,C0=I) from pure random noise -- pyitensor's own
    // _random_initial_state.
    static VumpsInit
    vumps_random_init(int D, int d_g, std::mt19937_64& rng)
        {
        auto A0 = vx_random_raw_tensor(D,d_g,rng);
        auto AL0 = vx_gram_schmidt_columns(A0,D*d_g,D);
        auto A0_rev = vx_transpose_lr(A0,D,d_g);
        auto AL0_rev = vx_gram_schmidt_columns(A0_rev,D*d_g,D);
        auto AR0 = vx_transpose_lr(AL0_rev,D,d_g);
        VumpsInit out; out.AL=AL0; out.AR=AR0; out.C=vx_eye(D);
        return out;
        }

    // (AL0,AR0,C0=I) at bond dimension D, warm-started by embedding a
    // smaller, already-converged D_old<D solution -- pyitensor's own
    // _grow_initial_state. See vumps_ground_state's own docstring for why
    // this D-ramp warm start matters (single-shot VUMPS from a purely
    // random D>1 start is documented, in pyitensor/vumps.py, to reliably
    // land at a worse-than-smaller-D energy).
    static VumpsInit
    vumps_grow_init(int D, int d_g, int D_old,
                     std::vector<Cplx> const& AL_old, std::vector<Cplx> const& AR_old,
                     std::mt19937_64& rng)
        {
        auto A0 = vx_random_raw_tensor(D,d_g,rng,&AL_old,D_old);
        auto AL0 = vx_gram_schmidt_columns(A0,D*d_g,D);
        auto AR_old_rev = vx_transpose_lr(AR_old,D_old,d_g);
        auto A0rev = vx_random_raw_tensor(D,d_g,rng,&AR_old_rev,D_old);
        auto AL0_rev = vx_gram_schmidt_columns(A0rev,D*d_g,D);
        auto AR0 = vx_transpose_lr(AL0_rev,D,d_g);
        VumpsInit out; out.AL=AL0; out.AR=AR0; out.C=vx_eye(D);
        return out;
        }

    // AC = AL@C -- einsum('lpm,mr->lpr',AL,C).
    static std::vector<Cplx>
    vumps_compose_AL_C(std::vector<Cplx> const& AL, std::vector<Cplx> const& C, int D, int d_g)
        {
        std::vector<Cplx> AC((size_t)D*d_g*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int p=0;p<d_g;++p)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int m=0;m<D;++m) acc += AL[(l*d_g+p)*D+m]*C[m*D+r];
            AC[(l*d_g+p)*D+r] = acc;
            }
        return AC;
        }

    // The d_g x d_g grouped-supersite operator matrix (M[i,o] convention,
    // same as idmrg_op_dense/vumps_onsite_matrix) obtained by placing each
    // ops_by_pos[p] (a dense d_p x d_p M[i,o] matrix) at sub-site p of the
    // n_uc-site unit cell and identity everywhere else, combined via a
    // Kronecker product in the SAME sequential (site 0 slowest-varying,
    // site n_uc-1 fastest-varying) order vumps_group_automaton itself uses
    // to build the grouped physical index (its own `pi = si0*d1+si1`) --
    // C++ analogue of pyitensor/vumps.py's own `_embed_group_operator`
    // (there built via `np.kron`, which produces exactly this row-major
    // layout for a chain of per-site (in,out) matrices). ops_by_pos may
    // name a given position at most once -- composing two operators at one
    // physical site (e.g. two_point_correlator's own r=0 same-site case)
    // is the caller's job (an ordinary d_p x d_p matrix product via
    // idmrg_matmul before calling this), not something this function does.
    std::vector<Cplx>
    vumps_embed_group_operator(std::map<int,std::vector<Cplx>> const& ops_by_pos) const
        {
        int n_uc = sites_.length();
        std::vector<Cplx> M;
        int dm = 1;
        for (int p=0;p<n_uc;++p)
            {
            int dn = dim(sites_.si(p+1));
            auto it = ops_by_pos.find(p);
            std::vector<Cplx> mat = (it != ops_by_pos.end()) ? it->second : vx_eye(dn);
            if (p==0) { M = std::move(mat); dm = dn; continue; }
            std::vector<Cplx> Mnew((size_t)dm*dn*dm*dn,Cplx(0,0));
            for (int i0=0;i0<dm;++i0)
            for (int j0=0;j0<dm;++j0)
            for (int i1=0;i1<dn;++i1)
            for (int j1=0;j1<dn;++j1)
                Mnew[(i0*dn+i1)*(dm*dn)+(j0*dn+j1)] = M[i0*dm+j0]*mat[i1*dn+j1];
            M = std::move(Mnew);
            dm = dm*dn;
            }
        return M;
        }

    // Groups n_uc<=2 per-sublattice dense automaton rows (idmrg_build_row's
    // own IdmrgAutomatonRow, S=0/F=1/pending=2.. channel convention) into
    // one single-supersite dense automaton W (Dw,Dw,d_g,d_g) -- C++
    // analogue of pyitensor/vumps.py's own _group_automaton
    // (np.einsum('Labm,mcdR->LacbdR',W,Wp) then reshape merging (a,c) and
    // (b,d)). row[p].flat's own (left,right,phys_in,phys_out) layout
    // differs from W's own (chan_l,chan_r,phys_in,phys_out) layout used
    // here -- both are internal-only conventions, chosen independently
    // for convenience; only W's own layout (defined by this function) is
    // used by every other vumps_*/vx_* helper below.
    void
    vumps_group_automaton(std::vector<IdmrgAutomatonRow> const& rows, int n_uc,
                           int& Dw, int& d_g, std::vector<Cplx>& W) const
        {
        if (n_uc==1)
            {
            Dw = rows[0].left_n; d_g = rows[0].d;
            if (rows[0].right_n != Dw)
                Error("Chain::vumps_ground_state: internal error -- n_uc=1 automaton's "
                      "own left/right channel counts differ");
            W.assign((size_t)Dw*Dw*d_g*d_g,Cplx(0,0));
            for (int l=0;l<Dw;++l)
            for (int r=0;r<Dw;++r)
            for (int si=0;si<d_g;++si)
            for (int so=0;so<d_g;++so)
                W[((l*Dw+r)*d_g+si)*d_g+so] =
                    rows[0].flat[((size_t)l*rows[0].right_n+r)*d_g*d_g + si*d_g+so];
            return;
            }
        // n_uc==2
        int Dleft = rows[0].left_n, Dmid = rows[0].right_n, Dright = rows[1].right_n;
        if (rows[1].left_n != Dmid || Dleft != Dright)
            Error("Chain::vumps_ground_state: internal error -- grouped automaton's "
                  "wraparound channel counts are inconsistent");
        Dw = Dleft;
        int d0 = rows[0].d, d1 = rows[1].d;
        d_g = d0*d1;
        W.assign((size_t)Dw*Dw*d_g*d_g,Cplx(0,0));
        for (int l=0;l<Dw;++l)
        for (int r=0;r<Dw;++r)
        for (int si0=0;si0<d0;++si0)
        for (int so0=0;so0<d0;++so0)
        for (int si1=0;si1<d1;++si1)
        for (int so1=0;so1<d1;++so1)
            {
            Cplx acc(0,0);
            for (int m=0;m<Dmid;++m)
                acc += rows[0].flat[((size_t)l*Dmid+m)*d0*d0+si0*d0+so0]
                     * rows[1].flat[((size_t)m*Dw+r)*d1*d1+si1*d1+so1];
            int pi = si0*d1+si1, po = so0*d1+so1;
            W[((l*Dw+r)*d_g+pi)*d_g+po] = acc;
            }
        }

    // Raises ITError if a grouped automaton has any nonzero transition
    // directly connecting two distinct pending channels -- the signature
    // of a bond with reach>1 supersite -- C++ analogue of pyitensor
    // idmrg_excitations._check_reach_one.
    static void
    vumps_check_reach_one(std::vector<Cplx> const& W, int Dw, int d_g)
        {
        for (int p=2;p<Dw;++p)
        for (int q=2;q<Dw;++q)
            for (int si=0;si<d_g;++si)
            for (int so=0;so<d_g;++so)
                if (std::abs(W[((p*Dw+q)*d_g+si)*d_g+so]) > 1e-12)
                    throw ITError("Chain::vumps_ground_state: a Hamiltonian term with "
                                   "reach>1 unit cell (a bond spanning more than 2 adjacent "
                                   "supersites) was detected -- VUMPS/the tangent-space "
                                   "excitation ansatz implemented here only support "
                                   "nearest-adjacent-unit-cell (reach<=1) couplings");
        }

    // W[S,:,:,F] -- the direct onsite Hamiltonian content, all-zero if
    // there is none -- C++ analogue of pyitensor idmrg_excitations._onsite_matrix.
    static std::vector<Cplx>
    vumps_onsite_matrix(std::vector<Cplx> const& W, int Dw, int d_g)
        {
        std::vector<Cplx> h1((size_t)d_g*d_g);
        for (int si=0;si<d_g;++si)
        for (int so=0;so<d_g;++so)
            h1[si*d_g+so] = W[((0*Dw+1)*d_g+si)*d_g+so];
        return h1;
        }

    // [(mat_a=W[S,:,:,p], mat_b=W[p,:,:,F])] for every pending channel p
    // -- C++ analogue of pyitensor idmrg_excitations._pending_channels.
    static std::vector<PendingChan>
    vumps_pending_channels(std::vector<Cplx> const& W, int Dw, int d_g)
        {
        std::vector<PendingChan> out;
        for (int p=2;p<Dw;++p)
            {
            PendingChan pc;
            pc.mat_a.assign((size_t)d_g*d_g,Cplx(0,0));
            pc.mat_b.assign((size_t)d_g*d_g,Cplx(0,0));
            for (int si=0;si<d_g;++si)
            for (int so=0;so<d_g;++so)
                {
                pc.mat_a[si*d_g+so] = W[((0*Dw+p)*d_g+si)*d_g+so];
                pc.mat_b[si*d_g+so] = W[((p*Dw+1)*d_g+si)*d_g+so];
                }
            out.push_back(std::move(pc));
            }
        return out;
        }

    // == apply_mpo / imps_sum private helpers (Chain::vumps_apply_mpo,
    // Chain::vumps_imps_sum, Chain::vumps_load_uniform_state, all public,
    // further below) -- C++ analogues of pyitensor/idmrg.py's own
    // _psd_sqrt_factor/_canonicalize_periodic and pyitensor/vumps.py's own
    // _complete_mixed_gauge/apply_mpo/imps_sum. _canonicalize_periodic is
    // specialized here to the trivial n_uc=1 periodic chain (no
    // per-sublattice propagation step needed) -- VUMPS's own apply_mpo/
    // imps_sum always reduce to exactly this n_uc=1 case even when the
    // ORIGINAL n_uc is 2, since VUMPS groups everything into one
    // grouped supersite up front (see pyitensor/vumps.py's own apply_mpo/
    // imps_sum module comments, and vumps_group_automaton's own doc
    // comment above, for why).

    // (D,D) transpose -- pyitensor's own A.T.
    static std::vector<Cplx>
    vx_transpose_square(std::vector<Cplx> const& A, int D)
        {
        std::vector<Cplx> B((size_t)D*D);
        for (int i=0;i<D;++i)
        for (int j=0;j<D;++j)
            B[j*D+i] = A[i*D+j];
        return B;
        }

    // Economy SVD that also returns the singular values (descending, as
    // LAPACK's zgesvd already produces them) -- vx_economy_svd's own
    // existing callers never needed `s`; the canonicalization below does,
    // to truncate.
    static void
    vx_economy_svd_full(std::vector<Cplx> const& A_row, int m, int n,
                         std::vector<Cplx>& U_row, std::vector<double>& s,
                         std::vector<Cplx>& Vt_row)
        {
        int k = std::min(m,n);
        std::vector<Cplx> Acol((size_t)m*n);
        for (int i=0;i<m;++i)
        for (int j=0;j<n;++j)
            Acol[i+(size_t)j*m] = A_row[i*n+j];
        std::vector<LAPACK_REAL> sv(k);
        std::vector<Cplx> Ucol((size_t)m*k), Vtcol((size_t)k*n);
        char jobz = 'S';
        LAPACK_INT mm=m, nn=n, info=0;
        zgesvd_wrapper(&jobz,&mm,&nn,Acol.data(),sv.data(),Ucol.data(),Vtcol.data(),&info);
        if (info!=0)
            throw ITError("Chain::vumps: zgesvd_wrapper failed (info="+std::to_string(info)+")");
        U_row.assign((size_t)m*k,Cplx(0,0));
        Vt_row.assign((size_t)k*n,Cplx(0,0));
        for (int i=0;i<m;++i) for (int j=0;j<k;++j) U_row[i*k+j] = Ucol[i+(size_t)j*m];
        for (int i=0;i<k;++i) for (int j=0;j<n;++j) Vt_row[i*n+j] = Vtcol[i+(size_t)j*k];
        s.assign(sv.begin(),sv.end());
        }

    // Mirrors pyitensor/svd.py's own _truncate: s descending, normalized
    // weights p=s^2/sum(s^2), drop the smallest first, stopping once
    // either mindim is reached or dropping the next would exceed `cutoff`
    // of discarded weight -- except maxdim (if >0) is a hard cap enforced
    // regardless of cutoff. maxdim<=0 means "no cap" (Python's maxdim=None).
    static int
    vx_truncate(std::vector<double> const& s, double cutoff, int maxdim, int mindim=1)
        {
        int n = (int)s.size();
        std::vector<double> p(n);
        double total = 0.0;
        for (int i=0;i<n;++i) { p[i] = s[i]*s[i]; total += p[i]; }
        if (total <= 0.0) return std::max(1,mindim);
        for (auto & v : p) v /= total;
        int keep = n;
        double discarded = 0.0;
        while (keep > mindim)
            {
            bool over_maxdim = maxdim>0 && keep>maxdim;
            if (!over_maxdim && discarded + p[keep-1] > cutoff) break;
            discarded += p[keep-1];
            keep -= 1;
            }
        return keep;
        }

    // Hermitian eigendecomposition, ALL eigenpairs (ascending eigenvalues;
    // eigenvectors as columns of `evecs_col`, itself kept COLUMN-major --
    // i.e. evecs_col[i+k*n] is the i-th component of the k-th eigenvector,
    // matching zheev_wrapper's own raw LAPACK output layout directly, no
    // transpose needed) -- generalizes vx_hermitian_ground_state (which
    // only keeps the lowest eigenpair) for the psd-sqrt-factor/
    // complete-mixed-gauge helpers below, both of which need the FULL
    // spectrum.
    static void
    vx_hermitian_eig_full(std::vector<Cplx> const& H_row, int n,
                           std::vector<double>& evals, std::vector<Cplx>& evecs_col)
        {
        evecs_col.assign((size_t)n*n,Cplx(0,0));
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            evecs_col[i+(size_t)j*n] = H_row[i*n+j];
        std::vector<LAPACK_REAL> d(n);
        auto info = zheev_wrapper(n,evecs_col.data(),d.data());
        if (info!=0)
            throw ITError("Chain::vumps: zheev_wrapper failed (info="+std::to_string(info)+")");
        evals.assign(d.begin(),d.end());
        }

    // X (keep x D row-major) such that rho ~= X^H X -- Hermitized square
    // root via full eigh, dropping eigenvalues at or below `rel_floor`
    // times the largest one (covers both clipping small-negative noise and
    // dropping small-but-positive/ill-conditioned directions) -- C++
    // analogue of pyitensor/idmrg.py's own _psd_sqrt_factor. `keep` (X's
    // own row count) is returned via the out-param.
    static std::vector<Cplx>
    vx_psd_sqrt_factor(std::vector<Cplx> const& rho, int D, int& keep, double rel_floor=1e-12)
        {
        auto herm = vx_hermitize(rho,D);
        std::vector<double> evals; std::vector<Cplx> evecs_col;
        vx_hermitian_eig_full(herm,D,evals,evecs_col);
        double top = evals.empty() ? 0.0 : evals.back();
        double floor = top > 0.0 ? rel_floor*top : 0.0;
        int start = 0;
        while (start<D && evals[start] <= std::max(floor,0.0)) ++start;
        keep = D - start;
        std::vector<Cplx> X((size_t)keep*D,Cplx(0,0));
        for (int k=0;k<keep;++k)
            {
            int idx = start+k;
            double sroot = std::sqrt(std::max(evals[idx],0.0));
            for (int i=0;i<D;++i)
                X[k*D+i] = sroot*std::conj(evecs_col[i+(size_t)idx*D]);
            }
        return X;
        }

    // Result of vx_canonicalize_n1 -- the n_uc=1-specialized C++ analogue
    // of pyitensor/idmrg.py's own _canonicalize_periodic (see this
    // section's own top comment for why n_uc=1 is the only case VUMPS's
    // own apply_mpo/imps_sum ever need).
    struct CanonResult { std::vector<Cplx> AL; int D; Cplx eta; };

    // B: raw (generally non-canonical) (D_in,d_g,D_in) periodic tensor
    // (already noPrime()d/physical-leg-resolved by the caller) -- the
    // trivial n_uc=1 case of pyitensor idmrg._canonicalize_periodic: factor
    // the dominant left/right transfer-matrix fixed points
    // rho_L_before=X^H X, rho_R_before=Y Y^H (vx_psd_sqrt_factor), SVD
    // M=X@Y=U S Vh, truncate (vx_truncate), and build the asymmetric gauge
    // G_left=U^H X (no S factor), G_right=Y Vh^H S^-1 (the full inverse) --
    // see pyitensor/idmrg.py's own _canonicalize_periodic docstring for the
    // full derivation of why this asymmetric split, rather than a symmetric
    // S^-1/2 one, is what reproduces a plain left-canonical tensor (no
    // separate bond-weight layer) directly.
    static CanonResult
    vx_canonicalize_n1(std::vector<Cplx> const& B, int D_in, int d_g,
                        double cutoff, int maxdim)
        {
        auto E = vx_op_transfer_matrix(B,D_in,d_g,B,false,{});
        auto [rho_R, eta_R] = vx_dominant_right_fixed_point(E,D_in);
        auto [rho_L, eta_L] = vx_dominant_left_fixed_point(E,D_in);
        if (std::abs(eta_L-eta_R) > 1e-6*std::max(1.0,std::abs(eta_R)))
            throw ITError("Chain::vumps: apply_mpo/imps_sum canonicalization -- dominant "
                           "left/right transfer eigenvalues disagree -- indicates a "
                           "near-degenerate transfer spectrum (e.g. summing two "
                           "ordinary, equally-normalized VUMPS ground states -- see "
                           "pyitensor/vumps.py's own imps_sum docstring)");
        // n_uc=1: pyitensor's own _all_left_fixed_points per-cut
        // renormalization is the identity here (scales=[1]), so its own
        // rho_L_before=(rho_L*scales[0]).T reduces to a plain transpose.
        auto rho_L_before = vx_transpose_square(rho_L,D_in);

        int kx=0, ky=0;
        auto X = vx_psd_sqrt_factor(rho_L_before,D_in,kx);       // (kx,D_in)
        auto Ydag = vx_psd_sqrt_factor(rho_R,D_in,ky);           // (ky,D_in) == Y^H
        auto Y = vx_dagger(Ydag,ky,D_in);                        // (D_in,ky)

        auto M = vx_matmul(X,kx,D_in,Y,D_in,ky);                 // (kx,ky)
        std::vector<Cplx> Usvd, Vt; std::vector<double> s;
        vx_economy_svd_full(M,kx,ky,Usvd,s,Vt);
        int r = (int)s.size();
        int keep = vx_truncate(s,cutoff,maxdim,1);

        std::vector<Cplx> Uk((size_t)kx*keep), Vtk((size_t)keep*ky);
        for (int i=0;i<kx;++i) for (int j=0;j<keep;++j) Uk[i*keep+j] = Usvd[i*r+j];
        for (int i=0;i<keep;++i) for (int j=0;j<ky;++j) Vtk[i*ky+j] = Vt[i*ky+j];

        auto Udag = vx_dagger(Uk,kx,keep);                       // (keep,kx)
        auto G_left = vx_matmul(Udag,keep,kx,X,kx,D_in);         // (keep,D_in)

        auto Vtkdag = vx_dagger(Vtk,keep,ky);                    // (ky,keep)
        auto YV = vx_matmul(Y,D_in,ky,Vtkdag,ky,keep);           // (D_in,keep)
        std::vector<Cplx> G_right((size_t)D_in*keep);
        for (int i=0;i<D_in;++i)
        for (int j=0;j<keep;++j)
            G_right[i*keep+j] = YV[i*keep+j]/s[j];

        // new_arr[a,p,d] = sum_{b,c} G_left[a,b]*B[b,p,c]*G_right[c,d]
        std::vector<Cplx> AL_new((size_t)keep*d_g*keep,Cplx(0,0));
        for (int a=0;a<keep;++a)
        for (int p=0;p<d_g;++p)
        for (int d=0;d<keep;++d)
            {
            Cplx acc(0,0);
            for (int b=0;b<D_in;++b)
                {
                Cplx inner(0,0);
                for (int c=0;c<D_in;++c) inner += B[(b*d_g+p)*D_in+c]*G_right[c*keep+d];
                acc += G_left[a*D_in+b]*inner;
                }
            AL_new[(a*d_g+p)*keep+d] = acc;
            }

        CanonResult out; out.AL = std::move(AL_new); out.D = keep; out.eta = eta_R;
        return out;
        }

    struct MixedGaugeResult { std::vector<Cplx> AR, C, AC; };

    // (AR,C,AC): completes a left-canonical grouped-supersite tensor AL
    // (D,d_g,D) to the full VUMPS mixed gauge -- C++ analogue of
    // pyitensor/vumps.py's own _complete_mixed_gauge (Vanderstraeten,
    // Haegeman, Verstraete, arXiv:1810.07006, Eq.(9)-(17), specialized to
    // an already-left-canonical input): factor AL's own dominant right
    // transfer-matrix fixed point r=C C^dagger (Hermitian PSD square
    // root), then AR:=C^-1 AL C, AC:=AL C. C's own pseudo-inverse is
    // obtained directly from the SAME eigendecomposition used to build C
    // (C is Hermitian PSD by construction: C = V diag(sqrt(evals)) V^H
    // implies pinv(C) = V diag(1/sqrt(evals) where >floor else 0) V^H),
    // rather than a separate general pinv routine.
    MixedGaugeResult
    vumps_complete_mixed_gauge(std::vector<Cplx> const& AL, int D, int d_g) const
        {
        auto E = vx_op_transfer_matrix(AL,D,d_g,AL,false,{});
        auto [r,eta] = vx_dominant_right_fixed_point(E,D);
        auto herm = vx_hermitize(r,D);
        std::vector<double> evals; std::vector<Cplx> evecs_col;
        vx_hermitian_eig_full(herm,D,evals,evecs_col);

        double top=0.0; for (double e : evals) top = std::max(top,e);
        double floor = top>0.0 ? 1e-12*top : 0.0;

        std::vector<Cplx> C((size_t)D*D,Cplx(0,0)), Cinv((size_t)D*D,Cplx(0,0));
        for (int i=0;i<D;++i)
        for (int j=0;j<D;++j)
            {
            Cplx accC(0,0), accCinv(0,0);
            for (int k=0;k<D;++k)
                {
                double ev = std::max(evals[k],0.0);
                double sroot = std::sqrt(ev);
                Cplx term = evecs_col[i+(size_t)k*D]*std::conj(evecs_col[j+(size_t)k*D]);
                accC += sroot*term;
                if (ev > floor) accCinv += (1.0/sroot)*term;
                }
            C[i*D+j] = accC;
            Cinv[i*D+j] = accCinv;
            }

        std::vector<Cplx> AR((size_t)D*d_g*D,Cplx(0,0));
        for (int a=0;a<D;++a)
        for (int p=0;p<d_g;++p)
        for (int d=0;d<D;++d)
            {
            Cplx acc(0,0);
            for (int b=0;b<D;++b)
                {
                Cplx inner(0,0);
                for (int c=0;c<D;++c) inner += AL[(b*d_g+p)*D+c]*C[c*D+d];
                acc += Cinv[a*D+b]*inner;
                }
            AR[(a*d_g+p)*D+d] = acc;
            }
        auto AC = vumps_compose_AL_C(AL,C,D,d_g);

        MixedGaugeResult out; out.AR=std::move(AR); out.C=std::move(C); out.AC=std::move(AC);
        return out;
        }

    // Contract a grouped-supersite MPO W (Dw,Dw,d_g,d_g, vumps_group_
    // automaton's own axis convention) against AL (D,d_g,D) sharing their
    // physical leg, Kronecker-merging (left_W,left_A) and (right_W,
    // right_A) into the grown tensor's own combined bond -- the n_uc=1
    // specialization of pyitensor idmrg.grow_by_mpo/_local_grow (VUMPS
    // always grows at the single-grouped-supersite level, see this
    // section's own top comment). Not yet canonicalized/truncated (see
    // vx_canonicalize_n1).
    static std::vector<Cplx>
    vx_grow_by_mpo_n1(std::vector<Cplx> const& W, int Dw, int d_g,
                       std::vector<Cplx> const& AL, int D)
        {
        int Dg = Dw*D;
        std::vector<Cplx> B((size_t)Dg*d_g*Dg,Cplx(0,0));
        for (int lw=0;lw<Dw;++lw)
        for (int la=0;la<D;++la)
        for (int so=0;so<d_g;++so)
        for (int rw=0;rw<Dw;++rw)
        for (int ra=0;ra<D;++ra)
            {
            Cplx acc(0,0);
            for (int si=0;si<d_g;++si)
                acc += W[((lw*Dw+rw)*d_g+si)*d_g+so]*AL[(la*d_g+si)*D+ra];
            int l = lw*D+la, rr = rw*D+ra;
            B[(size_t)(l*d_g+so)*Dg+rr] = acc;
            }
        return B;
        }

    // (e, source_l) -- C++ analogue of pyitensor/vumps.py's own
    // _energy_density_and_source_from_left, including the required
    // conj(r_AL) closing convention that function's own extensive
    // docstring derives and documents as load-bearing, not optional (see
    // that docstring for the full derivation and the D>1 bug it fixed).
    double
    vumps_energy_source_from_left(std::vector<Cplx> const& AL, int D, int d_g,
                                   std::vector<Cplx> const& h1,
                                   std::vector<PendingChan> const& pending,
                                   std::vector<Cplx> const& r_AL,
                                   std::vector<Cplx>& source_l) const
        {
        std::vector<Cplx> I = vx_eye(D);
        auto E_op = [&](std::vector<Cplx> const& M)
            { return vx_op_transfer_matrix(AL,D,d_g,AL,true,M); };
        source_l = vx_apply_transfer_from_left(E_op(h1),D,I);
        for (auto const& pc : pending)
            {
            auto inner = vx_apply_transfer_from_left(E_op(pc.mat_a),D,I);
            auto term = vx_apply_transfer_from_left(E_op(pc.mat_b),D,inner);
            for (int k=0;k<D*D;++k) source_l[k] += term[k];
            }
        return vx_trace_conjA_X(r_AL,source_l,D).real();
        }

    // Mirror of vumps_energy_source_from_left -- C++ analogue of
    // pyitensor/vumps.py's own _energy_density_and_source_from_right
    // (note the mat_a/mat_b order swap relative to the left version,
    // exactly matching that function's own term order).
    double
    vumps_energy_source_from_right(std::vector<Cplx> const& AR, int D, int d_g,
                                    std::vector<Cplx> const& h1,
                                    std::vector<PendingChan> const& pending,
                                    std::vector<Cplx> const& l_AR,
                                    std::vector<Cplx>& source_r) const
        {
        std::vector<Cplx> I = vx_eye(D);
        auto E_op = [&](std::vector<Cplx> const& M)
            { return vx_op_transfer_matrix(AR,D,d_g,AR,true,M); };
        source_r = vx_apply_transfer(E_op(h1),D,I);
        for (auto const& pc : pending)
            {
            auto inner = vx_apply_transfer(E_op(pc.mat_b),D,I);
            auto term = vx_apply_transfer(E_op(pc.mat_a),D,inner);
            for (int k=0;k<D*D;++k) source_r[k] += term[k];
            }
        return vx_trace_conjA_X(l_AR,source_r,D).real();
        }

    // GL -- C++ analogue of pyitensor/vumps.py's own _solve_left_environment.
    std::vector<Cplx>
    vumps_solve_left_environment(std::vector<Cplx> const& AL, int D, int d_g,
                                  std::vector<Cplx> const& r_AL,
                                  std::vector<Cplx> const& source_l, double e) const
        {
        auto E_id = vx_op_transfer_matrix(AL,D,d_g,AL,false,{});
        std::vector<Cplx> I = vx_eye(D);
        std::vector<Cplx> rhs((size_t)D*D);
        for (int k=0;k<D*D;++k) rhs[k] = source_l[k] - e*I[k];
        return vx_regularized_solve(D,Cplx(1.0,0.0),E_id,/*from_left=*/true,I,r_AL,rhs);
        }

    // GR -- mirror of vumps_solve_left_environment, C++ analogue of
    // pyitensor/vumps.py's own _solve_right_environment.
    std::vector<Cplx>
    vumps_solve_right_environment(std::vector<Cplx> const& AR, int D, int d_g,
                                   std::vector<Cplx> const& l_AR,
                                   std::vector<Cplx> const& source_r, double e) const
        {
        auto E_id = vx_op_transfer_matrix(AR,D,d_g,AR,false,{});
        std::vector<Cplx> I = vx_eye(D);
        std::vector<Cplx> rhs((size_t)D*D);
        for (int k=0;k<D*D;++k) rhs[k] = source_r[k] - e*I[k];
        return vx_regularized_solve(D,Cplx(1.0,0.0),E_id,/*from_left=*/false,I,l_AR,rhs);
        }

    // [(mat_a,mat_b,Lvec_a,Rvec_b)] -- one per pending channel, built once
    // per outer VUMPS iteration and reused across every H_AC/H_C matvec --
    // C++ analogue of pyitensor/vumps.py's own _precompute_bond_environments.
    std::vector<PendingBondEnv>
    vumps_precompute_bond_environments(std::vector<Cplx> const& AL, std::vector<Cplx> const& AR,
                                        int D, int d_g,
                                        std::vector<PendingChan> const& pending) const
        {
        std::vector<Cplx> I = vx_eye(D);
        std::vector<PendingBondEnv> out;
        out.reserve(pending.size());
        for (auto const& pc : pending)
            {
            PendingBondEnv be;
            be.mat_a = pc.mat_a; be.mat_b = pc.mat_b;
            be.Lvec_a = vx_apply_transfer_from_left(
                vx_op_transfer_matrix(AL,D,d_g,AL,true,pc.mat_a),D,I);
            be.Rvec_b = vx_apply_transfer(
                vx_op_transfer_matrix(AR,D,d_g,AR,true,pc.mat_b),D,I);
            out.push_back(std::move(be));
            }
        return out;
        }

    // H_AC[X] -- C++ analogue of pyitensor/vumps.py's own _h_ac_action.
    std::vector<Cplx>
    vumps_h_ac_action(std::vector<Cplx> const& X, int D, int d_g,
                       std::vector<Cplx> const& GL, std::vector<Cplx> const& GR,
                       std::vector<PendingBondEnv> const& bond_envs,
                       std::vector<Cplx> const& h1) const
        {
        auto Y = vx_apply_op_ket(h1,X,D,d_g);
        auto t1 = vx_cap_right(X,D,d_g,GR);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t1[i];
        auto t2 = vx_cap_left(GL,X,D,d_g);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t2[i];
        for (auto const& be : bond_envs)
            {
            auto t3 = vx_cap_right(vx_apply_op_ket(be.mat_a,X,D,d_g),D,d_g,be.Rvec_b);
            for (size_t i=0;i<Y.size();++i) Y[i]+=t3[i];
            auto t4 = vx_cap_left(be.Lvec_a,vx_apply_op_ket(be.mat_b,X,D,d_g),D,d_g);
            for (size_t i=0;i<Y.size();++i) Y[i]+=t4[i];
            }
        return Y;
        }

    // H_C[C] -- C++ analogue of pyitensor/vumps.py's own _h_c_action (C
    // treated as a (D,1,D) tensor, i.e. vx_cap_right/vx_cap_left called
    // with d_g=1 -- C's own flat (D,D) row-major layout already matches
    // that shape exactly, no reshaping needed).
    std::vector<Cplx>
    vumps_h_c_action(std::vector<Cplx> const& C, int D,
                      std::vector<Cplx> const& GL, std::vector<Cplx> const& GR,
                      std::vector<PendingBondEnv> const& bond_envs) const
        {
        auto Y = vx_cap_right(C,D,1,GR);
        auto t1 = vx_cap_left(GL,C,D,1);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t1[i];
        for (auto const& be : bond_envs)
            {
            auto t2 = vx_cap_right(vx_cap_left(be.Lvec_a,C,D,1),D,1,be.Rvec_b);
            for (size_t i=0;i<Y.size();++i) Y[i]+=t2[i];
            }
        return Y;
        }

    // Dense (D*d_g*D)x(D*d_g*D) matrix representing H_AC, built one basis
    // vector at a time.
    std::vector<Cplx>
    vumps_build_h_ac_dense(int D, int d_g, std::vector<Cplx> const& GL, std::vector<Cplx> const& GR,
                            std::vector<PendingBondEnv> const& bond_envs,
                            std::vector<Cplx> const& h1) const
        {
        int n = D*d_g*D;
        std::vector<Cplx> H((size_t)n*n,Cplx(0,0));
        std::vector<Cplx> e(n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            e[j]=Cplx(1,0);
            auto col = vumps_h_ac_action(e,D,d_g,GL,GR,bond_envs,h1);
            for (int i=0;i<n;++i) H[i*n+j]=col[i];
            e[j]=Cplx(0,0);
            }
        return H;
        }

    // Dense (D*D)x(D*D) matrix representing H_C.
    std::vector<Cplx>
    vumps_build_h_c_dense(int D, std::vector<Cplx> const& GL, std::vector<Cplx> const& GR,
                           std::vector<PendingBondEnv> const& bond_envs) const
        {
        int n = D*D;
        std::vector<Cplx> H((size_t)n*n,Cplx(0,0));
        std::vector<Cplx> e(n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            e[j]=Cplx(1,0);
            auto col = vumps_h_c_action(e,D,GL,GR,bond_envs);
            for (int i=0;i<n;++i) H[i*n+j]=col[i];
            e[j]=Cplx(0,0);
            }
        return H;
        }

    // New (AL,AR) from the just-solved (AC,C) via the orthogonal-
    // Procrustes least-squares update -- C++ analogue of pyitensor/
    // vumps.py's own _update_AL_AR. Uses the identity that, in ECONOMY
    // SVD form (A=U*diag(s)*Vt, k=min(m,n)), the nearest isometry to A
    // (in either the "left" or "right" scipy.linalg.polar sense) is
    // always U@Vt -- derived directly (see this file's own development
    // notes): scipy's polar(A,side) factors A=W*P ('right', W: m x n
    // isometry columns when m>=n) or A=P*W ('left', W: m x n isometry
    // rows when m<=n); in both cases W=U_econ@Vt_econ, only P differs
    // (V*diag(s)*Vt vs U*diag(s)*Ut) -- and _update_AL_AR only ever needs
    // W, never P, so both calls reduce to the same "economy SVD then
    // U@Vt" recipe regardless of which `side` pyitensor's own call used.
    static std::pair<std::vector<Cplx>,std::vector<Cplx>>
    vumps_update_AL_AR(std::vector<Cplx> const& AC, int D, int d_g, std::vector<Cplx> const& C)
        {
        // AC's own (l*d_g+p)*D+r row-major flat layout already IS AC
        // reshaped (D*d_g,D) row-major (row=(l,p), col=r) -- no data
        // movement needed for this first reshape.
        auto [U1,Vt1] = vx_economy_svd(AC,D*d_g,D);
        auto U_l = vx_matmul(U1,D*d_g,D,Vt1,D,D); // (D*d_g,D)
        auto [Uc,Vtc] = vx_economy_svd(C,D,D);
        auto U_c = vx_matmul(Uc,D,D,Vtc,D,D); // (D,D)
        auto U_c_dag = vx_dagger(U_c,D,D);
        auto AL_new = vx_matmul(U_l,D*d_g,D,U_c_dag,D,D); // (D*d_g,D) == (D,d_g,D) flat

        // AC's own flat layout, reinterpreted with row=l, col=(p*D+r), IS
        // exactly AC reshaped (D,d_g*D) row-major -- again no data movement.
        auto [U2,Vt2] = vx_economy_svd(AC,D,d_g*D);
        auto U_r = vx_matmul(U2,D,D,Vt2,D,d_g*D); // (D,d_g*D)
        auto AR_new = vx_matmul(U_c_dag,D,D,U_r,D,d_g*D); // (D,d_g*D) == (D,d_g,D) flat
        return {AL_new,AR_new};
        }

    // ||AC-AL@C||+||AC-C@AR||, normalized by ||AC|| -- C++ analogue of
    // pyitensor/vumps.py's own _gauge_mismatch (called with the SAME
    // (AL,AR) that were held fixed while solving for this AC/C -- see
    // that function's own docstring for why, not the freshly-refit
    // vumps_update_AL_AR(AC,C) output).
    static double
    vumps_gauge_mismatch(std::vector<Cplx> const& AC, std::vector<Cplx> const& C,
                          std::vector<Cplx> const& AL, std::vector<Cplx> const& AR,
                          int D, int d_g)
        {
        auto lhs1 = vumps_compose_AL_C(AL,C,D,d_g);
        std::vector<Cplx> lhs2((size_t)D*d_g*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int p=0;p<d_g;++p)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int m=0;m<D;++m) acc += C[l*D+m]*AR[(m*d_g+p)*D+r];
            lhs2[(l*d_g+p)*D+r] = acc;
            }
        double norm_ac=0.0;
        for (auto const& v : AC) norm_ac += std::norm(v);
        norm_ac = std::sqrt(norm_ac);
        if (norm_ac==0.0) return 0.0;
        double d1=0.0, d2=0.0;
        for (size_t i=0;i<AC.size();++i)
            {
            d1 += std::norm(AC[i]-lhs1[i]);
            d2 += std::norm(AC[i]-lhs2[i]);
            }
        return (std::sqrt(d1)+std::sqrt(d2))/norm_ac;
        }

    // GL/GR/e_cell/bond_envs from the current (AL,AR) -- one full
    // per-iteration environment build -- C++ analogue of pyitensor/
    // vumps.py's own _environments.
    VumpsEnv
    vumps_build_environments(std::vector<Cplx> const& AL, std::vector<Cplx> const& AR,
                              int D, int d_g, std::vector<Cplx> const& h1,
                              std::vector<PendingChan> const& pending) const
        {
        auto E_AL = vx_op_transfer_matrix(AL,D,d_g,AL,false,{});
        auto [r_AL_raw,eta_r] = vx_dominant_right_fixed_point(E_AL,D);
        (void)eta_r;
        auto r_AL = vx_hermitize(r_AL_raw,D);

        auto E_AR = vx_op_transfer_matrix(AR,D,d_g,AR,false,{});
        auto [l_AR_raw,eta_l] = vx_dominant_left_fixed_point(E_AR,D);
        (void)eta_l;
        auto l_AR = vx_hermitize(l_AR_raw,D);

        std::vector<Cplx> source_l, source_r;
        double e_L = vumps_energy_source_from_left(AL,D,d_g,h1,pending,r_AL,source_l);
        double e_R = vumps_energy_source_from_right(AR,D,d_g,h1,pending,l_AR,source_r);

        VumpsEnv env;
        env.GL = vumps_solve_left_environment(AL,D,d_g,r_AL,source_l,e_L);
        env.GR = vumps_solve_right_environment(AR,D,d_g,l_AR,source_r,e_R);
        env.bond_envs = vumps_precompute_bond_environments(AL,AR,D,d_g,pending);
        env.e_cell = 0.5*(e_L+e_R);
        return env;
        }

    // One full VUMPS attempt (random or warm-started start) run to
    // convergence or maxiter -- C++ analogue of pyitensor/vumps.py's own
    // _vumps_single_run.
    VumpsRunResult
    vumps_single_run(int D, int d_g, std::vector<Cplx> const& h1,
                      std::vector<PendingChan> const& pending,
                      double tol, int maxiter, VumpsInit const* init,
                      std::mt19937_64& rng) const
        {
        VumpsInit start = init ? *init : vumps_random_init(D,d_g,rng);
        auto AL = start.AL, AR = start.AR, C = start.C;
        auto AC = vumps_compose_AL_C(AL,C,D,d_g);

        bool converged=false; double mismatch=0.0; int it=0;
        VumpsEnv env;
        for (it=0; it<maxiter; ++it)
            {
            env = vumps_build_environments(AL,AR,D,d_g,h1,pending);

            auto HAC = vumps_build_h_ac_dense(D,d_g,env.GL,env.GR,env.bond_envs,h1);
            auto [eAC,AC_new] = vx_hermitian_ground_state(HAC,D*d_g*D);
            (void)eAC;
            auto HC = vumps_build_h_c_dense(D,env.GL,env.GR,env.bond_envs);
            auto [eC,C_new] = vx_hermitian_ground_state(HC,D*D);
            (void)eC;

            mismatch = vumps_gauge_mismatch(AC_new,C_new,AL,AR,D,d_g);

            AC = AC_new; C = C_new;
            auto [AL_new,AR_new] = vumps_update_AL_AR(AC,D,d_g,C);
            AL = AL_new; AR = AR_new;
            if (mismatch < tol) { converged=true; break; }
            }

        // GL/GR/e_cell above are one step stale (built from this
        // iteration's own INPUT AL/AR) -- refresh once more against the
        // FINAL AL/AR before returning, exactly like pyitensor's own
        // _vumps_single_run does at its own return point.
        env = vumps_build_environments(AL,AR,D,d_g,h1,pending);

        VumpsRunResult out;
        out.AL=AL; out.AR=AR; out.C=C; out.GL=env.GL; out.GR=env.GR;
        out.e_cell=env.e_cell; out.converged=converged;
        out.niter=std::min(it+1,maxiter); out.mismatch=mismatch;
        return out;
        }

    // Everything the excitation ansatz needs that does not depend on
    // momentum k -- built once per converged vumps_ground_state() and
    // cached (have_vumps_exc_env_) until the next one -- C++ analogue of
    // pyitensor/idmrg_excitations.py's own build_excitation_environment/
    // ExcitationEnvironment.
    void
    vumps_build_excitation_environment()
        {
        int D=vumps_D_, d_g=vumps_dg_, Dw=vumps_Dw_;
        auto const& AL=vumps_AL_; auto const& AR=vumps_AR_; auto const& C=vumps_C_;
        auto const& GL=vumps_GL_; auto const& GR=vumps_GR_; auto const& W=vumps_W_;
        auto h1 = vumps_onsite_matrix(W,Dw,d_g);
        auto pending = vumps_pending_channels(W,Dw,d_g);
        auto bond_envs = vumps_precompute_bond_environments(AL,AR,D,d_g,pending);
        std::vector<Cplx> I = vx_eye(D);

        vumps_GLfull_S_ = I; vumps_GLfull_F_ = GL;
        vumps_GRfull_F_ = I; vumps_GRfull_S_ = GR;
        vumps_GLfull_pending_.clear(); vumps_GRfull_pending_.clear();
        for (auto const& be : bond_envs)
            {
            vumps_GLfull_pending_.push_back(be.Lvec_a);
            vumps_GRfull_pending_.push_back(be.Rvec_b);
            }

        vumps_E_RL_ = vx_op_transfer_matrix(AR,D,d_g,AL,false,{}); // ket=AR,bra=AL
        vumps_E_LR_ = vx_op_transfer_matrix(AL,D,d_g,AR,false,{}); // ket=AL,bra=AR
        std::tie(vumps_r_RL_,vumps_l_RL_) = vx_mixed_fixed_points(vumps_E_RL_,D);
        std::tie(vumps_r_LR_,vumps_l_LR_) = vx_mixed_fixed_points(vumps_E_LR_,D);

        vumps_VL_ = vumps_null_space_left(AL,D,d_g);

        auto AC = vumps_compose_AL_C(AL,C,D,d_g);
        auto HAC_action = vumps_h_ac_action(AC,D,d_g,GL,GR,bond_envs,h1);
        Cplx num(0,0), den(0,0);
        for (size_t i=0;i<AC.size();++i)
            {
            num += std::conj(AC[i])*HAC_action[i];
            den += std::conj(AC[i])*AC[i];
            }
        vumps_lam_AC_ = (num/den).real();

        have_vumps_exc_env_ = true;
        }

    // V_L: (D*d_g, Dx) isometry spanning the null space of AL's own
    // (reshaped D*d_g x D) adjoint -- C++ analogue of pyitensor
    // idmrg_excitations._null_space_left, obtained here from the FULL
    // SVD's own U factor (columns D..D*d_g-1) rather than scipy's
    // dedicated null_space() -- mathematically identical (both are the
    // orthogonal complement of AL_mat's own column space, computed via
    // the same underlying SVD).
    static std::vector<Cplx>
    vumps_null_space_left(std::vector<Cplx> const& AL, int D, int d_g)
        {
        auto Ufull = vx_full_svd_U(AL,D*d_g,D); // (D*d_g,D*d_g) row-major
        int Dgd = D*d_g;
        int Dx = Dgd - D;
        std::vector<Cplx> V_L((size_t)Dgd*std::max(Dx,0));
        for (int i=0;i<Dgd;++i)
        for (int j=0;j<Dx;++j)
            V_L[i*Dx+j] = Ufull[i*Dgd + (D+j)];
        return V_L;
        }

    // The generic "left-channel-environment * ket * H[site] *
    // right-channel-environment" contraction H_eff(k) is built from three
    // times -- C++ analogue of pyitensor idmrg_excitations.
    // _full_channel_contraction.
    std::vector<Cplx>
    vumps_full_channel_contraction(std::vector<Cplx> const& left_S, std::vector<Cplx> const& left_F,
                                    std::vector<std::vector<Cplx>> const& left_pending,
                                    std::vector<Cplx> const& ket, int D, int d_g,
                                    std::vector<Cplx> const& right_S, std::vector<Cplx> const& right_F,
                                    std::vector<std::vector<Cplx>> const& right_pending,
                                    std::vector<Cplx> const& h1,
                                    std::vector<PendingChan> const& pending) const
        {
        auto Y = vx_cap_left(left_S, vx_cap_right(vx_apply_op_ket(h1,ket,D,d_g),D,d_g,right_F), D,d_g);
        auto t2 = vx_cap_left(left_F, vx_cap_right(ket,D,d_g,right_F), D,d_g);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t2[i];
        auto t3 = vx_cap_left(left_S, vx_cap_right(ket,D,d_g,right_S), D,d_g);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t3[i];
        for (size_t idx=0; idx<pending.size(); ++idx)
            {
            auto tA = vx_cap_left(left_S,
                vx_cap_right(vx_apply_op_ket(pending[idx].mat_a,ket,D,d_g),D,d_g,right_pending[idx]), D,d_g);
            for (size_t i=0;i<Y.size();++i) Y[i]+=tA[i];
            auto tB = vx_cap_left(left_pending[idx],
                vx_cap_right(vx_apply_op_ket(pending[idx].mat_b,ket,D,d_g),D,d_g,right_F), D,d_g);
            for (size_t i=0;i<Y.size();++i) Y[i]+=tB[i];
            }
        return Y;
        }

    // GBL(k) -- the channel-resolved "the excitation has already happened
    // somewhere to the left" environment (MPSKit's own lBs) -- C++
    // analogue of pyitensor idmrg_excitations._build_GBL. See that
    // function's own extensive docstring for the recursion derivation and
    // the two subtle points it documents (the momentum-phase direction,
    // and why channel F's own source needs a B inserted directly into its
    // identity self-loop, closed against GL_full[F] -- confirmed there to
    // be load-bearing for D>1 Hermiticity, not optional).
    VumpsChannelMap
    vumps_build_GBL(double k, std::vector<Cplx> const& B) const
        {
        int D=vumps_D_, d_g=vumps_dg_;
        auto const& AL=vumps_AL_; auto const& AR=vumps_AR_;
        auto h1 = vumps_onsite_matrix(vumps_W_,vumps_Dw_,d_g);
        auto pending = vumps_pending_channels(vumps_W_,vumps_Dw_,d_g);
        Cplx phase = std::exp(Cplx(0.0,1.0)*k);
        auto E_RL_id = vx_op_transfer_matrix(AR,D,d_g,AL,false,{});
        auto const& proj_out = vumps_l_RL_;
        auto const& proj_in = vumps_r_RL_;
        auto E_RL = [&](std::vector<Cplx> const& M){ return vx_op_transfer_matrix(AR,D,d_g,AL,true,M); };
        auto E_B  = [&](std::vector<Cplx> const& M){ return vx_op_transfer_matrix(B,D,d_g,AL,true,M); };
        auto E_B_none = vx_op_transfer_matrix(B,D,d_g,AL,false,{});

        auto src_S = vx_apply_transfer_from_left(E_B_none,D,vumps_GLfull_S_);
        std::vector<Cplx> src_S_scaled((size_t)D*D);
        for (int i=0;i<D*D;++i) src_S_scaled[i] = src_S[i]/phase;
        auto G_S = vx_regularized_solve(D,1.0/phase,E_RL_id,/*from_left=*/true,proj_out,proj_in,src_S_scaled);

        std::vector<std::vector<Cplx>> G_pending(pending.size());
        for (size_t idx=0; idx<pending.size(); ++idx)
            {
            auto term_bg = vx_apply_transfer_from_left(E_RL(pending[idx].mat_a),D,G_S);
            auto term_src = vx_apply_transfer_from_left(E_B(pending[idx].mat_a),D,vumps_GLfull_S_);
            std::vector<Cplx> sum((size_t)D*D);
            for (int i=0;i<D*D;++i) sum[i] = (term_bg[i]+term_src[i])/phase;
            G_pending[idx] = sum;
            }

        auto rhs_F1 = vx_apply_transfer_from_left(E_RL(h1),D,G_S);
        auto rhs_F2 = vx_apply_transfer_from_left(E_B(h1),D,vumps_GLfull_S_);
        auto rhs_F3 = vx_apply_transfer_from_left(E_B_none,D,vumps_GLfull_F_);
        std::vector<Cplx> rhs_F((size_t)D*D);
        for (int i=0;i<D*D;++i) rhs_F[i] = rhs_F1[i]+rhs_F2[i]+rhs_F3[i];
        for (size_t idx=0; idx<pending.size(); ++idx)
            {
            auto t1 = vx_apply_transfer_from_left(E_RL(pending[idx].mat_b),D,G_pending[idx]);
            auto t2 = vx_apply_transfer_from_left(E_B(pending[idx].mat_b),D,vumps_GLfull_pending_[idx]);
            for (int i=0;i<D*D;++i) rhs_F[i] += t1[i]+t2[i];
            }
        std::vector<Cplx> rhs_F_scaled((size_t)D*D);
        for (int i=0;i<D*D;++i) rhs_F_scaled[i] = rhs_F[i]/phase;
        auto G_F = vx_regularized_solve(D,1.0/phase,E_RL_id,/*from_left=*/true,proj_out,proj_in,rhs_F_scaled);

        VumpsChannelMap out; out.S=G_S; out.F=G_F; out.pending=G_pending;
        return out;
        }

    // GBR(k) -- mirror of vumps_build_GBL (MPSKit's own rBs) -- C++
    // analogue of pyitensor idmrg_excitations._build_GBR (note the
    // recursion MULTIPLIES by phase here, where GBL's own DIVIDES --
    // confirmed against MPSKit's own source not to be a typo, see that
    // function's own docstring).
    VumpsChannelMap
    vumps_build_GBR(double k, std::vector<Cplx> const& B) const
        {
        int D=vumps_D_, d_g=vumps_dg_;
        auto const& AL=vumps_AL_; auto const& AR=vumps_AR_;
        auto h1 = vumps_onsite_matrix(vumps_W_,vumps_Dw_,d_g);
        auto pending = vumps_pending_channels(vumps_W_,vumps_Dw_,d_g);
        Cplx phase = std::exp(Cplx(0.0,1.0)*k);
        auto E_LR_id = vx_op_transfer_matrix(AL,D,d_g,AR,false,{});
        auto const& proj_out = vumps_l_LR_;
        auto const& proj_in = vumps_r_LR_;
        auto E_LR = [&](std::vector<Cplx> const& M){ return vx_op_transfer_matrix(AL,D,d_g,AR,true,M); };
        auto E_B  = [&](std::vector<Cplx> const& M){ return vx_op_transfer_matrix(B,D,d_g,AR,true,M); };
        auto E_B_none = vx_op_transfer_matrix(B,D,d_g,AR,false,{});

        auto src_F = vx_apply_transfer(E_B_none,D,vumps_GRfull_F_);
        std::vector<Cplx> src_F_scaled((size_t)D*D);
        for (int i=0;i<D*D;++i) src_F_scaled[i] = src_F[i]*phase;
        auto G_F = vx_regularized_solve(D,phase,E_LR_id,/*from_left=*/false,proj_out,proj_in,src_F_scaled);

        std::vector<std::vector<Cplx>> G_pending(pending.size());
        for (size_t idx=0; idx<pending.size(); ++idx)
            {
            auto term_bg = vx_apply_transfer(E_LR(pending[idx].mat_b),D,G_F);
            auto term_src = vx_apply_transfer(E_B(pending[idx].mat_b),D,vumps_GRfull_F_);
            std::vector<Cplx> sum((size_t)D*D);
            for (int i=0;i<D*D;++i) sum[i] = (term_bg[i]+term_src[i])*phase;
            G_pending[idx] = sum;
            }

        auto rhs_S1 = vx_apply_transfer(E_LR(h1),D,G_F);
        auto rhs_S2 = vx_apply_transfer(E_B(h1),D,vumps_GRfull_F_);
        auto rhs_S3 = vx_apply_transfer(E_B_none,D,vumps_GRfull_S_);
        std::vector<Cplx> rhs_S((size_t)D*D);
        for (int i=0;i<D*D;++i) rhs_S[i] = rhs_S1[i]+rhs_S2[i]+rhs_S3[i];
        for (size_t idx=0; idx<pending.size(); ++idx)
            {
            auto t1 = vx_apply_transfer(E_LR(pending[idx].mat_a),D,G_pending[idx]);
            auto t2 = vx_apply_transfer(E_B(pending[idx].mat_a),D,vumps_GRfull_pending_[idx]);
            for (int i=0;i<D*D;++i) rhs_S[i] += t1[i]+t2[i];
            }
        std::vector<Cplx> rhs_S_scaled((size_t)D*D);
        for (int i=0;i<D*D;++i) rhs_S_scaled[i] = rhs_S[i]*phase;
        auto G_S = vx_regularized_solve(D,phase,E_LR_id,/*from_left=*/false,proj_out,proj_in,rhs_S_scaled);

        VumpsChannelMap out; out.F=G_F; out.S=G_S; out.pending=G_pending;
        return out;
        }

    // H_eff(k)[X] -- the momentum-dependent effective Hamiltonian acting
    // on a tangent-space parameter X ((Dx,D) flat, Dx=D*(d_g-1)) -- C++
    // analogue of pyitensor idmrg_excitations._h_eff_action (exactly the
    // 3 terms that function's own docstring/this file's own
    // vumps_excitation_energies doc comment describe).
    std::vector<Cplx>
    vumps_h_eff_action(double k, std::vector<Cplx> const& X) const
        {
        int D=vumps_D_, d_g=vumps_dg_;
        int Dgd = D*d_g, Dx = Dgd-D;
        std::vector<Cplx> B((size_t)Dgd*D,Cplx(0,0));
        for (int i=0;i<Dgd;++i)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int j=0;j<Dx;++j) acc += vumps_VL_[i*Dx+j]*X[j*D+r];
            B[i*D+r] = acc;
            }

        auto h1 = vumps_onsite_matrix(vumps_W_,vumps_Dw_,d_g);
        auto pending = vumps_pending_channels(vumps_W_,vumps_Dw_,d_g);

        auto Y = vumps_full_channel_contraction(
            vumps_GLfull_S_,vumps_GLfull_F_,vumps_GLfull_pending_, B,D,d_g,
            vumps_GRfull_S_,vumps_GRfull_F_,vumps_GRfull_pending_, h1,pending);

        auto GBL = vumps_build_GBL(k,B);
        auto t2 = vumps_full_channel_contraction(
            GBL.S,GBL.F,GBL.pending, AR_ref(),D,d_g,
            vumps_GRfull_S_,vumps_GRfull_F_,vumps_GRfull_pending_, h1,pending);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t2[i];

        auto GBR = vumps_build_GBR(k,B);
        auto t3 = vumps_full_channel_contraction(
            vumps_GLfull_S_,vumps_GLfull_F_,vumps_GLfull_pending_, AL_ref(),D,d_g,
            GBR.S,GBR.F,GBR.pending, h1,pending);
        for (size_t i=0;i<Y.size();++i) Y[i]+=t3[i];

        // Y is (D,d_g,D) flat, already == (D*d_g,D) row-major -- project
        // via V_L^dagger @ Y.
        std::vector<Cplx> out((size_t)Dx*D,Cplx(0,0));
        for (int j=0;j<Dx;++j)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int i=0;i<Dgd;++i) acc += std::conj(vumps_VL_[i*Dx+j])*Y[i*D+r];
            out[j*D+r] = acc;
            }
        return out;
        }

    std::vector<Cplx> const& AL_ref() const { return vumps_AL_; }
    std::vector<Cplx> const& AR_ref() const { return vumps_AR_; }

    // Dense (Dx*D)x(Dx*D) matrix representing H_eff(k).
    std::vector<Cplx>
    vumps_build_h_eff_dense(double k) const
        {
        int D=vumps_D_, d_g=vumps_dg_;
        int Dx = D*(d_g-1);
        int n = Dx*D;
        std::vector<Cplx> H((size_t)n*n,Cplx(0,0));
        std::vector<Cplx> e(n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            e[j]=Cplx(1,0);
            auto col = vumps_h_eff_action(k,e);
            for (int i=0;i<n;++i) H[i*n+j]=col[i];
            e[j]=Cplx(0,0);
            }
        return H;
        }

    // Flat, order-preserving (de)serialization of a rank-4 ITensor over
    // (i0,i1,i2,i3) -- used solely to warm-start idmrg_local_solve's
    // Arnoldi search from the *previous* macro-iteration's own converged
    // local ground vector at the same mstep position, mirroring
    // pyitensor's own idmrg.py _local_two_site_solve x0_warm mechanism
    // (see idmrg_local_solve's own comment for why this is not cosmetic).
    // Real ITensor v3 Index objects have no cross-macro-iteration identity
    // to reuse directly (HL_ket/HR_ket are freshly minted every
    // micro-step, see this method's own new_bond_u/new_bond_v), so --
    // exactly like pyitensor's own flat NumPy array, reused there purely
    // by size, not Index/object identity -- the warm-start vector is
    // threaded through as a plain std::vector<Cplx>, positionally aligned
    // to (i0,i1,i2,i3)'s own dims, and re-attached to whatever fresh
    // Index objects the *next* macro-iteration mints. eltC/. set() (not
    // raw storage access) are used deliberately, since an ITensor's
    // internal dense-array layout is not guaranteed to already match the
    // (i0,i1,i2,i3) order requested here (e.g. after the noPrime()/
    // replaceInds() calls in idmrg_local_solve's own matvec).
    static std::vector<Cplx>
    idmrg_tensor_to_flat4(ITensor const& T, Index const& i0, Index const& i1,
                           Index const& i2, Index const& i3)
        {
        int d0=dim(i0), d1=dim(i1), d2=dim(i2), d3=dim(i3);
        std::vector<Cplx> out((size_t)d0*d1*d2*d3);
        size_t idx=0;
        for (int a=1;a<=d0;++a)
        for (int b=1;b<=d1;++b)
        for (int c=1;c<=d2;++c)
        for (int e=1;e<=d3;++e)
            out[idx++] = eltC(T,i0(a),i1(b),i2(c),i3(e));
        return out;
        }

    static ITensor
    idmrg_flat4_to_tensor(std::vector<Cplx> const& flat, Index const& i0,
                           Index const& i1, Index const& i2, Index const& i3)
        {
        ITensor T(i0,i1,i2,i3);
        int d0=dim(i0), d1=dim(i1), d2=dim(i2), d3=dim(i3);
        size_t idx=0;
        for (int a=1;a<=d0;++a)
        for (int b=1;b<=d1;++b)
        for (int c=1;c<=d2;++c)
        for (int e=1;e<=d3;++e)
            {
            Cplx v = flat[idx++];
            if (v != Cplx(0,0)) T.set({i0(a),i1(b),i2(c),i3(e)},v);
            }
        return T;
        }

    // One micro-step's local ground-state solve: the effective 2-site
    // Hamiltonian sandwiched by (HL, W_pL, W_pR, HR), diagonalized via
    // arnoldi_smallest_real(..., Sel::SR) (the smallest-real-part Ritz
    // value of a Hermitian operator *is* its ground state, so the
    // existing non-Hermitian-capable Arnoldi engine already used by
    // nhdmrg_one_sweep above is reused directly rather than writing a
    // second, Hermitian-only Krylov solver) -- matches idmrg.py's own
    // _local_two_site_solve (kernels.make_matvec + dmrg._lanczos_ground_state),
    // but ITensor's own operator* already does all the index-matching
    // contraction that Python's kernels.py has to hand-roll. U and V
    // are returned *without* S (the singular values are discarded): HL/HR
    // are block operators built by a similarity transform through U (or
    // V) alone, exactly idmrg.py's own convention -- see that function's
    // own extensive comment on why absorbing sqrt(S) into both sides
    // instead is wrong (confirmed there against independent ED).
    //
    // `warm`: the previous macro-iteration's own converged local ground
    // vector at this same mstep position (idmrg_tensor_to_flat4's own
    // output, positionally aligned to (HL_ket,phys_L,phys_R,HR_ket)), or
    // empty for "no warm start available yet" (the very first-ever
    // micro-step, before bond dimension has saturated, or any macro-
    // iteration whose local Hilbert space dimension changed since the
    // stored vector was produced -- checked via a plain size comparison
    // below). Used as the Arnoldi start vector whenever both HL_ket/
    // HR_ket are present and the size matches; a fresh random vector
    // (randomITensorC) otherwise. This was, until now, ALWAYS a fresh
    // random vector regardless of `warm` -- a real, confirmed-elsewhere
    // gap relative to pyitensor/idmrg.py's own x0_warm mechanism (see that
    // function's own docstring): an always-random local solve lets
    // Arnoldi land on an arbitrary member of a (near-)degenerate local
    // ground manifold every macro-iteration (routine for gapless/
    // SU(2)-symmetric models) -- the reported *energy* still converges
    // fine (a degenerate manifold shares one eigenvalue), but the
    // converged idmrg_U_/idmrg_HL_/idmrg_HR_ snapshot this Chain stores
    // (consumed by td_dynamical_correlator_window's own IBC window
    // construction) keeps jumping between different members of that
    // manifold instead of settling into one self-consistent,
    // translationally-invariant state -- exactly the bug idmrg.py's own
    // x0_warm fix was written to avoid, not yet ported to this backend
    // despite this method's own surrounding comments describing this file
    // as "a line-for-line translation ... specifically to avoid
    // reintroducing bugs that were already found and fixed" in idmrg.py.
    std::tuple<double,ITensor,ITensor,Index,Index,std::vector<Cplx>>
    idmrg_local_solve(ITensor const& HL, ITensor const& W_pL, Index phys_L,
                       ITensor const& W_pR, Index phys_R, ITensor const& HR,
                       Index HL_bra, Index HL_ket, bool have_HL_ket,
                       Index HR_bra, Index HR_ket, bool have_HR_ket,
                       double cutoff, int maxdim, int krylovdim, int restarts,
                       std::vector<Cplx> const& warm) const
        {
        std::vector<Index> order_in;
        if (have_HL_ket) order_in.push_back(HL_ket);
        order_in.push_back(phys_L);
        order_in.push_back(phys_R);
        if (have_HR_ket) order_in.push_back(HR_ket);
        size_t dim_in = 1;
        for (auto const& ind : order_in) dim_in *= (size_t)dim(ind);
        ITensor x0;
        if (have_HL_ket && have_HR_ket && warm.size()==dim_in)
            x0 = idmrg_flat4_to_tensor(warm,HL_ket,phys_L,phys_R,HR_ket);
        else
            x0 = randomITensorC(IndexSet(order_in));

        // matvec must be an endomorphism on v's own index space (same
        // indices in and out) for arnoldi_smallest_real's own inner
        // products (dag(V.at(i))*w) to reduce to a scalar -- true for
        // nhdmrg_one_sweep's own apply_proj only because its psil/psir
        // *share* link-index identities by construction (bra vs ket
        // there is a *temporary* prime distinction within one
        // computation, undone by noPrime()); HL's/HR's own bra Index
        // (minted via sim(), not prime(), in idmrg_extend_HL/HR below)
        // is a *permanently* distinct object from their own ket Index,
        // so noPrime() alone does not reconcile the physical-leg
        // "out"->"in" relabeling with the link-leg "bra"->"ket" one --
        // do the latter explicitly.
        auto matvec = [&](ITensor v) -> ITensor
            {
            if (HL) v *= HL;
            v *= W_pL;
            v *= W_pR;
            if (HR) v *= HR;
            v.noPrime();
            if (have_HL_ket) v = replaceInds(v,{HL_bra},{HL_ket});
            if (have_HR_ket) v = replaceInds(v,{HR_bra},{HR_ket});
            return v;
            };

        // early_tol=1e-10 matches arnoldi_smallest_real's own hardcoded
        // between-restart residual tolerance exactly (see its early_tol
        // parameter comment) -- this doesn't loosen the convergence this
        // local solve accepts, it just detects that same criterion as
        // soon as it's met instead of only after paying for a full
        // krylovdim-size Krylov subspace every time.
        auto result = arnoldi_smallest_real(matvec,x0,krylovdim,restarts,
                                             Sel::SR,Cplx(0,0),1e-10);
        double energy = result.first.real();
        ITensor theta = result.second;

        std::vector<Cplx> theta_flat;
        if (have_HL_ket && have_HR_ket)
            theta_flat = idmrg_tensor_to_flat4(theta,HL_ket,phys_L,phys_R,HR_ket);

        std::vector<Index> left_inds;
        if (have_HL_ket) left_inds.push_back(HL_ket);
        left_inds.push_back(phys_L);
        auto [U,S,V] = svd(theta,IndexSet(left_inds),{"Cutoff",cutoff,"MaxDim",maxdim});
        Index new_bond_u = commonIndex(U,S);
        Index new_bond_v = commonIndex(S,V);
        return {energy,U,V,new_bond_u,new_bond_v,theta_flat};
        }

    // Absorb the newly-solved left-canonical site tensor U into HL, using
    // MPO tensor W_p for the sublattice just absorbed -- the C++ analogue
    // of idmrg.py's _extend_HL, but using ITensor's own operator*
    // (automatic index-matching contraction) instead of Python's
    // contract_many()+manual bra-Index substitution. left_ket_old is U's
    // own link toward the *existing* HL (only used if have_left_ket_old);
    // right_ket_new is U's own freshly-minted link away from HL (always
    // present). Returns (new_HL, new_HL_bra) where new_HL_bra is
    // right_ket_new's own bra-side counterpart, to be threaded into the
    // *next* call.
    std::pair<ITensor,Index>
    idmrg_extend_HL(ITensor const& HL, Index HL_bra, bool have_HL,
                     ITensor const& W_p, ITensor const& U,
                     Index left_ket_old, bool have_left_ket_old,
                     Index right_ket_new) const
        {
        // See arnoldi_smallest_real's own comment on this static.
        static const TagSet site_tag("Site");
        Index right_bra_new = sim(right_ket_new);
        ITensor bra_piece = dag(prime(U,site_tag));
        if (have_left_ket_old)
            bra_piece = replaceInds(bra_piece,{left_ket_old},{HL_bra});
        bra_piece = replaceInds(bra_piece,{right_ket_new},{right_bra_new});
        // Contraction order matters here even though the result doesn't:
        // U (and, before the replaceInds calls above, bra_piece) still
        // carries its own un-renamed copies of left_ket_old/right_ket_new,
        // which only cancel once HL itself is multiplied in (HL is the
        // only piece sharing HL_bra/HL_mpo/left_ket_old all at once). The
        // naive left-to-right "bra_piece*W_p*U" order defers that
        // cancellation to the very last step, so the (bra_piece*W_p)*U
        // contraction shares only the physical index (dim d) and produces
        // a huge near-outer-product intermediate carrying every one of
        // HL_bra/HL_ket_old/right_bra_new/right_ket_new's own maxdim-sized
        // legs at once (maxdim^4 elements) before immediately collapsing
        // back down to a small rank-3 result -- confirmed by profiling:
        // this was >90% of the C++ port's total iDMRG runtime before this
        // fix. Contracting HL in first keeps every intermediate rank<=4
        // and small throughout, matching the same, correct answer.
        ITensor new_HL;
        if (have_HL) new_HL = (HL * bra_piece) * W_p * U;
        else new_HL = bra_piece * W_p * U;
        return {new_HL,right_bra_new};
        }

    // Mirror of idmrg_extend_HL: V (the newly-solved right-canonical site
    // tensor) is absorbed into HR by prepending it on HR's left.
    std::pair<ITensor,Index>
    idmrg_extend_HR(ITensor const& HR, Index HR_bra, bool have_HR,
                     ITensor const& W_p, ITensor const& V,
                     Index right_ket_old, bool have_right_ket_old,
                     Index left_ket_new) const
        {
        // See idmrg_extend_HL's own comment on this static.
        static const TagSet site_tag("Site");
        Index left_bra_new = sim(left_ket_new);
        ITensor bra_piece = dag(prime(V,site_tag));
        if (have_right_ket_old)
            bra_piece = replaceInds(bra_piece,{right_ket_old},{HR_bra});
        bra_piece = replaceInds(bra_piece,{left_ket_new},{left_bra_new});
        // Same contraction-order fix as idmrg_extend_HL above (see its own
        // comment for the full derivation): contract HR in first so every
        // intermediate stays small instead of forming a maxdim^4-sized
        // near-outer-product that only collapses on the last multiply.
        ITensor new_HR;
        if (have_HR) new_HR = (HR * bra_piece) * W_p * V;
        else new_HR = bra_piece * W_p * V;
        return {new_HR,left_bra_new};
        }

    // -- IBC-window real-time dynamical correlator private helpers
    // (Chain::td_dynamical_correlator_window, public, above) --

    // idmrg_HL_/idmrg_HR_'s own "bra" leg is minted via sim() in
    // idmrg_extend_HL/HR (see this file's own idmrg_ground_state comment)
    // -- an independent, freshly-named Index unrelated by construction to
    // the "ket" leg. ITensor's own environment-propagation convention used
    // everywhere else in this codebase (LocalMPO::makeL/makeR's own
    // `dag(prime(psi(i)))`) instead expects a boundary tensor's own bra
    // leg to be *literally* prime() of its own ket leg -- this one-time
    // relabeling (pure Index substitution, no data movement) converts
    // between the two conventions so idmrg_HL_/idmrg_HR_ can be handed to
    // the ordinary tdvp()/LocalMPO machinery directly. See
    // td_dynamical_correlator_window's own top comment for the derivation
    // of why this is the only adaptation needed.
    static ITensor
    idmrg_relabel_bra_to_prime_ket(ITensor const& H, Index const& bra, Index const& ket)
        {
        if (!H) return H;
        return replaceInds(H,{bra},{prime(ket)});
        }

    // Fresh ITensor for named single-site operator `name` at 0-based
    // sublattice position p0based, built directly over the caller-supplied
    // physical Index `phys` (unprimed=in/ket, prime(phys)=out/bra) rather
    // than sites_.si(p0based+1) -- so it can be applied to a window
    // position's own (freshly minted, see idmrg_build_window) physical
    // Index instead of the original n_uc-site unit cell's. Reuses
    // idmrg_op_dense (below) for the dense matrix elements, same (in,out)
    // convention.
    ITensor
    idmrg_op_itensor(Index const& phys, int p0based, std::string const& name) const
        {
        int d = dim(phys);
        auto M = idmrg_op_dense(p0based,name);
        ITensor T(phys,prime(phys));
        for (int i=1;i<=d;++i)
        for (int j=1;j<=d;++j)
            {
            Cplx v = M[(size_t)(i-1)*d+(j-1)];
            if (v != Cplx(0,0)) T.set({phys(i),prime(phys)(j)},v);
            }
        return T;
        }

    // Doubled (ket (x) conj(bra)) transfer tensor for sublattice p, from
    // the converged, static idmrg_U_[p] -- C++ analogue of
    // pyitensor/idmrg.py's _transfer_matrices (one sublattice at a time).
    // Legs: (idmrg_U_left_[p], prime(idmrg_U_left_[p]), idmrg_U_right_[p],
    // prime(idmrg_U_right_[p])) -- consecutive p's own right/left legs
    // coincide by construction (idmrg_ground_state's own capture comment),
    // so composing several of these via ordinary ITensor multiply
    // contracts correctly with no further bookkeeping, exactly as
    // idmrg.py's own _compose does for plain NumPy arrays.
    ITensor
    idmrg_transfer_at(int p) const
        {
        ITensor const& U = idmrg_U_[p];
        Index l = idmrg_U_left_[p], r = idmrg_U_right_[p];
        ITensor bra = replaceInds(dag(U),{l,r},{prime(l),prime(r)});
        return U*bra; // shared (unprimed) Site leg contracts automatically
        }

    // Mirror of idmrg_transfer_at, with named operator `opname` applied to
    // the ket side first (idmrg.py's own _op_transfer).
    ITensor
    idmrg_transfer_at_with_op(int p, std::string const& opname) const
        {
        Index phys = sites_.si(p+1);
        ITensor OpT = idmrg_op_itensor(phys,p,opname);
        ITensor Uop = idmrg_U_[p]*OpT; Uop.noPrime("Site");
        Index l = idmrg_U_left_[p], r = idmrg_U_right_[p];
        ITensor bra = replaceInds(dag(idmrg_U_[p]),{l,r},{prime(l),prime(r)});
        return Uop*bra;
        }

    // Dominant right eigenvector of the full unit-cell transfer matrix
    // T_full=E_0*E_1*...*E_{n_uc-1} (viewed as a linear map on (l,prime(l))
    // "density matrices", l=idmrg_U_left_[0]), normalized to trace 1 --
    // C++ analogue of idmrg.py's own _dominant_right_fixed_point, but via
    // power iteration (repeatedly applying T_full and renormalizing)
    // rather than a dense eigensolve -- simpler to implement directly
    // against ITensor's own automatic index contraction (no flattening to
    // a dense chi^2 x chi^2 matrix needed) and robust regardless of chi.
    // Returned already relabeled onto (r,prime(r)) legs (r=
    // idmrg_U_right_[n_uc-1]) -- the natural convention every caller
    // needs, since idmrg_U_left_[0] and idmrg_U_right_[n_uc-1] are
    // generally *different* Index objects (only required to match in
    // *dimension*, not identity -- see idmrg_build_window's own comment).
    // Unlike idmrg.py's own dense-eigensolve version, this has no explicit
    // near-degeneracy check (see that function's own extensive comment on
    // why one matters) -- a genuinely (near-)degenerate dominant eigenvalue
    // will instead simply fail to converge within the iteration cap below,
    // surfacing as an explicit Error() rather than a silently-wrong single
    // arbitrary branch.
    ITensor
    idmrg_dominant_right_fixed_point() const
        {
        int n_uc = idmrg_n_uc_;
        ITensor T_full = idmrg_transfer_at(0);
        for (int p=1;p<n_uc;++p) T_full = T_full*idmrg_transfer_at(p);
        Index l = idmrg_U_left_[0], r = idmrg_U_right_[n_uc-1];
        if (dim(l) != dim(r))
            throw ITError("Chain::td_dynamical_correlator_window: the converged "
                  "unit cell's wraparound bond dimension is inconsistent "
                  "(idmrg_U_left_[0] and idmrg_U_right_[n_uc-1] differ) -- "
                  "try a different maxm/maxiter/etol combination for "
                  "gs_energy()/idmrg_ground_state()");
        Index lp = prime(l), rp = prime(r);
        ITensor rho = randomITensorC(IndexSet(l,lp));
        const int maxit = 2000; const double tol = 1e-12;
        for (int it=0; ; ++it)
            {
            if (it>=maxit)
                throw ITError("Chain::td_dynamical_correlator_window: transfer-"
                      "matrix dominant-eigenvector power iteration did not "
                      "converge in 2000 steps -- state may be gapless/"
                      "critical, poorly converged, or have a (near-)"
                      "degenerate dominant eigenvalue (see idmrg.py's own "
                      "_check_dominant_eigenvalue_nondegenerate for the "
                      "Python-side analogue of this failure mode)");
            ITensor new_rho = T_full*rho; // legs (r,rp)
            Cplx tr = eltC(new_rho*delta(r,rp));
            if (std::abs(tr)==0.0)
                throw ITError("Chain::td_dynamical_correlator_window: transfer-"
                      "matrix fixed-point iteration hit a zero-trace "
                      "density matrix");
            new_rho /= tr;
            ITensor rho_l = replaceInds(new_rho,{r,rp},{l,lp});
            double change = norm(rho_l-rho);
            rho = rho_l;
            if (change<tol) return replaceInds(rho,{l,lp},{r,rp});
            }
        }

    // rho_after[p] = the fixed-point "everything strictly after
    // sublattice p, wrapping back around" density matrix (legs
    // (idmrg_U_right_[p], prime(idmrg_U_right_[p]))), for every
    // p=0..n_uc-1 -- C++ analogue of idmrg.py's own
    // _all_right_fixed_points: one dominant-eigenvector computation
    // (p=n_uc-1) plus n_uc-1 cheap transfer-tensor applications.
    std::vector<ITensor>
    idmrg_all_right_fixed_points() const
        {
        int n_uc = idmrg_n_uc_;
        std::vector<ITensor> rho_after(n_uc);
        rho_after[n_uc-1] = idmrg_dominant_right_fixed_point();
        ITensor cur = rho_after[n_uc-1];
        for (int p=n_uc-1; p>0; --p)
            {
            ITensor step = idmrg_transfer_at(p)*cur; // legs (idmrg_U_left_[p],prime(...)) == (idmrg_U_right_[p-1],prime(...))
            Index l = idmrg_U_left_[p], lp = prime(l);
            Cplx tr = eltC(step*delta(l,lp));
            step /= tr;
            cur = step;
            rho_after[p-1] = cur;
            }
        return rho_after;
        }

    // <opname> at sublattice p (0..n_uc-1) of the converged, unperturbed
    // infinite chain -- C++ analogue of idmrg.py's own onsite_expectation,
    // used here only for td_dynamical_correlator_window's own connected-
    // background subtraction (<A><B>).
    Cplx
    idmrg_onsite_expectation(int p, std::string const& opname) const
        {
        auto rho_after = idmrg_all_right_fixed_points();
        ITensor E_op = idmrg_transfer_at_with_op(p,opname);
        ITensor val = E_op*rho_after[p]; // contracts (idmrg_U_right_[p],prime(...)), leaves (idmrg_U_left_[p],prime(...))
        Index l = idmrg_U_left_[p], lp = prime(l);
        return eltC(val*delta(l,lp));
        }

    // A 1-based window site index near the geometric middle whose own
    // sublattice position ((site-1)%n_uc) equals p_i -- C++ analogue of
    // idmrg_window.py's own _default_center.
    static int
    idmrg_window_center(int n_window, int n_uc, int p_i)
        {
        int n = n_window*n_uc;
        int mid = n/2+1;
        for (int offset=0; offset<n_uc; ++offset)
            {
            int site = mid+offset;
            if (1<=site && site<=n && (site-1)%n_uc==p_i) return site;
            }
        throw ITError("Chain::td_dynamical_correlator_window: idmrg_window_center "
              "found no matching site -- this should be unreachable");
        return -1;
        }

    struct IdmrgWindow
        {
        MPS psi;
        MPO mpo;
        int n = 0;
        std::vector<Index> phys; // phys[i-1] = window position i's own physical Index
        };

    // n_window*n_uc-site finite MPS/MPO, built by tiling idmrg_U_/
    // idmrg_rows_ n_window times and capping the two open ends with
    // idmrg_HL_ket_/idmrg_HR_ket_ (ket) and idmrg_HL_mpo_/idmrg_HR_mpo_
    // (mpo) -- C++ analogue of idmrg_window.py's own build_window/
    // _tile_periodic. Every window position gets its own genuinely fresh
    // physical Index (via idmrg_rows_[p].d, not sites_.si(p+1)) -- unlike
    // pyitensor's own _tile_periodic (which reuses sites_uc's own physical
    // Index across every copy of a given sublattice position, needing a
    // separate _refresh_physical_legs fix for n_uc==1, see that module's
    // own docstring), this sidesteps that whole bug class structurally:
    // no two window positions can ever collide on a shared physical
    // Index, regardless of n_uc.
    //
    // MPS::position()'s own initial canonicalization (called automatically
    // by TDVPWorker, see TDVP/tdvp.h, whenever `!isOrtho(psi) ||
    // psi.leftLim()!=0`) is relied on to establish a genuine canonical
    // form on first use, rather than this method asserting one itself --
    // psi.leftLim(0)/rightLim(n+1) below deliberately declares "nothing
    // orthogonal yet" (rightLim()-leftLim() != 1, so isOrtho() is false)
    // so that automatic canonicalization is guaranteed to trigger.
    IdmrgWindow
    idmrg_build_window(int n_window) const
        {
        int n_uc = idmrg_n_uc_;
        int n = n_window*n_uc;
        if (n_window>1 && dim(idmrg_U_right_[n_uc-1]) != dim(idmrg_U_left_[0]))
            throw ITError("Chain::td_dynamical_correlator_window: converged unit "
                  "cell's own wraparound bond dimension is inconsistent -- "
                  "a multi-copy window needs U's own left/right natural "
                  "bond dimensions to agree (see idmrg_window.py's own "
                  "build_window comment for the Python-side analogue) -- "
                  "try a different maxm/maxiter/etol combination for "
                  "gs_energy()/idmrg_ground_state()");

        static const TagSet link_tag("Link"), site_tag("Site");
        IdmrgWindow win;
        win.n = n;
        win.phys.resize(n);
        std::vector<ITensor> ket(n), mpoT(n);

        Index ket_left = idmrg_HL_ket_;
        Index mpo_left = idmrg_HL_mpo_;
        for (int c=0; c<n_window; ++c)
            {
            bool last_copy = (c==n_window-1);
            for (int p=0; p<n_uc; ++p)
                {
                int i = c*n_uc+p+1; // 1-based window position
                bool last_in_copy = (p==n_uc-1);
                bool at_end = last_copy && last_in_copy;

                win.phys[i-1] = Index(idmrg_rows_[p].d,site_tag);
                Index ket_right = at_end ? idmrg_HR_ket_ : sim(idmrg_U_right_[p]);
                ket[i-1] = replaceInds(idmrg_U_[p],
                    {idmrg_U_left_[p], sites_.si(p+1), idmrg_U_right_[p]},
                    {ket_left, win.phys[i-1], ket_right});

                Index mpo_right = at_end ? idmrg_HR_mpo_ : Index(idmrg_rows_[p].right_n,link_tag);
                mpoT[i-1] = idmrg_make_W(idmrg_rows_[p], mpo_left, mpo_right,
                                          win.phys[i-1], prime(win.phys[i-1]));

                ket_left = ket_right;
                mpo_left = mpo_right;
                }
            }

        MPS psi(n);
        MPO H(n);
        for (int i=1;i<=n;++i) { psi.set(i,ket[i-1]); H.set(i,mpoT[i-1]); }
        psi.leftLim(0); psi.rightLim(n+1);
        win.psi = psi;
        win.mpo = H;
        return win;
        }

    // Apply named single-site operator `opname` to window position `site`
    // (1-based), in place -- Sec. V.1 step 3 of arXiv:1804.09163
    // (A^dagger_0|psi>/B_0|psi>), C++ analogue of idmrg_window.py's own
    // apply_local_operator.
    void
    idmrg_window_apply_local_op(IdmrgWindow& win, int site, std::string const& opname) const
        {
        int n_uc = idmrg_n_uc_;
        int p = (site-1)%n_uc;
        Index phys = win.phys[site-1];
        ITensor OpT = idmrg_op_itensor(phys,p,opname);
        ITensor newT = win.psi.A(site)*OpT;
        newT.noPrime("Site");
        win.psi.set(site,newT);
        }

    // Dense (chi_l,d,chi_r) array (row-major) for ITensor T, given its own
    // (left,phys,right) Index triple -- C++ analogue of idmrg.py's own
    // _to_array_lpr, used only by idmrg_window_snapshot_correlator below
    // (see that method's own comment for why plain dense arrays, not
    // ITensor objects with their own Index bookkeeping, are the right tool
    // for this specific sub-computation).
    static std::vector<Cplx>
    idmrg_tensor_to_lpr_array(ITensor const& T, Index const& left, Index const& phys, Index const& right)
        {
        int chi_l = dim(left), d = dim(phys), chi_r = dim(right);
        std::vector<Cplx> out((size_t)chi_l*d*chi_r,Cplx(0,0));
        for (int l=1;l<=chi_l;++l)
        for (int s=1;s<=d;++s)
        for (int r=1;r<=chi_r;++r)
            out[((size_t)(l-1)*d+(s-1))*chi_r+(r-1)] = eltC(T,left(l),phys(s),right(r));
        return out;
        }

    // Dense (chi,chi) array (row-major) for a rank-2 (i1,i2) ITensor (e.g.
    // one of idmrg_all_right_fixed_points's own rho_after[p] tensors).
    static std::vector<Cplx>
    idmrg_matrix_to_array(ITensor const& M, Index const& i1, Index const& i2)
        {
        int n1 = dim(i1), n2 = dim(i2);
        std::vector<Cplx> out((size_t)n1*n2,Cplx(0,0));
        for (int a=1;a<=n1;++a)
        for (int b=1;b<=n2;++b)
            out[(size_t)(a-1)*n2+(b-1)] = eltC(M,i1(a),i2(b));
        return out;
        }

    // Sum over a chain of doubled (ket,conj(bra)) transfer steps, closed
    // on the left by a bare trace (correct for a left-canonical bra --
    // see td_dynamical_correlator_window's own top comment) and on the
    // right by rho_right (idmrg_all_right_fixed_points's own
    // rho_after[p_right], as a dense (chi,chi) array via
    // idmrg_matrix_to_array) -- C++ analogue of idmrg_window.py's own
    // _close_array_chain. bra_arrays/ket_arrays: same-length lists of flat
    // (chi_l,d,chi_r) arrays (idmrg_tensor_to_lpr_array's own layout),
    // aligned site by site; deliberately plain dense arrays rather than
    // ITensor objects, exactly as idmrg_window.py's own docstring explains
    // for the identical Python-side computation: ket_arrays come from the
    // window's own (generally non-uniform, TDVP-evolved) tensors, whose
    // *internal* bond Index identities are freely re-minted by every SVD
    // step and share no Index bookkeeping with bra_arrays' own fresh,
    // independently-built (static, periodic idmrg_U_-tiled) tensors at
    // all. Unlike an earlier version of this function, ket's own
    // (l,d,r) shape at each site is *not* assumed equal to bra's own --
    // ket's bond dimension generally *grows* under TDVP (up to maxdim) as
    // entanglement spreads from the perturbation, while bra's stays fixed
    // at the original converged unit cell's own (generally smaller)
    // natural bond dimension; only the shared physical dimension `d` is
    // required to match at each site (guaranteed by construction: both
    // sides are built over the same per-sublattice physical type).
    // Confirmed the hard way: an earlier version that assumed a shared
    // (l,r) indexed bra_arrays out of bounds using ket's own (generally
    // larger) dimension, segfaulting a few TDVP steps in once ket's bond
    // dimension had grown past bra's.
    //
    // The final right-edge contraction (`out=sum_{r,R} left_traced[r,R]*
    // rho_right[r,R]`) still requires ket's own accumulated right
    // dimension (Er) to equal rho_right's own dimension (chi_right,
    // idmrg_U_right_[p_right]'s natural dimension) -- a genuine
    // requirement of the method itself (idmrg_window.py's own
    // `_close_array_chain` has the identical constraint, via its final
    // `np.einsum('rR,rR->', left_traced, rho_after[p_right])`, which would
    // likewise raise a shape-mismatch error if this failed to hold), not
    // an artifact of this port: it holds once the window's own true
    // right-edge bond dimension (idmrg_HR_ket_, the accumulated growth
    // environment) has saturated to the *same* value as the per-unit-cell
    // natural SVD bond (idmrg_U_right_) -- typically true once maxm/niter
    // are large enough for both to hit the same maxm ceiling. If it
    // doesn't hold, this raises a clear Error() (see the check below)
    // rather than segfaulting or silently truncating.
    static Cplx
    idmrg_close_array_chain(std::vector<std::vector<Cplx>> const& bra_arrays,
                             std::vector<std::vector<Cplx>> const& ket_arrays,
                             std::vector<int> const& ket_l, std::vector<int> const& ket_r,
                             std::vector<int> const& bra_l, std::vector<int> const& bra_r,
                             std::vector<int> const& dims_d,
                             std::vector<Cplx> const& rho_right, int chi_right)
        {
        int nsite = (int)ket_arrays.size();
        // E: running (l,L,r,R) tensor (l,r: ket's own dims; L,R: bra's own),
        // flat [((l*EL+L)*Er+r)*ER+R]
        std::vector<Cplx> E; int El=1,EL=1,Er=1,ER=1;
        for (int site=0; site<nsite; ++site)
            {
            int l = ket_l[site], L = bra_l[site], d = dims_d[site], r = ket_r[site], R = bra_r[site];
            // step[l,L,r,R] = sum_s ket[l,s,r]*conj(bra[L,s,R])
            std::vector<Cplx> step((size_t)l*L*r*R,Cplx(0,0));
            auto const& K = ket_arrays[site]; auto const& B = bra_arrays[site];
            for (int li=0;li<l;++li)
            for (int Li=0;Li<L;++Li)
            for (int ri=0;ri<r;++ri)
            for (int Ri=0;Ri<R;++Ri)
                {
                Cplx acc(0,0);
                for (int s=0;s<d;++s)
                    acc += K[((size_t)li*d+s)*r+ri]*std::conj(B[((size_t)Li*d+s)*R+Ri]);
                step[((size_t)(li*L+Li)*r+ri)*R+Ri] = acc;
                }
            if (site==0) { E = step; El=l; EL=L; Er=r; ER=R; }
            else
                {
                // E[l,L,s,S] = sum_{r,R} E_old[l,L,r,R]*step[r,R,s,S]
                std::vector<Cplx> Enew((size_t)El*EL*r*R,Cplx(0,0));
                for (int li=0;li<El;++li)
                for (int Li=0;Li<EL;++Li)
                for (int si=0;si<r;++si)
                for (int Si=0;Si<R;++Si)
                    {
                    Cplx acc(0,0);
                    for (int ri=0;ri<Er;++ri)
                    for (int Ri=0;Ri<ER;++Ri)
                        acc += E[((size_t)(li*EL+Li)*Er+ri)*ER+Ri]*step[((size_t)(ri*R+Ri)*r+si)*R+Si];
                    Enew[((size_t)(li*EL+Li)*r+si)*R+Si] = acc;
                    }
                E = Enew; Er=r; ER=R;
                }
            }
        // left_traced[r,R] = sum_l E[l,l,r,R] (bare trace, l==L) -- valid
        // since El==EL by construction (site 0's own ket_l==bra_l, see
        // idmrg_window_snapshot_correlator's own comment: both equal
        // dim(idmrg_HL_ket_)==dim(idmrg_U_left_[0])).
        if (El!=EL)
            throw ITError("Chain::td_dynamical_correlator_window: window's own "
                  "left-edge bond dimension does not match the converged "
                  "unit cell's own natural left bond dimension needed to "
                  "close S(x,t) there -- should be unreachable given "
                  "idmrg_U_left_[0]==idmrg_HL_ket_ by construction");
        std::vector<Cplx> left_traced((size_t)Er*ER,Cplx(0,0));
        for (int li=0;li<El;++li)
        for (int ri=0;ri<Er;++ri)
        for (int Ri=0;Ri<ER;++Ri)
            left_traced[(size_t)ri*ER+Ri] += E[((size_t)(li*EL+li)*Er+ri)*ER+Ri];
        if (Er!=chi_right || ER!=chi_right)
            throw ITError("Chain::td_dynamical_correlator_window: window's own "
                  "right-edge bond dimension(s) do not match the converged "
                  "unit cell's own natural bond dimension needed to close "
                  "S(x,t) there -- see idmrg_close_array_chain's own "
                  "comment; try increasing maxm/niter so both saturate to "
                  "the same bond dimension, or a different maxm/maxiter/"
                  "etol combination for gs_energy()");
        Cplx out(0,0);
        for (int ri=0;ri<Er;++ri)
        for (int Ri=0;Ri<ER;++Ri)
            out += left_traced[(size_t)ri*ER+Ri]*rho_right[(size_t)ri*chi_right+Ri];
        return out;
        }

    // {x: <psi0|A_x|window>} for every x in x_values, at whatever time
    // `win` has already been evolved to -- C++ analogue of
    // idmrg_window.py's own snapshot_correlator, simplified to require
    // `center+x` stay within the window's own explicit range [1,win.n]
    // for every x (idmrg_window.py's own padding for x beyond the window
    // is not ported -- callers should increase n_window instead; this
    // Chain-level method raises a clear Error rather than silently
    // truncating x_values, see td_dynamical_correlator_window's own
    // caller-facing docstring in infinitechain.py). ket_arrays are read
    // directly off win's own (possibly TDVP-evolved) tensors; bra_arrays
    // are freshly built from the static, converged idmrg_U_/idmrg_rows_
    // (never touched by evolution), with opname_A inserted at position
    // center+x -- see idmrg_close_array_chain's own comment for why these
    // two sides deliberately don't share any ITensor Index identity.
    //
    // rho_flat/chi_right: idmrg_all_right_fixed_points()'s own
    // rho_after[p_right], already flattened by the caller (see
    // td_dynamical_correlator_window, which computes this once before its
    // own per-time-step loop) -- this depends only on the static,
    // unperturbed converged ground state, never on win's own (evolving)
    // state, so recomputing it on every one of this method's own nt calls
    // would just redo the same power iteration (capped at 2000 steps)
    // needlessly; confirmed directly via code review to be a real,
    // avoidable cost, not a correctness requirement.
    std::vector<Cplx>
    idmrg_window_snapshot_correlator(IdmrgWindow const& win, std::string const& opname_A,
                                      std::vector<int> const& x_values, int center,
                                      std::vector<Cplx> const& rho_flat, int chi_right) const
        {
        int n = win.n;
        int n_uc = idmrg_n_uc_;

        // ket_arrays: win's own tensors, positions 1..n. Left/right link
        // Index at each position is re-derived fresh via commonIndex (or
        // the window's own true edges) rather than cached, since TDVP's
        // own SVD steps are free to re-mint internal Link Index identities
        // at every call -- see idmrg_close_array_chain's own comment.
        std::vector<std::vector<Cplx>> ket_arrays(n);
        std::vector<int> ket_l(n), ket_r(n), dims_d(n);
        for (int i=1;i<=n;++i)
            {
            Index left = (i==1) ? idmrg_HL_ket_ : commonIndex(win.psi.A(i-1),win.psi.A(i));
            Index right = (i==n) ? idmrg_HR_ket_ : commonIndex(win.psi.A(i),win.psi.A(i+1));
            Index phys = win.phys[i-1];
            ket_l[i-1] = dim(left); dims_d[i-1] = dim(phys); ket_r[i-1] = dim(right);
            ket_arrays[i-1] = idmrg_tensor_to_lpr_array(win.psi.A(i),left,phys,right);
            }

        // bra_l/bra_r: the static, converged unit cell's own natural
        // per-position dimensions -- fixed for every x/t, unlike
        // ket_l/ket_r above (which change as win.psi evolves), so built
        // once here rather than inside the x loop below.
        std::vector<int> bra_l(n), bra_r(n);
        for (int i=1;i<=n;++i)
            {
            int p = (i-1)%n_uc;
            bra_l[i-1] = dim(idmrg_U_left_[p]);
            bra_r[i-1] = dim(idmrg_U_right_[p]);
            }

        std::vector<Cplx> out(x_values.size());
        for (size_t ix=0; ix<x_values.size(); ++ix)
            {
            int x = x_values[ix];
            int pos = center+x;
            if (pos<1 || pos>n)
                throw ITError("Chain::td_dynamical_correlator_window: x_values "
                      "must keep center+x within the window's own explicit "
                      "range [1,n] -- increase n_window (idmrg_window.py's "
                      "own padding for x beyond the window is not ported "
                      "here, see idmrg_window_snapshot_correlator's own "
                      "comment)");

            std::vector<std::vector<Cplx>> bra_arrays(n);
            for (int i=1;i<=n;++i)
                {
                int p = (i-1)%n_uc;
                Index l = idmrg_U_left_[p], s = sites_.si(p+1), r = idmrg_U_right_[p];
                if (i==pos)
                    {
                    ITensor OpT = idmrg_op_itensor(s,p,opname_A);
                    ITensor Uop = idmrg_U_[p]*OpT; Uop.noPrime("Site");
                    bra_arrays[i-1] = idmrg_tensor_to_lpr_array(Uop,l,s,r);
                    }
                else
                    bra_arrays[i-1] = idmrg_tensor_to_lpr_array(idmrg_U_[p],l,s,r);
                }

            out[ix] = idmrg_close_array_chain(bra_arrays,ket_arrays,ket_l,ket_r,
                                               bra_l,bra_r,dims_d,rho_flat,chi_right);
            }
        return out;
        }

    Args
    dmrg_args() const { return Args("Quiet",!verbose_,"Silent",!verbose_); }

    Sweeps
    make_sweeps(int ns, int maxdim) const
        {
        auto sweeps = Sweeps(ns);
        sweeps.maxdim() = maxdim;
        sweeps.cutoff() = cutoff_;
        sweeps.noise() = noise_;
        // noise only in the first half, mirrors get_sweeps.h
        for (int i=ns/2;i<ns;i++) sweeps.setnoise(i,0.0);
        return sweeps;
        }

    Sweeps
    make_sweeps() const { return make_sweeps(nsweeps_,maxm_); }

    double
    minimum_energy()
        {
        if (!have_bandwidth_min_)
            {
            bandwidth_emin_ = gs_energy(true);
            have_bandwidth_min_ = true;
            }
        return bandwidth_emin_;
        }

    double
    maximum_energy()
        {
        if (!have_bandwidth_max_)
            {
            // Reduced-effort DMRG on -H: this value is only ever consumed
            // as a spectral *bound* (scaled_hamiltonian()'s Chebyshev
            // window, excited_states()' overlap-penalty weight), never as
            // a physical result, so the full make_sweeps() ground-state
            // schedule (nsweeps_ sweeps at maxm_) is overkill. Safety:
            // DMRG is variational, so this always *under*estimates emax,
            // and scaled_hamiltonian()'s kpm_scale (default 0.7) already
            // maps the estimated spectrum only into
            // [-1/(2*kpm_scale),+1/(2*kpm_scale)] (~[-0.71,0.71]) of the
            // Chebyshev domain [-1,1] -- headroom that tolerates an emax
            // underestimate up to ~bandwidth/6 before the true spectrum
            // leaks outside [-1,1], orders of magnitude above what a few
            // sweeps at modest bond dimension miss by (the top edge of
            // -H is a low-entanglement, ferromagnet-like state). An
            // underestimate also *shrinks* the KPM moment count
            // (n ~ bandwidth/delta), so a looser solve can't backfire
            // into more moments; the residual risk of a too-tight bound
            // is caught loudly by kpm_moments_*'s divergence guard.
            auto psi = default_mps();
            auto sweeps = make_sweeps(std::min(nsweeps_,5),std::min(maxm_,20));
            auto negH = (-1.0)*H_;
            bandwidth_emax_ = -dmrg(psi,negH,sweeps,dmrg_args());
            have_bandwidth_max_ = true;
            }
        return bandwidth_emax_;
        }

    double
    bandwidth() { return maximum_energy()-minimum_energy(); }

    double
    energy_fluctuation(MPS psi1, MPO const& H) const
        {
        psi1.normalize();
        auto psi2 = apply_mpo(H,psi1,{"MaxDim",maxm_,"Cutoff",cutoff_});
        double de = innerC(psi1,psi2).real();
        de = innerC(psi2,psi2).real()-de*de;
        return de;
        }

    std::vector<MPS>
    gram_schmidt(std::vector<MPS> wfs) const
        {
        for (size_t i=1;i<wfs.size();i++)
            {
            for (size_t j=0;j<i;j++)
                {
                auto proj = innerC(wfs.at(j),wfs.at(i))*wfs.at(j);
                auto wf = sum_mps(wfs.at(i),(-1.0)*proj);
                wf.normalize();
                wfs.at(i) = wf;
                }
            }
        return wfs;
        }

    struct HamiltonianScale { MPO scaled_H; double emin, emax, scale; };

    HamiltonianScale
    scaled_hamiltonian(double kpm_scale)
        {
        double emin = minimum_energy();
        double emax = maximum_energy();
        double shift = -(emin+emax)/2.0;
        auto ampo = AutoMPO(sites_);
        ampo += shift,"Id",1;
        auto shift_mpo = toMPO(ampo);
        auto m = sum(H_,shift_mpo,{"MaxDim",mpomaxm_,"Cutoff",cutoff_});
        double scale = (emax-emin)*kpm_scale;
        scale = 1.0/scale;
        m = m*scale;
        HamiltonianScale out; out.scaled_H = m; out.emin = emin; out.emax = emax; out.scale = scale;
        return out;
        }

    bool
    same_mps(MPS const& vi, MPS const& vj, int maxm, double cutoff) const
        {
        auto d = sum(1.0*vi,-1.0*vj,{"MaxDim",maxm,"Cutoff",cutoff});
        double dd = sqrt(innerC(d,d).real());
        return dd<1e-10;
        }

    MPO
    sum_mpo(MPO const& A1, MPO const& A2) const
        {
        return sum(A1,A2,{"MaxDim",maxm_,"Cutoff",cutoff_});
        }

    // v3's nmultMPO() (unlike v2's) requires its two operands to have
    // genuinely distinct site indices -- since every MPO built on this
    // Chain shares the exact same physical Index objects (from sites_),
    // A2 needs an extra prime bump first, exactly as v3's own error
    // message for this case suggests ("You may have meant to call
    // nmultMPO(A,prime(B))") and as ITensor's own unittest/mpo_test.cc
    // ("nmultMPO") does. That leaves the *result*'s output leg at prime
    // level 2 (nmultMPO(A,prime(B)) composes A's level-0->1 map with B's
    // now-level-1->2 map, giving an overall level-0->2 operator) instead
    // of the standard level-0->1 single-application convention every other
    // MPO here uses -- e.g. custom_exp()/evoloperator() immediately
    // sum_mpo() this against plain (level-0->1) operators, which need
    // matching index structure. mapPrime(2,1) restores that convention;
    // link indices come back out at their standard level 0 regardless (see
    // the same unittest), so this only ever touches the site legs.
    MPO
    mult_mpo(MPO const& A1, MPO const& A2) const
        {
        auto ampo = AutoMPO(sites_);
        auto out = toMPO(ampo);
        nmultMPO(A1,prime(A2),out,{"MaxDim",mpomaxm_,"Cutoff",cutoff_});
        out.mapPrime(2,1);
        return out;
        }

    MPS
    bicstab(MPO const& A, MPS const& b, double tol, int max_it, Args const& args) const
        {
        MPS x = b;
        MPS r_old = sum(b,(-1.0)*apply_mpo(A,x,args),args);
        MPS r_ = r_old;
        MPS p = r_old;
        MPS s, Ap, As, r_new;
        Cplx alpha, beta, w;
        int k = 0;
        while (k<max_it)
            {
            Ap = apply_mpo(A,p,args);
            alpha = innerC(conjugate(r_old),r_) / innerC(conjugate(Ap),r_);
            s = sum(r_old,(-alpha)*Ap,args);
            As = apply_mpo(A,s,args);
            w = innerC(conjugate(As),s) / innerC(conjugate(As),As);
            x = sum(x,sum(alpha*p,w*s,args),args);
            r_new = sum(s,(-w)*As,args);
            double res = std::sqrt(std::abs(innerC(conjugate(r_new),r_new).real()));
            if (res<=tol) break;
            beta = (alpha/w) * innerC(conjugate(r_new),r_) / innerC(conjugate(r_old),r_);
            p = sum(r_new,beta*sum(p,(-w)*Ap,args),args);
            r_old = r_new;
            k++;
            }
        return x;
        }

    MPO
    custom_exp(MPO const& H, Cplx z) const
        {
        auto ampo = AutoMPO(sites_);
        ampo += 1.0,"Id",1;
        auto Iden = toMPO(ampo);
        auto out = sum_mpo(Iden,z*H);
        auto H2 = mult_mpo(H,H);
        out = sum_mpo(out,(0.5*z*z)*H2);
        return out;
        }

    // exp(-i*dt*H) Taylor-expanded to (nominally) 3rd order -- *verbatim*
    // port of mpscpp2/chain_session.h's evoloperator(), including the same
    // latent bug reproduced deliberately there: H3 (H^3) is computed but
    // the z^3/6 term multiplies H2 again instead of H3. See the v2 file's
    // comment for why this is preserved rather than fixed.
    MPO
    evoloperator(MPO const& H, Cplx dt) const
        {
        auto ampo = AutoMPO(sites_);
        ampo += 1.0,"Id",1;
        Cplx z = Cplx(0.0,-1.0)*dt;
        auto Iden = toMPO(ampo);
        auto out = sum_mpo(Iden,z*H);
        auto H2 = mult_mpo(H,H);
        auto H3 = mult_mpo(H,H2); // computed to match the original; unused below, see note
        (void)H3;
        out = sum_mpo(out,(0.5*z*z)*H2);
        out = sum_mpo(out,(z*z*z/6.0)*H2); // NOTE: original uses H2 here, not H3
        return out;
        }

    std::vector<std::complex<double>>
    kpm_moments_full(MPO const& m, MPS const& vi, MPS const& vj, int n,
                      int kpmmaxm, double kpmcutoff) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto v = 1.0*vi;
        auto am = 1.0*vi;
        auto a = apply_mpo(m,v,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
        auto ap = 1.0*a;
        // legitimate moments <vj|T_k|vi> are bounded by ||vi||*||vj||
        // (NOT by the zeroth moment <vj|vi>, which can be ~0 for a
        // near-orthogonal cross-correlator pair) -- see check_kpm_moment
        double bound = std::sqrt(innerC(vi,vi).real()*innerC(vj,vj).real());
        out.push_back(innerC(vj,v));
        out.push_back(innerC(vj,a));
        for (int i=0;i<n;i++)
            {
            ap = apply_mpo(m,a,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            out.push_back(innerC(vj,ap));
            check_kpm_moment(out,bound);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    std::vector<std::complex<double>>
    kpm_moments_accelerated(MPO const& m, MPS const& vi, int n,
                            int kpmmaxm, double kpmcutoff) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto a = apply_mpo(m,vi,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
        auto am = 1.0*vi;
        auto ap = 1.0*a;
        Cplx mu0 = innerC(vi,vi);
        Cplx mu1 = innerC(vi,a);
        // here vi==vj, so the moment bound ||vi||*||vj|| is just mu0
        double bound = std::abs(mu0);
        out.push_back(mu0);
        out.push_back(mu1);
        for (int i=0;i<n/2;i++)
            {
            ap = apply_mpo(m,a,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            Cplx bk = 2.0*innerC(a,a) - mu0;
            Cplx bk1 = 2.0*innerC(a,ap) - mu1;
            out.push_back(bk);
            out.push_back(bk1);
            check_kpm_moment(out,bound);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    // Chebyshev moments of a correctly scaled Hamiltonian (spectrum
    // inside [-1,1]) satisfy |<vj|T_k|vi>| <= ||vi||*||vj|| = `bound`
    // (passed by the caller; for the auto-correlator path that is just
    // the zeroth moment, for the cross path it is NOT -- <vj|vi> can be
    // ~0 for near-orthogonal pairs while the moments stay O(bound)).
    // Exponential growth beyond the bound means the scaled spectrum
    // leaked outside [-1,1] (band-edge estimate too tight for the
    // chosen kpm_scale) and every subsequent moment is garbage, so fail
    // loudly instead of returning a silently wrong correlator. The +1.0
    // keeps the threshold meaningful when both norms are tiny.
    //
    // Deliberately `throw ITError(...)` here, NOT the `Error(...)` macro
    // used everywhere else in this file for genuine internal invariant
    // violations: ITensor's own itensor::error() (itensor/util/error.h,
    // vendored -- not this file) prints the message and calls abort()
    // unconditionally, it never actually throws, despite ITError being a
    // real (std::runtime_error-derived) exception type. Confirmed
    // directly: an uncaught abort()'d process can't be told apart from a
    // real crash, and pybind11 has no chance to translate it into a
    // Python exception, so a user hitting this from Python (e.g. too
    // narrow a kpm_scale) gets a hard interpreter-killing SIGABRT instead
    // of a catchable RuntimeError -- confirmed to already affect the
    // pre-existing, untouched kpm_dynamical_correlator() path too, not
    // just kpm_dynamical_correlator_truncated() above. This is a genuine
    // recoverable-by-the-caller condition (the pyitensor backend's own
    // equivalent check, chain.py's _check_kpm_moment, already just
    // `raise`s a plain RuntimeError), unlike the many other Error(...)
    // call sites in this file that really do indicate a broken
    // precondition/internal bug -- so only this one is changed. ITError
    // derives from std::runtime_error, which pybind11 auto-translates to
    // Python's RuntimeError, matching pyitensor's own exception type and
    // message exactly (see test_kpm_energy_truncation_accuracy.py's
    // pytest.raises(RuntimeError, match="KPM moments diverging")).
    static void
    check_kpm_moment(std::vector<std::complex<double>> const& out, double bound)
        {
        if (std::abs(out.back()) > 1e3*(bound+1.0))
            throw ITError("KPM moments diverging: scaled spectrum outside [-1,1] "
                          "(band-edge estimate too tight; increasing kpm_scale "
                          "widens the safety margin)");
        }

    std::vector<std::complex<double>>
    kpm_moments(MPO const& m, MPS const& vi, MPS const& vj, int n,
                int kpmmaxm, double kpmcutoff, bool accelerate) const
        {
        if (accelerate && same_mps(vi,vj,maxm_,cutoff_))
            return kpm_moments_accelerated(m,vi,n,kpmmaxm,kpmcutoff);
        return kpm_moments_full(m,vi,vj,n,kpmmaxm,kpmcutoff);
        }

    // -- KPM energy truncation (Sec. III-B of the paper cited at
    // kpm_dynamical_correlator_truncated() above) -- everything below is
    // new, independent machinery; nothing in the KPM section above this
    // point (scaled_hamiltonian/kpm_moments_full/kpm_moments_accelerated/
    // check_kpm_moment/kpm_moments) is read or modified by it. Public
    // (scaled_hamiltonian_gs_anchored()/kpm_energy_truncate() are exposed
    // to Python -- see bindings.cc -- so this method's own algorithm can
    // be tested directly, independent of the Chebyshev recursion it
    // plugs into, mirroring pyitensor/kpm_energy_truncation.py's own
    // standalone-testable design).
    public:

    // Ground-state-anchored rescaling (paper's own Eq. 21b): places the
    // ground state energy E0 at -(1-safety/2) and E0+Ws at +(1-safety/2),
    // where Ws=(emax-emin)*kpm_scale is the *same* window-width formula
    // scaled_hamiltonian() uses (so kpm_scale keeps the same meaning/
    // units in both) -- but anchored at E0 rather than centered on the
    // full many-body bandwidth's midpoint. Used by
    // kpm_dynamical_correlator_truncated() instead of scaled_hamiltonian():
    // shrinking kpm_scale only delivers a resolution gain if the narrowed
    // window also sits where a typical correlator's spectral weight
    // actually lives -- just above the ground state, since its Chebyshev
    // vectors are built by acting operators on |0>, not around the
    // geometric middle of the entire many-body spectrum (confirmed
    // directly while building pyitensor's own version of this: with
    // scaled_hamiltonian()'s bandwidth-midpoint centering, any kpm_scale
    // below the safe floor clips the ground state itself out of the
    // window before ever reaching a genuinely useful, narrower regime,
    // since the ground state sits at the spectrum's own edge, not its
    // middle).
    //
    // Returns (m,emin,emax,scale) with the *same* meaning downstream
    // (kpm_dynamical_correlator_truncated()'s own moment-count, and the
    // Python energy-axis reconstruction it feeds, kpmdmrg.py) as
    // scaled_hamiltonian()'s own return value, even though emin/emax are
    // no longer literally the true many-body band edges: emin=E0 (still
    // correct -- the energy axis kpmdmrg.py reconstructs is "excitation
    // energy above the ground state", i.e. E_phys-E0, unchanged), and
    // emax=E0+Ws (the actual retained window's top, which is what should
    // set the polynomial order N, not the full, unused many-body
    // bandwidth).
    HamiltonianScale
    scaled_hamiltonian_gs_anchored(double kpm_scale, double safety=0.025)
        {
        double e0 = minimum_energy();
        double w = bandwidth();
        double ws = w*kpm_scale;
        double wp = 1.0-safety/2.0;
        double a = ws/(2.0*wp);
        double scale = 1.0/a;
        double shift = -e0-wp*a;
        auto ampo = AutoMPO(sites_);
        ampo += shift,"Id",1;
        auto shift_mpo = toMPO(ampo);
        auto m = sum(H_,shift_mpo,{"MaxDim",mpomaxm_,"Cutoff",cutoff_});
        m = m*scale;
        HamiltonianScale out; out.scaled_H = m; out.emin = e0; out.emax = e0+ws; out.scale = scale;
        return out;
        }

    // One site's worth of Eqs. (36)-(39): build an orthonormal Krylov
    // basis of dimension <= dK for the local effective Hamiltonian `PH`
    // (already positioned at this site by the caller), starting from this
    // site's current tensor `phi0`, diagonalize the dense projected
    // matrix, and return phi0 with any component whose Krylov-subspace
    // energy has |eigenvalue| >= threshold projected out, plus the
    // squared norm of the removed part (one site's contribution to
    // Eq. 40). Two deliberate departures from the paper's literal
    // Eqs. (36)-(41), matching pyitensor/kpm_energy_truncation.py's own
    // (see that module's docstring for the full rationale):
    //  - Krylov vectors are built via interleaved (modified) Gram-Schmidt
    //    -- same style as arnoldi_smallest_real() above -- rather than
    //    the paper's "compute all dK powers of H' first, then batch
    //    orthogonalize"; both build the same subspace mathematically, but
    //    interleaving is far better conditioned once vectors start
    //    converging toward H's locally dominant eigendirection.
    //  - Unlike arnoldi_smallest_real()'s Hessenberg-only bookkeeping
    //    (sufficient there since only one Ritz pair is ever extracted),
    //    every Krylov component below threshold must survive here, so
    //    the *full* dense Hermitian projected matrix is built -- but
    //    reusing the entries already computed while orthogonalizing
    //    below (see the comment on Hk_col), rather than a second,
    //    separate PH.product() pass over the whole basis: PH.product()
    //    (a full local-effective-Hamiltonian contraction) dominates this
    //    method's cost, so avoiding a redundant pass roughly halves it.
    //  - The projector keeps |eps_alpha| < threshold (both signs), not
    //    just eps_alpha < threshold as in Eq. (38): scaled_hamiltonian_
    //    gs_anchored() pins the ground state near -1 by construction, so
    //    in practice only the upper cut ever fires for a correlator's own
    //    Chebyshev vectors, but the symmetric form is used regardless
    //    since nothing here assumes which rescaling produced `PH`.
    std::pair<ITensor,double>
    kpm_local_krylov_projection(LocalMPO& PH, ITensor const& phi0, int dK, double threshold) const
        {
        double nrm = norm(phi0);
        if (nrm<1E-14) return {phi0,0.0};
        std::vector<ITensor> V;
        V.push_back(phi0/nrm);
        // Hk_col[j][i] = Hk(i,j) = <V_i|H|V_j> for i<=j: while extending
        // from V_j, w starts out as H*V_j, and each orthogonalization
        // pass below subtracts eltC(dag(Vi)*w)*Vi from it -- the *sum* of
        // those coefficients across both passes already equals
        // <V_i|H*V_j> to working precision (w's own residual component
        // along V_i is at machine-epsilon level after two passes), so
        // capturing them here needs no extra matvec. This only ever
        // yields entries for i<=j (V_j's own extension only ever
        // orthogonalizes against vectors already present, i.e. i<=j);
        // the missing i>j half is filled by Hermitian symmetry below,
        // and the one column this loop can never produce -- j=k-1, since
        // the loop stops one short of extending past the last accepted
        // vector -- gets exactly one extra PH.product() call instead of
        // the k calls the original full second pass used.
        std::vector<std::vector<Cplx>> Hk_col;
        for (int j=0;j<dK-1;++j)
            {
            ITensor w; PH.product(V.at(j),w);
            std::vector<Cplx> col(V.size(),0.0);
            for (int pass=0;pass<2;++pass)
                for (size_t i=0;i<V.size();++i)
                    {
                    Cplx c = eltC(dag(V.at(i))*w);
                    col[i] += c;
                    w -= c*V.at(i);
                    }
            Hk_col.push_back(col); // column j, entries i=0..j
            double nw = norm(w);
            if (nw<1E-12) break; // invariant subspace found; fewer than dK vectors is fine
            V.push_back(w/nw);
            }
        int k = (int)V.size();
        auto a = Index(k,TagSet("KPMEnergyTrunc,a"));
        auto Hk = ITensor(prime(a),a);
        for (int j=0;j<(int)Hk_col.size();++j)
            for (int i=0;i<(int)Hk_col[j].size();++i)
                {
                Hk.set(prime(a)(i+1),a(j+1),Hk_col[j][i]);
                if (i!=j) Hk.set(prime(a)(j+1),a(i+1),std::conj(Hk_col[j][i]));
                }
        if ((int)Hk_col.size()<k) // ran the full dK-1 extensions: last vector's own column was never formed above
            {
            ITensor w; PH.product(V.at(k-1),w);
            for (int i=0;i<k;++i)
                Hk.set(prime(a)(i+1),a(k),eltC(dag(V.at(i))*w));
            }
        Hk = 0.5*(Hk + dag(swapPrime(Hk,0,1))); // Hermitize away floating-point asymmetry
        ITensor U,D;
        diagHermitian(Hk,U,D);
        auto eig = uniqueIndex(U,Hk);
        ITensor cv(a);
        cv.set(a(1),nrm); // phi0 == nrm*V[0] == nrm*a(1) by construction
        auto ce = dag(U)*cv;
        double truncated_weight_sq = 0.0;
        for (int alpha=1;alpha<=k;++alpha)
            {
            double e = eltC(D,eig(alpha),prime(eig)(alpha)).real();
            if (std::abs(e)>=threshold)
                {
                truncated_weight_sq += std::norm(eltC(ce,eig(alpha)));
                ce.set(eig(alpha),0.0);
                }
            }
        auto new_cv = U*ce;
        ITensor new_phi;
        for (int i=0;i<k;++i)
            {
            auto coef = eltC(new_cv,a(i+1));
            if (i==0) new_phi = coef*V.at(0);
            else new_phi += coef*V.at(i);
            }
        return {new_phi,truncated_weight_sq};
        }

    // One directional pass over every site of psi (mutated in place),
    // projecting each site's local tensor via kpm_local_krylov_projection
    // and moving the orthogonality center between sites via plain
    // MPS::position() -- whose default (no MaxDim/Cutoff args, i.e. the
    // orthMPS/QR-based branch, not the DoSVDBond one) is already a
    // lossless gauge transformation, not an SVD-with-truncation one -- so
    // this sweep's own site-to-site moves never counteract the energy
    // truncation it just did. `PH` is shared across every half-sweep of
    // the whole kpm_energy_truncate() call (constructed once there),
    // exactly as ITensor's own two-site dmrg() reuses a single LocalMPO
    // across every half-sweep of a multi-sweep run: position() only
    // needs to *extend* its cached environment by the one site the bond
    // moved past, in whichever direction, not rebuild it from scratch.
    // Returns this pass's average truncated weight per site, Eq. (40).
    double
    kpm_truncate_half_sweep(MPS& psi, LocalMPO& PH, int dK, double threshold, bool forward) const
        {
        int N = psi.length();
        double total_truncated_sq = 0.0;
        int first = forward ? 1 : N;
        int last = forward ? N : 1;
        int step = forward ? 1 : -1;
        for (int b=first; ; b+=step)
            {
            PH.position(b,psi);
            auto res = kpm_local_krylov_projection(PH,psi.A(b),dK,threshold);
            total_truncated_sq += res.second;
            psi.set(b,res.first);
            if (b==last) break;
            psi.position(b+step);
            }
        return std::sqrt(total_truncated_sq/N);
        }

    // Project high-rescaled-energy components out of a Chebyshev vector
    // `psi` (a local copy; the projected result is returned, `psi` itself
    // is untouched at the call site). `H` must be the same rescaled/
    // shifted Hamiltonian MPO used to build psi's own Chebyshev recursion
    // (H' in the paper's notation) -- NOT the bare, unscaled Hamiltonian.
    // Runs n_sweeps directional passes over the chain, alternating
    // direction each pass (matching every other sweep-based algorithm in
    // this file, e.g. ITensor's own two-site dmrg()) so the orthogonality
    // center never needs an extra, wasted end-to-end reset between
    // passes. This is a line-for-line port of pyitensor/kpm_energy_
    // truncation.py's energy_truncate() algorithm, using ITensor's own
    // LocalMPO/diagHermitian instead of that module's hand-rolled
    // NumPy-array tensor primitives; see its docstring for the returned
    // diagnostics' meaning (Eqs. 40-41).
    std::pair<MPS,KPMTruncStats>
    kpm_energy_truncate(MPS psi, MPO const& H, int dK, int n_sweeps, double threshold) const
        {
        auto psi_before = psi;
        psi.position(1);
        LocalMPO PH(H,{"NumCenter",1});
        double avg_truncated_weight = 0.0;
        bool forward = true;
        for (int s=0;s<n_sweeps;++s)
            {
            avg_truncated_weight = kpm_truncate_half_sweep(psi,PH,dK,threshold,forward);
            forward = !forward;
            }
        double dchange = innerC(psi_before,psi_before).real()
                        + innerC(psi,psi).real()
                        - 2.0*innerC(psi_before,psi).real();
        KPMTruncStats stats; stats.avg_truncated_weight = avg_truncated_weight;
        stats.state_change_norm = std::max(dchange,0.0);
        return {psi,stats};
        }

    // Independent (not shared with kpm_moments_full/kpm_moments_accelerated
    // above) Chebyshev recursion with an energy-truncation call inserted
    // after every new vector is formed, before it is used for a moment or
    // fed into the next recursion step (mirrors the paper's own recursion
    // ordering, and pyitensor/chain.py's _kpm_moments_full/
    // _kpm_moments_accelerated wiring of the same call).
    std::vector<std::complex<double>>
    kpm_moments_truncated_full(MPO const& m, MPS const& vi, MPS const& vj, int n,
                      int kpmmaxm, double kpmcutoff, int dK, int n_sweeps, double threshold) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto v = 1.0*vi;
        auto am = 1.0*vi;
        auto a = apply_mpo(m,v,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
        auto ap = 1.0*a;
        double bound = std::sqrt(innerC(vi,vi).real()*innerC(vj,vj).real());
        out.push_back(innerC(vj,v));
        out.push_back(innerC(vj,a));
        for (int i=0;i<n;i++)
            {
            ap = apply_mpo(m,a,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = kpm_energy_truncate(ap,m,dK,n_sweeps,threshold).first;
            out.push_back(innerC(vj,ap));
            check_kpm_moment(out,bound);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    std::vector<std::complex<double>>
    kpm_moments_truncated_accelerated(MPO const& m, MPS const& vi, int n,
                            int kpmmaxm, double kpmcutoff, int dK, int n_sweeps, double threshold) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto a = apply_mpo(m,vi,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
        auto am = 1.0*vi;
        auto ap = 1.0*a;
        Cplx mu0 = innerC(vi,vi);
        Cplx mu1 = innerC(vi,a);
        double bound = std::abs(mu0);
        out.push_back(mu0);
        out.push_back(mu1);
        for (int i=0;i<n/2;i++)
            {
            ap = apply_mpo(m,a,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"MaxDim",kpmmaxm,"Cutoff",kpmcutoff});
            ap = kpm_energy_truncate(ap,m,dK,n_sweeps,threshold).first;
            Cplx bk = 2.0*innerC(a,a) - mu0;
            Cplx bk1 = 2.0*innerC(a,ap) - mu1;
            out.push_back(bk);
            out.push_back(bk1);
            check_kpm_moment(out,bound);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    std::vector<std::complex<double>>
    kpm_moments_truncated(MPO const& m, MPS const& vi, MPS const& vj, int n,
                int kpmmaxm, double kpmcutoff, bool accelerate,
                int dK, int n_sweeps, double threshold) const
        {
        if (accelerate && same_mps(vi,vj,maxm_,cutoff_))
            return kpm_moments_truncated_accelerated(m,vi,n,kpmmaxm,kpmcutoff,dK,n_sweeps,threshold);
        return kpm_moments_truncated_full(m,vi,vj,n,kpmmaxm,kpmcutoff,dK,n_sweeps,threshold);
        }

    private:

    SiteSet sites_;
    MPO H_; bool have_H_ = false;
    MPS wf0_; bool have_wf0_ = false;
    double wf0_energy_ = 0.0; bool have_wf0_energy_ = false;

    bool have_bandwidth_min_ = false, have_bandwidth_max_ = false;
    double bandwidth_emin_ = 0.0, bandwidth_emax_ = 0.0;

    int maxm_ = 30;
    int nsweeps_ = 15;
    double cutoff_ = 1e-12;
    double noise_ = 1e-1;
    int mpomaxm_ = 5000;
    bool verbose_ = false;

    // -- IBC-window real-time dynamical correlator (td_dynamical_correlator_window,
    // below) -- idmrg_ground_state's own converged-environment snapshot
    // (HL/HR "as they stood entering the last executed macro-iteration",
    // and the per-sublattice U tensors solved during that same
    // macro-iteration), cached here as private Chain state -- mirrors
    // pyitensor's own IDMRGResult (see idmrg_ground_state's own doc
    // comment), but kept private rather than added to the public
    // idmrg_ground_state() pybind11 binding's return value, so that
    // binding's existing (density,converged,niter_done) 3-tuple stays
    // unchanged (same pattern as wf0_/have_wf0_ above for gs_energy()).
    bool have_idmrg_snapshot_ = false;
    int idmrg_n_uc_ = 0;
    std::vector<IdmrgAutomatonRow> idmrg_rows_; // per-sublattice automaton row, dense data (see idmrg_make_W's own comment on why never a persistent ITensor)
    std::vector<ITensor> idmrg_U_; // idmrg_U_[p]: sublattice p's converged left-canonical ket tensor (legs: idmrg_U_left_[p], Site, idmrg_U_right_[p])
    std::vector<Index> idmrg_U_left_, idmrg_U_right_; // idmrg_U_[p]'s own natural (left,right) Link Indices, as solved (see idmrg_ground_state's own capture comment for why these must be tracked explicitly rather than re-derived from tags)
    ITensor idmrg_HL_, idmrg_HR_; // rank-3 (bra,mpo,ket) environment operators, entering the last executed macro-iteration
    Index idmrg_HL_bra_, idmrg_HL_ket_, idmrg_HL_mpo_;
    Index idmrg_HR_bra_, idmrg_HR_ket_, idmrg_HR_mpo_;
    bool idmrg_have_HL_ = false, idmrg_have_HR_ = false; // false only if idmrg_ground_state's own maxiter<2 guard were ever bypassed -- see that method's own comment

    // -- VUMPS ground state + tangent-space excitation ansatz snapshot
    // (Chain::vumps_ground_state/vumps_excitation_energies, public,
    // above) -- mirrors the idmrg_* snapshot fields immediately above,
    // but for the mixed-gauge {AL,AR,C,GL,GR} representation VUMPS
    // produces instead of idmrg_ground_state's own per-sublattice U_list.
    bool have_vumps_snapshot_ = false;
    int vumps_D_ = 0, vumps_dg_ = 0, vumps_Dw_ = 0;
    std::vector<Cplx> vumps_AL_, vumps_AR_, vumps_C_; // (D,d_g,D)/(D,d_g,D)/(D,D), row-major
    std::vector<Cplx> vumps_GL_, vumps_GR_;           // (D,D) row-major
    std::vector<Cplx> vumps_W_;                       // (Dw,Dw,d_g,d_g) row-major, see vumps_group_automaton

    // Excitation environment -- built lazily on the first
    // vumps_excitation_energies() call (have_vumps_exc_env_ false until
    // then), cached until the next vumps_ground_state() call invalidates
    // it. See Chain::vumps_build_excitation_environment.
    bool have_vumps_exc_env_ = false;
    std::vector<Cplx> vumps_VL_; // (D*d_g, Dx) isometry, Dx=D*(d_g-1)
    std::vector<Cplx> vumps_GLfull_S_, vumps_GLfull_F_; // (D,D) each
    std::vector<Cplx> vumps_GRfull_S_, vumps_GRfull_F_;
    std::vector<std::vector<Cplx>> vumps_GLfull_pending_, vumps_GRfull_pending_; // one (D,D) per pending channel
    std::vector<Cplx> vumps_E_RL_, vumps_E_LR_;       // (D,D,D,D) mixed transfer tensors (kept for reference; not reused after vx_mixed_fixed_points)
    std::vector<Cplx> vumps_r_RL_, vumps_l_RL_, vumps_r_LR_, vumps_l_LR_; // (D,D)
    double vumps_lam_AC_ = 0.0; // H_AC's own Rayleigh quotient on converged AC -- see vumps_build_excitation_environment's own comment
    };
