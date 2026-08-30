#include <tuple>
#include <map> // Chain::MettsEigCache (Chain::metts_vev)
#include <algorithm> // std::next_permutation (four_correlation_tensor_sweep)
#include <random> // std::mt19937_64/std::discrete_distribution (Chain::metts_vev)
#include <cmath> // std::isnan (Chain::gs_energy_generalized's lam0 sentinel)
#include <limits> // std::numeric_limits<double>::quiet_NaN() (ditto)
#include <stdexcept> // std::invalid_argument/runtime_error: ITensor's own Error()
                     // calls abort() (util/error.h), so the conserved-sector
                     // validation below throws instead, to reach Python as a
                     // catchable exception rather than killing the interpreter
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
// scope this port covers. Deliberately carries only the scalar summary:
// the converged *state* (unit cell, growth environments, final local
// superblock) lives on the Chain itself as a private snapshot, exactly
// like VumpsResult's, so idmrg_onsite_expectation/
// idmrg_two_point_correlator/idmrg_local_excitation_gap read it from
// there rather than through Python.
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
// carry_ferm: whether a Jordan-Wigner string is open between rel_a and
// rel_b (set by idmrg_classify_terms, consumed by idmrg_build_row).
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
        : sites_(SpinX(site_types)), site_types_(site_types), dense_sites_(sites_)
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

    // Ground-state bond-dimension ramp (see make_sweeps_ramped()). Kept
    // as its own setter rather than widening set_sweep_params(), which
    // has ~15 call sites across vev/kpm/excited/nhdmrg on the Python side
    // and in both bindings.cc files. enabled=false reproduces the
    // original flat-maxdim schedule exactly.
    void
    set_bond_ramp(bool enabled, int start, double fraction, double noise_decay)
        {
        ramp_ = enabled;
        ramp_start_ = start;
        ramp_fraction_ = fraction;
        ramp_noise_decay_ = noise_decay;
        }

    void
    set_verbose(bool verbose) { verbose_ = verbose; }

    // -- conserved sector: opt-in quantum-number conservation ----------
    //
    // `qns` is a list of (quantum-number name, target value) pairs -- e.g.
    // {{"Nf",6}} for "the ground state at exactly 6 particles", {{"Sz",0}}
    // for the total-Sz=0 sector (ITensor's own integer 2*Sz units, so a
    // value of 1 means one spin-1/2's worth), {{"Nf",6},{"Sz",0}} on a
    // Hubbard chain. An empty list turns sector mode back off, restoring
    // this session's default dense behavior exactly.
    //
    // This rebuilds sites_ with QN-carrying indices (get_sites.h's
    // SpinX(site_types,conserved)), which makes every MPS and MPO
    // previously built against this session's old Index objects
    // meaningless -- so the Hamiltonian, the stored ground state and every
    // cached band edge/iDMRG/VUMPS snapshot are dropped here, exactly as
    // if the chain had just been constructed. Python re-sends the
    // Hamiltonian on the next gs_energy() because the sector is part of
    // groundstate.py's send-cache key.
    //
    // Two consequences, both enforced rather than left to be discovered:
    // every operator built on a sector-mode chain must itself conserve the
    // sector's quantum numbers (sector_terms() rejects the ones that
    // don't, naming them, instead of letting ITensor abort the process),
    // and DMRG starts from a sum of random in-sector product states rather
    // than randomMPS -- this ITensor has no randomMPS(InitState,m>1) at
    // all (mps.cc's own Error), and randomMPS(SiteSet) refuses to guess a
    // sector.
    void
    set_conserved_sector(std::vector<std::pair<std::string,int>> const& qns)
        {
        std::set<std::string> names;
        for (auto const& q : qns)
            {
            if (q.first!="Nf" && q.first!="Sz" && q.first!="Nb")
                throw std::invalid_argument(tinyformat::format("Chain::set_conserved_sector: unknown quantum "
                      "number \"%s\" (known: Nf, Sz, Nb)",q.first));
            if (!names.insert(q.first).second)
                throw std::invalid_argument(tinyformat::format("Chain::set_conserved_sector: quantum number "
                      "\"%s\" given twice",q.first));
            }
        if (names.empty()) // sector mode off: back to plain dense sites
            {
            // dense_sites_, not a freshly minted SpinX(site_types_): it is
            // this chain's own original site set, so clearing the sector
            // restores exactly the indices the chain had before it was set
            // (and promote_to_dense() below can rebase a state onto them).
            sites_ = dense_sites_;
            sector_.clear();
            has_sector_ = false;
            forget_everything_built_on_sites();
            return;
            }
        // Every site must conserve something: ITensor has no SiteSet
        // mixing QN-carrying and dense indices, so a request that would
        // leave some site dense is rejected here rather than failing deep
        // inside a contraction later.
        std::set<std::string> offered;
        for (size_t k=0;k<site_types_.size();++k)
            {
            auto own = site_qn_names(site_types_.at(k));
            bool any = false;
            for (auto const& nm : own) if (names.count(nm)) { any = true; offered.insert(nm); }
            if (!any)
                throw std::invalid_argument(tinyformat::format("Chain::set_conserved_sector: site %d (type code "
                      "%d) conserves none of the requested quantum numbers; it offers %s",
                      int(k+1),site_types_.at(k),
                      own.empty() ? std::string("none at all (no sector support for this "
                                                "site type)") : join_strings(own)));
            }
        for (auto const& nm : names)
            if (!offered.count(nm))
                throw std::invalid_argument(tinyformat::format("Chain::set_conserved_sector: no site in this "
                      "chain carries a quantum number named \"%s\"",nm));
        // Switch, then check the target is actually reachable -- and roll
        // back if it is not, so a rejected request leaves the chain exactly
        // as it was rather than half-switched into a sector whose every
        // later call would fail. (The check needs the new sites_/sector_ in
        // place, hence the try/catch rather than a pre-check.)
        auto old_sites = sites_;
        auto old_sector = sector_;
        bool old_has = has_sector_;
        sites_ = SpinX(site_types_,names);
        sector_ = qns;
        has_sector_ = true;
        try { sector_state_plan(); }
        catch (...)
            {
            sites_ = old_sites;
            sector_ = old_sector;
            has_sector_ = old_has;
            throw; // nothing was invalidated, so the old H/wf0 stay valid
            }
        forget_everything_built_on_sites();
        }

    // The currently requested sector, empty when sector mode is off.
    std::vector<std::pair<std::string,int>>
    conserved_sector() const { return sector_; }

    // Leave sector mode *keeping* the state that was computed inside the
    // sector -- the one thing set_conserved_sector() with no quantum
    // numbers cannot do, since it drops everything built on the old sites.
    //
    // The point is the operations a sector forbids. While a sector is set,
    // sector_terms() rejects any operator that does not conserve it (a bare
    // C on an Nf-fixed chain, Sx on an Sz-fixed one, a pairing or transverse
    // field), and it has to: AutoMPO aborts the whole process over a
    // flux-violating term rather than raising anything catchable. After
    // promoting, the very same wavefunction lives on ordinary dense sites
    // and every one of those operators is legal again -- so the natural
    // workflow "solve for the ground state at exactly N particles, then
    // apply c_i to it / quench it with a pairing term / take a dynamical
    // correlator of C" works, with the sector used for the part that
    // benefits from it and dropped for the part that cannot tolerate it.
    //
    // The conversion is exact, not an approximation or a re-solve: a
    // QN-conserving MPS is the same wavefunction as its dense counterpart,
    // stored block-sparsely. ITensor's own removeQNs() (mps.cc) scatters
    // every block back into a full dense tensor and strips the QN store
    // from each index, and replaceSiteInds() rebases the physical legs onto
    // dense_sites_ -- the indices every operator built afterwards carries.
    // Nothing is truncated and no number changes; only the storage does,
    // together with the block-sparsity speedup, which is gone from here on.
    //
    // What does NOT survive: the Hamiltonian MPO and the band-edge/iDMRG/
    // VUMPS caches, all of which were built against the QN indices. Python
    // re-sends the Hamiltonian on the next gs_energy() anyway, because the
    // sector is part of groundstate.py's send-cache key.
    void
    promote_to_dense()
        {
        if (!has_sector_) return; // already dense: nothing to promote
        MPS promoted;
        bool had_wf0 = have_wf0_;
        if (had_wf0) promoted = to_dense_mps(wf0_);
        double energy = wf0_energy_;
        bool had_energy = have_wf0_energy_;
        sites_ = dense_sites_;
        sector_.clear();
        has_sector_ = false;
        forget_everything_built_on_sites();
        if (had_wf0) { wf0_ = promoted; have_wf0_ = true; }
        wf0_energy_ = energy;
        have_wf0_energy_ = had_energy;
        }

    // The same conversion applied to one MPS handle Python is already
    // holding (a gs_wavefunction() result, an apply_operator() output, an
    // excited state): those carry whichever site indices were current when
    // they were produced, so promoting the session alone would leave them
    // on stale QN indices and any later contraction against a freshly built
    // operator would not find a matching index. Safe to call on an already
    // dense MPS too -- removeQNs() returns it untouched and
    // replaceSiteInds() is a no-op once the indices already match.
    MPS
    promote_mps(MPS const& wf) const { return to_dense_mps(wf); }

    void
    set_hamiltonian(std::vector<MOTerm> const& terms)
        {
        set_hamiltonian_mpo(mpo_from_terms(terms));
        }

    // Set the Hamiltonian from an ALREADY-BUILT MPO, bypassing the symbolic
    // term list entirely -- the counterpart of set_hamiltonian(terms) for an
    // operator that was assembled with MPO algebra (build_operator +
    // multiply_operators/sum_operators/scale_operator, i.e. Python's
    // toMPO()/StaticOperator).
    //
    // Why this exists: some operators are cheap as an MPO and expensive as a
    // term list. The standing example is a total-spin penalty
    // g*S^2_total = g*((sum_i Sx_i)^2 + (sum_i Sy_i)^2 + (sum_i Sz_i)^2),
    // used to project onto a spin sector in a ground-state-only solver. As
    // symbols it is O(n^2) terms (measured: 2030 terms at 12 orbitals against
    // 302 for the Hamiltonian itself); as MPOs it is three squares of
    // extensive one-body operators, each of small bond dimension, so building
    // it as (Sx*Sx + Sy*Sy + Sz*Sz) is exact and stays compact regardless of
    // system size.
    //
    // Note this is a convenience-and-structure route, not automatically a
    // speed one: ITensor's own toMPO already SVD-compresses the term list,
    // and on that same measurement 6.7x the terms cost only 1.4x the sweep
    // time, so the term-count blow-up largely does NOT propagate into the
    // solve. What this buys is the ability to express an operator the way it
    // is naturally structured, and to avoid materializing a term list that
    // grows quadratically when the operator itself does not.
    void
    set_hamiltonian_mpo(MPO const& H)
        {
        H_ = H;
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
        // A warm start (an MPS handed in by set_wavefunction, or simply
        // the previous gs_energy() call's own solution) must not be
        // truncated by the ramp's early low-maxdim sweeps -- pass its
        // bond dimension as the ramp floor. A fresh start has no such
        // constraint, and is built directly at the ramp's starting bond
        // dimension rather than at the full maxm_ only to be truncated
        // on the first sweep.
        int floor_dim = 0;
        if (have_wf0_) floor_dim = maxLinkDim(wf0_);
        auto sweeps = make_sweeps_ramped(nsweeps_,maxm_,floor_dim);
        if (!have_wf0_)
            {
            wf0_ = default_mps(sweeps.maxdim(1));
            have_wf0_ = true;
            }
        double energy = dmrg(wf0_,H_,sweeps,dmrg_args());
        check_sector(wf0_,"gs_energy"); // no-op unless sector mode is on
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
        // An MPS built while a different sector setting was in force is
        // not merely in the wrong sector, its indices belong to a
        // different SiteSet entirely -- catch both here rather than in
        // the middle of a contraction.
        if (has_sector_ != hasQNs(wf(1)))
            throw std::invalid_argument("Chain::set_wavefunction: this wavefunction was built with a different "
                  "conserved-sector setting than the chain currently has");
        check_sector(wf,"set_wavefunction");
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
        auto A = mpo_from_terms(terms_a);
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
        auto H = mpo_from_terms(terms_h);
        auto HA = mpo_from_terms(terms_hadj);
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
        auto H = mpo_from_terms(terms_h);
        auto HA = mpo_from_terms(terms_hadj);
        auto A = mpo_from_terms(terms_a);
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
        auto A = mpo_from_terms(terms);
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
        auto A = mpo_from_terms(terms);
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
        auto A = mpo_from_terms(terms);
        return innerC(wf1,A,wf2);
        }

    MPS
    sum_mps(MPS const& wf1, MPS const& wf2) const
        {
        return sum(wf1,wf2,{"MaxDim",maxm_,"Cutoff",cutoff_});
        }

    // Rescale an MPS by a complex number. ITensor's own operator*=(MPS,Cplx)
    // multiplies a *single* site tensor (mps.cc: x.ref(x.leftLim()+1) *= z),
    // so this is O(chi^2 d) with no contraction, no bond growth and no
    // truncation. It exists because mps.py's MPS.__mul__ used to express
    // "multiply a wavefunction by a number" as a one-term "Id" AutoMPO sent
    // through apply_operator() -- i.e. a full truncating MPO sweep over the
    // whole chain per scalar. Measured on a 16-site Heisenberg CVM run
    // (cvm.py, cvm_maxm=40), that accounted for 15.4 s of 43.8 s: the CG
    // recurrence rescales an MPS five times per iteration, so a third of the
    // solver's wall time was spent applying identity operators.
    MPS
    scale_mps(MPS const& wf, Cplx z) const
        {
        auto out = wf;
        out *= z;
        return out;
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
        // Read the dim*dim matrix elements out explicitly, by index,
        // rather than with rho.visit(): visit() enumerates only the
        // *stored* elements, which for QN-carrying site indices
        // (set_conserved_sector) is a handful of 1x1 blocks rather than
        // the full dim*dim, so the caller in bindings.cc used to copy a
        // short vector into an uninitialized (dim,dim) array -- values
        // landing at the wrong offsets and the rest of the array left as
        // whatever was on the heap (confirmed by poisoning it: the poison
        // bytes came back inside the "density matrix"). Element (a,b) is
        // <a|rho|b>, matching the site-operator convention op_charge()
        // documents (unprimed = in, primed = out) and the ordering
        // pyitensor's own reduced_dm() produces.
        auto s = sites_.si(site);
        int d = dim(s);
        std::vector<std::complex<double>> out;
        out.reserve(d*d);
        for (int a=1;a<=d;a++)
        for (int b=1;b<=d;b++)
            out.push_back(eltC(rho,s(a),prime(s)(b)));
        return out;
        }

    MPS
    exponential_apply(std::vector<MOTerm> const& terms, MPS const& wf,
                       std::complex<double> tau, int nsteps) const
        {
        auto H = mpo_from_terms(terms);
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
        return mpo_from_terms(terms);
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
        auto H = mpo_from_terms(terms_h);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = ampo_from_terms(terms_h);
        ampo += -EGS,"Id",1;
        auto expH = evoloperator(toMPO(ampo),dt);
        auto A1 = mpo_from_terms(terms_i);
        auto A2 = mpo_from_terms(terms_j);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure *before* evolving, so correlator[k] is C(k*dt),
            // matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            // pairs it with -- see the "measure before evolving" comment
            // on evolution_dmrg_DC() there for why the old
            // evolve-then-measure order was a real (not cosmetic) bug.
            out.correlator.push_back(innerC(psi2,psi1));
            if (fit_td) psi1 = apply_mpo(expH,psi1,psi1,args);
            else psi1 = apply_mpo(expH,psi1,args);
            psi1.normalize();
            psi1 *= norm0;
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
        auto ampo = ampo_from_terms(terms_h);
        auto expH = evoloperator(toMPO(ampo),dt);
        auto A = mpo_from_terms(terms_op);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure *before* evolving, so correlator[k] is C(k*dt),
            // matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            // pairs it with -- see the "measure before evolving" comment
            // on evolution_dmrg_DC() there for why the old
            // evolve-then-measure order was a real (not cosmetic) bug.
            out.correlator.push_back(innerC(psi,A,psi));
            if (fit_td) psi = apply_mpo(expH,psi,psi,args);
            else psi = apply_mpo(expH,psi,args);
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
        auto H = mpo_from_terms(terms_h);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = ampo_from_terms(terms_h);
        ampo += -EGS,"Id",1;
        auto Hshift = toMPO(ampo);
        auto A1 = mpo_from_terms(terms_i);
        auto A2 = mpo_from_terms(terms_j);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure *before* evolving, so correlator[k] is C(k*dt),
            // matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            // pairs it with -- see the "measure before evolving" comment
            // on evolution_dmrg_DC() there for why the old
            // evolve-then-measure order was a real (not cosmetic) bug.
            out.correlator.push_back(innerC(psi2,psi1));
            psi1 = tdvp_step(Hshift,psi1,dt);
            psi1.normalize();
            psi1 *= norm0;
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure_tdvp(std::vector<MOTerm> const& terms_h,
                             std::vector<MOTerm> const& terms_op,
                             MPS const& wf, int nt, double dt)
        {
        auto H = mpo_from_terms(terms_h);
        auto A = mpo_from_terms(terms_op);
        auto psi = wf;
        // tdvp_step() passes {"DoNormalize",true}, so without restoring the
        // input norm after every step this silently rescaled the whole
        // trajectory by 1/||wf||^2 for any non-unit-norm start -- exactly
        // what evolution_ABA() hands in (wfA = A*wf, for a non-unitary A).
        // C(0) was still right (the loop measures before evolving), so the
        // error looked like a discontinuity at the first step rather than
        // an overall factor: on a 5-site XXZ chain with A=Sx[0]
        // (||A|gs>||^2 = 0.25) the t>0 values came back 4x too large,
        // recovering ED to ~1e-10 once multiplied by 0.25. quench_tdvp(),
        // quench_tdvp_gse() and quench_tebd() have always done this; these
        // two were the ones that did not. The same defect was found and
        // fixed on the Julia backend (mpsjulialive/tdvp.jl), see
        // docs/documentation.md.
        auto norm0 = std::sqrt(std::real(innerC(psi,psi)));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure *before* evolving, so correlator[k] is C(k*dt),
            // matching the ts=[0,dt,...,(nt-1)*dt] grid timedependent.py
            // pairs it with -- see the "measure before evolving" comment
            // on evolution_dmrg_DC() there for why the old
            // evolve-then-measure order was a real (not cosmetic) bug.
            out.correlator.push_back(innerC(psi,A,psi));
            psi = tdvp_step(H,psi,dt);
            psi.normalize();
            psi *= norm0;
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
        auto H = mpo_from_terms(terms_h);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto ampo = ampo_from_terms(terms_h);
        ampo += -EGS,"Id",1;
        auto Hshift = toMPO(ampo);
        auto A1 = mpo_from_terms(terms_i);
        auto A2 = mpo_from_terms(terms_j);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure before evolving -- see the matching comment in
            // quench() above.
            out.correlator.push_back(innerC(psi2,psi1));
            if (it<gse_sweeps)
                psi1 = global_subspace_expand(Hshift,psi1,krylov_order,gse_cutoff,0);
            psi1 = tdvp_step(Hshift,psi1,dt,1);
            psi1.normalize();
            psi1 *= norm0;
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
        auto H = mpo_from_terms(terms_h);
        auto A = mpo_from_terms(terms_op);
        auto psi = wf;
        // see evolve_and_measure_tdvp() above for why the input norm has
        // to be restored after every normalizing tdvp_step()
        auto norm0 = std::sqrt(std::real(innerC(psi,psi)));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure before evolving -- see the matching comment in
            // quench() above.
            out.correlator.push_back(innerC(psi,A,psi));
            if (it<gse_sweeps)
                psi = global_subspace_expand(H,psi,krylov_order,gse_cutoff,0);
            psi = tdvp_step(H,psi,dt,1);
            psi.normalize();
            psi *= norm0;
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
        reject_sector("quench_tebd"); // see tebd.h: its gate assembly is dense-only
        if (!have_wf0_) gs_energy();
        auto H = mpo_from_terms(terms_h);
        auto args = Args("Cutoff",cutoff_,"MaxDim",maxm_);
        double EGS = innerC(wf0_,H,wf0_).real()/innerC(wf0_,wf0_).real();
        auto terms_shifted = terms_h;
        MOTerm shift_term;
        shift_term.coef = Cplx(-EGS,0.0);
        shift_term.factors = {{"Id",1}};
        terms_shifted.push_back(shift_term);
        auto h_bonds = bond_hamiltonians(sites_,terms_shifted);
        auto gates = build_tebd_gates(sites_,h_bonds,dt);
        auto A1 = mpo_from_terms(terms_i);
        auto A2 = mpo_from_terms(terms_j);
        auto psi1 = apply_mpo(A1,wf0_,args);
        auto psi2 = apply_mpo(A2,wf0_,args);
        Cplx norm0 = std::sqrt(innerC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure before evolving -- see the matching comment in
            // quench() above.
            out.correlator.push_back(innerC(psi2,psi1));
            tebd_step(psi1,gates,cutoff_,maxm_);
            psi1.normalize();
            psi1 *= norm0;
            }
        out.final_wf = psi1;
        return out;
        }

    TimeEvolutionResult
    evolve_and_measure_tebd(std::vector<MOTerm> const& terms_h,
                             std::vector<MOTerm> const& terms_op,
                             MPS const& wf, int nt, double dt)
        {
        reject_sector("evolve_and_measure_tebd"); // see tebd.h: its gate assembly is dense-only
        auto h_bonds = bond_hamiltonians(sites_,terms_h);
        auto gates = build_tebd_gates(sites_,h_bonds,dt);
        auto A = mpo_from_terms(terms_op);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            // Measure before evolving -- see the matching comment in
            // quench() above.
            out.correlator.push_back(innerC(psi,A,psi));
            tebd_step(psi,gates,cutoff_,maxm_);
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
        reject_sector("metts_vev");
        if (!have_H_) Error("Chain::metts_vev called before set_hamiltonian");
        if (T<=0) Error("Chain::metts_vev: T must be > 0");
        if (basis_ops.empty()) Error("Chain::metts_vev: basis_ops must be non-empty");
        if (nsamples<1) Error("Chain::metts_vev: nsamples must be >= 1");
        if (terms_ops.empty()) Error("Chain::metts_vev: terms_ops must be non-empty");
        std::vector<MPO> ops;
        ops.reserve(terms_ops.size());
        for (auto const& terms_op : terms_ops) ops.push_back(mpo_from_terms(terms_op));
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
        reject_sector("metts_dynamical_correlator");
        if (!have_H_) Error("Chain::metts_dynamical_correlator called before set_hamiltonian");
        if (T<=0) Error("Chain::metts_dynamical_correlator: T must be > 0");
        if (basis_ops.empty()) Error("Chain::metts_dynamical_correlator: basis_ops must be non-empty");
        if (nsamples<1) Error("Chain::metts_dynamical_correlator: nsamples must be >= 1");
        if (nt<1) Error("Chain::metts_dynamical_correlator: nt must be >= 1");

        auto A = mpo_from_terms(terms_a);
        auto B = mpo_from_terms(terms_b);
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
        auto S1 = mpo_from_terms(terms_i);
        auto S2 = mpo_from_terms(terms_j);
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
        auto A = mpo_from_terms(terms);
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
    // distinct* go through the fast sweep below; the O(N^3) entries with a
    // repeated index are handled separately at the tail of this function.
    // Repeats need same-site multi-operator products (e.g. Cdag_i C_i = N_i,
    // or Cdag_i Cdag_i = 0) that the sweep below -- keyed on four *strictly
    // increasing* site positions -- isn't built to produce. They used to
    // fall back to the per-tuple AutoMPO method of
    // four_correlation_tensor(), on the reasoning that a lower order of the
    // total work can afford a slower method; measured, that was backwards
    // (see the tail), and they are folded locally now.
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
    // `accelerate` here only gates the O(N^3) repeated-
    // index entries below, unlike four_correlation_tensor()/
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

        // Repeated-index entries (not covered by the sweep above). These
        // used to go through the same per-tuple AutoMPO+toMPO+innerC method
        // as four_correlation_tensor(), on the reasoning that there are only
        // O(N^3) of them against the sweep's O(N^4) -- see the "subdominant"
        // wording this comment replaces. Measured, that is backwards: fewer
        // tuples times a far more expensive per-tuple cost dominates. On the
        // pure-Python twin of this method (pyitensor/chain.py, same
        // algorithm) the fallback was 96% of total runtime at N=12, and the
        // C++ ratio is the same shape -- an AutoMPO build plus a full-chain
        // innerC per tuple, against a handful of local contractions for a
        // whole (a,b,c,d) leaf.
        //
        // Every one of these operators is a product of four Cdag/C factors
        // on at most 3 distinct sites, i.e. *local*: there is no need to
        // compile an MPO over the whole chain and sweep it. Resolve it to
        // per-site operators once (Jordan-Wigner threading included) and
        // fold it over [mn,mx] only, reusing this method's own `fold`
        // pattern and the same mixed-canonical shortcut at both ends.
        //
        // The sign is the parity of the permutation that sorts the four
        // operators into site order (every pair at distinct sites
        // anticommutes; same-site pairs keep their relative order and
        // contribute nothing) -- NOT the `perms` table above, whose extra
        // mask-3/6/9/12 correction is specific to that table's own
        // strictly-increasing-site construction. Threading F by running
        // parity, as below, needs only the plain permutation parity;
        // verified directly against AutoMPO+toMPO+innerC on the Python side
        // over 67 tuples before being ported here.
        {
        struct Quad { int i,j,k,l; };
        std::vector<std::vector<Quad>> buckets(N);
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
            int mn0 = std::min(std::min(i,j),std::min(k,l));
            buckets[mn0].push_back({i,j,k,l});
            }
        // Grouped by minimum site so one position() call serves each bucket.
        for (int mn0=0; mn0<N; ++mn0)
            {
            if (buckets[mn0].empty()) continue;
            int mn = mn0+1;
            psi.position(mn);
            for (auto const& q : buckets[mn0])
                {
                int i=q.i, j=q.j, k=q.k, l=q.l;
                int mx = std::max(std::max(i,j),std::max(k,l))+1;
                std::pair<std::string,int> factors[4] =
                    {{"Cdag",i+1},{"C",j+1},{"Cdag",k+1},{"C",l+1}};
                ITensor E;
                if (mn>1)
                    {
                    auto ln = commonIndex(psi.A(mn-1),psi.A(mn));
                    E = delta(dag(prime(ln)),ln);
                    }
                bool carry = false;
                for (int s=mn; s<=mx; ++s)
                    {
                    // Compose this site's factors in slot order, then the
                    // JW F if the running parity demands one. Both are
                    // appended on the RIGHT of the accumulated operator (F
                    // acts on the ket first, the earliest-slot factor last),
                    // matching what autompo.HTerm.resolve() produces.
                    //
                    // The ITensor idiom for that is `A*prime(M)`, NOT
                    // `M*prime(A)`: with op() returning (out=s', in=s),
                    // contracting M's *out* against prime(A)'s *in* yields
                    // A applied last, i.e. the reverse of what is wanted.
                    // Getting this backwards is not a subtle numerical
                    // effect -- it moved entries by ~0.86 against the
                    // (independently validated) pure-Python twin.
                    ITensor M; // default-constructed = identity, skip apply
                    int nhere = 0;
                    for (auto const& f : factors) if (f.second==s)
                        {
                        ++nhere;
                        auto A = op(sites_,f.first,s);
                        if (!M) M = A;
                        else { M = A*prime(M); M.mapPrime(2,1); }
                        }
                    bool odd = (nhere%2==1); // every Cdag/C factor is fermionic
                    if (carry != odd)
                        {
                        auto F = op(sites_,"F",s);
                        if (!M) M = F;
                        else { M = F*prime(M); M.mapPrime(2,1); }
                        }
                    carry = (carry != odd);
                    ITensor T = psi.A(s);
                    if (M) { T = T*M; T.noPrime(TagSet("Site")); }
                    E = (E ? E*T : T);
                    E = E*dag(prime(psi.A(s),TagSet("Link")));
                    }
                if (mx<N)
                    {
                    auto rn = commonIndex(psi.A(mx),psi.A(mx+1));
                    E = E*delta(dag(prime(rn)),rn);
                    }
                int ss[4] = {i,j,k,l};
                int inv = 0;
                for (int x=0;x<4;++x) for (int y=x+1;y<4;++y) if (ss[x]>ss[y]) ++inv;
                Cplx c = ((inv%2) ? -1.0 : 1.0)*eltC(E);
                out[idx(i,j,k,l)] = c;
                std::tuple<int,int,int,int> current{i,j,k,l};
                std::tuple<int,int,int,int> conjugate{l,k,j,i};
                if (!accelerate || current!=conjugate) out[idx(l,k,j,i)] = std::conj(c);
                }
            }
        }
        return out;
        }

    // <Cdag_i C_j Cdag_k C_l> over flat fermionic *modes*, by local operator
    // folds -- no AutoMPO/MPO is built for any tuple.
    //
    // `cdag_names[m]`/`cdag_sites[m]` (and the C pair) give the operator
    // name and 1-based site of mode m, so this covers a spinless chain
    // (mode == site, "Cdag"/"C") and a native spinful one (two modes per
    // ElectronSite, "Cdagup"/"Cup"/"Cdagdn"/"Cdn") with no flavor-aware
    // logic here: Jordan-Wigner threading comes from the site's own "F" and
    // running parity, exactly as AutoMPO's rewriteFermionic derives it.
    //
    // This is the port of pyitensor/chain.py's four_correlation_tensor_fold,
    // added because for native spinful sites the compiled backend only had
    // four_correlation_tensor_spinful -- a per-tuple AutoMPO build -- and
    // was measured 3.3x SLOWER than the pure-Python fold (5 sites: 2.56s vs
    // 0.77s). Algorithm, not language.
    //
    // The sign is the parity of the permutation sorting the four operators
    // into site order (see four_correlation_tensor_sweep's own note on why
    // that is the right rule here and its `perms` table is not).
    std::vector<std::complex<double>>
    four_correlation_tensor_fold(MPS const& wf,
                                 std::vector<std::string> const& cdag_names,
                                 std::vector<int> const& cdag_sites,
                                 std::vector<std::string> const& c_names,
                                 std::vector<int> const& c_sites,
                                 bool accelerate=true) const
        {
        int nm = int(c_names.size());
        if (int(cdag_names.size())!=nm || int(cdag_sites.size())!=nm
            || int(c_sites.size())!=nm)
            Error("four_correlation_tensor_fold: mode lists have different lengths");
        int N = sites_.length();
        std::vector<std::complex<double>> out(size_t(nm)*nm*nm*nm, 0.0);
        auto idx = [nm](int i,int j,int k,int l) { return ((size_t(i)*nm+j)*nm+k)*nm+l; };

        auto psi = wf;
        struct Quad { int i,j,k,l; };
        std::vector<std::vector<Quad>> buckets(N+1);
        for (int i=0;i<nm;i++)
        for (int j=0;j<nm;j++)
        for (int k=0;k<nm;k++)
        for (int l=0;l<nm;l++)
            {
            std::tuple<int,int,int,int> current{i,j,k,l}, conjugate{l,k,j,i};
            if (accelerate && current>conjugate) continue;
            int mn = std::min(std::min(cdag_sites[i],c_sites[j]),
                              std::min(cdag_sites[k],c_sites[l]));
            buckets[mn].push_back({i,j,k,l});
            }

        for (int mn=1; mn<=N; ++mn)
            {
            if (buckets[mn].empty()) continue;
            psi.position(mn);
            for (auto const& q : buckets[mn])
                {
                int i=q.i, j=q.j, k=q.k, l=q.l;
                std::pair<std::string,int> factors[4] = {
                    {cdag_names[i], cdag_sites[i]}, {c_names[j], c_sites[j]},
                    {cdag_names[k], cdag_sites[k]}, {c_names[l], c_sites[l]}};
                int mx = std::max(std::max(factors[0].second,factors[1].second),
                                  std::max(factors[2].second,factors[3].second));
                ITensor E;
                if (mn>1)
                    {
                    auto ln = commonIndex(psi.A(mn-1),psi.A(mn));
                    E = delta(dag(prime(ln)),ln);
                    }
                bool carry = false;
                for (int s=mn; s<=mx; ++s)
                    {
                    ITensor M;
                    int nhere = 0;
                    for (auto const& f : factors) if (f.second==s)
                        {
                        ++nhere;
                        auto A = op(sites_,f.first,s);
                        // append on the RIGHT: A*prime(M), see this file's
                        // four_correlation_tensor_sweep tail for why the
                        // opposite order silently applies them backwards
                        if (!M) M = A;
                        else { M = A*prime(M); M.mapPrime(2,1); }
                        }
                    bool odd = (nhere%2==1);
                    if (carry != odd)
                        {
                        auto F = op(sites_,"F",s);
                        if (!M) M = F;
                        else { M = F*prime(M); M.mapPrime(2,1); }
                        }
                    carry = (carry != odd);
                    ITensor T = psi.A(s);
                    if (M) { T = T*M; T.noPrime(TagSet("Site")); }
                    E = (E ? E*T : T);
                    E = E*dag(prime(psi.A(s),TagSet("Link")));
                    }
                if (mx<N)
                    {
                    auto rn = commonIndex(psi.A(mx),psi.A(mx+1));
                    E = E*delta(dag(prime(rn)),rn);
                    }
                int ss[4] = {factors[0].second,factors[1].second,
                             factors[2].second,factors[3].second};
                int inv=0;
                for (int x=0;x<4;++x) for (int y=x+1;y<4;++y) if (ss[x]>ss[y]) ++inv;
                Cplx c = ((inv%2)?-1.0:1.0)*eltC(E);
                out[idx(i,j,k,l)] = c;
                std::tuple<int,int,int,int> current{i,j,k,l}, conjugate{l,k,j,i};
                if (!accelerate || current!=conjugate) out[idx(l,k,j,i)] = std::conj(c);
                }
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
        auto m1 = mpo_from_terms(terms_i);
        auto m2 = mpo_from_terms(terms_j);
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
        auto m = mpo_from_terms(terms_x);
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
        auto hs = mpo_from_terms(terms_hs);
        auto hsd = mpo_from_terms(terms_hs_dag);

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
        auto m1 = mpo_from_terms(terms_i);
        auto m2 = mpo_from_terms(terms_j);
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
    // infinite-size algorithm (White 1992, with an n_uc-site unit cell)
    // -- see that module's own docstring for the full derivation and bug
    // history. It began as a line-for-line translation of its logic
    // (automaton construction, growth loop, local 2-site solve, HL/HR
    // extension), not an independent re-derivation, specifically to avoid
    // reintroducing bugs already found and fixed there.
    //
    // The port is kept deliberately in step with the Python side. Three
    // things pyitensor/idmrg.py gained after the first translation --
    // McCulloch's wavefunction prediction (_wavefunction_prediction) as
    // each micro-step's Krylov start, extraction of the converged unit
    // cell from a single micro-step's theta (_theta_cell/
    // _canonical_theta_cell, IDMRGResult.cell_list) rather than chaining
    // the per-micro-step factors this port snapshots as idmrg_U_, and a
    // per-site energy baseline subtracted from both growing environments
    // (_subtract_energy_baseline, the reference idmrg.h's
    // HL += -energy*IL, with the density adding the shift back) -- were
    // for a while absent here, which is why this comment used to declare
    // the two implementations non-equivalent. All three are ported now:
    // idmrg_wavefunction_prediction, idmrg_theta_cell +
    // ic_canonicalize_cell (idmrg_cell_), and
    // idmrg_subtract_energy_baseline respectively. With the unit cell in
    // hand, the static observables built on it are ported too --
    // idmrg_onsite_expectation/idmrg_two_point_correlator, plus
    // idmrg_local_excitation_gap over the stored final superblock -- so
    // this backend is no longer energy-density-only.
    //
    // What remains itensor_version="python"-only: local_excitation_gap's
    // `window>0` variant (pyitensor's own local_excitation_gap_windowed is
    // an explicit prototype rather than stable API) and the iMPS algebra
    // built on IDMRGResult (imps_overlap/apply_mpo/imps_sum/grow_by_mpo),
    // none of which infinitechain.py exposes for gs_method="idmrg"; the
    // VUMPS path above already covers the latter on this backend
    // (vumps_apply_mpo/vumps_imps_sum). See CLAUDE.md's mpscpp3 section.
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
    // pairing). Fermionic terms ARE supported: idmrg_classify_terms
    // threads the Jordan-Wigner string itself (locally, between each
    // term's own two endpoints -- the terms arrive untransformed, see that
    // function's own comment), and idmrg_build_row propagates it along the
    // pending channel; a term with odd total fermion parity, whose string
    // could never close within the unit cell, is rejected there.
    //
    // On return, this Chain's own private snapshots are set: the growth
    // environments (have_idmrg_snapshot_, consumed by
    // td_dynamical_correlator_window's IBC window), the gauge-consistent
    // unit cell (have_idmrg_cell_, consumed by idmrg_onsite_expectation/
    // idmrg_two_point_correlator), and the final local superblock
    // (have_idmrg_superblock_, consumed by idmrg_local_excitation_gap) --
    // mirroring vumps_ground_state's own have_vumps_snapshot_ pattern.
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
    // White noise: how strong, and for how many macro-iterations -- see
    // idmrg_noisy_isometry and idmrg_local_solve's own start-vector
    // perturbation for what the two halves do, and
    // docs/known_issue_idmrg_product_state_collapse.md for the failure they
    // exist to break (a product state the growth loop cannot leave, returned
    // with converged=True in a fraction of a second).
    //
    // The schedule matters as much as the noise. Noise is OFF for the tail
    // of the run, for two independent reasons: the returned unit cell and
    // its Schmidt values are re-seeded every macro-iteration, so a clean
    // tail is what keeps <H_uc> = n_uc*e0 and the vev/correlator gauge
    // exact; and a perturbed state must never be what etol declares
    // converged, so the convergence test is suppressed entirely while noise
    // is active. That suppression also removes the "converged=True in 0 s"
    // pathology by construction: no run can now finish in fewer than
    // noise_iters+2 macro-iterations. Same values as pyitensor/idmrg.py's
    // own _NOISE_DEFAULT/_NOISE_ITERS.
    static constexpr double idmrg_noise_default_ = 1e-4;
    static constexpr int idmrg_noise_iters_ = 40;
    static constexpr int idmrg_noise_tail_ = 6;

    IdmrgResult
    idmrg_ground_state(std::vector<MOTerm> const& terms_intra,
                        std::vector<MOTerm> const& terms_inter,
                        int maxm, double cutoff, int maxiter, double etol,
                        int krylovdim, int restarts,
                        double noise = idmrg_noise_default_,
                        int noise_iters = idmrg_noise_iters_)
        {
        reject_sector("idmrg_ground_state"); // dense environments/unit cell below
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

        // Per-site energy baseline subtracted from BOTH growing
        // environments every micro-step (idmrg_subtract_energy_baseline),
        // latched once as soon as a density estimate exists and never
        // revised afterwards -- see the density computation at the bottom
        // of the macro-iteration loop for why constancy matters, and
        // idmrg_subtract_energy_baseline's own comment for why the
        // baseline is needed at all. 0.0 means "not latched yet".
        double eshift = 0.0;

        // The last density actually computed, kept separately from
        // prev_density (which the baseline latch below deliberately
        // discards along with prev_energy, since both were measured
        // against un-shifted environments). pyitensor/idmrg.py returns its
        // own prev_density directly and therefore reports None in the one
        // edge case where maxiter runs out on exactly the latching
        // iteration; this port reports the best estimate it actually has
        // instead, which agrees with Python in every other case.
        double report_density = 0.0;
        bool have_report_density = false;

        // The *immediately preceding* micro-step's own SVD factors, as
        // dense (chi_l,d,chi_r) arrays plus their singular values -- the
        // raw material of McCulloch's wavefunction prediction (see
        // idmrg_wavefunction_prediction). `s_prev` lags one step further
        // behind still, since it is the bond matrix the prediction has to
        // invert. Deliberately NOT indexed by mstep: the prediction reads
        // the step whose bond bases the current local problem actually
        // lives in, which is always the one just executed.
        std::vector<Cplx> last_U, last_V;
        int last_Ul=0, last_Ud=0, last_Ur=0, last_Vl=0, last_Vd=0, last_Vr=0;
        std::vector<double> last_S, prev_S;
        bool have_last_S=false, have_prev_S=false, have_last_UV=false;

        // (U, S, V, lambda_outer) of the last executed macro-iteration's
        // own FIRST micro-step, kept as dense arrays -- turned into the
        // gauge-consistent 2-site unit cell (idmrg_theta_cell) once, after
        // the loop, exactly as pyitensor/idmrg.py's own cell_seed is (its
        // re-gauging step is an eigendecomposition of a chi^2 x chi^2
        // transfer matrix, far too expensive to pay every macro-iteration).
        std::vector<Cplx> seed_U, seed_V;
        int seed_Ul=0, seed_Ud=0, seed_Ur=0, seed_Vl=0, seed_Vd=0, seed_Vr=0;
        std::vector<double> seed_S, seed_lam_outer;
        bool have_cell_seed=false;

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

        // Capped well short of maxiter so a short run still gets a clean
        // tail: a caller asking for maxiter=10 must not spend 8 of them
        // perturbed and then report whatever the last one happened to be.
        int noise_iters_eff = std::min(noise_iters,std::max(0,maxiter-2));
        int noise_left = 0, noise_spent = 0;
        for (macro_iter=0; macro_iter<maxiter; ++macro_iter)
            {
            double noise_now = (noise_left>0 && noise_spent<noise_iters_eff) ? noise : 0.0;
            if (noise_now > 0.0) ++noise_spent;
            bool product_this_iter = true;
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

                // McCulloch's wavefunction prediction: the previous
                // micro-step's converged state translated into *this*
                // step's enlarged bond bases (idmrg_wavefunction_
                // prediction). Falls back to prev_local[mstep] -- the
                // plain same-position reuse this port had before the
                // prediction existed -- and, inside idmrg_local_solve,
                // to a fresh random vector, for the early micro-steps
                // where the previous step's shapes don't line up yet.
                std::vector<Cplx> x0_warm;
                if (have_last_UV && have_last_S && have_prev_S)
                    x0_warm = idmrg_wavefunction_prediction(
                        last_S,last_V,last_Vl,last_Vd,last_Vr,
                        prev_S,last_U,last_Ul,last_Ud,last_Ur,
                        have_HL_ket ? dim(HL_ket) : 0, dim(phys_L_in),
                        dim(phys_R_in), have_HR_ket ? dim(HR_ket) : 0);
                if (x0_warm.empty()) x0_warm = prev_local[mstep];

                auto [energy_here,U,V,new_bond_u,new_bond_v,theta_flat,svals,theta,purity] =
                    idmrg_local_solve(
                        HL,W_pL,phys_L_in,W_pR,phys_R_in,HR,
                        HL_bra,HL_ket,have_HL_ket,HR_bra,HR_ket,have_HR_ket,
                        cutoff,maxm,krylovdim,restarts,x0_warm,noise_now);
                energy = energy_here;
                if (!theta_flat.empty()) prev_local[mstep] = std::move(theta_flat);
                // A product state at ANY micro-step is enough to keep the
                // noise armed: that is the trap (see idmrg_noisy_isometry).
                if (!(purity > 1.0-1e-10)) product_this_iter = false;

                // The raw ingredients of *this* micro-step's own 2-site
                // local eigenproblem, captured before the extend step
                // below mutates HL/HR -- overwritten every micro-step, so
                // what survives the loop is the final one actually run.
                // Kept so idmrg_local_excitation_gap() can re-diagonalize
                // that same effective Hamiltonian for its second-lowest
                // eigenpair without rerunning the growing algorithm; C++
                // analogue of pyitensor/idmrg.py's own last_superblock/
                // IDMRGResult.local_superblock.
                idmrg_sb_HL_ = HL; idmrg_sb_HL_bra_ = HL_bra; idmrg_sb_HL_ket_ = HL_ket;
                idmrg_sb_have_HL_ket_ = have_HL_ket;
                idmrg_sb_HR_ = HR; idmrg_sb_HR_bra_ = HR_bra; idmrg_sb_HR_ket_ = HR_ket;
                idmrg_sb_have_HR_ket_ = have_HR_ket;
                idmrg_sb_W_pL_ = W_pL; idmrg_sb_phys_L_ = phys_L_in;
                idmrg_sb_W_pR_ = W_pR; idmrg_sb_phys_R_ = phys_R_in;
                idmrg_sb_theta_ = theta; idmrg_sb_energy_ = energy_here;
                have_idmrg_superblock_ = true;

                // Carry this step's SVD factors forward for the next
                // step's prediction, and -- at mstep 0 -- stash the raw
                // ingredients of the gauge-consistent unit cell. `last_S`
                // is still the *previous* micro-step's center matrix at
                // this point, which is exactly the lambda_o on this
                // theta's own outer bond (see idmrg_theta_cell).
                auto U_lpr = idmrg_svd_factor_lpr(U,HL_ket,have_HL_ket,
                                                   phys_L_in,new_bond_u,true);
                auto V_lpr = idmrg_svd_factor_lpr(V,new_bond_v,true,
                                                   phys_R_in,HR_ket,have_HR_ket);
                int Ul = have_HL_ket ? dim(HL_ket) : 1, Ud = dim(phys_L_in), Ur = dim(new_bond_u);
                int Vl = dim(new_bond_v), Vd = dim(phys_R_in), Vr = have_HR_ket ? dim(HR_ket) : 1;
                if (mstep==0 && have_last_S)
                    {
                    seed_U = U_lpr; seed_Ul=Ul; seed_Ud=Ud; seed_Ur=Ur;
                    seed_V = V_lpr; seed_Vl=Vl; seed_Vd=Vd; seed_Vr=Vr;
                    seed_S = svals; seed_lam_outer = last_S;
                    have_cell_seed = true;
                    }
                prev_S = last_S; have_prev_S = have_last_S;
                last_S = svals; have_last_S = true;
                last_U = std::move(U_lpr); last_Ul=Ul; last_Ud=Ud; last_Ur=Ur;
                last_V = std::move(V_lpr); last_Vl=Vl; last_Vd=Vd; last_Vr=Vr;
                have_last_UV = true;

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

                if (eshift != 0.0)
                    {
                    // Channel 0 is "S" and channel 1 is "F" in every
                    // chans[p] (idmrg_channels_at builds them in that
                    // order), so these indices are the same for both
                    // environments -- only their roles swap.
                    HL = idmrg_subtract_energy_baseline(HL,HL_mpo,1,0,eshift);
                    HR = idmrg_subtract_energy_baseline(HR,HR_mpo,0,1,eshift);
                    }

                snap_U[p_L] = U; snap_U_left[p_L] = left_ket_old; snap_U_right[p_L] = new_bond_u;
                }

            // Each micro-step subtracted `eshift` from both environments,
            // i.e. 2*n_uc*eshift per macro-iteration, so add it back to
            // recover the true density. Exact as long as `eshift` was
            // constant across the two macro-iterations being differenced,
            // which is why it is latched once and never revised (below).
            bool have_density = have_prev_energy;
            double density = have_density ? (energy-prev_energy)/(2.0*n_uc)+eshift : 0.0;
            if (verbose_)
                println("idmrg macro-iter ",macro_iter,": E=",energy,
                        " density=",(have_density?density:0.0));
            // Demand-driven schedule: arm the noise whenever the state is
            // still a product state, and let it run for a few further
            // macro-iterations after that clears, so it does not switch off
            // the instant the basis opens and let the state fall straight
            // back. A healthy (already entangled) model never arms it at
            // all, and so is bit-for-bit unaffected.
            if (product_this_iter) noise_left = idmrg_noise_tail_;
            else if (noise_left>0) --noise_left;

            if (have_density && have_prev_density && noise_now==0.0 &&
                std::abs(density-prev_density)<etol)
                converged = true;
            prev_energy = energy; have_prev_energy = true;
            if (have_density)
                {
                prev_density = density; have_prev_density = true;
                report_density = density; have_report_density = true;
                }
            if (converged) break;
            if (eshift == 0.0 && have_density)
                {
                // Latch the baseline once, as soon as there is an estimate
                // to use. prev_energy is dropped in the same breath: it was
                // measured against un-shifted environments, so differencing
                // the next iteration's energy against it would report the
                // transient rather than a density. One macro-iteration is
                // given up doing so.
                eshift = density;
                have_prev_energy = false; have_prev_density = false;
                prev_energy = 0.0; prev_density = 0.0;
                }
            }

        IdmrgResult out;
        out.density = have_report_density ? report_density : prev_density;
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

        // The gauge-consistent 2-site unit cell every static observable on
        // this backend tiles (idmrg_onsite_expectation/
        // idmrg_two_point_correlator below, and
        // td_dynamical_correlator_window's own connected background) --
        // extracted once here, from the last executed macro-iteration's own
        // first micro-step, exactly as pyitensor/idmrg.py does. The cell is
        // always 2 sites: for n_uc=1 that is two repeats of the one-site
        // cell, for n_uc=2 exactly one unit cell (theta's own two sites are
        // sublattices p_L=0 and p_R=1 at micro-step 0). Either way cell
        // position k carries sublattice k%n_uc.
        have_idmrg_cell_ = false;
        ic_cache_valid_ = false;
        ic_Es_.clear(); ic_chis_.clear(); ic_rho_after_.clear(); ic_l_before_.clear();
        idmrg_cell_.clear(); idmrg_cell_l_.clear();
        idmrg_cell_d_.clear(); idmrg_cell_r_.clear();
        have_idmrg_cell_raw_ = false;
        iw_cache_valid_ = false;
        iw_Es_.clear(); iw_chis_.clear(); iw_rho_after_.clear(); iw_l_before_.clear();
        idmrg_cell_raw_.clear(); idmrg_cell_raw_l_.clear();
        idmrg_cell_raw_d_.clear(); idmrg_cell_raw_r_.clear();
        if (have_cell_seed)
            {
            std::vector<Cplx> a1, a2;
            if (idmrg_theta_cell(seed_U,seed_Ul,seed_Ud,seed_Ur,seed_S,
                                  seed_V,seed_Vl,seed_Vd,seed_Vr,
                                  seed_lam_outer,a1,a2))
                {
                idmrg_cell_ = {std::move(a1),std::move(a2)};
                idmrg_cell_l_ = {seed_Ul,seed_Vl};
                idmrg_cell_d_ = {seed_Ud,seed_Vd};
                idmrg_cell_r_ = {seed_Ur,seed_Vr};
                // Keep the raw cell BEFORE ic_canonicalize_cell mutates it
                // in place -- the window tiles this one, see
                // idmrg_cell_raw_'s own declaration.
                idmrg_cell_raw_ = idmrg_cell_;
                idmrg_cell_raw_l_ = idmrg_cell_l_;
                idmrg_cell_raw_d_ = idmrg_cell_d_;
                idmrg_cell_raw_r_ = idmrg_cell_r_;
                have_idmrg_cell_raw_ = true;
                ic_canonicalize_cell(idmrg_cell_,idmrg_cell_l_,
                                      idmrg_cell_d_,idmrg_cell_r_);
                have_idmrg_cell_ = true;
                }
            }

        return out;
        }

    // <opname> at sublattice p (0..n_uc-1) of a converged
    // idmrg_ground_state()'s infinite chain -- C++ port of
    // pyitensor/idmrg.py's own onsite_expectation, tiling the
    // gauge-consistent unit cell idmrg_theta_cell extracted (see that
    // function's docstring for why the per-micro-step idmrg_U_ chain
    // cannot be tiled instead: doing so put <Sz> at -0.13 against an exact
    // 0 on the XX chain).
    Cplx
    idmrg_onsite_expectation(std::string const& opname, int p) const
        {
        ic_require_cell("idmrg_onsite_expectation");
        int n_uc = idmrg_n_uc_;
        if (p<0 || p>=n_uc)
            throw ITError("Chain::idmrg_onsite_expectation: p must be in 0.."+
                           std::to_string(n_uc-1)+" (n_uc-1), got "+std::to_string(p));
        ic_build_cache();
        auto const& chis = ic_chis_;
        auto const& A = idmrg_cell_[p];
        int cl = idmrg_cell_l_[p], cd = idmrg_cell_d_[p], cr = idmrg_cell_r_[p];
        auto x_num = ic_apply_site_transfer(A,cl,cd,cr,
                                             idmrg_op_dense(p%n_uc,opname),
                                             ic_rho_after_[p]);
        auto x_den = ic_apply_site_transfer(A,cl,cd,cr,{},ic_rho_after_[p]);
        return ic_close_expectation(ic_l_before_[p%(int)ic_Es_.size()],
                                     chis[p]*chis[p],x_num,x_den);
        }

    // <opname_i(site p_i) opname_j(site p_i + r)> of a converged
    // idmrg_ground_state()'s infinite chain, r measured in physical sites
    // (r>=0) -- C++ port of pyitensor/idmrg.py's own two_point_correlator,
    // same r=0 same-site convention (M_j @ M_i, see
    // idmrg_correlator_endpoints) and the same Jordan-Wigner treatment for
    // fermionic operators (the string is threaded across every site
    // strictly between the two endpoints; a pair of odd total fermion
    // parity is rejected, since its string can never close).
    //
    // Positions advance through the *tiled 2-site cell*, while the site
    // type at each position is set by its sublattice (position % n_uc) --
    // the two coincide unless the cell is longer than the unit cell, which
    // is exactly the n_uc=1 case.
    //
    // One deliberate departure from a literal transcription of
    // pyitensor/idmrg.py: that version *composes* the per-site transfer
    // matrices along the string into one (chi^2, chi^2) product and only
    // then applies it to the right fixed point, which costs O(chi^6) per
    // intervening site (a chi^2-by-chi^2 matrix product). Here the same
    // chain is instead applied to the fixed point one site at a time,
    // right to left -- O(chi^4) per site. This is a pure re-association of
    // the same matrix chain, not a different formula, so it returns
    // exactly the same number; it just stops the cost of a correlator
    // sweep from being dominated by matrix-matrix products. Measured on a
    // maxm=24 chain: ~4.1s for a single r=3 correlator before, ~0.01s
    // after.
    Cplx
    idmrg_two_point_correlator(std::string const& opname_i, int p_i,
                                std::string const& opname_j, int r) const
        {
        ic_require_cell("idmrg_two_point_correlator");
        if (r<0)
            throw ITError("Chain::idmrg_two_point_correlator: r must be >= 0");
        int n_uc = idmrg_n_uc_;
        if (p_i<0 || p_i>=n_uc)
            throw ITError("Chain::idmrg_two_point_correlator: p_i must be in 0.."+
                           std::to_string(n_uc-1)+" (n_uc-1), got "+std::to_string(p_i));
        int n_cell = (int)idmrg_cell_.size();
        ic_build_cache();
        auto const& Es = ic_Es_;
        auto const& chis = ic_chis_;
        auto const& rho_after = ic_rho_after_;

        std::vector<Cplx> mat_i, mat_j;
        bool strung=false;
        idmrg_correlator_endpoints(opname_i,p_i,opname_j,(p_i+r)%n_uc,r==0,
                                    mat_i,mat_j,strung);

        auto cell_transfer = [&](int k, std::vector<Cplx> const& M)
            { return ic_transfer(idmrg_cell_[k],idmrg_cell_l_[k],idmrg_cell_d_[k],
                                  idmrg_cell_r_[k],M); };

        if (r==0)
            {
            auto const& A = idmrg_cell_[p_i];
            auto x_num = ic_apply_site_transfer(A,idmrg_cell_l_[p_i],idmrg_cell_d_[p_i],
                                                 idmrg_cell_r_[p_i],mat_i,rho_after[p_i]);
            auto x_den = ic_apply_site_transfer(A,idmrg_cell_l_[p_i],idmrg_cell_d_[p_i],
                                                 idmrg_cell_r_[p_i],{},rho_after[p_i]);
            return ic_close_expectation(ic_l_before_[p_i%n_cell],
                                         chis[p_i]*chis[p_i],x_num,x_den);
            }

        int k_j = (p_i+r)%n_cell;
        std::vector<Cplx> x_num = rho_after[k_j], x_den = rho_after[k_j];
        for (int step=r; step>=0; --step)
            {
            int k = (p_i+step)%n_cell;
            std::vector<Cplx> M;
            if (step==r)      M = mat_j;
            else if (step==0) M = mat_i;
            // An open Jordan-Wigner string puts F on each intervening
            // site; the denominator (x_den) stays the plain transfer
            // throughout either way -- it is the norm, and carries no
            // operator content at all.
            else if (strung)  M = idmrg_op_dense(k%n_uc,"F");
            auto const& A = idmrg_cell_[k];
            int cl = idmrg_cell_l_[k], cd = idmrg_cell_d_[k], cr = idmrg_cell_r_[k];
            x_num = ic_apply_site_transfer(A,cl,cd,cr,M,x_num);
            x_den = ic_apply_site_transfer(A,cl,cd,cr,{},x_den);
            }
        return ic_close_expectation(ic_l_before_[p_i % n_cell],
                                     chis[p_i]*chis[p_i],x_num,x_den);
        }

    // The "local superblock gap": re-diagonalize the *same* 2-site
    // effective Hamiltonian the growing algorithm already solved for its
    // ground state at the very last micro-step, but for its second-lowest
    // eigenvalue instead, and return the difference -- C++ port of
    // pyitensor/idmrg.py's own local_excitation_gap. See that function's
    // own docstring for the full rationale (in particular why the
    // orthogonality constraint is enforced here as an exact projector
    // P = I - |psi0><psi0| rather than via finite-DMRG's soft
    // overlap-penalty weight) and for its accuracy caveat: this has no
    // momentum label and reuses HL/HR exactly as they converged for the
    // *ground* state, so it is a cheap cross-check, not a variationally
    // optimal excited state.
    //
    // The stored superblock energy carries idmrg_ground_state's own
    // per-site baseline shift, but this returns a *difference* of two
    // eigenvalues of the same shifted operator, so the shift cancels
    // exactly.
    // (gap, e0 as re-solved here, e0 as the growth loop reported it) --
    // see idmrg_local_excitation_gap_detail's own comment below for why
    // the last two are worth handing back separately.
    struct LocalGap { double gap, e0_fresh, e0_stored; };

    double
    idmrg_local_excitation_gap(int niter) const
        { return idmrg_local_excitation_gap_detail(niter).gap; }

    // Same calculation, also returning the ground eigenvalue this method
    // re-solved for and the one idmrg_ground_state's own final micro-step
    // reported. They agree whenever the growth loop actually found its own
    // local effective Hamiltonian's ground state, which is the ordinary
    // case; a large disagreement means it did not -- the growth loop
    // warm-starts every micro-step from the previous one and so stays on
    // whatever branch it started on, while this method also searches from
    // a random start, and with ConserveQNs=false (mpscpp3/get_sites.h)
    // nothing confines that search to the converged state's own particle-
    // number sector. infinitechain.py surfaces the disagreement as a
    // Python warning rather than swallowing it, since in that situation
    // the returned gap is a genuine spectral gap of the stored operator
    // but no longer a gap "above the state whose observables vev/
    // correlator report".
    LocalGap
    idmrg_local_excitation_gap_detail(int niter) const
        {
        if (!have_idmrg_superblock_)
            throw ITError("Chain::idmrg_local_excitation_gap: called before "
                           "idmrg_ground_state (no stored local superblock)");
        Index phys_L = idmrg_sb_phys_L_, phys_R = idmrg_sb_phys_R_;
        Index HL_ket = idmrg_sb_HL_ket_, HR_ket = idmrg_sb_HR_ket_;
        Index HL_bra = idmrg_sb_HL_bra_, HR_bra = idmrg_sb_HR_bra_;
        bool have_HL_ket = idmrg_sb_have_HL_ket_, have_HR_ket = idmrg_sb_have_HR_ket_;
        ITensor const& HL = idmrg_sb_HL_;
        ITensor const& HR = idmrg_sb_HR_;

        size_t dim_in = (size_t)dim(phys_L)*dim(phys_R);
        if (have_HL_ket) dim_in *= (size_t)dim(HL_ket);
        if (have_HR_ket) dim_in *= (size_t)dim(HR_ket);
        if (dim_in < 2)
            throw ITError("Chain::idmrg_local_excitation_gap: the local Hilbert space "
                           "has dimension "+std::to_string(dim_in)+" -- too small to hold "
                           "a state orthogonal to the ground state");

        // Same effective Hamiltonian idmrg_local_solve built, rebuilt here
        // verbatim from the stored pieces.
        auto matvec = [&](ITensor v) -> ITensor
            {
            if (HL) v *= HL;
            v *= idmrg_sb_W_pL_;
            v *= idmrg_sb_W_pR_;
            if (HR) v *= HR;
            v.noPrime();
            if (have_HL_ket) v = replaceInds(v,{HL_bra},{HL_ket});
            if (have_HR_ket) v = replaceInds(v,{HR_bra},{HR_ket});
            return v;
            };

        ITensor psi0 = idmrg_sb_theta_;
        double n0 = norm(psi0);
        if (n0 == 0.0)
            throw ITError("Chain::idmrg_local_excitation_gap: stored local ground "
                           "vector has zero norm");
        psi0 /= n0;

        int kd = std::max(niter,2);
        kd = std::min<int>(kd,(int)std::min<size_t>(dim_in,200));

        // BOTH eigenvalues below come from a fresh solve of this same
        // stored operator, at the same solver strength -- the ground one
        // is NOT taken from idmrg_sb_energy_ (the growth loop's own
        // value), even though that value is by construction the Ritz value
        // of idmrg_sb_theta_ under exactly this operator. The two need not
        // be solved equally hard: the growth loop runs at
        // krylovdim=Infinite_Many_Body_Chain.niter (30 by default) with 2
        // restarts, warm-started from McCulloch's wavefunction prediction,
        // while this method runs at krylovdim=niter (200 by default) with
        // 4 restarts from a random start. Arnoldi's own convergence test
        // bounds the distance to *some* eigenvalue, not to the lowest one,
        // so the growth loop can stop on a low-lying excited eigenvalue
        // that this method's stronger solve would step past -- and mixing
        // the two is then not a difference of two eigenvalues of one
        // operator at all. (This is a separate matter from the shifted
        // deflation below, which is what actually caused the reported
        // negative gaps; solving both here is what makes the returned
        // number a well-defined function of the stored superblock alone,
        // independent of how hard the growth loop happened to solve, and
        // it is what idmrg_local_excitation_gap_detail reports back so
        // infinitechain.py can warn when the two disagree.)
        //
        // The warm start is still tried (it is nearly free -- a Krylov
        // space grown from a near-exact eigenvector hits happy breakdown
        // almost immediately) and kept only if it wins, since it converges
        // to full accuracy in the ordinary case where the growth loop did
        // find the local ground state.
        auto ground_solve = [&](ITensor const& start)
            {
            return arnoldi_smallest_real(matvec,start,kd,4,Sel::SR,Cplx(0,0),1e-10);
            };
        auto best = ground_solve(psi0);
        auto from_random = ground_solve(randomITensorC(psi0.inds()));
        if (from_random.first.real() < best.first.real()) best = from_random;
        double e0 = best.first.real();
        psi0 = best.second;
        psi0 /= norm(psi0);

        // The excited solve looks for the lowest eigenvalue of this same
        // operator restricted to psi0's orthogonal complement. The obvious
        // way to write that -- the bare projector P = I - |psi0><psi0| and
        // the operator P H P -- is WRONG here, and was the actual cause of
        // the large negative gaps this method used to return. P H P is
        // Hermitian and does agree with H on the complement, but it also
        // carries psi0 itself as an eigenvector with eigenvalue exactly
        // ZERO. That extra eigenvalue is harmless only while the rest of
        // the spectrum lies below it; the moment the operator's own ground
        // eigenvalue is POSITIVE, 0 becomes the smallest eigenvalue of
        // P H P and is exactly what a smallest-eigenvalue solver is asked
        // to find. Deflation keeps psi0 out of the Krylov space in exact
        // arithmetic only -- rounding regrows the component every
        // iteration, and a restarted solver targeting the bottom of the
        // spectrum then locks onto it. The returned "gap" is 0 - e0, i.e.
        // exactly minus the stored superblock energy: confirmed on a
        // spinful (c,f) Kondo chain at maxm=16, where this method returned
        // -358.003 meV against a stored energy of +0.358002587681, digit
        // for digit, and reproduced in isolation on a random Hermitian
        // matrix shifted to a positive ground energy (the deflated solve
        // returns 0.0 with |<psi0|v>|^2 = 1, i.e. psi0 itself).
        //
        // Whether the stored superblock energy is positive is not under
        // anyone's control: idmrg_subtract_energy_baseline leaves it as a
        // small residual boundary term of either sign, which is why this
        // struck one model at maxm=16 and not the same model at maxm=8
        // (energy -1.05, negative, correct answer).
        //
        // The fix is to SHIFT rather than merely project: send psi0's own
        // eigenvalue to sigma, far above the bottom of the spectrum, so it
        // can never be the answer, while leaving the complement untouched.
        // A = P H P + sigma |psi0><psi0| is still Hermitian, still agrees
        // with H on the complement, and its smallest eigenvalue now IS the
        // first excited state.
        double sigma = e0 + 10.0*(1.0+std::abs(e0));
        auto deflate = [&](ITensor v) -> ITensor
            { return v - psi0*eltC(dag(psi0)*v); };
        auto deflated = [&](ITensor v) -> ITensor
            {
            Cplx c = eltC(dag(psi0)*v);
            ITensor out = deflate(matvec(v - psi0*c));
            return out + (sigma*c)*psi0;
            };

        // Two independent safety nets on top of the shift, both no-ops in
        // the ordinary case. (1) If the solve still comes back sitting on
        // psi0, sigma was not far enough above the spectrum -- push it up
        // and retry. (2) If it comes back BELOW the ground eigenvalue
        // accepted above, then that ground solve, not this one, was the
        // one that missed: promote the lower state and re-deflate against
        // it. Each promotion strictly lowers e0, so this terminates; the
        // spurious psi0 mode can no longer be promoted, since the shift
        // put its eigenvalue at sigma.
        double e1 = 0.0;
        for (int round=0; round<4; ++round)
            {
            auto result = arnoldi_smallest_real(deflated,
                    deflate(randomITensorC(psi0.inds())),kd,4,Sel::SR,Cplx(0,0),1e-10);
            e1 = result.first.real();
            ITensor v1 = result.second;
            double nv = norm(v1);
            // A zero-norm eigenvector means the solver gave up without
            // producing one; there is nothing to promote or test against,
            // so keep the value it did return rather than deflating
            // against a zero vector on the next round.
            if (nv <= 0.0) break;
            v1 /= nv;
            if (std::norm(eltC(dag(psi0)*v1)) > 0.5)
                { sigma = e0 + 100.0*(sigma-e0); continue; }
            if (e1 >= e0 - 1e-10*(1.0+std::abs(e0))) break;
            e0 = e1;
            psi0 = v1;
            sigma = e0 + 10.0*(1.0+std::abs(e0));
            }
        // e1 and e0 are eigenvalues of the same Hermitian operator, the
        // second one restricted to the first's orthogonal complement, so
        // the difference is non-negative by construction -- but only up to
        // the solver's own precision when the two are the same eigenvalue
        // of a degenerate level. Flatten a rounding-level negative to
        // exactly zero; anything genuinely negative is left visible.
        double gap = e1-e0;
        if (gap<0.0 && gap>-1e-10*(1.0+std::abs(e0))) gap = 0.0;
        return {gap,e0,idmrg_sb_energy_};
        }

    // <opname> at site p of a converged multi-site (n_uc>2) VUMPS state.
    // AC[p] is already the exactly normalized single-site reduced state by
    // construction of the mixed canonical gauge, so this is a plain
    // contraction with no eigenproblem -- same formula the grouped
    // vumps_onsite_expectation uses, with the relevant AC picked out of the
    // cell instead of embedded into a supersite operator.
    Cplx
    vms_onsite_expectation(std::string const& opname, int p) const
        {
        if (!have_vms_snapshot_)
            throw ITError("Chain::vms_onsite_expectation: called before a "
                           "multi-site vumps_ground_state");
        int n_uc = (int)vms_rows_.size();
        if (p<0 || p>=n_uc)
            throw ITError("Chain::vms_onsite_expectation: p out of range");
        int D = vms_D_, d = vms_rows_[p].d;
        auto M = idmrg_op_dense(p,opname);
        auto const& AC = vms_AC_[p];
        auto AC_op = vx_apply_op_ket(M,AC,D,d);
        Cplx num(0,0), den(0,0);
        for (size_t k=0;k<AC.size();++k)
            {
            num += std::conj(AC[k])*AC_op[k];
            den += std::conj(AC[k])*AC[k];
            }
        return num/den;
        }

    // <opname_i(p_i) opname_j(p_i+r)> of a converged multi-site VUMPS
    // state, r in physical sites. Puts the mixed-gauge centre at p_i (so
    // AC[p_i] carries the exact normalization and everything to its left is
    // already the identity), then walks right through AR, which is exactly
    // right-orthonormal so the tail beyond the second operator closes to
    // the identity with no eigenproblem either. Jordan-Wigner strings are
    // threaded exactly as idmrg_two_point_correlator does.
    Cplx
    vms_two_point_correlator(std::string const& opname_i, int p_i,
                              std::string const& opname_j, int r) const
        {
        if (!have_vms_snapshot_)
            throw ITError("Chain::vms_two_point_correlator: called before a "
                           "multi-site vumps_ground_state");
        int n_uc = (int)vms_rows_.size();
        if (p_i<0 || p_i>=n_uc)
            throw ITError("Chain::vms_two_point_correlator: p_i out of range");
        if (r<0)
            throw ITError("Chain::vms_two_point_correlator: r must be >= 0");
        int D = vms_D_;
        bool strung=false;
        std::vector<Cplx> mat_i, mat_j;
        int p_j = (p_i+r) % n_uc;
        idmrg_correlator_endpoints(opname_i,p_i,opname_j,p_j,/*same_site=*/r==0,
                                    mat_i,mat_j,strung);
        if (r==0)
            {
            // idmrg_correlator_endpoints already composed the two matrices
            // into mat_i for the same-site case (its own r=0 convention).
            auto M = mat_i;
            int d = vms_rows_[p_i].d;
            auto const& AC = vms_AC_[p_i];
            auto AC_op = vx_apply_op_ket(M,AC,D,d);
            Cplx num(0,0), den(0,0);
            for (size_t k=0;k<AC.size();++k)
                {
                num += std::conj(AC[k])*AC_op[k];
                den += std::conj(AC[k])*AC[k];
                }
            return num/den;
            }
        auto const& AC = vms_AC_[p_i];
        int d_i = vms_rows_[p_i].d;
        Cplx den(0,0);
        for (size_t k=0;k<AC.size();++k) den += std::conj(AC[k])*AC[k];
        auto X = vms_push_left(vx_eye(D),AC,D,d_i,mat_i,true);
        for (int step=1; step<=r; ++step)
            {
            int pos = (p_i+step) % n_uc;
            int d = vms_rows_[pos].d;
            if (step==r)
                X = vms_push_left(X,vms_AR_[pos],D,d,mat_j,true);
            else if (strung)
                X = vms_push_left(X,vms_AR_[pos],D,d,idmrg_op_dense(pos,"F"),true);
            else
                X = vms_push_left(X,vms_AR_[pos],D,d,{},false);
            }
        Cplx tr(0,0);
        for (int i=0;i<D;++i) tr += X[i*D+i];
        return tr/den;
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
                        int D, double tol, int maxiter, int nrestarts,
                        int niter_lanczos = vx_lanczos_niter_)
        {
        reject_sector("vumps_ground_state"); // dense uniform tensors below
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
        // n_uc > 2 goes to the sequential multi-site algorithm, which never
        // groups the cell and so costs linearly rather than exponentially
        // in n_uc -- see the "Sequential multi-site VUMPS" block above and
        // pyitensor/vumps_ms.py's own module docstring. n_uc <= 2 stays on
        // the grouped path purely so its already-validated values do not
        // move; the two agree to machine precision where both apply.
        if (n_uc>2)
            {
            std::mt19937_64 rng_ms(std::random_device{}());
            auto run = vms_ground_state(rows,D,tol,maxiter,niter_lanczos,
                                         nrestarts,rng_ms);
            vms_AL_ = run.AL; vms_AR_ = run.AR; vms_C_ = run.C; vms_AC_ = run.AC;
            vms_rows_ = rows; vms_D_ = D;
            have_vms_snapshot_ = true;
            have_vumps_snapshot_ = false;   // the grouped consumers cannot read this
            VumpsResult out_ms;
            out_ms.e0 = run.e_cell / n_uc;
            out_ms.converged = run.converged;
            out_ms.niter_done = run.niter;
            out_ms.gauge_mismatch = run.mismatch;
            return out_ms;
            }

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
        // The bond dimensions actually solved at on the way to D: 1,2,4,
        // 8,...,D (doubling, always ending exactly at D) -- pyitensor/
        // vumps.py's own _d_ramp, same sequence. Every step is a full
        // multi-restart VUMPS solve costing ~O(D^3), so ramping one
        // integer at a time makes the driver O(D^4) and the ramp, not the
        // solve at the target D, dominates (~100 complete solves to reach
        // D=30, versus 6 here). What the ramp is for is unchanged: each
        // step warm-starts from the previous one via vumps_grow_init,
        // which embeds the old solution in the new tensors' leading block
        // plus noise and assumes nothing about how big the jump is.
        std::vector<int> d_ramp;
        for (int d_cur=1; d_cur<D; d_cur*=2) d_ramp.push_back(d_cur);
        d_ramp.push_back(D);
        for (int D_cur : d_ramp)
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
                                               niter_lanczos,has_init?&init:nullptr,rng);
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

    // Endpoint operator matrices, ordering coefficient, and "is a
    // Jordan-Wigner string open between them" for a two-point correlator
    // <opname_i(p_i) opname_j(p_j)>, p_i <= p_j -- the two-factor case of
    // pyitensor/idmrg.py's _term_site_matrices, and deliberately the same
    // rules idmrg_classify_terms applies to a Hamiltonian's own 2-site
    // terms, so a correlator and a Hamiltonian term written with the same
    // operator names cannot disagree about the convention (which endpoint
    // carries the extra F, and whether the sites in between carry one).
    // Throws for a pair of odd total fermion parity: its string can never
    // close, so the correlator is not a well-defined object here.
    // p_i/p_j are SUBLATTICE indices (0..n_uc-1, what idmrg_op_dense and
    // sites_.si expect); `same_site` says whether the two operators are on
    // one and the same physical site (p_i==p_j AND no whole cells between
    // them), which sublattice indices alone cannot express.
    void
    idmrg_correlator_endpoints(std::string const& opname_i, int p_i,
                                std::string const& opname_j, int p_j,
                                bool same_site,
                                std::vector<Cplx>& mat_i,
                                std::vector<Cplx>& mat_j,
                                bool& strung) const
        {
        bool fi = idmrg_is_fermionic(opname_i), fj = idmrg_is_fermionic(opname_j);
        if (same_site)
            {
            int d = dim(sites_.si(p_i+1));
            mat_i = idmrg_matmul(idmrg_op_dense(p_i,opname_j),
                                  idmrg_op_dense(p_i,opname_i),d);
            bool site_ferm = (fi != fj);   // odd number of fermionic factors
            if (site_ferm)
                throw ITError("Chain: two_point_correlator: the operator pair has "
                               "odd total fermion parity -- its Jordan-Wigner string "
                               "cannot be closed");
            mat_j.clear();
            strung = false;
            return;
            }
        if (fi != fj)
            throw ITError("Chain: two_point_correlator: the operator pair has "
                           "odd total fermion parity -- its Jordan-Wigner string "
                           "cannot be closed");
        int di = dim(sites_.si(p_i+1));
        mat_i = idmrg_op_dense(p_i,opname_i);
        mat_j = idmrg_op_dense(p_j,opname_j);
        if (fi) // carry enters false, so a fermionic first endpoint picks up F
            mat_i = idmrg_matmul(idmrg_op_dense(p_i,"F"),mat_i,di);
        // the second endpoint's own parity matches the carry, so it gets none
        strung = fi;
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

        std::vector<Cplx> mat_i, mat_j;
        bool strung = false;
        idmrg_correlator_endpoints(opname_i,p_i,opname_j,p_j,
                                    (cell_offset==0 && p_j==p_i),
                                    mat_i,mat_j,strung);
        // {p: F} for lo <= p < hi -- the piece of the string inside one cell
        auto string_ops = [&](int lo, int hi)
            {
            std::map<int,std::vector<Cplx>> out;
            if (strung)
                for (int p=lo;p<hi;++p) out[p] = idmrg_op_dense(p,"F");
            return out;
            };

        if (cell_offset==0)
            {
            std::vector<Cplx> M;
            if (p_j==p_i)
                {
                M = vumps_embed_group_operator({{p_i, mat_i}});
                }
            else
                {
                auto ops = string_ops(p_i+1,p_j);   // string strictly between
                ops[p_i] = mat_i;
                ops[p_j] = mat_j;
                M = vumps_embed_group_operator(ops);
                }
            auto AC_op = vx_apply_op_ket(M,AC,D,d_g);
            return vx_dot_conj(AC,AC_op)/norm;
            }

        // The string, when open, runs from just right of p_i to the end of
        // the AC cell, across every fully-crossed cell, and from the start
        // of the final cell up to just left of p_j.
        auto ops_i = string_ops(p_i+1,n_uc);
        ops_i[p_i] = mat_i;
        auto Mi_embed = vumps_embed_group_operator(ops_i);
        auto AC_op = vx_apply_op_ket(Mi_embed,AC,D,d_g);
        // Open right-bond object: bra/ket both AC, operator on the ket side
        // only -- this already IS the full left closure (AC's own left leg
        // is summed away here), leaving just the (ket-bond, bra-bond) legs
        // open.
        auto X = vx_left_close(AC_op,AC,D,d_g);
        for (auto& x : X) x /= norm;

        if (cell_offset>1)
            {
            // A fully-crossed cell carries F on ALL of its sub-sites when
            // the string is open, and the plain transfer otherwise.
            auto AR_cross = vumps_AR_;
            if (strung)
                AR_cross = vx_apply_op_ket(
                    vumps_embed_group_operator(string_ops(0,n_uc)),vumps_AR_,D,d_g);
            auto E_AR = vx_op_transfer_matrix(AR_cross,D,d_g,vumps_AR_,false,{});
            for (int i=0;i<cell_offset-1;++i)
                X = vx_apply_transfer_from_left(E_AR,D,X);
            }

        auto ops_j = string_ops(0,p_j);              // string up to p_j
        ops_j[p_j] = mat_j;
        auto Mj_embed = vumps_embed_group_operator(ops_j);
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
        if (idmrg_is_fermionic(opname_A)!=idmrg_is_fermionic(opname_B))
            throw std::invalid_argument(
                  "Chain::td_dynamical_correlator_window: the operator pair ("
                  +opname_A+", "+opname_B+") has odd total fermion parity -- "
                  "its Jordan-Wigner string can never close, so <A_x(t) B_0> "
                  "is not a well-defined object on an infinite chain");
        int n_uc = idmrg_n_uc_;
        if (!(0<=p_i && p_i<n_uc))
            throw ITError("Chain::td_dynamical_correlator_window: p_i out of range");
        if (n_window<1)
            throw ITError("Chain::td_dynamical_correlator_window: n_window must be >= 1");

        auto win = idmrg_build_window(n_window);
        int n = win.n;
        // win.n_window_uc, not the requested n_window: the tiling unit is
        // the 2-site cell, so the realized window can be one unit cell
        // longer than asked for (see IdmrgWindow::n_window_uc).
        int center = idmrg_window_center(win.n_window_uc,n_uc,p_i);

        // ITensor's own environment-propagation convention (see this
        // method's own top comment) -- convert once, outside the time
        // loop. Computed *before* perturbing win.psi below, since eshift
        // (just underneath) needs it evaluated on the unperturbed ground
        // window too.
        ITensor LH = idmrg_relabel_bra_to_prime_ket(idmrg_HL_,idmrg_HL_bra_,idmrg_HL_ket_);
        ITensor RH = idmrg_relabel_bra_to_prime_ket(idmrg_HR_,idmrg_HR_bra_,idmrg_HR_ket_);

        // NO eshift/global-phase correction here, unlike the pre-2026-08-29
        // version of this method and unlike an ordinary finite-chain
        // quench. The window's own effective H is built from idmrg_HL_/
        // idmrg_HR_, which are NOT energy-baseline-subtracted, so evolving
        // under it multiplies the state by a large, run-dependent
        // exp(-i*E_win*t). This used to be undone by measuring E_win once
        // (a null tdvp() step on a throwaway copy, whose Rayleigh quotient
        // traces the window's two dangling boundary legs) and applying
        // exp(+i*E_win*t) to every snapshot. That was only ever consistent
        // while the snapshot itself closed those same legs by a bare
        // trace; it does not survive closing them with the cell's own
        // transfer-matrix fixed points (idmrg_close_array_chain), which is
        // what the gauge fix requires -- the two closings genuinely see
        // different energies, exactly as pyitensor/idmrg_window.py's own
        // window_total_energy docstring (point 3) records. Measured after
        // the gauge fix but before this one: v3 tracked the python backend
        // to 2e-6 at t=0 and then drifted linearly, 5.2e-2 by t=0.15.
        //
        // The vacuum branch below removes the factor exactly instead of
        // estimating it: an unperturbed copy of the same window, evolved
        // with the identical eshift-free tdvp() calls and measured every
        // step through the *identical* contraction with the identity
        // operator, carries the same exp(-i*E_win*t); dividing cancels it
        // whatever E_win is and whichever closure is used. Direct port of
        // idmrg_window.py's own dynamical_correlator_td vacuum
        // normalization, including its shallow-copy trick: the two
        // branches need independent tensor lists but must keep the same
        // Index identities the environments LH/RH are keyed on, which an
        // MPS copy gives (every operation below returns a new ITensor
        // rather than mutating one).

        // Sec. V.1 step 3: perturb (B_0|psi>), in place -- keeping an
        // unperturbed copy first, for the vacuum branch (see above).
        IdmrgWindow win_ground = win;
        idmrg_window_apply_local_op(win,center,opname_B);
        // Establish a genuine orthogonality center now (via ITensor's own
        // robust QR-sweep-based position(), which -- unlike inner()/
        // innerC() below -- correctly handles win.psi's own two extra
        // dangling boundary legs at sites 1/n): TDVPWorker (TDVP/tdvp.h)
        // would otherwise do this on its own first call, but the
        // norm-save/restore fix just below needs it established *before*
        // that first call too (norm(MPS) requires isOrtho()).
        win.psi.position(1);
        IdmrgWindow win_vac = win_ground; // unperturbed, evolved alongside
        win_vac.psi.position(1);

        // The disconnected background <A><B>, from the same
        // gauge-consistent unit cell every other static observable on this
        // backend tiles -- this used to tile the raw per-micro-step
        // idmrg_U_ chain instead, i.e. it carried exactly the gauge error
        // idmrg_theta_cell exists to remove (see idmrg_onsite_expectation).
        // One deliberate divergence from pyitensor/idmrg_window.py here:
        // that side's _correlator_cell silently falls back to tiling
        // U_list when no cell was extracted, whereas
        // idmrg_onsite_expectation throws (ic_require_cell). Only a run so
        // short that no macro-iteration ever produced a cell seed can
        // reach it, and a clear error beats a silently gauge-wrong number.
        Cplx mean_B = Cplx(0,0);
        std::map<int,Cplx> background;
        if (connected)
            {
            if (idmrg_is_fermionic(opname_B))
                {
                // A parity-odd operator has <O> = 0 identically (its own
                // Jordan-Wigner string cannot close on a single site), so
                // the disconnected background is zero by symmetry rather
                // than by measurement -- and idmrg_onsite_expectation would
                // return a stringless artifact for it. Mirrors
                // idmrg_window.py's own dynamical_correlator_td.
                for (int x : x_values) background[x] = Cplx(0,0);
                }
            else
                {
                int p_B = (center-1)%n_uc;
                mean_B = idmrg_onsite_expectation(opname_B,p_B);
                for (int x : x_values)
                    {
                    int p_A = ((center+x-1)%n_uc+n_uc)%n_uc;
                    background[x] = mean_B*idmrg_onsite_expectation(opname_A,p_A);
                    }
                }
            }

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
            auto snap = idmrg_window_snapshot_correlator(win,opname_A,x_values,center);
            // <psi|psi(t)> through the identical contraction (same
            // closure, same calibration) -- divides out the spurious
            // exp(-i*E_win*t) both branches carry, see the block above.
            // Applied to the *raw* snapshot value, before background
            // subtraction, since `background` is computed from the static,
            // un-evolved ground state and never carries that factor.
            Cplx vac = idmrg_window_snapshot_correlator(
                            win_vac,"Id",std::vector<int>{0},center)[0];
            if (std::abs(vac) < 1e-300)
                throw std::runtime_error("Chain::td_dynamical_correlator_window: "
                      "the vacuum amplitude <psi|psi(t)> vanished -- the "
                      "window's own evolution has lost its norm entirely");
            for (size_t ix=0; ix<x_values.size(); ++ix)
                {
                Cplx val = snap[ix]/vac;
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
                // The vacuum branch, stepped identically -- same dt, same
                // sweeps/truncation, same norm save/restore -- so that
                // whatever global factor and whatever truncation drift the
                // perturbed branch picks up, this one picks up too.
                Cplx norm0_vac = Cplx(norm(win_vac.psi),0.0);
                tdvp(win_vac.psi,win_vac.mpo,t,LH,RH,sweeps,args);
                win_vac.psi.normalize();
                win_vac.psi *= norm0_vac;
                }
            }
        return out;
        }

    private:

    // -- conserved sector (QN mode) ------------------------------------
    //
    // Everything below is inert while sector mode is off (has_sector_ ==
    // false): the checks return immediately and mpo_from_terms()/
    // default_mps() take exactly the paths they always did.

    // dense_sites_'s physical indices in site order: the target
    // replaceSiteInds() rebases a de-QN'd MPS onto in to_dense_mps().
    IndexSet
    dense_site_inds() const
        {
        std::vector<Index> v;
        v.reserve(site_types_.size());
        for (size_t k=0;k<site_types_.size();++k) v.push_back(dense_sites_(int(k+1)));
        return IndexSet(std::move(v));
        }

    // Block-sparse (QN) MPS -> the identical wavefunction stored densely on
    // dense_sites_. See promote_to_dense() for why this is exact.
    //
    // The one assumption is that the QN-carrying and dense site indices
    // enumerate their local basis states in the same order, so that basis
    // state k means the same thing on both. That holds because both come
    // from the same ITensor site header with only ConserveQNs flipped
    // (spinhalf.h declares Up then Dn either way, electron.h Emp/Up/Dn/UpDn,
    // ...), and op_charge() right below already depends on exactly this
    // correspondence -- it labels dense_sites_'s matrix elements with
    // sites_'s quantum numbers.
    MPS
    to_dense_mps(MPS const& wf) const
        {
        return replaceSiteInds(removeQNs(wf),dense_site_inds());
        }

    // Drop every piece of session state that was built against the old
    // sites_ -- which, after a SiteSet swap, is all of it.
    void
    forget_everything_built_on_sites()
        {
        have_H_ = false;
        have_wf0_ = false;
        have_wf0_energy_ = false;
        have_bandwidth_min_ = false;
        have_bandwidth_max_ = false;
        have_idmrg_snapshot_ = false;
        have_idmrg_cell_ = false;
        ic_cache_valid_ = false;
        have_idmrg_cell_raw_ = false;
        iw_cache_valid_ = false;
        have_idmrg_superblock_ = false;
        have_vms_snapshot_ = false;
        have_vumps_snapshot_ = false;
        have_vumps_exc_env_ = false;
        }

    static std::string
    join_strings(std::vector<std::string> const& v)
        {
        std::string out;
        for (size_t k=0;k<v.size();++k) { if (k) out += ", "; out += v.at(k); }
        return out;
        }

    // A QN as "Nf=3, Sz=0", listing only its active components.
    static std::string
    qn_string(QN const& q)
        {
        std::string out;
        for (auto const& qv : q.store())
            {
            if (!isActive(qv)) continue;
            if (!out.empty()) out += ", ";
            out += tinyformat::format("%s=%d",std::string(qv.name().c_str()),qv.val());
            }
        return out.empty() ? std::string("(no charge)") : out;
        }

    std::string
    sector_string() const
        {
        std::string out;
        for (auto const& sq : sector_)
            {
            if (!out.empty()) out += ", ";
            out += tinyformat::format("%s=%d",sq.first,sq.second);
            }
        return out.empty() ? std::string("(no sector)") : out;
        }

    static std::string
    term_string(MOTerm const& t)
        {
        std::string out;
        for (auto const& f : t.factors)
            out += tinyformat::format("%s%s(%d)",out.empty()?"":"*",f.name,f.site);
        return out.empty() ? std::string("(empty term)") : out;
        }

    // True if every active component of q is zero. Spelled out rather than
    // compared against QN(), whose own operator== is not a "same charge"
    // test between differently-populated QNs.
    static bool
    qn_is_zero(QN const& q)
        {
        for (auto const& qv : q.store()) if (isActive(qv) && qv.val()!=0) return false;
        return true;
        }

    // The charge a named single-site operator carries, inferred from its
    // matrix elements on the *dense* site set plus the QN labels the sector
    // grading gives each local basis state. Returns false if the operator
    // has no definite charge at all -- Sx on Sz-conserving spin sites is
    // the standing example: its nonzero elements raise *and* lower.
    //
    // The reason this reads dense matrix elements instead of just calling
    // flux(sites_.op(name,i)): for exactly those non-definite operators,
    // building them on QN-carrying indices aborts the process from inside
    // ITensor::set ("cannot set element with flux different from ITensor
    // flux") rather than throwing anything catchable. Confirmed directly.
    bool
    op_charge(std::string const& name, int site, QN& out) const
        {
        auto O = dense_sites_.op(name,site);
        auto sd = dense_sites_(site);
        int d = dim(sd);
        bool first = true;
        QN f;
        for (int a=1;a<=d;a++)
        for (int b=1;b<=d;b++)
            {
            // Site operators are built as ITensor(dag(s),prime(s)), i.e.
            // unprimed = in, primed = out (see sites/*.h), so element
            // (a,b) maps basis state a to basis state b.
            if (std::abs(eltC(O,sd(a),prime(sd)(b))) < 1e-12) continue;
            QN delta = qn(sites_(site)(b)) - qn(sites_(site)(a));
            if (first) { f = delta; first = false; }
            else if (!qn_is_zero(delta-f)) return false;
            }
        out = f; // an identically-zero operator carries no charge at all
        return true;
        }

    // Sector-mode normalization of a term list: expand Sx/Sy into ladder
    // operators and cancel the strings that only conserve the sector as a
    // sum (the textbook Sx.Sx+Sy.Sy Heisenberg case -- see mo_terms.h's
    // own comment for why the cancellation has to happen here rather than
    // being left to AutoMPO), then reject whatever is left that still does
    // not conserve it, by name, before ITensor can abort over it.
    // Dense d x d matrix of a single-site operator, in op_charge()'s own
    // convention: element (a,b) is the amplitude taking basis state a to
    // basis state b.
    std::vector<Cplx>
    site_matrix(std::string const& name, int site) const
        {
        auto O = dense_sites_.op(name,site);
        auto sd = dense_sites_(site);
        int d = dim(sd);
        std::vector<Cplx> m(d*d,Cplx(0.,0.));
        for (int a=1;a<=d;a++)
        for (int b=1;b<=d;b++)
            m[(a-1)*d+(b-1)] = eltC(O,sd(a),prime(sd)(b));
        return m;
        }

    // Can a single site of this SiteSet carry a charge change of `delta`
    // at all? (i.e. is there a pair of basis states whose QNs differ by
    // it) A fermionic site holds 0 or 1 particle, so it can do -1, 0, +1
    // and nothing else.
    bool
    site_delta_representable(int site, QN const& delta) const
        {
        auto s = sites_(site);
        int d = dim(s);
        for (int a=1;a<=d;a++)
        for (int b=1;b<=d;b++)
            if (qn_is_zero(qn(s(b))-qn(s(a))-delta)) return true;
        return false;
        }

    std::vector<MOTerm>
    sector_terms(std::vector<MOTerm> const& terms) const
        {
        auto raw = combine_terms(expand_xy_terms(dense_sites_,terms));
        // Same-site factors, before the net-charge check below. That check
        // sums each term's charge over ALL its factors, so a term that is
        // net-neutral could still pile an impossible charge onto ONE site
        // -- Cdag(0) C(1) Cdag(0) C(1), which the four-point correlator
        // produces as a matter of course, is net zero but +2 on site 0.
        // AutoMPO then aborted the whole process ("Index does not contain
        // given QN block"), uncatchably, since ITensor's Error() calls
        // abort(). Such a product is usually identically zero by Pauli
        // exclusion (mode="ED" answers 0), so drop those terms outright
        // and raise a catchable error only for the ones that are not.
        std::vector<MOTerm> out;
        for (auto const& t : raw)
            {
            std::map<int,std::vector<std::string>> by_site;
            for (auto const& f : t.factors) by_site[f.site].push_back(f.name);
            bool vanishes = false;
            for (auto const& kv : by_site)
                {
                int site = kv.first;
                auto const& names = kv.second;
                if (names.size()<2) continue;
                int d = dim(dense_sites_(site));
                // composite of the same-site factors: written left to
                // right, the rightmost acts first, so with the (in,out)
                // element convention above the matrices multiply in
                // reverse order
                std::vector<Cplx> comp(d*d,Cplx(0.,0.));
                for (int a=0;a<d;a++) comp[a*d+a] = Cplx(1.,0.);
                for (auto it=names.rbegin();it!=names.rend();++it)
                    {
                    auto m = site_matrix(*it,site);
                    std::vector<Cplx> prod(d*d,Cplx(0.,0.));
                    for (int a=0;a<d;a++)
                    for (int c=0;c<d;c++)
                        {
                        Cplx acc(0.,0.);
                        for (int b=0;b<d;b++) acc += comp[a*d+b]*m[b*d+c];
                        prod[a*d+c] = acc;
                        }
                    comp = prod;
                    }
                double nrm = 0.;
                for (auto const& z : comp) nrm += std::norm(z);
                if (nrm < 1e-24) { vanishes = true; break; }
                // not zero: is the charge it accumulates on this site one
                // the site can actually carry?
                QN acc;
                bool ok = true;
                for (auto const& nm : names)
                    {
                    QN q;
                    if (!op_charge(nm,site,q)) { ok = false; break; }
                    acc = acc + q;
                    }
                if (ok && !site_delta_representable(site,acc))
                    throw std::invalid_argument(tinyformat::format(
                        "Conserved-sector mode (%s): the term %s applies "
                        "several operators to site %d that together change "
                        "its charge by %s, which a single site cannot carry. "
                        "(Without a sector this is built and evaluated "
                        "normally; here it would abort inside ITensor.)",
                        sector_string(),term_string(t),site,qn_string(acc)));
                }
            if (!vanishes) out.push_back(t);
            }
        for (auto const& t : out)
            {
            QN total;
            for (auto const& f : t.factors)
                {
                QN q;
                if (!op_charge(f.name,f.site,q))
                    throw std::invalid_argument(tinyformat::format("Conserved-sector mode (%s): operator \"%s\" "
                          "at site %d carries no definite charge, so no operator "
                          "containing it can be built while a sector is set. Write the "
                          "model with S+/S- (or Cdag/C) instead, or clear the sector.",
                          sector_string(),f.name,f.site));
                total = total + q;
                }
            if (!qn_is_zero(total))
                throw std::invalid_argument(tinyformat::format("Conserved-sector mode (%s): the term %s changes "
                      "the conserved quantum numbers by %s -- every operator built on a "
                      "chain with a sector set must conserve them. (Terms are named after "
                      "the exact Sx/Sy -> S+/S- expansion and the cancellation of strings "
                      "that only conserve the sector as a sum, so this may not be spelled "
                      "the way it was written.)",
                      sector_string(),term_string(t),qn_string(total)));
            }
        return out;
        }

    MPO
    mpo_from_terms(std::vector<MOTerm> const& terms) const
        {
        if (!has_sector_) return build_mpo(sites_,terms,mpomaxm_);
        return build_mpo(sites_,sector_terms(terms),mpomaxm_);
        }

    AutoMPO
    ampo_from_terms(std::vector<MOTerm> const& terms) const
        {
        if (!has_sector_) return build_ampo(sites_,terms);
        return build_ampo(sites_,sector_terms(terms));
        }

    // The charge basis state `st` of site `i` carries, as a plain integer
    // per conserved quantity (same order as sector_).
    std::vector<int>
    state_charge(int i, int st) const
        {
        auto q = qn(sites_(i)(st));
        std::vector<int> out;
        out.reserve(sector_.size());
        for (auto const& sq : sector_) out.push_back(q.val(sq.first));
        return out;
        }

    // One per-site basis-state assignment whose charges add up to exactly
    // the requested sector: both the reachability check
    // set_conserved_sector() runs (so an impossible target errors at
    // request time rather than mid-sweep) and the arrangement sector_mps()
    // starts DMRG from.
    //
    // A small dynamic program over reachable partial charges rather than a
    // greedy fill, because with more than one conserved quantity -- a
    // Hubbard chain at fixed Nf *and* Sz -- greedy choices dead-end. Each
    // layer keeps one representative per distinct partial charge, so its
    // size is bounded by the number of reachable charge tuples, not by the
    // number of assignments.
    std::vector<int>
    sector_state_plan() const
        {
        int N = sites_.length();
        int nq = int(sector_.size());
        // per-site reachable charge range, used to prune partial sums that
        // can no longer reach the target
        std::vector<std::vector<int>> smin(N+2,std::vector<int>(nq,0)),
                                      smax(N+2,std::vector<int>(nq,0));
        for (int i=N;i>=1;--i)
            {
            std::vector<int> lo(nq,0), hi(nq,0);
            for (int st=1;st<=dim(sites_(i));++st)
                {
                auto c = state_charge(i,st);
                for (int k=0;k<nq;k++)
                    {
                    if (st==1) { lo[k] = hi[k] = c[k]; }
                    else { lo[k] = std::min(lo[k],c[k]); hi[k] = std::max(hi[k],c[k]); }
                    }
                }
            for (int k=0;k<nq;k++)
                {
                smin[i][k] = lo[k] + smin[i+1][k];
                smax[i][k] = hi[k] + smax[i+1][k];
                }
            }
        std::vector<int> target;
        target.reserve(nq);
        for (auto const& sq : sector_) target.push_back(sq.second);

        struct Node { std::vector<int> charge; int prev; int state; };
        std::vector<std::vector<Node>> layers(N+1);
        layers[0].push_back(Node{std::vector<int>(nq,0),-1,0});
        for (int i=1;i<=N;i++)
            {
            std::map<std::vector<int>,int> seen;
            for (size_t pi=0;pi<layers[i-1].size();++pi)
                {
                auto const& node = layers[i-1].at(pi);
                for (int st=1;st<=dim(sites_(i));++st)
                    {
                    auto c = state_charge(i,st);
                    std::vector<int> partial(nq);
                    bool ok = true;
                    for (int k=0;k<nq;k++)
                        {
                        partial[k] = node.charge.at(k) + c.at(k);
                        // can the remaining sites still close the gap?
                        if (partial[k]+smin[i+1][k] > target[k]) { ok = false; break; }
                        if (partial[k]+smax[i+1][k] < target[k]) { ok = false; break; }
                        }
                    if (!ok) continue;
                    if (seen.count(partial)) continue; // one representative is enough
                    seen.emplace(partial,int(layers[i].size()));
                    layers[i].push_back(Node{partial,int(pi),st});
                    }
                }
            if (layers[i].empty())
                throw std::invalid_argument(tinyformat::format("Chain::set_conserved_sector: no state of this "
                      "chain has %s -- the requested sector is empty",sector_string()));
            }
        // the final layer was pruned to exactly the target
        std::vector<int> plan(N,1);
        int at = 0;
        for (int i=N;i>=1;--i)
            {
            auto const& node = layers[i].at(at);
            plan[i-1] = node.state;
            at = node.prev;
            }
        return plan;
        }

    // `k` in-sector product states: the same multiset of local basis
    // states as sector_state_plan()'s, shuffled among the sites that can
    // host them. Permuting within one site-type code keeps every state
    // legal on the site it lands on, and any permutation is charge-neutral
    // by construction, so all of these live in the requested sector.
    std::vector<std::vector<int>>
    sector_state_arrangements(int k) const
        {
        auto base = sector_state_plan();
        int N = int(base.size());
        std::map<int,std::vector<int>> positions; // type code -> site positions
        for (int i=0;i<N;i++) positions[site_types_.at(i)].push_back(i);
        std::mt19937 rng(1234u + unsigned(sector_draws_++));
        std::vector<std::vector<int>> out;
        out.push_back(base); // the deterministic plan is always one of them
        for (int r=1;r<k;r++)
            {
            auto plan = base;
            for (auto const& kv : positions)
                {
                auto pos = kv.second;
                std::vector<int> vals;
                for (auto i : pos) vals.push_back(base.at(i));
                std::shuffle(vals.begin(),vals.end(),rng);
                for (size_t j=0;j<pos.size();++j) plan[pos.at(j)] = vals.at(j);
                }
            out.push_back(plan);
            }
        return out;
        }

    // A bond-dimension-1 product MPS in the requested sector, from a
    // per-site basis-state assignment. Mirrors ITensor's own
    // init_tensors() (mps.cc, backing MPS(InitState)) -- the QN-carrying
    // Link indices are the reason this cannot be assembled from plain
    // dim-1 links the way build_product_state() does: an ITensor's indices
    // are either all QN-carrying or all dense, never mixed. Taking the
    // basis states as IndexVals rather than names (which is all InitState
    // itself stores) avoids needing a per-site-type table of state names.
    MPS
    sector_product_mps(std::vector<int> const& states) const
        {
        int N = sites_.length();
        auto a = std::vector<Index>(N+1);
        auto qa = std::vector<QN>(N+1); // qa[i] = QN on the i-th bond
        for (int i=1;i<=N;i++) qa[0] -= qn(sites_(i)(states.at(i-1)))*In;
        for (int i=1;i<=N;i++) // zero total divergence, solve bond by bond
            qa[i] = Out*(-qa[i-1]*In - qn(sites_(i)(states.at(i-1))));
        for (int i=1;i<N;i++) a[i] = Index(qa[i],1,tinyformat::format("Link,l=%d",i));
        MPS psi(N);
        if (N==1) { psi.set(1,setElt(sites_(1)(states.at(0)))); }
        else
            {
            psi.set(1,setElt(sites_(1)(states.at(0)),a[1](1)));
            for (int i=2;i<N;i++)
                psi.set(i,setElt(dag(a[i-1])(1),sites_(i)(states.at(i-1)),a[i](1)));
            psi.set(N,setElt(dag(a[N-1])(1),sites_(N)(states.at(N-1))));
            }
        psi.leftLim(0); // MPS(InitState)'s own convention: OC at site 1
        psi.rightLim(2);
        return psi;
        }

    // The sector-mode counterpart of randomMPS(sites_,m): a normalized sum
    // of random in-sector product states. Every summand has the same
    // (correct) flux, so the sum does too; bond dimension is bounded by
    // the number of summands, and two-site DMRG grows it from there the
    // same way it does from a randomMPS start. Measured to reach the same
    // energy as the dense/randomMPS path to twelve digits on a Heisenberg
    // chain (see get_sites.h's comment).
    MPS
    sector_mps(int m) const
        {
        int k = std::max(1,std::min(m,16));
        auto plans = sector_state_arrangements(k);
        std::vector<MPS> parts;
        parts.reserve(plans.size());
        for (auto const& plan : plans) parts.push_back(sector_product_mps(plan));
        MPS psi = (parts.size()==1) ? parts.front()
                                    : sum(parts,{"Cutoff",1e-14,"MaxDim",std::max(1,m)});
        psi.normalize();
        return psi;
        }

    // Refuse a code path that assembles tensors by hand -- dense Link
    // indices (METTS product states, iDMRG/VUMPS environments) or bond
    // gates summed before their fluxes agree (TEBD). An ITensor's indices
    // are either all QN-carrying or all dense, and tensors of different
    // flux cannot be added, so such a path does not merely give a wrong
    // answer under sector mode: it aborts the process inside ITensor
    // (Error() -> abort(), nothing to catch). Better to say so.
    //
    // Not a general statement about time evolution or correlators -- both
    // TDVP variants and every AutoMPO-based path work under a sector, as
    // long as the operators involved conserve it.
    void
    reject_sector(char const* what) const
        {
        if (!has_sector_) return;
        throw std::invalid_argument(tinyformat::format("Chain::%s does not support "
              "conserved-sector mode (%s): this code path assembles tensors that are "
              "dense or of mixed flux by construction. Clear the sector with "
              "set_conserved_sector() first (or, for time evolution, use the default "
              "tevol_method=\"TDVP\", which does support it).",
              what,sector_string()));
        }

    // Check an MPS actually came out in the requested sector. Cheap (the
    // flux is read off the tensors, not measured), and it is the one thing
    // that would silently produce a wrong answer if any of the above were
    // subtly wrong.
    void
    check_sector(MPS const& psi, char const* where) const
        {
        if (!has_sector_) return;
        auto q = totalQN(psi);
        for (auto const& sq : sector_)
            if (q.val(sq.first) != sq.second)
                throw std::runtime_error(tinyformat::format("Chain::%s: the state is in sector %s, not the "
                      "requested %s",where,qn_string(q),sector_string()));
        }

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
    default_mps() const { return default_mps(maxm_); }

    // Same, at an explicit bond dimension -- used by gs_energy() so that
    // a ramped schedule builds its random starting state directly at the
    // ramp's first maxdim instead of at the full maxm_.
    MPS
    default_mps(int m) const
        {
        // Sector mode cannot use randomMPS at all: randomMPS(SiteSet)
        // refuses to guess a sector under QNs, and this ITensor has no
        // randomMPS(InitState,m>1) (both are hard Errors in mps.cc).
        if (has_sector_) return sector_mps(m<1?1:m);
        return randomMPS(sites_,m<1?1:m);
        }

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
        reject_sector("build_product_state"); // dense links below (METTS)
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

    // Name-based fermion classification, mirroring ITensor's own
    // isFermionic() (ITensor/itensor/mps/autompo.cc) and pyitensor's
    // sites/base.py::is_fermionic: any operator name starting with 'C' is
    // fermionic for Jordan-Wigner purposes, on every site type.
    static bool
    idmrg_is_fermionic(std::string const& name)
        { return !name.empty() && name[0]=='C'; }

    // Classifies every term (terms_intra ++ terms_inter, MOTerm's own
    // 1-based sites converted to 0-based here to match idmrg.py's own
    // arithmetic exactly) into 1-site (onsite) and 2-site (bonds) pieces,
    // composing multiple same-site operator factors via reversed matrix
    // product exactly like _term_site_matrices does (mats[-1] first,
    // then left-multiplied by each earlier factor in reverse -- i.e. for
    // factors written [A,B,C] at one site, the composed matrix is C@B@A).
    //
    // The terms arrive *without* any Jordan-Wigner transform applied
    // (infinitechain.py passes to_terms(jordan_wigner_transform=False) --
    // the finite-chain transform's strings run from site 1 of the chain,
    // which an infinite chain does not have). This function threads the
    // string itself, locally and translation-invariantly, exactly like
    // pyitensor/idmrg.py's _term_site_matrices:
    //   * fermionic reordering sign when sorting factors into site order
    //     (by_site below does that sorting, so the sign has to be supplied
    //     here -- autompo.cc's HTerm::add() applies the same sign while
    //     insertion-sorting its own factors);
    //   * an endpoint F factor on any touched site whose own parity differs
    //     from the parity carried in from its left (HTerm::resolve()'s own
    //     need_F);
    //   * carry_ferm on the resulting 2-site bond, so idmrg_build_row
    //     propagates F rather than Id along the pending channel spanning
    //     the term's two endpoints;
    //   * rejection of odd total fermion parity, whose string could never
    //     be closed within the unit cell.
    // See _term_site_matrices' own (much longer) docstring for the
    // derivation of each of these and for what silently breaks without
    // them (in particular F@C == -C flips the sign of a C-leading hopping
    // term, i.e. a silently non-Hermitian Hamiltonian).
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
            Cplx coef = term.coef;
            // Sorting the factors into site order is a *fermionic*
            // reordering: every pair of fermionic factors swapped past each
            // other contributes a -1.
            for (size_t x=0;x<term.factors.size();++x)
            for (size_t y=x+1;y<term.factors.size();++y)
                if (term.factors[x].site > term.factors[y].site
                        && idmrg_is_fermionic(term.factors[x].name)
                        && idmrg_is_fermionic(term.factors[y].name))
                    coef = -coef;
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
            std::vector<bool> ferm;
            bool carry = false; // parity carried in from the touched sites to the left
            for (int rel : rel_sites)
                {
                int p = ((rel % n_uc)+n_uc) % n_uc;
                auto const& names = by_site[rel];
                int d = dim(sites_.si(p+1));
                auto combined = idmrg_op_dense(p,names.back());
                for (int k=(int)names.size()-2;k>=0;--k)
                    combined = idmrg_matmul(combined,idmrg_op_dense(p,names[k]),d);
                int nferm = 0;
                for (auto const& nm : names) if (idmrg_is_fermionic(nm)) ++nferm;
                bool site_ferm = (nferm%2)==1;
                if (carry != site_ferm)
                    combined = idmrg_matmul(idmrg_op_dense(p,"F"),combined,d);
                carry = (carry != site_ferm);
                mats.push_back(std::move(combined));
                ferm.push_back(site_ferm);
                }
            if (carry)
                Error("Chain::idmrg_ground_state: a term has odd total "
                      "fermion parity -- its Jordan-Wigner string cannot be "
                      "closed within the unit cell");
            if (rel_sites.size()==1)
                onsite.push_back(IdmrgOnsite{rel_sites[0],coef,mats[0]});
            else if (rel_sites.size()==2)
                {
                int a = rel_sites[0], b = rel_sites[1];
                if (a>=n_uc)
                    Error("Chain::idmrg_ground_state: internal error -- "
                          "inter-cell term does not touch the central cell");
                // carry_ferm: the string is open between the two endpoints
                // exactly when the first one has odd parity of its own.
                bonds.push_back(IdmrgBond{a,b,mats[0],mats[1],ferm[0],coef});
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
    // function's docstring for the derivation of each case), including its
    // Jordan-Wigner string handling: the last branch (a pending channel
    // propagating one more site) carries this site's own "F" whenever the
    // bond it belongs to has carry_ferm set, and plain Id otherwise.
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
                    setmat(li,ri,bonds[lch.bond_index].carry_ferm
                                    ? idmrg_op_dense(p,"F") : eye());
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
        // vx_build_linear_map applies `action` D*D times purely to
        // MATERIALIZE the matrix, then vx_solve is another O(D^6) LU on it
        // -- the single largest cost in the grouped VUMPS iteration. Above
        // the threshold, solve the same operator matrix-free instead,
        // exactly as the multi-site path's own solve_reg already does.
        if (D*D > vms_dense_solve_max_)
            {
            bool ok=false;
            auto sol = vx_bicgstab(action,rhs,D*D,1e-12,20*D*D,ok);
            // Fall through to the dense path rather than trust a
            // non-converged iterate -- never less robust than before.
            if (ok) return sol;
            }
        auto Mat = vx_build_linear_map(D,action);
        return vx_solve(Mat,D*D,rhs);
        }

    // Relative tolerances for the transfer-matrix dominant-eigenvalue
    // selection, mirroring pyitensor/idmrg.py's _DEGENERACY_RTOL and
    // _PERRON_TIE_RTOL. See _check_dominant_eigenvalue_nondegenerate's own
    // docstring there for the full reasoning. The short version:
    //
    // A magnitude-only test rejects legitimate period-p physics. A
    // transfer matrix's peripheral spectrum is rho*exp(2 pi i k / p) for a
    // period-p state, so a period-2 state (any half-filled/2k_F=pi chain,
    // and a gapless Heisenberg chain along some trajectories) carries an
    // eigenvalue at -rho whose magnitude approaches the leader's as the
    // correlation length diverges. Those are DISTINCT eigenvalues with
    // distinct eigenvectors -- nothing about them is ill-posed. Only a
    // repeated EIGENVALUE (the "cat state", both copies at +rho) is
    // genuinely ambiguous.
    //
    // Measured before this fix: the grouped VUMPS path here could not
    // complete a D=8 Heisenberg ground state at all ("every attempt at
    // D=4 failed (degenerate...)"), and on the pyitensor side the same
    // magnitude test killed ~42% of individual VUMPS attempts and ~16% of
    // whole solves on a half-filled free-fermion chain at D=16, where
    // every firing had the signature |lambda|=(1,0.999999999),
    // arg=(0,+-pi).
    static constexpr double vx_degeneracy_rtol_ = 1e-9;
    // Wider than vx_degeneracy_rtol_ on purpose: the observed tie sits at
    // exactly 1e-9, i.e. on that constant's own boundary, so reusing it
    // would leave the Perron selection decided by rounding.
    static constexpr double vx_perron_tie_rtol_ = 1e-6;

    // Reorder `order` (already sorted by descending |evals|) so that among
    // the entries tied in MAGNITUDE with the leader, the one with the
    // largest real part comes first -- the Perron root, the only member of
    // a peripheral spectrum whose eigenvector is the positive
    // (semidefinite) matrix every fixed-point caller treats as a density
    // matrix and divides by the trace of. Without this, which member comes
    // back is whatever the sort happened to do with a tie, so a period-2
    // state can hand out the -rho eigenvector as a "density matrix".
    static void
    vx_perron_reorder(std::vector<int>& order, std::vector<Cplx> const& evals)
        {
        if (order.size() < 2) return;
        double top = std::abs(evals[order[0]]);
        if (!(top > 0.0)) return;
        size_t ntied = 0;
        while (ntied < order.size()
               && std::abs(evals[order[ntied]]) >= (1.0-vx_perron_tie_rtol_)*top)
            ++ntied;
        if (ntied < 2) return;
        std::stable_sort(order.begin(), order.begin()+ntied,
                          [&](int a,int b)
                          { return evals[a].real() > evals[b].real(); });
        }

    // Degenerate means the same EIGENVALUE twice, not merely the same
    // magnitude -- see vx_degeneracy_rtol_ above.
    static void
    vx_check_perron_nondegenerate(Cplx a, Cplx b, std::string const& where)
        {
        if (std::abs(b-a) <= vx_degeneracy_rtol_*std::abs(a))
            throw ITError(where+": the transfer matrix's dominant eigenvalue is "
                           "genuinely degenerate (two copies of the same eigenvalue) "
                           "-- a single dominant fixed point is not well-defined here. "
                           "See pyitensor/idmrg.py's own "
                           "_check_dominant_eigenvalue_nondegenerate for the physical "
                           "reason (a \"cat state\" superposition of two branches with "
                           "matched per-site norm)");
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
        vx_perron_reorder(order,evals);
        if (n>1) vx_check_perron_nondegenerate(evals[order[0]],evals[order[1]],
                                                "Chain::vumps");
        int idx = order[0];
        std::vector<Cplx> vec(n);
        for (int i=0;i<n;++i) vec[i] = vr[i + (size_t)idx*n];
        return {evals[idx], vec};
        }

    // (eta, v): the dominant eigenpair of the transfer tensor E acting
    // either from the right (from_left=false, apply_transfer) or from the
    // left (from_left=true, apply_transfer_from_left), dense below
    // vx_dense_eig_max_ and matrix-free restarted Arnoldi above it.
    //
    // The dense route is O(D^6) -- zgeev on the whole (D^2,D^2) matrix --
    // and it runs TWICE per VUMPS iteration (once per side), which is what
    // made the grouped n_uc<=2 path far slower than the pure-Python one it
    // ports. The multi-site (vms_*) path already went matrix-free for
    // exactly this reason; this brings the grouped path in line, reusing
    // the same ic_arnoldi_dominant rather than a second implementation.
    //
    // The left case uses the from-left ACTION rather than materializing
    // E's transpose: transposing is itself O(D^4) storage and defeats the
    // point of never forming the matrix.
    //
    // Below the threshold nothing changes -- the small cases every v3
    // VUMPS test was validated against stay bit-identical to the dense
    // path they were validated on.
    static std::pair<Cplx,std::vector<Cplx>>
    vx_dominant_pair(std::vector<Cplx> const& E, int D, bool from_left)
        {
        int n = D*D;
        if (n <= vx_dense_eig_max_)
            {
            if (!from_left) return vx_dominant_eig(E,n);
            std::vector<Cplx> Et((size_t)n*n);
            for (int i=0;i<n;++i)
            for (int j=0;j<n;++j)
                Et[(size_t)i*n+j] = E[(size_t)j*n+i];
            return vx_dominant_eig(Et,n);
            }
        auto act = [&](std::vector<Cplx> const& x)
            {
            return from_left ? vx_apply_transfer_from_left(E,D,x)
                             : vx_apply_transfer(E,D,x);
            };
        Cplx eta(0,0), eta1(0,0);
        std::vector<Cplx> vec;
        ic_arnoldi_dominant(act,n,eta,eta1,vec);
        vx_check_perron_nondegenerate(eta,eta1,"Chain::vumps");
        return {eta,vec};
        }

    // Dominant RIGHT fixed point rho of transfer tensor E (apply_transfer(E,rho)=eta*rho),
    // normalized to trace(rho)=1 -- C++ analogue of pyitensor idmrg.py's
    // _dominant_right_fixed_point. E's own row-major (l,L,r,R) flattening
    // already IS the (D*D,D*D) matrix vx_dominant_eig needs (see this
    // section's own top comment) -- no reshaping.
    static std::pair<std::vector<Cplx>,Cplx>
    vx_dominant_right_fixed_point(std::vector<Cplx> const& E, int D)
        {
        auto [eta,vec] = vx_dominant_pair(E,D,false);
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
        auto [eta,vec] = vx_dominant_pair(E,D,true);
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

    // Above this, the multi-site environment's regularized solve goes
    // matrix-free (vx_bicgstab); at or below it the dense build+LU is
    // cheaper and is kept. Same shape of split, and the same reason, as
    // ic_dense_eig_max_ and vx_dense_eig_max_ above -- and the same reason
    // pyitensor's own _solve_linear_map switches to GMRES past
    // _DENSE_SOLVE_MAX. vx_regularized_solve's dense route costs D^2
    // applications of the action to BUILD the matrix plus an O(D^6) LU to
    // solve it, twice per outer iteration, on every rung of the D-ramp:
    // measured on a 3-site TFIM cell, D=24 took 64.7s against the
    // pure-Python backend's 9.4s with everything else already fixed.
    static constexpr int vms_dense_solve_max_ = 256;   // i.e. D <= 16

    // Matrix-free BiCGSTAB for a general (non-Hermitian) linear map given
    // only by its action on a flat length-n complex vector. Returns the
    // solution and sets `ok` -- a caller that gets ok=false is expected to
    // fall back to the dense path rather than trust the result, exactly as
    // pyitensor's _solve_linear_map falls back when GMRES does not
    // converge, so this is never LESS robust than the solve it replaces.
    template <typename Fn>
    static std::vector<Cplx>
    vx_bicgstab(Fn&& A, std::vector<Cplx> const& b, int n, double tol,
                 int maxit, bool& ok)
        {
        auto dot = [](std::vector<Cplx> const& x, std::vector<Cplx> const& y)
            {
            Cplx s(0,0);
            for (size_t i=0;i<x.size();++i) s += std::conj(x[i])*y[i];
            return s;
            };
        auto nrm = [&](std::vector<Cplx> const& x)
            { return std::sqrt(std::abs(dot(x,x))); };

        double bnorm = nrm(b);
        ok = false;
        std::vector<Cplx> x(n,Cplx(0,0));
        if (bnorm <= 0.0) { ok = true; return x; }
        std::vector<Cplx> r = b, rhat = b;
        std::vector<Cplx> v(n,Cplx(0,0)), p(n,Cplx(0,0));
        Cplx rho(1,0), alpha(1,0), omega(1,0);
        for (int it=0; it<maxit; ++it)
            {
            Cplx rho_new = dot(rhat,r);
            if (std::abs(rho_new) < 1e-300) return x;         // breakdown
            if (it==0)
                p = r;
            else
                {
                Cplx beta = (rho_new/rho)*(alpha/omega);
                for (int i=0;i<n;++i) p[i] = r[i] + beta*(p[i] - omega*v[i]);
                }
            v = A(p);
            Cplx den = dot(rhat,v);
            if (std::abs(den) < 1e-300) return x;             // breakdown
            alpha = rho_new/den;
            std::vector<Cplx> sv(n);
            for (int i=0;i<n;++i) sv[i] = r[i] - alpha*v[i];
            if (nrm(sv) < tol*bnorm)
                {
                for (int i=0;i<n;++i) x[i] += alpha*p[i];
                ok = true;
                return x;
                }
            auto t = A(sv);
            Cplx tt = dot(t,t);
            if (std::abs(tt) < 1e-300) return x;              // breakdown
            omega = dot(t,sv)/tt;
            for (int i=0;i<n;++i) x[i] += alpha*p[i] + omega*sv[i];
            for (int i=0;i<n;++i) r[i] = sv[i] - omega*t[i];
            if (nrm(r) < tol*bnorm) { ok = true; return x; }
            if (std::abs(omega) < 1e-300) return x;           // breakdown
            rho = rho_new;
            }
        return x;
        }

    static constexpr int _vms_S = 0;   // automaton channel indices,
    static constexpr int _vms_F = 1;   // idmrg_build_row's convention

    // ==== Sequential multi-site VUMPS (any n_uc) =========================
    //
    // C++ port of pyitensor/vumps_ms.py -- see that module's own docstring
    // for the algorithm and the literature. In one line: the grouped route
    // right below folds the whole cell into one supersite of dimension
    // prod(d_p), which is exponential in the cell size and is why this
    // backend capped n_uc at 2; the sequential algorithm instead keeps a
    // LIST of per-site tensors and sweeps the cell, solving one H_AC[n] per
    // site and one H_C[n] per bond, at a cost linear in n_uc.
    //
    // The environments are fully channel-resolved here (one (D,D) matrix
    // per automaton channel, per bond) rather than the grouped path's
    // reach-1 triple of {GL, GR, bond_envs}: that specialization is exactly
    // what grouping bought (any coupling inside the cell is on-site once
    // the cell is one supersite), and ungrouped it no longer holds.
    //
    // The per-site automaton is idmrg_build_row's own IdmrgAutomatonRow,
    // used directly with no grouping step: row.flat[(l*right_n+r)*d*d +
    // si*d + so], channels S=0/F=1/pending=2.., and left_n/right_n are the
    // channel counts on this site's two bonds (they may differ, so every
    // loop below takes them from the row rather than assuming one Dw).

    // X[l*D+L] -> out[r*D+R], one site forward. Same conventions as
    // vx_cap_left/vx_cap_right (ket index first, bra second).
    static std::vector<Cplx>
    vms_push_left(std::vector<Cplx> const& X, std::vector<Cplx> const& A,
                   int D, int d, std::vector<Cplx> const& M, bool hasM)
        {
        std::vector<Cplx> B = hasM ? vx_apply_op_ket(M,A,D,d) : A;
        // T[L,p,r] = sum_l X[l,L] B[l,p,r]
        std::vector<Cplx> T((size_t)D*d*D,Cplx(0,0));
        for (int Lb=0;Lb<D;++Lb)
        for (int p=0;p<d;++p)
        for (int r=0;r<D;++r)
            {
            Cplx acc(0,0);
            for (int l=0;l<D;++l) acc += X[l*D+Lb]*B[(l*d+p)*D+r];
            T[(Lb*d+p)*D+r] = acc;
            }
        // out[r,R] = sum_{L,p} T[L,p,r] conj(A[L,p,R])
        std::vector<Cplx> out((size_t)D*D,Cplx(0,0));
        for (int r=0;r<D;++r)
        for (int R=0;R<D;++R)
            {
            Cplx acc(0,0);
            for (int Lb=0;Lb<D;++Lb)
            for (int p=0;p<d;++p)
                acc += T[(Lb*d+p)*D+r]*std::conj(A[(Lb*d+p)*D+R]);
            out[r*D+R] = acc;
            }
        return out;
        }

    // Y[r*D+R] -> out[l*D+L], one site backward.
    static std::vector<Cplx>
    vms_push_right(std::vector<Cplx> const& Y, std::vector<Cplx> const& A,
                    int D, int d, std::vector<Cplx> const& M, bool hasM)
        {
        std::vector<Cplx> B = hasM ? vx_apply_op_ket(M,A,D,d) : A;
        // T[l,p,R] = sum_r B[l,p,r] Y[r,R]
        std::vector<Cplx> T((size_t)D*d*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int p=0;p<d;++p)
        for (int R=0;R<D;++R)
            {
            Cplx acc(0,0);
            for (int r=0;r<D;++r) acc += B[(l*d+p)*D+r]*Y[r*D+R];
            T[(l*d+p)*D+R] = acc;
            }
        // out[l,L] = sum_{p,R} T[l,p,R] conj(A[L,p,R])
        std::vector<Cplx> out((size_t)D*D,Cplx(0,0));
        for (int l=0;l<D;++l)
        for (int Lb=0;Lb<D;++Lb)
            {
            Cplx acc(0,0);
            for (int p=0;p<d;++p)
            for (int R=0;R<D;++R)
                acc += T[(l*d+p)*D+R]*std::conj(A[(Lb*d+p)*D+R]);
            out[l*D+Lb] = acc;
            }
        return out;
        }

    // One site forward for a channel-resolved left environment:
    // G'[b] = sum_a push(G[a], W[a,:,:,b]). Zero automaton slabs are
    // skipped -- the automaton is sparse, and skipping is what keeps the
    // per-site cost proportional to the number of live channels.
    static std::vector<std::vector<Cplx>>
    vms_push_site_left(std::vector<std::vector<Cplx>> const& G,
                        std::vector<Cplx> const& A,
                        IdmrgAutomatonRow const& row, int D)
        {
        int nl = row.left_n, nr = row.right_n, d = row.d;
        std::vector<std::vector<Cplx>> out(nr, std::vector<Cplx>((size_t)D*D,Cplx(0,0)));
        std::vector<Cplx> slab((size_t)d*d);
        for (int a=0;a<nl;++a)
            {
            bool any_a = false;
            for (auto const& z : G[a]) if (z!=Cplx(0,0)) { any_a=true; break; }
            if (!any_a) continue;
            for (int b=0;b<nr;++b)
                {
                bool any_w = false;
                for (int si=0;si<d;++si)
                for (int so=0;so<d;++so)
                    {
                    slab[si*d+so] = row.flat[((size_t)a*nr+b)*d*d + si*d+so];
                    if (slab[si*d+so]!=Cplx(0,0)) any_w = true;
                    }
                if (!any_w) continue;
                auto piece = vms_push_left(G[a],A,D,d,slab,true);
                for (size_t k=0;k<piece.size();++k) out[b][k] += piece[k];
                }
            }
        return out;
        }

    // Mirror: one site backward, G'[a] = sum_b push(G[b], W[a,:,:,b]).
    static std::vector<std::vector<Cplx>>
    vms_push_site_right(std::vector<std::vector<Cplx>> const& G,
                         std::vector<Cplx> const& A,
                         IdmrgAutomatonRow const& row, int D)
        {
        int nl = row.left_n, nr = row.right_n, d = row.d;
        std::vector<std::vector<Cplx>> out(nl, std::vector<Cplx>((size_t)D*D,Cplx(0,0)));
        std::vector<Cplx> slab((size_t)d*d);
        for (int b=0;b<nr;++b)
            {
            bool any_b = false;
            for (auto const& z : G[b]) if (z!=Cplx(0,0)) { any_b=true; break; }
            if (!any_b) continue;
            for (int a=0;a<nl;++a)
                {
                bool any_w = false;
                for (int si=0;si<d;++si)
                for (int so=0;so<d;++so)
                    {
                    slab[si*d+so] = row.flat[((size_t)a*nr+b)*d*d + si*d+so];
                    if (slab[si*d+so]!=Cplx(0,0)) any_w = true;
                    }
                if (!any_w) continue;
                auto piece = vms_push_right(G[b],A,D,d,slab,true);
                for (size_t k=0;k<piece.size();++k) out[a][k] += piece[k];
                }
            }
        return out;
        }

    // The whole cell's plain (operator-free) transfer, as an explicit
    // E4[l,L,r,R] -- built by pushing each of the D^2 basis matrices
    // through the cell.
    //
    // Materializing it costs O(n_uc * D^5 * d) once per outer iteration,
    // where the matrix-free form would be O(D^3) per application. That is
    // a deliberate trade: it lets the two things built on it -- the
    // dominant fixed point and the regularized environment solve -- reuse
    // vx_dominant_*_fixed_point and vx_regularized_solve unchanged, which
    // are the validated single-site routines this backend already trusts,
    // rather than growing matrix-free variants of both. The grouped path
    // already builds a single-site E4 at O(D^4 d) for exactly the same
    // reason, so this stays within this backend's existing design point;
    // if D ever grows past what that supports, both paths need the same
    // matrix-free treatment, not just this one.
    static std::vector<Cplx>
    vms_cell_transfer(std::vector<std::vector<Cplx>> const& A_list,
                       std::vector<IdmrgAutomatonRow> const& rows, int D,
                       bool from_left)
        {
        int n_uc = (int)A_list.size();
        std::vector<Cplx> E((size_t)D*D*D*D,Cplx(0,0));
        std::vector<Cplx> X((size_t)D*D,Cplx(0,0));
        for (int i=0;i<D;++i)
        for (int j=0;j<D;++j)
            {
            std::fill(X.begin(),X.end(),Cplx(0,0));
            X[i*D+j] = Cplx(1,0);
            std::vector<Cplx> Y = X;
            if (from_left)
                for (int n=0;n<n_uc;++n)
                    Y = vms_push_left(Y,A_list[n],D,rows[n].d,{},false);
            else
                for (int n=n_uc-1;n>=0;--n)
                    Y = vms_push_right(Y,A_list[n],D,rows[n].d,{},false);
            for (int r=0;r<D;++r)
            for (int R=0;R<D;++R)
                {
                // from_left: Y[r,R] = sum_{l,L} X[l,L] E[l,L,r,R], so this
                // basis element fills E[i,j,:,:]. from_right: the mirror,
                // E[:,:,i,j] via vx_apply_transfer's own convention.
                if (from_left) E[((i*D+j)*D+r)*D+R] = Y[r*D+R];
                else           E[((r*D+R)*D+i)*D+j] = Y[r*D+R];
                }
            }
        return E;
        }

    // (G, source): the channel-resolved left environment at the cell's left
    // edge for every channel EXCEPT F, plus the F content one cell
    // traversal produces from it.
    //
    // G[S]=I by definition; the pending block is nilpotent (a Hamiltonian
    // term has finite reach, so a partially-applied term must close within
    // a bounded number of sites), so plain iteration reaches it exactly
    // rather than needing a solve. F alone needs the regularized solve --
    // it is the one channel that feeds back into itself (W[F,:,:,F]=I), so
    // summing it is a geometric series over the half-infinite chain.
    static void
    vms_nilpotent_left(std::vector<std::vector<Cplx>> const& AL,
                        std::vector<IdmrgAutomatonRow> const& rows, int D,
                        std::vector<std::vector<Cplx>>& G,
                        std::vector<Cplx>& source)
        {
        int n_uc = (int)AL.size();
        int n0 = rows[0].left_n;
        G.assign(n0, std::vector<Cplx>((size_t)D*D,Cplx(0,0)));
        G[_vms_S] = vx_eye(D);
        auto push_cell = [&](std::vector<std::vector<Cplx>> const& Gin)
            {
            auto cur = Gin;
            for (int n=0;n<n_uc;++n) cur = vms_push_site_left(cur,AL[n],rows[n],D);
            return cur;
            };
        for (int sweep=0; sweep<n0+1; ++sweep)
            {
            auto next = push_cell(G);
            next[_vms_S] = vx_eye(D);
            if ((int)next.size() > _vms_F)
                std::fill(next[_vms_F].begin(),next[_vms_F].end(),Cplx(0,0));
            double diff = 0.0;
            for (size_t a=0;a<next.size() && a<G.size();++a)
            for (size_t k=0;k<next[a].size();++k)
                diff = std::max(diff,std::abs(next[a][k]-G[a][k]));
            G = next;
            if (diff < 1e-14) break;
            }
        source = push_cell(G)[_vms_F];
        }

    // Mirror of vms_nilpotent_left, with S and F swapping roles: on the
    // right it is F that means "nothing applied yet" (G[F]=I) and S that
    // accumulates.
    static void
    vms_nilpotent_right(std::vector<std::vector<Cplx>> const& AR,
                         std::vector<IdmrgAutomatonRow> const& rows, int D,
                         std::vector<std::vector<Cplx>>& G,
                         std::vector<Cplx>& source)
        {
        int n_uc = (int)AR.size();
        int nlast = rows[n_uc-1].right_n;
        G.assign(nlast, std::vector<Cplx>((size_t)D*D,Cplx(0,0)));
        G[_vms_F] = vx_eye(D);
        auto push_cell = [&](std::vector<std::vector<Cplx>> const& Gin)
            {
            auto cur = Gin;
            for (int n=n_uc-1;n>=0;--n) cur = vms_push_site_right(cur,AR[n],rows[n],D);
            return cur;
            };
        for (int sweep=0; sweep<nlast+1; ++sweep)
            {
            auto next = push_cell(G);
            next[_vms_F] = vx_eye(D);
            if ((int)next.size() > _vms_S)
                std::fill(next[_vms_S].begin(),next[_vms_S].end(),Cplx(0,0));
            double diff = 0.0;
            for (size_t a=0;a<next.size() && a<G.size();++a)
            for (size_t k=0;k<next[a].size();++k)
                diff = std::max(diff,std::abs(next[a][k]-G[a][k]));
            G = next;
            if (diff < 1e-14) break;
            }
        source = push_cell(G)[_vms_S];
        }

    // Per-bond channel-resolved environments for the whole cell, and the
    // energy density per unit cell. GL[n] sits immediately LEFT of site n,
    // GR[n] immediately right. Built once at the cell edges and pushed
    // across -- one traversal gives every site's environment, which is the
    // point of the sequential formulation.
    struct VmsEnv
        {
        std::vector<std::vector<std::vector<Cplx>>> GL, GR;
        double e_cell = 0.0;
        };

    VmsEnv
    vms_environments(std::vector<std::vector<Cplx>> const& AL,
                      std::vector<std::vector<Cplx>> const& AR,
                      std::vector<IdmrgAutomatonRow> const& rows, int D) const
        {
        int n_uc = (int)AL.size();
        // The two dominant fixed points are found MATRIX-FREE above
        // ic_dense_eig_max_, by the same restarted Arnoldi the iDMRG
        // correlator path uses. The dense route (vx_dominant_eig) runs a
        // full zgeev on the (D^2,D^2) transfer -- O(D^6), twice per outer
        // iteration, across every rung of the D-ramp -- which is exactly
        // the cost ic_dense_eig_max_'s own comment records as having
        // "dominated an entire iDMRG solve" at chi=24. Measured here on a
        // 3-site TFIM cell before this split: D=16 took 24.5s against the
        // pure-Python backend's 2.0s, i.e. the C++ port was 12x SLOWER
        // than the reference it was meant to accelerate.
        auto dominant = [&](bool from_left, std::vector<std::vector<Cplx>> const& A)
            {
            int n = D*D;
            std::vector<Cplx> vec;
            if (n <= ic_dense_eig_max_)
                {
                auto E = vms_cell_transfer(A,rows,D,from_left);
                auto [v,eta] = from_left ? vx_dominant_left_fixed_point(E,D)
                                          : vx_dominant_right_fixed_point(E,D);
                (void)eta;
                return v;
                }
            auto act = [&](std::vector<Cplx> const& X)
                {
                auto Y = X;
                if (from_left)
                    for (int m=0;m<n_uc;++m)
                        Y = vms_push_left(Y,A[m],D,rows[m].d,{},false);
                else
                    for (int m=n_uc-1;m>=0;--m)
                        Y = vms_push_right(Y,A[m],D,rows[m].d,{},false);
                return Y;
                };
            Cplx e0c(0,0), e1c(0,0);
            ic_arnoldi_dominant(act,n,e0c,e1c,vec);
            vx_check_perron_nondegenerate(e0c,e1c,"Chain::vms_environments");
            Cplx tr(0,0);
            for (int i=0;i<D;++i) tr += vec[i*D+i];
            if (std::abs(tr) < 1e-13)
                throw ITError("Chain::vms_environments: dominant fixed point has "
                               "~zero trace -- degenerate/ill-defined normalization");
            for (auto& z : vec) z /= tr;
            return vec;
            };
        auto r_AL = vx_hermitize(dominant(false,AL),D);
        auto l_AR = vx_hermitize(dominant(true,AR),D);
        // The eigensolver leaves the scale free; these close a normalized
        // state against the other side's exact (identity) fixed point, so
        // their trace must be 1.
        Cplx tr_r(0,0), tr_l(0,0);
        for (int i=0;i<D;++i) { tr_r += r_AL[i*D+i]; tr_l += l_AR[i*D+i]; }
        if (std::abs(tr_r)>1e-300) for (auto& z : r_AL) z /= tr_r;
        if (std::abs(tr_l)>1e-300) for (auto& z : l_AR) z /= tr_l;

        std::vector<std::vector<Cplx>> GL0, GRlast;
        std::vector<Cplx> src_l, src_r;
        vms_nilpotent_left(AL,rows,D,GL0,src_l);
        vms_nilpotent_right(AR,rows,D,GRlast,src_r);

        // trace(conj(r_AL) @ source_l), NOT an elementwise sum -- the
        // transpose in the second factor is the whole difference between
        // the two, and getting it wrong leaves the STATE correct (vev and
        // correlators still match to 1e-9) while the reported energy is
        // silently off, which is exactly how this was caught. Uses the
        // same helper the grouped path does so the two cannot diverge.
        double e_L = vx_trace_conjA_X(r_AL,src_l,D).real();
        double e_R = vx_trace_conjA_X(l_AR,src_r,D).real();

        // GL[F] / GR[S]: the same regularized (I - T_cell + P) solve the
        // grouped path does, with the single-site transfer replaced by the
        // whole cell's. Each side is regularized against ITS OWN energy --
        // never a shared average, which vumps_solve_left_environment's own
        // Python counterpart documents as a confirmed period-2 limit cycle.
        std::vector<Cplx> I = vx_eye(D);
        std::vector<Cplx> rhs_l((size_t)D*D), rhs_r((size_t)D*D);
        for (int k=0;k<D*D;++k)
            {
            rhs_l[k] = src_l[k] - e_L*I[k];
            rhs_r[k] = src_r[k] - e_R*I[k];
            }
        // (I - T_cell + P)[G] = source - e*I, with P[X] = I*trace(conj(fp)@X).
        // Matrix-free above vms_dense_solve_max_: the action is a cell
        // traversal (O(n_uc D^3 d)), where the dense route needs D^2 of
        // those just to BUILD the matrix plus an O(D^6) LU on top.
        auto make_reg = [&](std::vector<std::vector<Cplx>> const& A, bool from_left,
                             std::vector<Cplx> const& fp)
            {
            return [&,from_left](std::vector<Cplx> const& x)
                {
                auto Tx = x;
                if (from_left)
                    for (int m=0;m<n_uc;++m)
                        Tx = vms_push_left(Tx,A[m],D,rows[m].d,{},false);
                else
                    for (int m=n_uc-1;m>=0;--m)
                        Tx = vms_push_right(Tx,A[m],D,rows[m].d,{},false);
                std::vector<Cplx> out((size_t)D*D);
                Cplx sc = vx_trace_conjA_X(fp,x,D);
                for (int k=0;k<D*D;++k) out[k] = x[k] - Tx[k] + I[k]*sc;
                return out;
                };
            };
        auto solve_reg = [&](std::vector<std::vector<Cplx>> const& A, bool from_left,
                              std::vector<Cplx> const& fp,
                              std::vector<Cplx> const& rhs)
            {
            auto act = make_reg(A,from_left,fp);
            if (D*D > vms_dense_solve_max_)
                {
                bool ok=false;
                auto sol = vx_bicgstab(act,rhs,D*D,1e-12,20*D*D,ok);
                if (ok) return sol;
                // Fall through to the dense path rather than trust a
                // non-converged iterate -- never less robust than before.
                }
            auto Mat = vx_build_linear_map(D,act);
            return vx_solve(Mat,D*D,rhs);
            };
        GL0[_vms_F] = solve_reg(AL,true,r_AL,rhs_l);
        GRlast[_vms_S] = solve_reg(AR,false,l_AR,rhs_r);

        VmsEnv env;
        env.GL.resize(n_uc);
        env.GR.resize(n_uc);
        env.GL[0] = GL0;
        for (int n=0;n+1<n_uc;++n)
            env.GL[n+1] = vms_push_site_left(env.GL[n],AL[n],rows[n],D);
        env.GR[n_uc-1] = GRlast;
        for (int n=n_uc-1;n>0;--n)
            env.GR[n-1] = vms_push_site_right(env.GR[n],AR[n],rows[n],D);
        env.e_cell = 0.5*(e_L+e_R);
        return env;
        }

    // H_AC[n] applied to X (D,d,D):
    //   Y[L,o,R] = sum GL[a][l,L] X[l,i,r] W[a,i,o,b] GR[b][r,R]
    // -- every term (on-site, the two half-open bond terms, the
    // accumulated background on either side) in one contraction, because
    // the environments are channel-resolved.
    static std::vector<Cplx>
    vms_h_ac_action(std::vector<Cplx> const& X,
                     std::vector<std::vector<Cplx>> const& GLn,
                     std::vector<std::vector<Cplx>> const& GRn,
                     IdmrgAutomatonRow const& row, int D)
        {
        int nl = row.left_n, nr = row.right_n, d = row.d;
        std::vector<Cplx> Y((size_t)D*d*D,Cplx(0,0));
        std::vector<Cplx> slab((size_t)d*d);
        for (int a=0;a<nl;++a)
            {
            bool any_a = false;
            for (auto const& z : GLn[a]) if (z!=Cplx(0,0)) { any_a=true; break; }
            if (!any_a) continue;
            // XL[L,(i,r)] = sum_l GL[a][l,L] X[l,(i,r)], i.e. GL[a]^T @ X
            // with X read as a (D, d*D) matrix -- its own flat layout
            // already IS that, so no data movement beyond the transpose.
            std::vector<Cplx> GLt((size_t)D*D);
            for (int l=0;l<D;++l)
            for (int Lb=0;Lb<D;++Lb)
                GLt[Lb*D+l] = GLn[a][l*D+Lb];
            auto XL = vx_matmul(GLt,D,D,X,D,d*D);
            for (int b=0;b<nr;++b)
                {
                bool any_b = false;
                for (auto const& z : GRn[b]) if (z!=Cplx(0,0)) { any_b=true; break; }
                if (!any_b) continue;
                bool any_w = false;
                for (int si=0;si<d;++si)
                for (int so=0;so<d;++so)
                    {
                    slab[si*d+so] = row.flat[((size_t)a*nr+b)*d*d + si*d+so];
                    if (slab[si*d+so]!=Cplx(0,0)) any_w = true;
                    }
                if (!any_w) continue;
                auto XLW = vx_apply_op_ket(slab,XL,D,d);
                auto term = vx_cap_right(XLW,D,d,GRn[b]);
                for (size_t k=0;k<Y.size();++k) Y[k] += term[k];
                }
            }
        return Y;
        }

    // H_C[n] applied to the bond matrix C (D,D):
    //   Y[L,R] = sum_a GL[a][l,L] C[l,r] GR[a][r,R]
    // -- the zero-site effective Hamiltonian, a plain sum of outer
    // products because the same channel closes on both sides.
    static std::vector<Cplx>
    vms_h_c_action(std::vector<Cplx> const& C,
                    std::vector<std::vector<Cplx>> const& GL_bond,
                    std::vector<std::vector<Cplx>> const& GR_n, int D)
        {
        int nch = (int)std::min(GL_bond.size(),GR_n.size());
        std::vector<Cplx> Y((size_t)D*D,Cplx(0,0));
        for (int a=0;a<nch;++a)
            {
            bool any_l=false, any_r=false;
            for (auto const& z : GL_bond[a]) if (z!=Cplx(0,0)) { any_l=true; break; }
            if (!any_l) continue;
            for (auto const& z : GR_n[a]) if (z!=Cplx(0,0)) { any_r=true; break; }
            if (!any_r) continue;
            // Y += GL[a]^T C GR[a], as TWO matrix products rather than one
            // four-deep loop: the naive form is O(D^4) per channel and was
            // measured to dominate the whole solve (D=16 took 24s against
            // the pure-Python backend's 2s until this was staged).
            std::vector<Cplx> GLt((size_t)D*D);
            for (int l=0;l<D;++l)
            for (int Lb=0;Lb<D;++Lb)
                GLt[Lb*D+l] = GL_bond[a][l*D+Lb];
            auto T = vx_matmul(GLt,D,D,C,D,D);
            auto term = vx_matmul(T,D,D,GR_n[a],D,D);
            for (size_t k=0;k<Y.size();++k) Y[k] += term[k];
            }
        return Y;
        }

    // New (AL,AR) for ONE site from its AC and the bond matrices on either
    // side. Unlike the grouped single-site case -- where one C sits on both
    // sides of the only site, so the same polar factor serves twice -- the
    // two bonds are different matrices here, and reusing one for both
    // silently gauges the state wrong.
    static std::pair<std::vector<Cplx>,std::vector<Cplx>>
    vms_update_AL_AR(std::vector<Cplx> const& AC, int D, int d,
                      std::vector<Cplx> const& C_left,
                      std::vector<Cplx> const& C_right)
        {
        auto [U1,Vt1] = vx_economy_svd(AC,D*d,D);
        auto U_l = vx_matmul(U1,D*d,D,Vt1,D,D);
        auto [Ucr,Vtcr] = vx_economy_svd(C_right,D,D);
        auto U_cr = vx_matmul(Ucr,D,D,Vtcr,D,D);
        auto AL_new = vx_matmul(U_l,D*d,D,vx_dagger(U_cr,D,D),D,D);

        auto [U2,Vt2] = vx_economy_svd(AC,D,d*D);
        auto U_r = vx_matmul(U2,D,D,Vt2,D,d*D);
        auto [Ucl,Vtcl] = vx_economy_svd(C_left,D,D);
        auto U_cl = vx_matmul(Ucl,D,D,Vtcl,D,D);
        auto AR_new = vx_matmul(vx_dagger(U_cl,D,D),D,D,U_r,D,d*D);
        return {AL_new,AR_new};
        }

    // One sequential multi-site VUMPS attempt: build the environments once
    // for the whole cell, then sweep it solving H_AC[n] at every site and
    // H_C[n] at every bond, updating (AL[n],AR[n]) from the results.
    struct VmsRun
        {
        std::vector<std::vector<Cplx>> AL, AR, C, AC;
        VmsEnv env;
        double e_cell = 0.0, mismatch = 0.0;
        bool converged = false;
        int niter = 0;
        };

    VmsRun
    vms_single_run(std::vector<IdmrgAutomatonRow> const& rows, int D,
                    double tol, int maxiter, int niter_lanczos,
                    VmsRun const* init, std::mt19937_64& rng) const
        {
        int n_uc = (int)rows.size();
        std::vector<std::vector<Cplx>> AL(n_uc), AR(n_uc), C(n_uc), AC(n_uc);
        if (init)
            { AL = init->AL; AR = init->AR; C = init->C; }
        else
            for (int n=0;n<n_uc;++n)
                {
                auto st = vumps_random_init(D,rows[n].d,rng);
                AL[n] = st.AL; AR[n] = st.AR; C[n] = st.C;
                }
        for (int n=0;n<n_uc;++n)
            AC[n] = vumps_compose_AL_C(AL[n],C[n],D,rows[n].d);

        bool converged = false;
        double mismatch = 0.0;
        int it = 0;
        VmsEnv env;
        for (it=0; it<maxiter; ++it)
            {
            env = vms_environments(AL,AR,rows,D);
            std::vector<std::vector<Cplx>> AC_new(n_uc), C_new(n_uc);
            for (int n=0;n<n_uc;++n)
                {
                int d = rows[n].d;
                auto act_ac = [&](std::vector<Cplx> const& X)
                    { return vms_h_ac_action(X,env.GL[n],env.GR[n],rows[n],D); };
                AC_new[n] = vx_lanczos_ground_state(act_ac,AC[n],D*d*D,
                                                     niter_lanczos).second;
                // H_C[n] lives on the bond to the RIGHT of site n: its left
                // environment is the one left of site n+1, its right one is
                // GR[n].
                auto GL_bond = (n+1<n_uc) ? env.GL[n+1]
                    : vms_push_site_left(env.GL[n],AL[n],rows[n],D);
                auto act_c = [&](std::vector<Cplx> const& X)
                    { return vms_h_c_action(X,GL_bond,env.GR[n],D); };
                C_new[n] = vx_lanczos_ground_state(act_c,C[n],D*D,
                                                    niter_lanczos).second;
                }
            mismatch = 0.0;
            for (int n=0;n<n_uc;++n)
                mismatch = std::max(mismatch,
                    vumps_gauge_mismatch(AC_new[n],C_new[n],AL[n],AR[n],D,rows[n].d));
            AC = AC_new; C = C_new;
            for (int n=0;n<n_uc;++n)
                {
                auto [alnew,arnew] = vms_update_AL_AR(
                    AC[n],D,rows[n].d,C[(n-1+n_uc)%n_uc],C[n]);
                AL[n] = alnew; AR[n] = arnew;
                }
            if (mismatch < tol) { converged = true; break; }
            }

        // Refresh against the FINAL AL/AR -- what was built above reflects
        // this iteration's INPUT tensors, one update behind what is being
        // returned. Same reason the grouped single run does this.
        env = vms_environments(AL,AR,rows,D);
        VmsRun out;
        out.AL=AL; out.AR=AR; out.C=C; out.AC=AC; out.env=env;
        out.e_cell=env.e_cell; out.mismatch=mismatch; out.converged=converged;
        out.niter=std::min(it+1,maxiter);
        return out;
        }

    // Sequential multi-site VUMPS with the same restart/D-ramp discipline
    // the grouped driver uses, for the same measured reason (a pure random
    // start at large D lands in a bad basin often enough to matter).
    VmsRun
    vms_ground_state(std::vector<IdmrgAutomatonRow> const& rows, int D,
                      double tol, int maxiter, int niter_lanczos,
                      int nrestarts, std::mt19937_64& rng) const
        {
        std::vector<int> ramp;
        for (int d_cur=1; d_cur<D; d_cur*=2) ramp.push_back(d_cur);
        ramp.push_back(D);

        auto better = [](VmsRun const& a, VmsRun const& b)
            {
            if (a.converged != b.converged) return a.converged;
            return a.e_cell < b.e_cell;
            };
        VmsRun best; bool have_best=false;
        VmsRun prev; bool have_prev=false;
        for (int D_cur : ramp)
            {
            int n_here = (D_cur==D) ? nrestarts : std::min(nrestarts,3);
            VmsRun local; bool have_local=false;
            for (int attempt=0; attempt<n_here; ++attempt)
                {
                try
                    {
                    // Warm-starting across the ramp needs the previous
                    // rung's tensors re-embedded at the new D; only the
                    // same-D case is reused directly here, and a fresh
                    // random start is used otherwise -- the ramp still
                    // helps through the "beat the smaller D" check below.
                    bool reuse = (attempt==0 && have_prev
                                   && (int)prev.AL.size()==(int)rows.size()
                                   && prev.AL[0].size()==(size_t)D_cur*rows[0].d*D_cur);
                    auto r = vms_single_run(rows,D_cur,tol,maxiter,niter_lanczos,
                                             reuse ? &prev : nullptr, rng);
                    if (verbose_)
                        println("vumps_ms D=",D_cur," attempt=",attempt,": e_cell=",
                                r.e_cell," converged=",r.converged);
                    if (!have_local || better(r,local)) { local=r; have_local=true; }
                    }
                catch (ITError const& e)
                    {
                    if (verbose_)
                        println("vumps_ms D=",D_cur," attempt=",attempt,
                                ": failed (",e.what(),")");
                    }
                }
            if (!have_local)
                throw ITError("Chain::vms_ground_state: every attempt at D="+
                               std::to_string(D_cur)+" failed -- try increasing "
                               "nrestarts");
            prev = local; have_prev = true;
            best = local; have_best = true;
            }
        if (!have_best)
            throw ITError("Chain::vms_ground_state: no attempt succeeded");
        return best;
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

    // Above this, VUMPS's own H_AC/H_C ground-state solves go matrix-free
    // (vx_lanczos_ground_state); at or below it they are diagonalized
    // densely, exactly as this port always did. Same shape of split, and
    // the same reason for it, as ic_dense_eig_max_ above.
    //
    // The dense route costs O(n^3) for the zheev PLUS n calls to the
    // action to build the matrix at all, with n = D*d_g*D for H_AC. That
    // was written when this port's only validated models were TFIM/
    // Heisenberg at D<=3 with a ONE-site unit cell (d_g = 2 or 4, so
    // n <= 36 -- genuinely negligible, as vx_hermitian_ground_state's own
    // comment says). It is not negligible for a two-site unit cell of
    // native spinful sites: d_g = 4*4 = 16, so D=8 gives n = 1024 and
    // D=16 gives n = 4096, per outer iteration, per attempt, at every
    // step of the D-ramp. Measured on a (c,f) Kondo chain at D=8 in its
    // free limit, this backend had not returned an energy density after
    // 10 minutes at 100% CPU (the ramp had reached D=4), against 22.7 s
    // for the identical model through pyitensor's own VUMPS -- which
    // never forms these matrices, using _lanczos_ground_state on
    // _h_ac_action/_h_c_action instead. See
    // docs/known_issue_v3_vumps_no_return.md.
    //
    // Where to put the threshold follows from the two costs: dense is n
    // calls to the action (to build the matrix) plus O(n^3) for zheev,
    // while Lanczos is at most niter_lanczos calls plus an O(m^3) on a
    // tridiagonal of size m <= niter_lanczos. Dense therefore only wins
    // while n is of order niter_lanczos or below, which puts the natural
    // crossover in the tens, not the hundreds.
    //
    // 64 also keeps every case this port was originally validated on
    // (tests/test_vumps_v3.py, tests/test_vumps_correlator_v3.py) on the
    // exact dense path they were validated against -- those are TFIM/
    // Heisenberg with d_g = 2 or 4 at D <= 3, so n = D*d_g*D <= 36 -- and
    // so this is a pure extension to the regime that did not work at all
    // rather than a change to one that did.
    static constexpr int vx_dense_eig_max_ = 64;
    static constexpr int vx_lanczos_niter_ = 60;

    // Lowest eigenpair of a Hermitian operator given ONLY by its action on
    // a flat length-n complex vector -- a direct C++ port of pyitensor/
    // dmrg.py's own _lanczos_ground_state (Lanczos with full
    // reorthogonalization, stopping early once the lowest Ritz value
    // stabilizes to `tol` relative), which is what pyitensor/vumps.py uses
    // for exactly these two solves. The small tridiagonal eigenproblem is
    // handed to the existing dense solver: it is niter x niter at most, so
    // its own O(m^3) is irrelevant next to one call to the action.
    template <typename Fn>
    static std::pair<double,std::vector<Cplx>>
    vx_lanczos_ground_state(Fn&& action, std::vector<Cplx> v0, int n,
                             int niter, double tol=1e-12)
        {
        auto dot = [](std::vector<Cplx> const& a, std::vector<Cplx> const& b)
            {
            Cplx s(0,0);
            for (size_t i=0;i<a.size();++i) s += std::conj(a[i])*b[i];
            return s;
            };
        auto nrm = [&](std::vector<Cplx> const& a)
            { return std::sqrt(std::abs(dot(a,a))); };

        double beta0 = nrm(v0);
        if (beta0 <= 0.0)
            throw ITError("Chain::vx_lanczos_ground_state: zero start vector");
        for (auto& z : v0) z /= beta0;

        std::vector<std::vector<Cplx>> qs;
        qs.push_back(v0);
        std::vector<double> alphas, betas;

        auto w = action(qs[0]);
        double alpha = dot(qs[0],w).real();
        alphas.push_back(alpha);
        for (size_t i=0;i<w.size();++i) w[i] -= alpha*qs[0][i];

        // Lowest eigenvalue (and, on the iteration that returns, the
        // eigenvector) of the real symmetric tridiagonal built so far --
        // pyitensor's _tridiag_ground_value/_tridiag_ground_ritz.
        auto tridiag_ground = [&](bool want_vec, std::vector<Cplx>& vec)
            {
            int m = (int)alphas.size();
            std::vector<Cplx> T((size_t)m*m,Cplx(0,0));
            for (int i=0;i<m;++i) T[(size_t)i*m+i] = Cplx(alphas[i],0);
            for (int i=0;i+1<m;++i)
                {
                T[(size_t)i*m+(i+1)] = Cplx(betas[i],0);
                T[(size_t)(i+1)*m+i] = Cplx(betas[i],0);
                }
            auto [ev,evec] = vx_hermitian_ground_state(T,m);
            if (want_vec) vec = evec;
            return ev;
            };

        std::vector<Cplx> dummy;
        double prev_eval = tridiag_ground(false,dummy);
        int m = std::min(niter,n);
        for (int step=1; step<m; ++step)
            {
            double beta = nrm(w);
            if (beta < tol) break;
            betas.push_back(beta);
            std::vector<Cplx> q_new(w.size());
            for (size_t i=0;i<w.size();++i) q_new[i] = w[i]/beta;
            qs.push_back(q_new);
            w = action(q_new);
            alpha = dot(q_new,w).real();
            alphas.push_back(alpha);
            auto const& q_prev = qs[qs.size()-2];
            for (size_t i=0;i<w.size();++i)
                w[i] -= alpha*q_new[i] + beta*q_prev[i];
            // Full reorthogonalization against every earlier Lanczos
            // vector, exactly as the Python original does -- the short
            // recurrence alone loses orthogonality fast enough to matter.
            for (size_t k=0;k+1<qs.size();++k)
                {
                Cplx c = dot(qs[k],w);
                for (size_t i=0;i<w.size();++i) w[i] -= c*qs[k][i];
                }
            double cur_eval = tridiag_ground(false,dummy);
            if (std::abs(cur_eval-prev_eval) < tol*std::max(1.0,std::abs(cur_eval)))
                break;
            prev_eval = cur_eval;
            }

        std::vector<Cplx> yvec;
        double eval = tridiag_ground(true,yvec);
        std::vector<Cplx> out(n,Cplx(0,0));
        for (size_t k=0;k<qs.size();++k)
            for (int i=0;i<n;++i)
                out[i] += yvec[k]*qs[k][i];
        return {eval,out};
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
                      double tol, int maxiter, int niter_lanczos,
                      VumpsInit const* init, std::mt19937_64& rng) const
        {
        VumpsInit start = init ? *init : vumps_random_init(D,d_g,rng);
        auto AL = start.AL, AR = start.AR, C = start.C;
        auto AC = vumps_compose_AL_C(AL,C,D,d_g);

        bool converged=false; double mismatch=0.0; int it=0;
        VumpsEnv env;
        for (it=0; it<maxiter; ++it)
            {
            env = vumps_build_environments(AL,AR,D,d_g,h1,pending);

            // Dense diagonalization below vx_dense_eig_max_, matrix-free
            // Lanczos above it -- see that constant's own comment for why
            // the split exists and why the threshold sits where it does.
            // The Lanczos start vector is the CURRENT AC (resp. C), which
            // is what pyitensor/vumps.py passes too: successive outer
            // iterations move the solution only a little, so a warm start
            // converges in a handful of Krylov vectors.
            int n_ac = D*d_g*D;
            std::vector<Cplx> AC_new;
            if (n_ac <= vx_dense_eig_max_)
                {
                auto HAC = vumps_build_h_ac_dense(D,d_g,env.GL,env.GR,env.bond_envs,h1);
                AC_new = vx_hermitian_ground_state(HAC,n_ac).second;
                }
            else
                {
                auto act = [&](std::vector<Cplx> const& X)
                    { return vumps_h_ac_action(X,D,d_g,env.GL,env.GR,env.bond_envs,h1); };
                AC_new = vx_lanczos_ground_state(act,AC,n_ac,niter_lanczos).second;
                }
            int n_c = D*D;
            std::vector<Cplx> C_new;
            if (n_c <= vx_dense_eig_max_)
                {
                auto HC = vumps_build_h_c_dense(D,env.GL,env.GR,env.bond_envs);
                C_new = vx_hermitian_ground_state(HC,n_c).second;
                }
            else
                {
                auto act = [&](std::vector<Cplx> const& X)
                    { return vumps_h_c_action(X,D,env.GL,env.GR,env.bond_envs); };
                C_new = vx_lanczos_ground_state(act,C,n_c,niter_lanczos).second;
                }

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

    // (chi_l,d,chi_r) dense row-major array for one of svd()'s own factors,
    // whose left (for U) or right (for V) bond may be absent entirely at
    // the very first micro-steps -- a missing bond becomes a length-1 axis,
    // exactly as pyitensor/idmrg.py's own _u_array_lpr/_v_array_lpr insert
    // one by reshape. Kept separate from idmrg_tensor_to_lpr_array (which
    // requires all three Index objects to exist) for that reason alone.
    static std::vector<Cplx>
    idmrg_svd_factor_lpr(ITensor const& T, Index const& left, bool have_left,
                          Index const& phys, Index const& right, bool have_right)
        {
        int d = dim(phys);
        int chi_l = have_left ? dim(left) : 1;
        int chi_r = have_right ? dim(right) : 1;
        std::vector<Cplx> out((size_t)chi_l*d*chi_r,Cplx(0,0));
        for (int l=1;l<=chi_l;++l)
        for (int s=1;s<=d;++s)
        for (int r=1;r<=chi_r;++r)
            {
            Cplx v(0,0);
            if (have_left && have_right) v = eltC(T,left(l),phys(s),right(r));
            else if (have_left)          v = eltC(T,left(l),phys(s));
            else if (have_right)         v = eltC(T,phys(s),right(r));
            else                         v = eltC(T,phys(s));
            out[((size_t)(l-1)*d+(s-1))*chi_r+(r-1)] = v;
            }
        return out;
        }

    // Relative floor below which a previous micro-step's singular value is
    // treated as zero (its inverse set to 0) when building the wavefunction
    // prediction or the theta unit cell -- the direct analogue of ITensor's
    // own iDMRG `InverseCut` argument (detail::PseudoInvert, default 1e-8,
    // applied to the pseudo-inverted center matrix in
    // mpscpp2/ITensor/itensor/mps/idmrg.h), and the same constant
    // pyitensor/idmrg.py carries as _PREDICTION_INVERSE_CUT. Taken relative
    // to the largest singular value, since svd() here returns raw
    // (un-normalized) singular values.
    static constexpr double idmrg_inverse_cut_ = 1e-8;

    // lam^-1 with everything at or below idmrg_inverse_cut_*max(lam) sent
    // to 0 -- shared by idmrg_wavefunction_prediction and idmrg_theta_cell.
    static std::vector<double>
    idmrg_pseudo_inverse(std::vector<double> const& lam)
        {
        double top = 0.0;
        for (double v : lam) top = std::max(top,v);
        double floor = idmrg_inverse_cut_*top;
        std::vector<double> out(lam.size(),0.0);
        for (size_t i=0;i<lam.size();++i)
            if (lam[i] > floor) out[i] = 1.0/lam[i];
        return out;
        }

    // McCulloch's iDMRG wavefunction prediction: the trial two-site tensor
    // for the *next* micro-step, built from the previous micro-step's own
    // SVD factors, in that step's own bond bases -- C++ port of
    // pyitensor/idmrg.py's _wavefunction_prediction, see that function's
    // own docstring for the derivation:
    //
    //     theta = lambda_k . B_k . lambda_{k-1}^-1 . A_k . lambda_k
    //
    // with A_k=U, B_k=V, lambda_k=S from micro-step k's own SVD. This is
    // not a heuristic warm start but the exact translation of the converged
    // state into the new step's enlarged bases, and it is what makes
    // successive macro-iterations gauge-compatible (the ingredient this
    // port was missing entirely -- see idmrg_ground_state's own comment).
    // Note the swap: the previous step's right-canonical V supplies the new
    // LEFT physical slot and its left-canonical U the new RIGHT one (the
    // reference implementation's own swapUnitCells, idmrg.h), which lines
    // the sublattice labels up for n_uc<=2 and would not for n_uc>=3.
    //
    // Returns the flat vector positionally aligned to
    // (HL_ket,phys_L,phys_R,HR_ket) -- idmrg_tensor_to_flat4's own layout,
    // hence directly usable as idmrg_local_solve's `warm` argument -- or an
    // empty vector when the previous step's shapes don't line up with the
    // current local problem (the early micro-steps, and while the bond
    // dimension is still growing).
    static std::vector<Cplx>
    idmrg_wavefunction_prediction(std::vector<double> const& s_last,
                                   std::vector<Cplx> const& V_last, int Vl, int Vd, int Vr,
                                   std::vector<double> const& s_prev,
                                   std::vector<Cplx> const& U_last, int Ul, int Ud, int Ur,
                                   int shape0, int shape1, int shape2, int shape3)
        {
        int ns_last = (int)s_last.size(), ns_prev = (int)s_prev.size();
        if (Vr != ns_prev || Ul != ns_prev) return {};
        if (Vl != ns_last || Ur != ns_last) return {};
        if (shape0 != ns_last || shape1 != Vd || shape2 != Ud || shape3 != ns_last)
            return {};
        auto s_inv = idmrg_pseudo_inverse(s_prev);
        std::vector<Cplx> theta((size_t)ns_last*Vd*Ud*ns_last,Cplx(0,0));
        double norm2 = 0.0;
        for (int i=0;i<ns_last;++i)
        for (int s=0;s<Vd;++s)
        for (int r=0;r<Ud;++r)
        for (int j=0;j<ns_last;++j)
            {
            Cplx acc(0,0);
            for (int b=0;b<ns_prev;++b)
                acc += V_last[((size_t)i*Vd+s)*Vr+b]*s_inv[b]*U_last[((size_t)b*Ud+r)*Ur+j];
            acc *= s_last[i]*s_last[j];
            theta[(((size_t)i*Vd+s)*Ud+r)*ns_last+j] = acc;
            norm2 += std::norm(acc);
            }
        double norm = std::sqrt(norm2);
        if (!std::isfinite(norm) || norm == 0.0) return {};
        for (auto & v : theta) v /= norm;
        return theta;
        }

    // The gauge-consistent two-site unit cell, extracted from ONE
    // micro-step's own solved theta -- C++ port of pyitensor/idmrg.py's
    // _theta_cell. See that function's own docstring for why the
    // per-micro-step idmrg_U_ chain (whose two ends live in bond bases
    // minted by *different* micro-steps) cannot be tiled instead, and for
    // the measured size of the error that mistake produces.
    //
    // Writing theta in Vidal form
    // theta = lambda_o . Gamma_L . lambda_c . Gamma_R . lambda_o, its SVD
    // gives U = lambda_o.Gamma_L and S.V = lambda_c.Gamma_R.lambda_o, so
    // the two left-canonical site tensors are A_1 = U and
    // A_2 = S . V . lambda_o^-1, with lambda_o (`lam_outer`) the Schmidt
    // values on theta's own outer bond -- i.e. the *previous* micro-step's
    // center matrix, the same one idmrg_wavefunction_prediction inverts.
    // A_2's right bond is then back in A_1's left basis, so the cell tiles.
    //
    // Returns false (leaving a1/a2 untouched) if the shapes don't line up.
    static bool
    idmrg_theta_cell(std::vector<Cplx> const& U_lpr, int Ul, int Ud, int Ur,
                      std::vector<double> const& svals,
                      std::vector<Cplx> const& V_lpr, int Vl, int Vd, int Vr,
                      std::vector<double> const& lam_outer,
                      std::vector<Cplx>& a1, std::vector<Cplx>& a2)
        {
        if (Vr != (int)lam_outer.size() || Vr != Ul) return false;
        if (Vl != (int)svals.size()) return false;
        auto lam_inv = idmrg_pseudo_inverse(lam_outer);
        a1 = U_lpr;
        a2.assign((size_t)Vl*Vd*Vr,Cplx(0,0));
        for (int l=0;l<Vl;++l)
        for (int p=0;p<Vd;++p)
        for (int r=0;r<Vr;++r)
            a2[((size_t)l*Vd+p)*Vr+r] = V_lpr[((size_t)l*Vd+p)*Vr+r]*svals[l]*lam_inv[r];
        (void)Ud; (void)Ur;
        return true;
        }

    // env[dst_chan] -= shift * env[src_chan] along env's own mpo axis --
    // C++ port of pyitensor/idmrg.py's _subtract_energy_baseline, itself
    // this algorithm's equivalent of the reference implementation's
    // HL += -energy*IL (mpscpp2/ITensor/itensor/mps/idmrg.h).
    //
    // No separate identity environment has to be accumulated for it: the
    // automaton's channel space already carries one, since chans[p] is
    // [S, F] + pending and a growing environment's S/F components *are* the
    // identity and energy accumulations respectively (hence HL[F] -= e*HL[S]
    // on the left, and the mirrored HR[S] -= e*HR[F] on the right). Channel
    // indices here are 0-based, ITensor's own are 1-based, hence the +1.
    //
    // Subtracting a per-site baseline keeps the superblock eigenvalue
    // bounded instead of growing like 2*n_uc*k*e0 -- which matters because
    // the Arnoldi solve stops on a *relative* criterion, so the absolute
    // error in the eigenvalue (and so in the finite-difference energy
    // density) otherwise grows linearly with the macro-iteration count.
    static ITensor
    idmrg_subtract_energy_baseline(ITensor const& env, Index const& mpo_idx,
                                    int dst_chan, int src_chan, double shift)
        {
        ITensor src = env*setElt(mpo_idx(src_chan+1)); // mpo axis projected onto src_chan
        return env - shift*(src*setElt(mpo_idx(dst_chan+1)));
        }

    // White's density-matrix perturbation for one side of this micro-step's
    // two-site wavefunction -- C++ port of pyitensor/idmrg.py's own
    // _noise_perturbed_split, itself modelled on ITensor's LocalOp::deltaRho
    // (mpscpp3/ITensor/itensor/mps/localop.h), including its Hermitization
    // and trace normalization.
    //
    // Returns the eigenvectors of  rho + noise * drho  as an isometry from
    // `group` onto a fresh bond Index, truncated by (cutoff, maxdim), where
    //     rho  = tr_other |theta><theta|
    //     drho = (env . theta . W)(env . theta . W)^dag
    // with the MPO's own crossing index left OPEN in drho, so contracting it
    // is precisely White's sum over channels. See _noise_perturbed_split's
    // docstring for why this is needed at all (a product state is an exact
    // fixed point of the growth loop that no amount of solver refinement can
    // escape) and why the channels, not H itself, are what opens new
    // directions.
    struct NoisyIso { ITensor U; Index bond; std::vector<double> svals; double purity; };

    NoisyIso
    idmrg_noisy_isometry(ITensor const& theta, std::vector<Index> const& group,
                          ITensor const& env, Index env_bra, Index env_ket,
                          bool have_env, ITensor const& W,
                          double noise, double cutoff, int maxdim) const
        {
        auto [cmb,ci] = combiner(IndexSet(group),{"Tags","cmb"});
        ITensor thetac = cmb*theta;
        ITensor rho = thetac*dag(prime(thetac,ci));
        // Trace-normalize before adding, so `noise` is a pure relative weight
        // regardless of theta's own normalization.
        Cplx trc = eltC(rho*delta(dag(ci),prime(ci)));
        if (std::abs(trc) > 1e-300) rho /= trc.real();
        // Purity of the NOISE-FREE reduced state: exactly 1 for a product
        // state, strictly less for any entangled one. This, not the enriched
        // basis, is what the noise schedule keys on -- see
        // idmrg_ground_state's own comment.
        // tr(rho^2): swapPrime turns rho into its own transpose, so the
        // product contracts BOTH legs at once and lands on the trace.
        double purity = eltC(rho*swapPrime(rho,0,1)).real();
        if (noise > 0.0 && W)
            {
            ITensor t = theta;
            if (have_env && env)
                {
                t *= env;
                t = replaceInds(t,{env_bra},{env_ket});
                }
            t *= W;
            t.noPrime();
            ITensor drho = cmb*t;
            drho = drho*dag(prime(drho,ci));
            // "Expedient to ensure drho is Hermitian", as ITensor's own
            // deltaRho puts it.
            drho = drho + dag(swapPrime(drho,0,1));
            drho /= 2.0;
            Cplx trd = eltC(drho*delta(dag(ci),prime(ci)));
            if (std::abs(trd) > 1e-16) { drho /= trd.real(); rho += noise*drho; }
            }
        ITensor Uc,D;
        diagHermitian(rho,Uc,D);
        Index dc = commonIndex(Uc,D);
        // diagHermitian does not truncate: sort descending and keep by the
        // same (cutoff, maxdim) rule svd() would have applied, with the
        // cutoff taken on the discarded WEIGHT so numerically-zero
        // directions are never retained.
        int n = dim(dc);
        std::vector<std::pair<double,int>> ev;
        for (int k=1;k<=n;++k) ev.push_back({eltC(D,dc(k),prime(dc)(k)).real(),k});
        std::sort(ev.begin(),ev.end(),
                   [](auto const& a, auto const& b){ return a.first > b.first; });
        double total = 0.0;
        for (auto const& e : ev) total += std::max(e.first,0.0);
        int keep = 0; double discarded = 0.0;
        for (int k=n-1;k>=0;--k)
            {
            double w = std::max(ev[k].first,0.0);
            if (total > 0.0 && (discarded+w)/total > cutoff) { keep = k+1; break; }
            discarded += w;
            }
        if (keep < 1) keep = 1;
        if (maxdim > 0 && keep > maxdim) keep = maxdim;
        Index nb(keep,TagSet("Link"));
        ITensor Uk(ci,nb);
        std::vector<double> svals(keep,0.0);
        for (int k=0;k<keep;++k)
            {
            for (int i=1;i<=dim(ci);++i)
                Uk.set(ci(i),nb(k+1),eltC(Uc,ci(i),dc(ev[k].second)));
            svals[k] = std::sqrt(std::max(ev[k].first,0.0));
            }
        return {dag(cmb)*Uk,nb,svals,purity};
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
    std::tuple<double,ITensor,ITensor,Index,Index,std::vector<Cplx>,
                std::vector<double>,ITensor,double>
    idmrg_local_solve(ITensor const& HL, ITensor const& W_pL, Index phys_L,
                       ITensor const& W_pR, Index phys_R, ITensor const& HR,
                       Index HL_bra, Index HL_ket, bool have_HL_ket,
                       Index HR_bra, Index HR_ket, bool have_HR_ket,
                       double cutoff, int maxdim, int krylovdim, int restarts,
                       std::vector<Cplx> const& warm, double noise) const
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
            {
            x0 = idmrg_flat4_to_tensor(warm,HL_ket,phys_L,phys_R,HR_ket);
            if (noise > 0.0)
                {
                // Perturbing the START VECTOR is the other half of the
                // noise, and on the product-state trap it is the half that
                // does the work. Enlarging the basis (the density-matrix
                // term in idmrg_noisy_isometry) is necessary but NOT
                // sufficient: a product state like the vacuum is an exact
                // eigenstate of the FULL Hamiltonian, so it stays an exact
                // eigenvector of the local effective Hamiltonian in any
                // basis, however enriched -- and a Krylov solve started
                // exactly on an eigenvector breaks down on its first step
                // and hands the same vector straight back. Confirmed on the
                // Python side, where the density-matrix term alone grew chi
                // away from 1 and still left the energy identically 0 for
                // every macro-iteration. sqrt(noise) so the perturbation's
                // WEIGHT is O(noise), matching the density-matrix term.
                double nx = norm(x0);
                if (nx > 0.0)
                    {
                    ITensor kick = randomITensorC(IndexSet(order_in));
                    double nk = norm(kick);
                    if (nk > 0.0) x0 = x0/nx + (std::sqrt(noise)/nk)*kick;
                    }
                }
            }
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
        ITensor U,S,V;
        Index new_bond_u, new_bond_v;
        std::vector<double> noisy_svals;
        double purity = 1.0;
        if (noise > 0.0)
            {
            // White's density-matrix perturbation instead of a plain SVD --
            // see idmrg_noisy_isometry's own comment. U and V come from two
            // INDEPENDENTLY perturbed density matrices (this loop extends
            // HL and HR in the same micro-step, unlike a finite sweep which
            // only ever needs one direction), so they are no longer an SVD
            // pair of theta -- which nothing here requires: extend_HL and
            // extend_HR each transform through one factor alone, and the two
            // bonds have always been independent.
            std::vector<Index> right_inds;
            right_inds.push_back(phys_R);
            if (have_HR_ket) right_inds.push_back(HR_ket);
            auto L = idmrg_noisy_isometry(theta,left_inds,HL,HL_bra,HL_ket,
                                           have_HL_ket,W_pL,noise,cutoff,maxdim);
            auto R = idmrg_noisy_isometry(theta,right_inds,HR,HR_bra,HR_ket,
                                           have_HR_ket,W_pR,noise,cutoff,maxdim);
            // V is the RIGHT-orthonormal factor, i.e. the conjugate of the
            // isometry that maps the right group onto its own new bond.
            U = L.U; new_bond_u = L.bond;
            V = dag(R.U); new_bond_v = R.bond;
            noisy_svals = L.svals;
            purity = L.purity;
            }
        else
            {
            std::tie(U,S,V) = svd(theta,IndexSet(left_inds),{"Cutoff",cutoff,"MaxDim",maxdim});
            new_bond_u = commonIndex(U,S);
            new_bond_v = commonIndex(S,V);
            }
        // The singular values themselves -- discarded by the growth loop's
        // own HL/HR extension (see this method's own comment above), but
        // needed verbatim by McCulloch's wavefunction prediction and by the
        // gauge-consistent unit-cell extraction (idmrg_wavefunction_prediction/
        // idmrg_theta_cell), both of which work in Vidal's Gamma/lambda form.
        // C++ analogue of pyitensor/idmrg.py's own _svd_singular_values.
        std::vector<double> svals;
        if (noise > 0.0)
            svals = noisy_svals;   // sqrt of rho's retained eigenvalues
        else
            {
            svals.assign((size_t)dim(new_bond_u),0.0);
            for (int k=1;k<=dim(new_bond_u);++k)
                svals[k-1] = eltC(S,new_bond_u(k),new_bond_v(k)).real();
            }
        if (noise <= 0.0)
            {
            // Same purity the noisy branch reports, read off the Schmidt
            // values the SVD already produced.
            double tot=0.0, sq=0.0;
            for (double sv : svals) { double w=sv*sv; tot+=w; sq+=w*w; }
            purity = (tot>0.0) ? sq/(tot*tot) : 1.0;
            }
        return {energy,U,V,new_bond_u,new_bond_v,theta_flat,svals,theta,purity};
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

    // -- iDMRG static-correlator private helpers (the "ic_" prefix) --
    //
    // Dense, LAPACK-backed port of pyitensor/idmrg.py's own transfer-matrix
    // machinery (_transfer_matrices/_compose/_apply_transfer/
    // _dominant_*_fixed_point/_all_*_fixed_points/_expectation), operating
    // on the gauge-consistent 2-site unit cell idmrg_theta_cell extracts
    // (idmrg_cell_ below) rather than on ITensor objects -- the same
    // design choice, and for the same reasons, as the VUMPS "vx_" helpers
    // above (see IdmrgResult/VumpsResult's own struct comments): the
    // objects here are (chi,d,chi) arrays with chi at most maxm, every
    // contraction is a plain dense matmul, and ITensor's Index identity
    // machinery buys nothing while costing exactly the kind of duplicate-
    // index bookkeeping idmrg_ground_state's own top comment describes.
    //
    // The one thing these cannot borrow from the vx_ helpers is squareness:
    // VUMPS works at a single fixed bond dimension D throughout, whereas
    // the theta cell's two bonds genuinely differ (the outer bond is the
    // environment bond theta was solved against, chi_o; the inner one is
    // its own SVD bond, chi_c = min(chi_o*d_L, d_R*chi_o, maxm)). Every
    // helper below is therefore rectangular-aware, indexed by a `chis`
    // array with chis[k] = cell position k's own LEFT bond dimension and
    // chis[n_cell] = chis[0] (the cell wraps).
    //
    // A transfer tensor for cell position k is stored as the
    // (chis[k]^2, chis[k+1]^2) row-major matrix E[(l,L),(r,R)] --
    // idmrg.py's own einsum('lpr,LpR->lLrR', A, conj(A)) flattened
    // exactly as vx_op_transfer_matrix already flattens its square case,
    // so "compose two sites in a row" is an ordinary matrix product.

    // E for cell tensor A (chi_l,d,chi_r), with operator matrix M
    // ((d,d) row-major, (in,out) convention as idmrg_op_dense returns)
    // optionally applied to the ket side first -- idmrg.py's
    // _op_transfer_mat (M empty = the plain transfer tensor).
    static std::vector<Cplx>
    ic_transfer(std::vector<Cplx> const& A, int chi_l, int d, int chi_r,
                 std::vector<Cplx> const& M)
        {
        std::vector<Cplx> ket = A;
        if (!M.empty())
            {
            ket.assign((size_t)chi_l*d*chi_r,Cplx(0,0));
            for (int l=0;l<chi_l;++l)
            for (int o=0;o<d;++o)
            for (int r=0;r<chi_r;++r)
                {
                Cplx acc(0,0);
                for (int i=0;i<d;++i) acc += M[(size_t)i*d+o]*A[((size_t)l*d+i)*chi_r+r];
                ket[((size_t)l*d+o)*chi_r+r] = acc;
                }
            }
        std::vector<Cplx> E((size_t)chi_l*chi_l*chi_r*chi_r,Cplx(0,0));
        for (int l=0;l<chi_l;++l)
        for (int Lb=0;Lb<chi_l;++Lb)
        for (int r=0;r<chi_r;++r)
        for (int R=0;R<chi_r;++R)
            {
            Cplx acc(0,0);
            for (int p=0;p<d;++p)
                acc += ket[((size_t)l*d+p)*chi_r+r]*std::conj(A[((size_t)Lb*d+p)*chi_r+R]);
            E[((size_t)l*chi_l+Lb)*((size_t)chi_r*chi_r)+((size_t)r*chi_r+R)] = acc;
            }
        return E;
        }

    // One site's transfer tensor applied to a boundary matrix, WITHOUT
    // ever forming the (chi_l^2, chi_r^2) transfer matrix itself:
    //
    //   out[l,L] = sum_{p,r,R} (M^T A)[l,p,r] rho[r,R] conj(A)[L,p,R]
    //
    // A is the cell tensor (chi_l,d,chi_r) row-major; M its (d,d)
    // (in,out)-convention operator matrix (empty = the plain transfer).
    // rho is (chi_r,chi_r), the result (chi_l,chi_l), both row-major.
    //
    // Same quantity ic_matvec(ic_transfer(A,...,M), ..., rho) computes,
    // re-associated: building the transfer tensor costs O(chi^4 d) and
    // applying it another O(chi^4), whereas contracting rho in first makes
    // both halves O(chi^3 d). Exact, not approximate -- only the
    // contraction order changes. Mirrors pyitensor/idmrg.py's own
    // _apply_site_transfer, which is where this was measured first: with
    // the transfer matrix formed per site, this backend's r=1..7 sweep on
    // a maxm=32 chain cost 0.19s against the Python side's 0.0022s, i.e.
    // the "fast" backend was 85x slower at the one thing both had just
    // been optimized for.
    static std::vector<Cplx>
    ic_apply_site_transfer(std::vector<Cplx> const& A, int chi_l, int d, int chi_r,
                            std::vector<Cplx> const& M, std::vector<Cplx> const& rho)
        {
        // Aop[l,o,r] = sum_i M[i,o] A[l,i,r]   (M empty -> Aop = A)
        std::vector<Cplx> Aop;
        if (M.empty()) Aop = A;
        else
            {
            Aop.assign((size_t)chi_l*d*chi_r,Cplx(0,0));
            for (int l=0;l<chi_l;++l)
            for (int o=0;o<d;++o)
            for (int r=0;r<chi_r;++r)
                {
                Cplx acc(0,0);
                for (int i=0;i<d;++i) acc += M[(size_t)i*d+o]*A[((size_t)l*d+i)*chi_r+r];
                Aop[((size_t)l*d+o)*chi_r+r] = acc;
                }
            }
        // tmp[l,p,R] = sum_r Aop[l,p,r] rho[r,R]
        std::vector<Cplx> tmp((size_t)chi_l*d*chi_r,Cplx(0,0));
        for (int l=0;l<chi_l;++l)
        for (int p=0;p<d;++p)
        for (int R=0;R<chi_r;++R)
            {
            Cplx acc(0,0);
            for (int r=0;r<chi_r;++r)
                acc += Aop[((size_t)l*d+p)*chi_r+r]*rho[(size_t)r*chi_r+R];
            tmp[((size_t)l*d+p)*chi_r+R] = acc;
            }
        // out[l,L] = sum_{p,R} tmp[l,p,R] conj(A[L,p,R])
        std::vector<Cplx> out((size_t)chi_l*chi_l,Cplx(0,0));
        for (int l=0;l<chi_l;++l)
        for (int Lb=0;Lb<chi_l;++Lb)
            {
            Cplx acc(0,0);
            for (int p=0;p<d;++p)
            for (int R=0;R<chi_r;++R)
                acc += tmp[((size_t)l*d+p)*chi_r+R]
                        *std::conj(A[((size_t)Lb*d+p)*chi_r+R]);
            out[(size_t)l*chi_l+Lb] = acc;
            }
        return out;
        }

    // M (m,n) row-major times v (n) -- rectangular counterpart of
    // vx_matvec_row, i.e. idmrg.py's own _apply_transfer.
    static std::vector<Cplx>
    ic_matvec(std::vector<Cplx> const& M, int m, int n, std::vector<Cplx> const& v)
        {
        std::vector<Cplx> out((size_t)m,Cplx(0,0));
        for (int i=0;i<m;++i)
            {
            Cplx acc(0,0);
            for (int j=0;j<n;++j) acc += M[(size_t)i*n+j]*v[j];
            out[i] = acc;
            }
        return out;
        }

    // v (m) times M (m,n) row-major -- rectangular counterpart of
    // vx_vecmat_row, i.e. idmrg.py's own _apply_transfer_from_left.
    static std::vector<Cplx>
    ic_vecmat(std::vector<Cplx> const& v, std::vector<Cplx> const& M, int m, int n)
        {
        std::vector<Cplx> out((size_t)n,Cplx(0,0));
        for (int j=0;j<n;++j)
            {
            Cplx acc(0,0);
            for (int i=0;i<m;++i) acc += v[i]*M[(size_t)i*n+j];
            out[j] = acc;
            }
        return out;
        }

    // trace of a flattened (chi,chi) row-major matrix.
    static Cplx
    ic_trace(std::vector<Cplx> const& rho, int chi)
        {
        Cplx tr(0,0);
        for (int i=0;i<chi;++i) tr += rho[(size_t)i*chi+i];
        return tr;
        }

    // The full unit-cell transfer matrix T = E_0 @ E_1 @ ... @ E_{n-1},
    // a (chis[0]^2, chis[0]^2) matrix since the cell wraps.
    static std::vector<Cplx>
    ic_full_transfer(std::vector<std::vector<Cplx>> const& Es,
                      std::vector<int> const& chis)
        {
        int n = (int)Es.size();
        std::vector<Cplx> T = Es[0];
        int rows = chis[0]*chis[0], cols = chis[1]*chis[1];
        for (int p=1;p<n;++p)
            {
            int next = chis[p+1]*chis[p+1];
            T = vx_matmul(T,rows,cols,Es[p],cols,next);
            cols = next;
            }
        return T;
        }

    // Apply the whole unit cell's transfer map to a flat (chis[0]^2)
    // vector WITHOUT ever forming the (chis[0]^2, chis[0]^2) product --
    // one ic_matvec/ic_vecmat per cell position instead, i.e. O(n*chi^4)
    // rather than the O(chi^6) ic_full_transfer costs just to build T.
    // `right`: the direction ic_dominant_fixed_point's own "right" case
    // needs (E_{n-1} applied first, mirroring idmrg.py's own
    // reversed(Es) + _apply_transfer); otherwise the left/forward one.
    static std::vector<Cplx>
    ic_apply_chain(std::vector<std::vector<Cplx>> const& Es,
                    std::vector<int> const& chis, std::vector<Cplx> const& x,
                    bool right)
        {
        int n = (int)Es.size();
        std::vector<Cplx> out = x;
        if (right)
            for (int p=n-1;p>=0;--p)
                out = ic_matvec(Es[p],chis[p]*chis[p],chis[p+1]*chis[p+1],out);
        else
            for (int p=0;p<n;++p)
                out = ic_vecmat(out,Es[p],chis[p]*chis[p],chis[p+1]*chis[p+1]);
        return out;
        }

    // All eigenvalues and right eigenvectors of a general (non-Hermitian)
    // complex n x n row-major matrix, via itensor::zgeev_wrapper --
    // vx_dominant_eig's own decomposition step, factored out because the
    // small Hessenberg eigenproblem inside ic_arnoldi_dominant below needs
    // the FULL spectrum (evecs_col stays column-major, matching LAPACK's
    // own output directly).
    static void
    ic_geev_all(std::vector<Cplx> const& A_row, int n,
                 std::vector<Cplx>& evals, std::vector<Cplx>& evecs_col)
        {
        std::vector<Cplx> Acol((size_t)n*n);
        for (int i=0;i<n;++i)
        for (int j=0;j<n;++j)
            Acol[i+(size_t)j*n] = A_row[(size_t)i*n+j];
        evals.assign(n,Cplx(0,0));
        evecs_col.assign((size_t)n*n,Cplx(0,0));
        std::vector<Cplx> vl(1,Cplx(0,0));
        auto info = zgeev_wrapper('N','V',n,Acol.data(),evals.data(),
                                   vl.data(),evecs_col.data());
        if (info!=0)
            throw ITError("Chain::idmrg static correlators: zgeev_wrapper failed (info="+
                           std::to_string(info)+")");
        }

    // Transfer-matrix eigenproblems (chi^2 x chi^2) at or below this size
    // are solved densely; above it, iteratively (ic_arnoldi_dominant).
    // Same threshold, and the same reason for it, as pyitensor/idmrg.py's
    // own _DENSE_EIG_MAX: the dense route costs O(chi^6) TWICE over
    // (ic_full_transfer builds T at chi^6, then zgeev on it is another
    // chi^6), which at chi=24 already dominated an entire iDMRG solve --
    // measured directly, ~5s per vev()/correlator() call before this split
    // existed, against ~0.4s at chi=16.
    static constexpr int ic_dense_eig_max_ = 64;
    static constexpr int ic_arnoldi_krylov_ = 40;
    static constexpr int ic_arnoldi_restarts_ = 20;
    // Magnitude window within which two of Arnoldi's RITZ values count as
    // peripheral-tied. Much looser than vx_perron_tie_rtol_ on purpose:
    // that constant is applied to exact zgeev eigenvalues, these are
    // Krylov approximations, and a subdominant branch's Ritz magnitude is
    // routinely off by far more than 1e-6 while still being the same
    // peripheral eigenvalue physically.
    static constexpr double ic_arnoldi_tie_rtol_ = 1e-3;

    // The two largest-|lambda| eigenvalues, plus the dominant eigenvector,
    // of a linear map given only by its action on a flat length-n complex
    // vector -- a restarted Arnoldi, the matrix-free counterpart of
    // vx_dominant_eig (which forms and decomposes the whole matrix). The
    // runner-up eigenvalue is returned because the caller's
    // near-degeneracy guard needs it, exactly as pyitensor/idmrg.py asks
    // ARPACK for k=2 rather than k=1 for the same reason.
    //
    // The start vector is the identity, not a random one -- deterministic
    // (the dense route this replaces is, and every downstream
    // normalization diagnostic would otherwise vary run to run), and close
    // to the fixed point of an approximately-canonical cell anyway. Same
    // choice pyitensor makes for its own ARPACK v0.
    static void
    ic_arnoldi_dominant(std::function<std::vector<Cplx>(std::vector<Cplx> const&)> const& act,
                         int n, Cplx& eval0, Cplx& eval1, std::vector<Cplx>& vec0)
        {
        int m = std::min(ic_arnoldi_krylov_,n);
        if (m < 2) m = std::min(2,n);
        int chi = (int)std::lround(std::sqrt((double)n));
        std::vector<Cplx> v((size_t)n,Cplx(0,0));
        if (chi*chi == n) for (int i=0;i<chi;++i) v[(size_t)i*chi+i] = Cplx(1,0);
        else v.assign((size_t)n,Cplx(1,0));
        auto normalize = [](std::vector<Cplx>& x)
            {
            double s = 0.0;
            for (auto const& e : x) s += std::norm(e);
            s = std::sqrt(s);
            if (s > 0.0) for (auto & e : x) e /= s;
            return s;
            };
        normalize(v);
        eval0 = eval1 = Cplx(0,0);

        for (int rs=0; rs<=ic_arnoldi_restarts_; ++rs)
            {
            std::vector<std::vector<Cplx>> V;
            V.reserve(m);
            V.push_back(v);
            std::vector<Cplx> H((size_t)m*m,Cplx(0,0)); // row-major
            int mused = m;
            for (int j=0;j<m;++j)
                {
                auto w = act(V[j]);
                // Two Gram-Schmidt passes: one is not enough to keep the
                // Krylov basis orthogonal for a transfer matrix, whose
                // spectrum is strongly graded (the dominant eigenvalue sits
                // at ~1 and the rest decay), exactly the regime classical
                // Arnoldi loses orthogonality in.
                for (int pass=0;pass<2;++pass)
                for (int i=0;i<=j;++i)
                    {
                    Cplx c(0,0);
                    for (int k=0;k<n;++k) c += std::conj(V[i][k])*w[k];
                    H[(size_t)i*m+j] += c;
                    for (int k=0;k<n;++k) w[k] -= c*V[i][k];
                    }
                double beta = 0.0;
                for (auto const& e : w) beta += std::norm(e);
                beta = std::sqrt(beta);
                if (j+1 < m)
                    {
                    if (beta < 1e-13) { mused = j+1; break; } // invariant subspace
                    H[(size_t)(j+1)*m+j] = Cplx(beta,0);
                    for (auto & e : w) e /= beta;
                    V.push_back(std::move(w));
                    }
                }

            std::vector<Cplx> Hs((size_t)mused*mused);
            for (int i=0;i<mused;++i)
            for (int j=0;j<mused;++j)
                Hs[(size_t)i*mused+j] = H[(size_t)i*m+j];
            std::vector<Cplx> evals, evecs_col;
            ic_geev_all(Hs,mused,evals,evecs_col);
            std::vector<int> order(mused);
            for (int i=0;i<mused;++i) order[i]=i;
            std::sort(order.begin(),order.end(),
                       [&](int a,int b){ return std::abs(evals[a]) > std::abs(evals[b]); });
            // Perron selection happens HERE, not at the call site: only one
            // eigenvector is reconstructed and returned, so a caller handed
            // the magnitude-dominant Ritz pair has no way to recover the
            // +rho one if -rho edged ahead. See vx_perron_reorder.
            vx_perron_reorder(order,evals);
            // Among the peripheral (magnitude-tied) candidates, take the one
            // whose eigenvector actually is a FIXED POINT. Every caller of
            // this function wants a transfer-map fixed point and divides by
            // its trace, and that is exactly what separates the branches: a
            // period-2 state's -rho eigenvector carries opposite signs on
            // the two sublattices, so its trace cancels to ~0, while the
            // +rho (Perron) one is a positive matrix with an O(1) trace.
            //
            // Selecting on the eigenVALUE alone is not enough here. Ritz
            // values are approximate, so which of a +rho/-rho pair comes
            // out with the larger magnitude is partly numerical accident --
            // measured directly: on a half-filled/critical chain at D=16
            // this returned the -rho branch and the caller then died on
            // "dominant left fixed point has ~zero trace". The trace test
            // asks the question the caller actually cares about instead of
            // inferring it from an approximate eigenvalue.
            int chi_sq = (chi*chi == n) ? chi : 0;
            if (chi_sq > 0 && mused > 1)
                {
                double top = std::abs(evals[order[0]]);
                int ntied = 0;
                while (ntied < mused
                       && std::abs(evals[order[ntied]]) >= (1.0-ic_arnoldi_tie_rtol_)*top)
                    ++ntied;
                if (ntied > 1)
                    {
                    int pick = 0; double best = -1.0;
                    for (int c=0;c<ntied;++c)
                        {
                        std::vector<Cplx> xc((size_t)n,Cplx(0,0));
                        for (int i=0;i<mused;++i)
                            {
                            Cplx y = evecs_col[(size_t)i+(size_t)order[c]*mused];
                            for (int k=0;k<n;++k) xc[k] += y*V[i][k];
                            }
                        double nrm2 = 0.0;
                        for (auto const& e : xc) nrm2 += std::norm(e);
                        if (!(nrm2 > 0.0)) continue;
                        Cplx tr(0,0);
                        for (int i=0;i<chi_sq;++i) tr += xc[(size_t)i*chi_sq+i];
                        double score = std::abs(tr)/std::sqrt(nrm2);
                        if (score > best) { best = score; pick = c; }
                        }
                    if (pick != 0) std::swap(order[0],order[pick]);
                    }
                }
            eval0 = evals[order[0]];
            eval1 = mused>1 ? evals[order[1]] : Cplx(0,0);

            std::vector<Cplx> x((size_t)n,Cplx(0,0));
            for (int i=0;i<mused;++i)
                {
                Cplx y = evecs_col[(size_t)i+(size_t)order[0]*mused];
                for (int k=0;k<n;++k) x[k] += y*V[i][k];
                }
            if (normalize(x) == 0.0) break;
            v = std::move(x);

            auto Av = act(v);
            double res = 0.0;
            for (int k=0;k<n;++k) res += std::norm(Av[k]-eval0*v[k]);
            res = std::sqrt(res);
            if (res < 1e-11*std::max(1.0,std::abs(eval0))) break;
            }
        vec0 = std::move(v);
        }

    // (rho, eta): the dominant right (right=true) or left (right=false)
    // fixed point of the unit-cell transfer map, as a (chis[0],chis[0])
    // matrix normalized to trace 1 -- idmrg.py's own
    // _dominant_fixed_point, including its near-degeneracy guard (a single
    // dominant fixed point is not well defined when the leading eigenvalue
    // is (near-)degenerate; see _check_dominant_eigenvalue_nondegenerate's
    // own docstring for what that means physically). Dense below
    // ic_dense_eig_max_, matrix-free Arnoldi above it.
    static std::pair<std::vector<Cplx>,Cplx>
    ic_dominant_fixed_point(std::vector<std::vector<Cplx>> const& Es,
                             std::vector<int> const& chis, bool right)
        {
        int chi = chis[0], n = chi*chi;
        if (n <= ic_dense_eig_max_)
            {
            auto T = ic_full_transfer(Es,chis);
            return right ? vx_dominant_right_fixed_point(T,chi)
                          : vx_dominant_left_fixed_point(T,chi);
            }
        auto act = [&](std::vector<Cplx> const& x)
            { return ic_apply_chain(Es,chis,x,right); };
        Cplx eta(0,0), eta1(0,0);
        std::vector<Cplx> vec;
        ic_arnoldi_dominant(act,n,eta,eta1,vec);
        vx_check_perron_nondegenerate(eta,eta1,"Chain::idmrg static correlators");
        Cplx tr = ic_trace(vec,chi);
        if (std::abs(tr) < 1e-13)
            throw ITError("Chain::idmrg static correlators: dominant fixed point has "
                           "~zero trace -- degenerate/ill-defined normalization");
        for (auto & e : vec) e /= tr;
        return {vec,eta};
        }

    // rho_after[p] = the fixed-point "everything strictly after cell
    // position p, wrapping back around" density matrix ((chis[p+1],
    // chis[p+1]), trace 1) -- idmrg.py's own _all_right_fixed_points.
    static std::pair<std::vector<std::vector<Cplx>>,Cplx>
    ic_all_right_fixed_points(std::vector<std::vector<Cplx>> const& Es,
                               std::vector<int> const& chis)
        {
        int n = (int)Es.size();
        auto [rho_full,eta] = ic_dominant_fixed_point(Es,chis,true);
        std::vector<std::vector<Cplx>> rho_after(n);
        rho_after[n-1] = rho_full;
        auto cur = rho_full;
        for (int p=n-1;p>0;--p)
            {
            cur = ic_matvec(Es[p],chis[p]*chis[p],chis[p+1]*chis[p+1],cur);
            Cplx tr = ic_trace(cur,chis[p]);
            if (std::abs(tr) < 1e-13)
                throw ITError("Chain::idmrg static correlators: a right fixed-point "
                               "propagation step produced a ~zero-trace density matrix");
            for (auto & v : cur) v /= tr;
            rho_after[p-1] = cur;
            }
        return {rho_after,eta};
        }

    // rho_before[p] = the fixed point for "everything strictly before cell
    // position p, wrapping back around" ((chis[p],chis[p])), plus the
    // accumulated per-step trace factors divided out along the way
    // (`scales`, needed only by ic_canonicalize_cell -- see idmrg.py's own
    // _all_left_fixed_points for why G_left, unlike G_right, is not
    // invariant to that renormalization).
    static std::tuple<std::vector<std::vector<Cplx>>,Cplx,std::vector<double>>
    ic_all_left_fixed_points(std::vector<std::vector<Cplx>> const& Es,
                              std::vector<int> const& chis)
        {
        int n = (int)Es.size();
        auto [rho_full,eta] = ic_dominant_fixed_point(Es,chis,false);
        std::vector<std::vector<Cplx>> rho_before(n);
        std::vector<double> scales(n,1.0);
        rho_before[0] = rho_full;
        auto cur = rho_full;
        double scale = 1.0;
        for (int p=0;p<n-1;++p)
            {
            cur = ic_vecmat(cur,Es[p],chis[p]*chis[p],chis[p+1]*chis[p+1]);
            Cplx tr = ic_trace(cur,chis[p+1]);
            if (std::abs(tr) < 1e-13)
                throw ITError("Chain::idmrg static correlators: a left fixed-point "
                               "propagation step produced a ~zero-trace density matrix");
            for (auto & v : cur) v /= tr;
            scale *= tr.real();
            scales[p+1] = scale;
            rho_before[p+1] = cur;
            }
        return {rho_before,eta,scales};
        }

    // <...> for an operator string already composed into `running`, closed
    // against the left fixed point *at the string's own starting position*
    // rather than an implicit identity, and divided by the same
    // contraction with the operators dropped (`running_id`) -- direct port
    // of idmrg.py's own _expectation, including both of the details that
    // function's docstring flags as load-bearing (which left fixed point,
    // and that the denominator is the same contraction rather than a bare
    // trace).
    static Cplx
    ic_expectation(std::vector<std::vector<Cplx>> const& l_before,
                    int n_cell,
                    std::vector<Cplx> const& running,
                    std::vector<Cplx> const& running_id,
                    int nl, int nr,
                    std::vector<Cplx> const& rho_after_j, int k_start)
        {
        auto const& l = l_before[k_start % n_cell];
        auto num_v = ic_matvec(running,nl,nr,rho_after_j);
        auto den_v = ic_matvec(running_id,nl,nr,rho_after_j);
        Cplx num(0,0), den(0,0);
        for (int i=0;i<nl;++i) { num += l[i]*num_v[i]; den += l[i]*den_v[i]; }
        if (std::abs(den) < 1e-300)
            throw ITError("Chain::idmrg static correlators: zero normalization");
        return num/den;
        }

    // ic_expectation's own final closure, given the left fixed point at the
    // string's starting position and the two already-propagated boundary
    // matrices -- split out so callers that propagate site by site
    // (ic_apply_site_transfer) never build a composed transfer matrix just
    // to apply it once. Mirrors pyitensor/idmrg.py's _close_expectation.
    static Cplx
    ic_close_expectation(std::vector<Cplx> const& l, int nl,
                          std::vector<Cplx> const& x_num,
                          std::vector<Cplx> const& x_den)
        {
        Cplx num(0,0), den(0,0);
        for (int i=0;i<nl;++i) { num += l[i]*x_num[i]; den += l[i]*x_den[i]; }
        if (std::abs(den) < 1e-300)
            throw ITError("Chain::idmrg static correlators: zero normalization");
        return num/den;
        }

    // Build (once) the converged cell's transfer tensors and both families
    // of fixed points, and cache them on this Chain until the next
    // idmrg_ground_state() run invalidates them.
    //
    // pyitensor/idmrg.py recomputes all of this inside every single
    // onsite_expectation/two_point_correlator call (and _expectation
    // recomputes the left family a second time within the same call). That
    // is pure repetition -- none of it depends on which operator is being
    // measured, and the converged cell never changes between measurements
    // -- and it is not cheap: these are the transfer-matrix eigenproblems,
    // by far the dominant cost of a measurement. Caching is why a full
    // correlator sweep on this backend costs about one measurement's worth
    // of eigensolve rather than one per r.
    void
    ic_build_cache() const
        {
        if (ic_cache_valid_) return;
        int n = (int)idmrg_cell_.size();
        ic_chis_.assign(n+1,0);
        for (int k=0;k<n;++k) ic_chis_[k] = idmrg_cell_l_[k];
        ic_chis_[n] = idmrg_cell_l_[0];
        ic_Es_.resize(n);
        for (int k=0;k<n;++k)
            ic_Es_[k] = ic_transfer(idmrg_cell_[k],idmrg_cell_l_[k],
                                     idmrg_cell_d_[k],idmrg_cell_r_[k],{});
        auto rr = ic_all_right_fixed_points(ic_Es_,ic_chis_);
        ic_rho_after_ = std::move(rr.first);
        auto ll = ic_all_left_fixed_points(ic_Es_,ic_chis_);
        ic_l_before_ = std::move(std::get<0>(ll));
        ic_cache_valid_ = true;
        }

    // The window's own counterpart of ic_build_cache, on the RAW (un-re-
    // gauged) cell instead of the canonical one -- see idmrg_cell_raw_'s
    // own declaration for why the window cannot use the canonical cell,
    // and pyitensor/idmrg_window.py's _close_array_chain for the
    // Python-side original. Same contents (transfer tensors, both
    // families of fixed points), same invalidation, and the same reason
    // for caching: none of it depends on which operator is measured, on
    // x, or on the time step, so recomputing the transfer-matrix
    // eigenproblems inside every one of a run's nt*len(x_values)
    // snapshots would be pure repetition.
    void
    iw_build_cache() const
        {
        if (iw_cache_valid_) return;
        iw_require_cell("iw_build_cache");
        int n = (int)idmrg_cell_raw_.size();
        // chis[n] = chis[0]: the cell wraps. That is only meaningful if
        // the last position's own right bond really has the first
        // position's own left dimension -- i.e. dim(HR_ket)==dim(HL_ket)
        // for the raw cell, whose two outer legs are exactly those. A
        // mismatch means the growing algorithm has not settled into a
        // translationally-invariant state yet.
        if (idmrg_cell_raw_r_[n-1] != idmrg_cell_raw_l_[0])
            throw std::runtime_error(
                  "Chain::td_dynamical_correlator_window: the converged unit "
                  "cell's own wraparound bond dimension is inconsistent (its "
                  "right edge is " + std::to_string(idmrg_cell_raw_r_[n-1]) +
                  ", its left edge " + std::to_string(idmrg_cell_raw_l_[0]) +
                  ") -- a tiled window needs these to agree; try a different "
                  "maxm/maxiter/etol combination for gs_energy()/"
                  "idmrg_ground_state()");
        iw_chis_.assign(n+1,0);
        for (int k=0;k<n;++k) iw_chis_[k] = idmrg_cell_raw_l_[k];
        iw_chis_[n] = idmrg_cell_raw_l_[0];
        iw_Es_.resize(n);
        for (int k=0;k<n;++k)
            iw_Es_[k] = ic_transfer(idmrg_cell_raw_[k],idmrg_cell_raw_l_[k],
                                     idmrg_cell_raw_d_[k],idmrg_cell_raw_r_[k],{});
        auto rr = ic_all_right_fixed_points(iw_Es_,iw_chis_);
        iw_rho_after_ = std::move(rr.first);
        auto ll = ic_all_left_fixed_points(iw_Es_,iw_chis_);
        iw_l_before_ = std::move(std::get<0>(ll));
        iw_cache_valid_ = true;
        }

    void
    iw_require_cell(const char* who) const
        {
        if (!have_idmrg_cell_raw_)
            throw std::runtime_error(std::string("Chain::")+who+": called before "
                           "idmrg_ground_state, or the growing algorithm never "
                           "produced a gauge-consistent unit cell (maxiter must be "
                           "large enough for at least two macro-iterations)");
        }

    void
    ic_require_cell(const char* who) const
        {
        if (!have_idmrg_cell_)
            throw ITError(std::string("Chain::")+who+": called before "
                           "idmrg_ground_state, or the growing algorithm never "
                           "produced a gauge-consistent unit cell (maxiter must be "
                           "large enough for at least two macro-iterations)");
        }

    // The gauge-consistent unit cell re-gauged into left-canonical form and
    // normalized to a per-cell transfer eigenvalue of exactly 1 -- C++ port
    // of pyitensor/idmrg.py's _canonical_theta_cell + the two-position case
    // of _canonicalize_periodic (vx_canonicalize_n1 above is the same
    // construction specialized to VUMPS's single grouped supersite; this
    // one has to stay rectangular, see this section's own top comment).
    // A_1 is already an exact isometry (it is svd()'s own U), but
    // A_2 = S.V.lambda_o^-1 is only approximately one, which shows up as a
    // transfer eigenvalue of 1.00001 rather than 1 if left alone.
    // Returns false (leaving the cell raw) when the transfer spectrum is
    // (near-)degenerate -- expectation values still come out right via
    // ic_expectation, which normalizes by its own denominator; only the
    // norm diagnostic is then slightly off. Same fallback pyitensor takes.
    bool
    ic_canonicalize_cell(std::vector<std::vector<Cplx>>& cell,
                          std::vector<int>& cl, std::vector<int> const& cd,
                          std::vector<int>& cr) const
        {
        int n = (int)cell.size();
        std::vector<std::vector<Cplx>> Es(n);
        std::vector<int> chis(n+1);
        for (int k=0;k<n;++k) chis[k] = cl[k];
        chis[n] = cl[0];
        for (int k=0;k<n;++k) Es[k] = ic_transfer(cell[k],cl[k],cd[k],cr[k],{});

        std::vector<std::vector<Cplx>> rho_after, rho_L;
        Cplx eta_R(0,0), eta_L(0,0);
        std::vector<double> scales;
        try
            {
            auto rr = ic_all_right_fixed_points(Es,chis);
            rho_after = rr.first; eta_R = rr.second;
            auto ll = ic_all_left_fixed_points(Es,chis);
            rho_L = std::get<0>(ll); eta_L = std::get<1>(ll); scales = std::get<2>(ll);
            }
        catch (ITError const&)
            {
            return false; // degenerate transfer spectrum -- leave the cell raw
            }
        if (std::abs(eta_L-eta_R) > 1e-6*std::max(1.0,std::abs(eta_R)))
            return false;

        // rho_R_before[p] = "everything after position p-1" = rho_after[p-1].
        std::vector<std::vector<Cplx>> rho_R_before(n);
        for (int p=0;p<n;++p) rho_R_before[p] = rho_after[(p-1+n)%n];
        // Undo ic_all_left_fixed_points' own per-cut renormalization, and
        // transpose into the (bra,ket) convention the isometry identity
        // below needs -- see idmrg.py's own _canonicalize_periodic for the
        // full derivation of both steps (each was a separately-diagnosed
        // bug there).
        for (int p=0;p<n;++p)
            {
            int D = chis[p];
            std::vector<Cplx> t((size_t)D*D);
            for (int i=0;i<D;++i)
            for (int j=0;j<D;++j)
                t[(size_t)j*D+i] = rho_L[p][(size_t)i*D+j]*scales[p];
            rho_L[p] = std::move(t);
            }

        std::vector<std::vector<Cplx>> G_left(n), G_right(n);
        std::vector<int> keep_n(n,0);
        for (int p=0;p<n;++p)
            {
            int D = chis[p];
            int kx=0, ky=0;
            auto X = vx_psd_sqrt_factor(rho_L[p],D,kx);        // (kx,D)
            auto Ydag = vx_psd_sqrt_factor(rho_R_before[p],D,ky); // (ky,D) == Y^H
            auto Y = vx_dagger(Ydag,ky,D);                     // (D,ky)
            auto M = vx_matmul(X,kx,D,Y,D,ky);                 // (kx,ky)
            std::vector<Cplx> Usvd, Vt; std::vector<double> s;
            vx_economy_svd_full(M,kx,ky,Usvd,s,Vt);
            int rk = (int)s.size();
            int keep = vx_truncate(s,0.0,0,1);
            keep_n[p] = keep;
            std::vector<Cplx> Uk((size_t)kx*keep), Vtk((size_t)keep*ky);
            for (int i=0;i<kx;++i) for (int j=0;j<keep;++j) Uk[(size_t)i*keep+j] = Usvd[(size_t)i*rk+j];
            for (int i=0;i<keep;++i) for (int j=0;j<ky;++j) Vtk[(size_t)i*ky+j] = Vt[(size_t)i*ky+j];
            auto Udag = vx_dagger(Uk,kx,keep);                 // (keep,kx)
            G_left[p] = vx_matmul(Udag,keep,kx,X,kx,D);        // (keep,D)
            auto Vtkdag = vx_dagger(Vtk,keep,ky);              // (ky,keep)
            auto YV = vx_matmul(Y,D,ky,Vtkdag,ky,keep);        // (D,keep)
            G_right[p].assign((size_t)D*keep,Cplx(0,0));
            for (int i=0;i<D;++i)
            for (int j=0;j<keep;++j)
                G_right[p][(size_t)i*keep+j] = YV[(size_t)i*keep+j]/s[j];
            }

        std::vector<std::vector<Cplx>> new_cell(n);
        for (int p=0;p<n;++p)
            {
            int a_n = keep_n[p], d = cd[p], b_n = chis[p], c_n = chis[p+1];
            int e_n = keep_n[(p+1)%n];
            std::vector<Cplx> arr((size_t)a_n*d*e_n,Cplx(0,0));
            for (int a=0;a<a_n;++a)
            for (int q=0;q<d;++q)
            for (int e=0;e<e_n;++e)
                {
                Cplx acc(0,0);
                for (int b=0;b<b_n;++b)
                    {
                    Cplx inner(0,0);
                    for (int c=0;c<c_n;++c)
                        inner += cell[p][((size_t)b*d+q)*c_n+c]*G_right[(p+1)%n][(size_t)c*e_n+e];
                    acc += G_left[p][(size_t)a*b_n+b]*inner;
                    }
                arr[((size_t)a*d+q)*e_n+e] = acc;
                }
            new_cell[p] = std::move(arr);
            }

        // Left-canonicality check, at the same 1e-4 tolerance (and for the
        // same conditioning reason) idmrg.py's own _canonicalize_periodic
        // uses. Failing it abandons the re-gauging and keeps the raw cell,
        // rather than failing the whole ground-state run -- exactly what
        // pyitensor does, since _canonical_theta_cell wraps its
        // _canonicalize_periodic call in `except RuntimeError: return
        // cell` and that function raises RuntimeError for this very check.
        // Not hypothetical: the critical (gapless) Heisenberg chain at
        // maxm=30 lands at ~4e-4 here, because the growing algorithm's own
        // outer environment bond is still visibly short of a true Schmidt
        // basis of the converged state at that point. Expectation values
        // survive it -- ic_expectation closes against the left fixed point
        // and divides by the same contraction with the operators dropped,
        // so a non-canonical cell is normalized away rather than
        // mis-measured (that is precisely why _expectation stopped assuming
        // an identity left fixed point, see its own docstring); only the
        // eta-based norm diagnostic would be off, and nothing here consumes
        // one.
        for (int p=0;p<n;++p)
            {
            int a_n = keep_n[p], d = cd[p], e_n = keep_n[(p+1)%n];
            double worst = 0.0;
            for (int c=0;c<e_n;++c)
            for (int e=0;e<e_n;++e)
                {
                Cplx acc(0,0);
                for (int a=0;a<a_n;++a)
                for (int q=0;q<d;++q)
                    acc += std::conj(new_cell[p][((size_t)a*d+q)*e_n+c])
                            *new_cell[p][((size_t)a*d+q)*e_n+e];
                if (c==e) acc -= Cplx(1,0);
                worst = std::max(worst,std::abs(acc));
                }
            if (worst > 1e-4)
                {
                if (verbose_)
                    println("idmrg: unit-cell canonicalization left a ",worst,
                            " deviation from left-canonical form -- keeping the raw cell");
                return false;
                }
            }

        // ...and normalize to a per-cell transfer eigenvalue of exactly 1
        // (re-gauging alone does not get there -- see idmrg.py's own
        // _canonical_theta_cell). A pure normalization: ic_expectation
        // divides by the same fixed-point contraction, so it cannot change
        // an expectation value.
        double eta = eta_R.real();
        double scale = std::abs(eta) > 0.0 ? std::pow(eta,-1.0/(2.0*n)) : 1.0;
        for (int p=0;p<n;++p)
            for (auto & v : new_cell[p]) v *= scale;

        cell = std::move(new_cell);
        for (int p=0;p<n;++p) { cl[p] = keep_n[p]; cr[p] = keep_n[(p+1)%n]; }
        return true;
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
        // Realized window size in UNIT CELLS. Generally >= the n_window
        // the caller asked for: the tiling unit is the 2-site cell
        // (idmrg_cell_raw_), which is n_cell/n_uc unit cells long -- 1 for
        // n_uc=2 (nothing changes), 2 for n_uc=1, where an odd n_window is
        // rounded up so the realized window is never smaller than asked.
        // Every downstream site/sublattice computation (idmrg_window_center
        // above all) must use this, not the requested value.
        int n_window_uc = 0;
        std::vector<Index> phys; // phys[i-1] = window position i's own physical Index
        };

    // Dense (chi_l,d,chi_r) row-major array back into an ITensor over the
    // given (left,phys,right) Index triple -- inverse of
    // idmrg_tensor_to_lpr_array, used to lift the raw unit cell (kept as
    // plain arrays, see idmrg_cell_raw_) into the ITensor objects a real
    // MPS needs.
    static ITensor
    idmrg_lpr_array_to_tensor(std::vector<Cplx> const& A, Index const& left,
                               Index const& phys, Index const& right)
        {
        int chi_l = dim(left), d = dim(phys), chi_r = dim(right);
        if ((size_t)chi_l*d*chi_r != A.size())
            throw std::runtime_error("Chain::idmrg_lpr_array_to_tensor: shape mismatch");
        ITensor T(left,phys,right);
        for (int l=1;l<=chi_l;++l)
        for (int sI=1;sI<=d;++sI)
        for (int r=1;r<=chi_r;++r)
            T.set(left(l),phys(sI),right(r),
                   A[((size_t)(l-1)*d+(sI-1))*chi_r+(r-1)]);
        return T;
        }

    // Finite MPS/MPO covering `n_window` unit cells (rounded up to a whole
    // number of 2-site cells, see IdmrgWindow::n_window_uc), built by
    // tiling the gauge-consistent RAW unit cell idmrg_cell_raw_ and
    // idmrg_rows_ and capping the two open ends with idmrg_HL_ket_/
    // idmrg_HR_ket_ (ket) and idmrg_HL_mpo_/idmrg_HR_mpo_ (mpo) -- C++
    // analogue of idmrg_window.py's own build_window/_tile_periodic.
    //
    // Tiles idmrg_cell_raw_, NOT the raw per-micro-step idmrg_U_ factors
    // this used to tile. idmrg_U_'s two ends live in bond bases minted by
    // *different* micro-steps, so every copy boundary in a multi-copy
    // window silently identifies two different bases: the window's energy
    // stays right while its static observables do not. That is exactly the
    // failure idmrg_theta_cell exists to remove, and it was measured here
    // as an S(x,t=0) missing the exact static two_point_correlator by up
    // to 1.7e-1 on a plain Heisenberg chain (no fermions, no strings, any
    // operator; x=0 exact, error growing with |x|). See idmrg_cell_raw_'s
    // own declaration for why the RAW cell rather than the canonicalized
    // idmrg_cell_ static observables use: only the raw cell's outer legs
    // are still literally idmrg_HL_ket_/idmrg_HR_ket_, i.e. only it
    // attaches to the environment caps below with no assumption at all.
    //
    // Every window position gets its own genuinely fresh physical Index
    // (via idmrg_rows_[p].d, not sites_.si(p+1)) -- unlike pyitensor's own
    // _tile_periodic (which reuses sites_uc's own physical Index across
    // every copy of a given sublattice position, needing a separate
    // _refresh_physical_legs fix for n_uc==1, see that module's own
    // docstring), this sidesteps that whole bug class structurally: no two
    // window positions can ever collide on a shared physical Index,
    // regardless of n_uc.
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
        iw_require_cell("idmrg_build_window");
        iw_build_cache(); // also runs the wraparound-dimension check
        int n_uc = idmrg_n_uc_;
        int n_cell = (int)idmrg_cell_raw_.size();
        // n_window counts unit cells but the tiling unit is the cell.
        // Round the copy count up so the realized window is never smaller
        // than the caller asked for, and record the realized size in unit
        // cells so every downstream site/sublattice computation stays
        // exact (mirrors build_window's own uc_per_cell/n_copies).
        int uc_per_cell = n_cell/n_uc;
        int n_copies = (n_window + uc_per_cell - 1)/uc_per_cell;
        int n = n_copies*n_cell;

        static const TagSet link_tag("Link"), site_tag("Site");
        IdmrgWindow win;
        win.n = n;
        win.n_window_uc = n_copies*uc_per_cell;
        win.phys.resize(n);
        std::vector<ITensor> ket(n), mpoT(n);

        Index ket_left = idmrg_HL_ket_;
        Index mpo_left = idmrg_HL_mpo_;
        if (dim(ket_left) != idmrg_cell_raw_l_[0])
            throw std::runtime_error(
                  "Chain::idmrg_build_window: the converged cell's own left "
                  "bond dimension does not match the stored environment "
                  "snapshot's -- should be unreachable, since the cell is "
                  "seeded at the very micro-step that snapshot is taken at");
        for (int c=0; c<n_copies; ++c)
            {
            bool last_copy = (c==n_copies-1);
            for (int k=0; k<n_cell; ++k)
                {
                int i = c*n_cell+k+1; // 1-based window position
                int p = k%n_uc;       // sublattice at this cell position
                bool at_end = last_copy && (k==n_cell-1);

                win.phys[i-1] = Index(idmrg_rows_[p].d,site_tag);
                Index ket_right = at_end ? idmrg_HR_ket_
                                          : Index(idmrg_cell_raw_r_[k],link_tag);
                ket[i-1] = idmrg_lpr_array_to_tensor(idmrg_cell_raw_[k],ket_left,
                                                      win.phys[i-1],ket_right);

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
    //
    // Fermionic (parity-odd) `opname`: a bare C/Cdag matrix is not the
    // physical fermionic operator -- under Jordan-Wigner the operator at
    // site s is (prod_{j<s} F_j) O_s, so the string is applied explicitly
    // on every window site left of `site`, truncated at the window's own
    // left edge exactly as idmrg_window.py's own apply_local_operator does
    // (see its docstring for why truncating there is the same boundary
    // approximation the IBC window already makes, and why it is absent
    // altogether at t=0).
    void
    idmrg_window_apply_local_op(IdmrgWindow& win, int site, std::string const& opname) const
        {
        int n_uc = idmrg_n_uc_;
        if (idmrg_is_fermionic(opname))
            for (int j=1;j<site;++j)
                {
                int pj = (j-1)%n_uc;
                ITensor FT = idmrg_op_itensor(win.phys[j-1],pj,"F");
                ITensor Tj = win.psi.A(j)*FT;
                Tj.noPrime("Site");
                win.psi.set(j,Tj);
                }
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

    // Sum over a chain of doubled (ket,conj(bra)) transfer steps, closed
    // on the left by the tiled cell's own left transfer fixed point
    // (l_left, iw_build_cache's own iw_l_before_[0]) and on the right by
    // its right one (rho_right, iw_rho_after_[n_cell-1]) -- C++ analogue
    // of idmrg_window.py's own _close_array_chain, whose own comment
    // records why a bare trace on the left is not good enough here. The
    // caller divides by the same contraction over the unperturbed ground
    // state, see idmrg_window_snapshot_correlator.
    // bra_arrays/ket_arrays: same-length lists of flat
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
    // dimension (Er) to equal rho_right's own dimension (chi_right, the
    // converged cell's own wraparound bond dimension) -- a genuine
    // requirement of the method itself (idmrg_window.py's own
    // `_close_array_chain` has the identical constraint, via its final
    // `np.einsum('rR,rR->', left_traced, rho_after[p_right])`, which would
    // likewise raise a shape-mismatch error if this failed to hold), not
    // an artifact of this port: it holds once the window's own true
    // right-edge bond dimension (idmrg_HR_ket_, the accumulated growth
    // environment) has saturated to the *same* value as the cell's own
    // wraparound bond -- typically true once maxm/niter
    // are large enough for both to hit the same maxm ceiling. If it
    // doesn't hold, this raises a clear Error() (see the check below)
    // rather than segfaulting or silently truncating.
    static Cplx
    idmrg_close_array_chain(std::vector<std::vector<Cplx>> const& bra_arrays,
                             std::vector<std::vector<Cplx>> const& ket_arrays,
                             std::vector<int> const& ket_l, std::vector<int> const& ket_r,
                             std::vector<int> const& bra_l, std::vector<int> const& bra_r,
                             std::vector<int> const& dims_d,
                             std::vector<Cplx> const& l_left,
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
                // E[l,L,s,S] = sum_{r,R} E_old[l,L,r,R]*step[r,R,s,S].
                // step's own strides are (l:L*r*R, L:r*R, r:R, R:1), so
                // its second axis has extent L (the BRA's left dimension),
                // not R. Indexing it with R here instead was a real bug,
                // silent for as long as every tiled tensor happened to be
                // square (L==R, true of every cell with a uniform bond
                // dimension) and wrong the moment it is not: on a
                // dimerized spinless chain at maxm=20 the converged cell
                // comes out (18,20)/(20,18) and this put S(x,0) at 0.83
                // against an exact 0.51, for every operator including a
                // plain N.
                std::vector<Cplx> Enew((size_t)El*EL*r*R,Cplx(0,0));
                for (int li=0;li<El;++li)
                for (int Li=0;Li<EL;++Li)
                for (int si=0;si<r;++si)
                for (int Si=0;Si<R;++Si)
                    {
                    Cplx acc(0,0);
                    for (int ri=0;ri<Er;++ri)
                    for (int Ri=0;Ri<ER;++Ri)
                        acc += E[((size_t)(li*EL+Li)*Er+ri)*ER+Ri]*step[((size_t)(ri*L+Ri)*r+si)*R+Si];
                    Enew[((size_t)(li*EL+Li)*r+si)*R+Si] = acc;
                    }
                E = Enew; Er=r; ER=R;
                }
            }
        // left_traced[r,R] = sum_{l,L} l_left[l,L] E[l,L,r,R] -- closed
        // against the tiled cell's own LEFT transfer fixed point, not a
        // bare trace. The bare trace is the special case where every tiled
        // tensor is exactly left-canonical, which held for the old
        // idmrg_U_ tiling but does not for idmrg_cell_raw_ (its second
        // tensor S.V.lambda_o^-1 is an isometry only to ~1e-3). Left as a
        // bare trace it breaks the exact checks outright -- confirmed on
        // the Python side first (idmrg_window.py's own _close_array_chain
        // comment): S(x=0,t=0), which must equal <Sz Sz> = 0.25 for
        // spin-1/2, came out at -0.0776. El==EL is still required (the
        // fixed point is square): site 0's own ket_l==bra_l, both equal to
        // dim(idmrg_HL_ket_)==idmrg_cell_raw_l_[0].
        if (El!=EL || (size_t)El*EL != l_left.size())
            throw std::runtime_error("Chain::td_dynamical_correlator_window: window's own "
                  "left-edge bond dimension does not match the converged "
                  "unit cell's own natural left bond dimension needed to "
                  "close S(x,t) there -- should be unreachable given "
                  "idmrg_cell_raw_l_[0]==dim(idmrg_HL_ket_) by construction");
        std::vector<Cplx> left_traced((size_t)Er*ER,Cplx(0,0));
        for (int li=0;li<El;++li)
        for (int Li=0;Li<EL;++Li)
        for (int ri=0;ri<Er;++ri)
        for (int Ri=0;Ri<ER;++Ri)
            left_traced[(size_t)ri*ER+Ri] +=
                l_left[(size_t)li*EL+Li]*E[((size_t)(li*EL+Li)*Er+ri)*ER+Ri];
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
    // are freshly built from the static, converged gauge-consistent unit
    // cell idmrg_cell_raw_ (never touched by evolution -- the SAME cell
    // idmrg_build_window tiles the ket from, which is what makes the
    // t=0 identity S(x,0)==two_point_correlator exact; building the bra
    // from the raw per-micro-step idmrg_U_ chain instead is the gauge bug
    // idmrg_build_window's own comment records), with opname_A inserted at
    // position center+x -- see idmrg_close_array_chain's own comment for
    // why these two sides deliberately don't share any ITensor Index
    // identity.
    //
    // Both closures come from iw_build_cache's own fixed points of that
    // same cell: the left one (iw_l_before_[0], the window's first
    // position is cell position 0) because the raw cell is not exactly
    // left-canonical, and the right one (iw_rho_after_[n_cell-1], since
    // win.n is always a whole number of cells) as the correct weighting
    // for "everything beyond the window". The result is divided by the
    // same contraction run over the ground state itself, so a chain
    // carrying no operator returns exactly 1: both caps are independently
    // trace-1-normalized and the left one is not the identity, so without
    // this calibration the raw overlaps are off by an arbitrary constant
    // -- measured on the Python side (idmrg_window.py's own
    // _close_array_chain), S(x=0,t=0) came out at 0.025 instead of 0.25
    // with a chi-only rescaling and -0.078 with none. The denominator
    // depends only on the cell and on win.n, never on x or on t, so it is
    // computed once here rather than per x.
    std::vector<Cplx>
    idmrg_window_snapshot_correlator(IdmrgWindow const& win, std::string const& opname_A,
                                      std::vector<int> const& x_values, int center) const
        {
        int n = win.n;
        int n_uc = idmrg_n_uc_;
        int n_cell = (int)idmrg_cell_raw_.size();
        bool ferm_A = idmrg_is_fermionic(opname_A);
        iw_build_cache();
        auto const& l_left = iw_l_before_[0];
        auto const& rho_right = iw_rho_after_[n_cell-1];
        int chi_right = iw_chis_[n_cell];

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

        // bra_l/bra_r: the converged cell's own natural per-position
        // dimensions -- fixed for every x/t, unlike ket_l/ket_r above
        // (which change as win.psi evolves), so built once here rather
        // than inside the x loop below.
        std::vector<int> bra_l(n), bra_r(n);
        for (int i=1;i<=n;++i)
            {
            int k = (i-1)%n_cell;
            bra_l[i-1] = idmrg_cell_raw_l_[k];
            bra_r[i-1] = idmrg_cell_raw_r_[k];
            }

        // The calibration denominator: the identical closure run over the
        // unperturbed ground state on both sides, i.e. the plain product
        // of this window's own tiled transfer tensors, closed with the
        // same two caps. Applied right to left so no (chi^2,chi^2) matrix
        // product is ever formed.
        Cplx den(0,0);
        {
        std::vector<Cplx> v = rho_right;
        for (int i=n-1;i>=0;--i)
            {
            int k = i%n_cell;
            v = ic_matvec(iw_Es_[k],iw_chis_[k]*iw_chis_[k],
                           iw_chis_[k+1]*iw_chis_[k+1],v);
            }
        for (size_t idx=0; idx<l_left.size(); ++idx) den += l_left[idx]*v[idx];
        if (std::abs(den) < 1e-300)
            throw std::runtime_error("Chain::td_dynamical_correlator_window: the "
                  "window's own normalization is zero -- the converged unit "
                  "cell's transfer fixed points are degenerate or the state "
                  "is not converged");
        }

        std::vector<Cplx> out(x_values.size());
        for (size_t ix=0; ix<x_values.size(); ++ix)
            {
            int x = x_values[ix];
            int pos = center+x;
            if (pos<1 || pos>n)
                throw std::runtime_error("Chain::td_dynamical_correlator_window: x_values "
                      "must keep center+x within the window's own explicit "
                      "range [1,n] -- increase n_window (idmrg_window.py's "
                      "own padding for x beyond the window is not ported "
                      "here, see idmrg_window_snapshot_correlator's own "
                      "comment)");

            std::vector<std::vector<Cplx>> bra_arrays(n);
            for (int i=1;i<=n;++i)
                {
                // bulk tensor by *cell* position, operator matrix by
                // sublattice (n_cell is a multiple of n_uc, so the two
                // agree modulo n_uc).
                int k = (i-1)%n_cell, p = (i-1)%n_uc;
                int chi_l = idmrg_cell_raw_l_[k], d = idmrg_cell_raw_d_[k];
                int chi_r = idmrg_cell_raw_r_[k];
                std::vector<Cplx> M;
                if (i==pos)
                    {
                    // idmrg_close_array_chain conjugates the bra, so
                    // whatever is applied here acts as its own ADJOINT on
                    // the result: applying A^dagger is what makes this the
                    // documented <psi|A_x ...|psi>. (Applying A itself --
                    // what this did before, mirroring idmrg_window.py's own
                    // pre-fix bra -- silently computed <psi|A^dagger_x ...>
                    // instead: invisible for every Hermitian name, exactly
                    // wrong for C/Cdag/Sp/Sm.) idmrg_op_dense's own layout
                    // is [in,out] row-major, so the adjoint is
                    // Mdag[i,o] = conj(M[o,i]).
                    auto Mop = idmrg_op_dense(p,opname_A);
                    M.assign((size_t)d*d,Cplx(0,0));
                    for (int a=0;a<d;++a)
                    for (int b=0;b<d;++b)
                        M[(size_t)a*d+b] = std::conj(Mop[(size_t)b*d+a]);
                    }
                else if (ferm_A && i<pos)
                    {
                    // opname_A's own Jordan-Wigner string, truncated at the
                    // same left reference (window site 1) that
                    // idmrg_window_apply_local_op truncates the ket's at.
                    // F is Hermitian, so applying it to the bra inserts F.
                    M = idmrg_op_dense(p,"F");
                    }
                if (M.empty())
                    bra_arrays[i-1] = idmrg_cell_raw_[k];
                else
                    {
                    // Bop[l,o,r] = sum_i M[i,o] A[l,i,r] -- ic_transfer's
                    // own operator convention.
                    std::vector<Cplx> Bop((size_t)chi_l*d*chi_r,Cplx(0,0));
                    auto const& A = idmrg_cell_raw_[k];
                    for (int lI=0;lI<chi_l;++lI)
                    for (int o=0;o<d;++o)
                    for (int rI=0;rI<chi_r;++rI)
                        {
                        Cplx acc(0,0);
                        for (int iI=0;iI<d;++iI)
                            acc += M[(size_t)iI*d+o]*A[((size_t)lI*d+iI)*chi_r+rI];
                        Bop[((size_t)lI*d+o)*chi_r+rI] = acc;
                        }
                    bra_arrays[i-1] = std::move(Bop);
                    }
                }

            out[ix] = idmrg_close_array_chain(bra_arrays,ket_arrays,ket_l,ket_r,
                                               bra_l,bra_r,dims_d,l_left,
                                               rho_right,chi_right)/den;
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
        // Noise only in the first half of the schedule (mirrors the old
        // file-based backend's get_sweeps.h). Sweeps is 1-based -- its
        // arrays are sized nsweep+1 and ITensor's dmrg() loops
        // `for(sw=1; sw<=nsweeps; ++sw)` -- so the range that actually
        // covers the second half is [ns/2+1, ns]. The obvious-looking
        // 0-based `for (i=ns/2;i<ns;i++)` used here before silently left
        // the *final* sweep running with the full noise term, i.e. the
        // returned state and energy were always the output of a noisy
        // sweep. That is a real off-by-one, not a reproduced-on-purpose
        // quirk, and it is fixed here.
        for (int i=ns/2+1;i<=ns;i++) sweeps.setnoise(i,0.0);
        return sweeps;
        }

    Sweeps
    make_sweeps() const { return make_sweeps(nsweeps_,maxm_); }

    // Ground-state sweep schedule with a bond-dimension ramp: instead of
    // running every sweep at the full target maxdim, grow the bond
    // dimension geometrically over the first few sweeps and then hold it
    // at maxdim for the rest. The early sweeps are then much cheaper
    // (two-site DMRG cost is ~O(m^3)) while still doing the bulk of the
    // work of finding the right variational subspace, so the expensive
    // full-maxdim sweeps start from an already-good state -- the standard
    // ITensor sweep-table idiom (see sweeps.h's own example table).
    //
    // floor_dim is the bond dimension of an already-present starting
    // wavefunction, or 0 when starting fresh. It matters: ramping down
    // from a *converged* warm state (Chain::set_wavefunction, or
    // groundstate.py's maxde reconverge path, which re-runs at doubled
    // maxm from the previous solution) would truncate it and throw away
    // exactly the information the caller wanted kept. The ramp therefore
    // never starts below what the incoming state already carries.
    //
    // Noise: the noise term matters while the state is still being
    // built, so it starts at noise_, decays geometrically by
    // ramp_noise_decay_ per ramping sweep, and is switched off entirely
    // once the schedule reaches the full maxdim -- the final, converged
    // sweeps are noise-free. This is the standard ITensor sweep-table
    // noise ladder (see sweeps.h's own example table).
    Sweeps
    make_sweeps_ramped(int ns, int maxdim, int floor_dim) const
        {
        if (!ramp_ || ns<2 || maxdim<=1) return make_sweeps(ns,maxdim);
        int start = ramp_start_>0 ? ramp_start_ : 10;
        if (floor_dim>start) start = floor_dim;
        if (start>=maxdim) return make_sweeps(ns,maxdim);
        // The ramp spends a *fixed fraction* of the schedule growing, and
        // interpolates geometrically from start to maxdim across it, so
        // the bond dimension really does go up at (nearly) every one of
        // those sweeps. Growing as fast as possible instead -- doubling
        // each sweep until maxdim is reached, which was the first thing
        // tried here -- minimizes the number of *cheap* sweeps and so
        // gives away almost all of the speedup: on a 30-site
        // inhomogeneous Heisenberg-Hubbard chain at nsweeps=20 doubling
        // leaves only 4 of 20 sweeps below a maxdim of 150, where filling
        // half the schedule leaves 10 cheap sweeps costing a few percent
        // of a full one apiece. Measured end to end by
        // examples/groundstate/bond_dimension_ramp (BLAS pinned to one
        // thread, two dedicated cores): 2.15x at maxdim=90 and 1.7-2.0x
        // at maxdim=60 across runs, for the same energy to ~1e-8. The
        // speedup grows with maxdim, since that is what makes the ramped
        // sweeps cheap relative to the full ones.
        int nr = (int)((double)ns*ramp_fraction_);
        if (nr<1) nr = 1;
        if (nr>ns-1) nr = ns-1; // always at least one sweep at full maxdim
        auto sweeps = Sweeps(ns);
        sweeps.cutoff() = cutoff_;
        for (int i=0;i<ns;i++)
            {
            int m = maxdim;
            if (i<nr)
                {
                double x = std::pow((double)maxdim/(double)start,
                                    (double)i/(double)nr);
                m = (int)std::llround((double)start*x);
                if (m<1) m = 1;
                if (m>maxdim) m = maxdim;
                }
            // The noise term decays geometrically across the ramp and is
            // off entirely once the schedule reaches full maxdim, so the
            // final, converged sweeps are noise-free -- the standard
            // ITensor sweep-table noise ladder (sweeps.h's own example).
            double nz = i<nr ? noise_*std::pow(ramp_noise_decay_,i) : 0.0;
            sweeps.setmaxdim(i+1,m); // Sweeps is 1-based, see make_sweeps()
            sweeps.setnoise(i+1,nz);
            }
        return sweeps;
        }

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
            {
            double nrm_before = innerC(ap,ap).real();
            auto tr = kpm_energy_truncate(ap,m,dK,n_sweeps,threshold);
            check_kpm_truncation(nrm_before,tr.second.state_change_norm);
            ap = tr.first;
            }
            out.push_back(innerC(vj,ap));
            check_kpm_moment(out,bound);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    // Guard on Holzner's own Eq. (41) diagnostic, which both truncated
    // moment loops below used to compute and throw away (`.first`).
    // Energy truncation is only valid while the correlator's spectral
    // weight actually fits inside the ground-state-anchored window
    // [E0,E0+Ws] that kpm_scale sets. Push kpm_scale far enough down and
    // it does not: the discarded weight is not removed but piled onto the
    // window edges, and since the resulting moments stay bounded by mu0,
    // check_kpm_moment() cannot see it either -- the spectrum comes back
    // peaked at ~0.96x the window top with several times the correct
    // amplitude and no warning at all. The sharp, unambiguous end of that
    // spectrum is the truncation annihilating the whole Chebyshev vector
    // (relative state change ~1); healthy runs on the same chain sit
    // around 0.01-0.17. Only that regime is rejected here -- see
    // docs/known_issue_kpm_energy_truncation_window.md for what is still
    // silent between the two.
    void
    check_kpm_truncation(double norm_before, double state_change_norm) const
        {
        if (norm_before <= 0.) return;
        double rel = state_change_norm/norm_before;
        if (rel > 0.5)
            throw std::runtime_error(tinyformat::format(
                "kpm_energy_truncate: the energy truncation removed %.1f%% "
                "of a Chebyshev vector, i.e. essentially all of it -- the "
                "retained energy window is far too small for this "
                "correlator, and the moments it produces are meaningless. "
                "Raise kpm_scale, or set kpm_energy_truncate=False.",
                100.*rel));
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
            {
            double nrm_before = innerC(ap,ap).real();
            auto tr = kpm_energy_truncate(ap,m,dK,n_sweeps,threshold);
            check_kpm_truncation(nrm_before,tr.second.state_change_norm);
            ap = tr.first;
            }
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
    // -- conserved sector (opt-in QN mode), see set_conserved_sector() --
    // site_types_ is what SpinX needs in order to rebuild sites_ either
    // way. dense_sites_ is the chain's own original ConserveQNs=false site
    // set, kept from construction and never rebuilt: it is what an
    // operator's matrix elements are read from without ever building that
    // operator on QN-carrying indices (see op_charge(), and why that
    // matters), what clearing the sector puts back into sites_, and what
    // promote_to_dense() rebases a sector state onto. sector_ is the
    // requested (QN name, value) list, empty whenever sector mode is off.
    std::vector<int> site_types_;
    SiteSet dense_sites_;
    std::vector<std::pair<std::string,int>> sector_;
    bool has_sector_ = false;
    mutable int sector_draws_ = 0; // varies sector_mps()'s arrangement per call
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
    // Bond-dimension ramp for the ground-state sweep schedule, see
    // set_bond_ramp()/make_sweeps_ramped(). Defaults mirror
    // Many_Body_Chain's own (manybodychain.py), which is what actually
    // drives them in practice -- these only apply if Python never calls
    // set_bond_ramp() at all.
    bool ramp_ = true;
    int ramp_start_ = 10;
    double ramp_fraction_ = 0.5;
    double ramp_noise_decay_ = 0.1;

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

    // The gauge-consistent 2-site unit cell (idmrg_theta_cell, re-gauged by
    // ic_canonicalize_cell) and its own per-position bond/physical
    // dimensions -- what every static observable on this backend tiles.
    // idmrg_cell_[k] is a (idmrg_cell_l_[k], idmrg_cell_d_[k],
    // idmrg_cell_r_[k]) row-major array, with idmrg_cell_r_[k] ==
    // idmrg_cell_l_[(k+1)%2] by construction (the cell wraps). Cell
    // position k carries sublattice k%idmrg_n_uc_.
    bool have_idmrg_cell_ = false;
    std::vector<std::vector<Cplx>> idmrg_cell_;
    std::vector<int> idmrg_cell_l_, idmrg_cell_d_, idmrg_cell_r_;
    // The SAME cell before ic_canonicalize_cell re-gauges it -- what the
    // IBC window (idmrg_build_window / idmrg_window_snapshot_correlator)
    // must tile, and deliberately NOT idmrg_cell_ above. theta's own two
    // outer legs literally ARE idmrg_HL_ket_/idmrg_HR_ket_ (the
    // environment snapshot is taken entering the very micro-step theta is
    // solved at, see the snap_* block in idmrg_ground_state), so the raw
    // cell attaches to the window's own environment caps exactly, with no
    // assumption at all; ic_canonicalize_cell's re-gauging is precisely
    // what destroys that correspondence. Same choice, and for the same
    // reason, as pyitensor/idmrg_window.py's own _window_cell picking
    // IDMRGResult.cell_raw over cell_list. The two are the same state, so
    // nothing is given up -- the price is only that the raw cell is not
    // exactly left-canonical (its second tensor S.V.lambda_o^-1 is an
    // isometry only to ~1e-3), which is why the window's own closures use
    // genuine left/right transfer fixed points plus a calibration
    // denominator (iw_build_cache / idmrg_close_array_chain) rather than
    // the bare trace a left-canonical chain would allow.
    bool have_idmrg_cell_raw_ = false;
    std::vector<std::vector<Cplx>> idmrg_cell_raw_;
    std::vector<int> idmrg_cell_raw_l_, idmrg_cell_raw_d_, idmrg_cell_raw_r_;
    // The raw cell's own transfer tensors and both families of fixed
    // points, cached exactly as ic_build_cache caches the canonical
    // cell's (see iw_build_cache).
    mutable bool iw_cache_valid_ = false;
    mutable std::vector<std::vector<Cplx>> iw_Es_;
    mutable std::vector<std::vector<Cplx>> iw_rho_after_, iw_l_before_;
    mutable std::vector<int> iw_chis_;

    // Lazily-built, measurement-independent part of every static
    // observable on that cell (ic_build_cache): the per-position transfer
    // tensors and both families of fixed points. Mutable because the
    // observables themselves are const and this is a pure memoization;
    // invalidated at the top of every idmrg_ground_state() run.
    mutable bool ic_cache_valid_ = false;
    mutable std::vector<std::vector<Cplx>> ic_Es_, ic_rho_after_, ic_l_before_;
    mutable std::vector<int> ic_chis_;

    // The very last micro-step's own 2-site local eigenproblem, kept so
    // idmrg_local_excitation_gap() can re-diagonalize that same effective
    // Hamiltonian for a second eigenpair -- C++ analogue of
    // pyitensor/idmrg.py's own IDMRGResult.local_superblock.
    bool have_idmrg_superblock_ = false;
    // Sequential multi-site VUMPS snapshot (n_uc>2) -- the per-site
    // tensors vms_onsite_expectation/vms_two_point_correlator read. Kept
    // separate from the grouped vumps_* snapshot rather than overloading
    // it: they are different representations, and a grouped consumer
    // reading a per-site list as one tensor is exactly the failure this
    // separation prevents.
    std::vector<std::vector<Cplx>> vms_AL_, vms_AR_, vms_C_, vms_AC_;
    std::vector<IdmrgAutomatonRow> vms_rows_;
    int vms_D_ = 0;
    bool have_vms_snapshot_ = false;

    ITensor idmrg_sb_HL_, idmrg_sb_HR_, idmrg_sb_W_pL_, idmrg_sb_W_pR_, idmrg_sb_theta_;
    Index idmrg_sb_HL_bra_, idmrg_sb_HL_ket_, idmrg_sb_HR_bra_, idmrg_sb_HR_ket_;
    Index idmrg_sb_phys_L_, idmrg_sb_phys_R_;
    bool idmrg_sb_have_HL_ket_ = false, idmrg_sb_have_HR_ket_ = false;
    double idmrg_sb_energy_ = 0.0;
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
