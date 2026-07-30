#include <tuple>
#include <map> // Chain::MettsEigCache (Chain::metts_vev)
#include <algorithm> // std::next_permutation (four_correlation_tensor_sweep)
#include <random> // std::mt19937_64/std::discrete_distribution (Chain::metts_vev)
#include <cmath> // std::isnan (Chain::gs_energy_generalized's lam0 sentinel)
#include <limits> // std::numeric_limits<double>::quiet_NaN() (ditto)

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

        for (macro_iter=0; macro_iter<maxiter; ++macro_iter)
            {
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

                auto [energy_here,U,V,new_bond_u,new_bond_v] = idmrg_local_solve(
                    HL,W_pL,phys_L_in,W_pR,phys_R_in,HR,
                    HL_bra,HL_ket,have_HL_ket,HR_bra,HR_ket,have_HR_ket,
                    cutoff,maxm,krylovdim,restarts);
                energy = energy_here;

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
    MPS
    apply_mpo(MPO const& K, MPS const& x, Args const& args) const
        {
        auto out = applyMPO(K,x,args);
        out.noPrime(TagSet("Site"));
        return out;
        }

    MPS
    apply_mpo(MPO const& K, MPS const& x, MPS const& x0, Args const& args) const
        {
        auto out = applyMPO(K,x,x0,args);
        out.noPrime(TagSet("Site"));
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
                auto mat = eye();
                if (!onsite_mat.empty())
                    for (int k=0;k<d*d;++k) mat[k] += onsite_mat[k];
                setmat(li,ri,mat);
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

    // One micro-step's local ground-state solve: the effective 2-site
    // Hamiltonian sandwiched by (HL, W_pL, W_pR, HR), diagonalized via
    // arnoldi_smallest_real(..., Sel::SR) (the smallest-real-part Ritz
    // value of a Hermitian operator *is* its ground state, so the
    // existing non-Hermitian-capable Arnoldi engine already used by
    // nhdmrg_one_sweep above is reused directly rather than writing a
    // second, Hermitian-only Krylov solver) -- matches idmrg.py's own
    // _local_two_site_solve (kernels.make_matvec + dmrg._lanczos_ground_state),
    // but ITensor's own operator* already does all the index-matching
    // contraction that Python's kernels.py has to hand-roll. No
    // pre-existing ket tensor to seed the local guess from (every
    // micro-step inserts brand-new sites), so the Arnoldi start vector is
    // always a fresh random one (randomITensorC), matching
    // idmrg.py's own rationale (see that function's own docstring on why
    // mpscpp3/pyitensor never seed DMRG from a product state). U and V
    // are returned *without* S (the singular values are discarded): HL/HR
    // are block operators built by a similarity transform through U (or
    // V) alone, exactly idmrg.py's own convention -- see that function's
    // own extensive comment on why absorbing sqrt(S) into both sides
    // instead is wrong (confirmed there against independent ED).
    //
    // A dedicated Hermitian Lanczos solver with early-exit convergence
    // checking (mirroring pyitensor's own _lanczos_ground_state) was
    // tried here to cut the remaining Krylov-solve cost further, but was
    // reverted after it produced measurably wrong ground-state energies
    // once bond dimension grew large (traced to some accumulated-error or
    // logic issue in that from-scratch routine that this port's own
    // regression tests didn't catch quickly enough to trust shipping it)
    // -- arnoldi_smallest_real's own double-reorthogonalization,
    // proven-correct machinery is kept here instead. See this codebase's
    // iDMRG optimization history for the (much larger) fix that *is* kept:
    // idmrg_extend_HL/HR's contraction-order reordering below.
    std::tuple<double,ITensor,ITensor,Index,Index>
    idmrg_local_solve(ITensor const& HL, ITensor const& W_pL, Index phys_L,
                       ITensor const& W_pR, Index phys_R, ITensor const& HR,
                       Index HL_bra, Index HL_ket, bool have_HL_ket,
                       Index HR_bra, Index HR_ket, bool have_HR_ket,
                       double cutoff, int maxdim, int krylovdim, int restarts) const
        {
        std::vector<Index> order_in;
        if (have_HL_ket) order_in.push_back(HL_ket);
        order_in.push_back(phys_L);
        order_in.push_back(phys_R);
        if (have_HR_ket) order_in.push_back(HR_ket);
        auto x0 = randomITensorC(IndexSet(order_in));

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

        std::vector<Index> left_inds;
        if (have_HL_ket) left_inds.push_back(HL_ket);
        left_inds.push_back(phys_L);
        auto [U,S,V] = svd(theta,IndexSet(left_inds),{"Cutoff",cutoff,"MaxDim",maxdim});
        Index new_bond_u = commonIndex(U,S);
        Index new_bond_v = commonIndex(S,V);
        return {energy,U,V,new_bond_u,new_bond_v};
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
    //    the *full* dense Hermitian projected matrix is built (one extra
    //    matvec pass over the already-orthonormalized basis) rather than
    //    reusing only the entries encountered during Gram-Schmidt.
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
        for (int j=0;j<dK-1;++j)
            {
            ITensor w; PH.product(V.at(j),w);
            for (int pass=0;pass<2;++pass)
                for (auto const& Vi : V)
                    w -= eltC(dag(Vi)*w)*Vi;
            double nw = norm(w);
            if (nw<1E-12) break; // invariant subspace found; fewer than dK vectors is fine
            V.push_back(w/nw);
            }
        int k = (int)V.size();
        auto a = Index(k,TagSet("KPMEnergyTrunc,a"));
        auto Hk = ITensor(prime(a),a);
        for (int j=0;j<k;++j)
            {
            ITensor w; PH.product(V.at(j),w);
            for (int i=0;i<k;++i)
                Hk.set(prime(a)(i+1),a(j+1),eltC(dag(V.at(i))*w));
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
    };
