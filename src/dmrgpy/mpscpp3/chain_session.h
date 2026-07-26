#include <tuple>
#include <map> // Chain::MettsEigCache (Chain::metts_vev)
#include <algorithm> // std::next_permutation (four_correlation_tensor_sweep)
#include <random> // std::mt19937_64/std::discrete_distribution (Chain::metts_vev)

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
    NHDMRGResult
    nhdmrg(std::vector<MOTerm> const& terms_h,
           std::vector<MOTerm> const& terms_hadj,
           int krylovdim=20, int restarts=2)
        {
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto HA = build_mpo(sites_,terms_hadj,mpomaxm_);
        int N = sites_.length();
        // always a fresh random start (never wf0_): the non-Hermitian
        // energy is not a variational bound, so a rare stalled run can
        // only be detected by the caller's eigen-residual check and cured
        // by re-running -- which requires every run to draw its own
        // random initial state (see nhdmrg.py's retry loop)
        MPS psir = default_mps();
        psir.position(1);
        psir.normalize();
        MPS psil = psir;
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
        Cplx energy = 0;
        bool have_energy = false;
        for (int sw=1;sw<=nsweeps_;++sw)
            {
            // mirrors make_sweeps(): noise only in the first half-schedule
            double noise = (sw<=nsweeps_/2) ? noise_ : 0.0;
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
    template <typename Fn>
    std::pair<Cplx,ITensor>
    arnoldi_smallest_real(Fn&& A, ITensor x0, int krylovdim, int restarts,
                          Sel sel=Sel::SR, Cplx target=0) const
        {
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
                if (j+1<krylovdim) V.push_back(w/nw);
                }
            auto a = Index(m,"a");
            auto Hm = ITensor(prime(a),a);
            for (int i=0;i<m;++i)
            for (int j=0;j<m;++j)
                if (i<=j+1) Hm.set(prime(a)(i+1),a(j+1),h.at(i).at(j));
            ITensor W,D;
            eigen(Hm,W,D);
            auto c = commonIndex(W,D);
            std::vector<Cplx> ev(m);
            for (int k=1;k<=m;++k) ev[k-1] = eltC(D,c(k),prime(c)(k));
            // Ritz-value selection; the global-min-then-candidates
            // SRTieBreak formulation is deliberately identical to
            // pyitensor/nhdmrg.py's _select_ritz
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
    static void
    check_kpm_moment(std::vector<std::complex<double>> const& out, double bound)
        {
        if (std::abs(out.back()) > 1e3*(bound+1.0))
            Error("KPM moments diverging: scaled spectrum outside [-1,1] "
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
