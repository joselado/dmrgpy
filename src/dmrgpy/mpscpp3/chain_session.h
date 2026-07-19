#include <tuple>

// TDVP/tdvp.h is a quoted include, so it resolves relative to this file's
// own directory (mpscpp3/) with no Makefile/include-path change needed.
// mpscpp3-only: TDVP/ has no v2 counterpart (mpscpp2's Chain has no TDVP
// methods at all), so tdvp_step()/quench_tdvp()/evolve_and_measure_tdvp()
// below only exist in this file. Only the two-site TDVP algorithm is used
// (NumCenter=2 in tdvp_step() below) -- two-site TDVP grows bond dimension
// via SVD the same way two-site DMRG does, which is enough for the
// quench/correlator use cases this class exposes, so the global subspace
// expansion machinery in TDVP/basisextension.h (needed mainly for one-site
// TDVP, or long-range Hamiltonians where two-site growth is too slow) is
// deliberately not wired in here. See TDVP/README.md for the algorithm.
#include "TDVP/tdvp.h"

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
    // has an equivalent one.
    MPS
    tdvp_step(MPO const& H, MPS psi, Cplx dt) const
        {
        Cplx t = Cplx(0.0,-1.0)*dt;
        auto sweeps = Sweeps(1);
        sweeps.maxdim() = maxm_;
        sweeps.cutoff() = cutoff_;
        sweeps.niter() = 50;
        tdvp(psi,H,t,sweeps,{"Quiet",!verbose_,"Silent",!verbose_,
                              "NumCenter",2,"DoNormalize",true});
        return psi;
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
