#include <tuple>

// Phase 3 session/handle model: one Chain instance per Python
// Many_Body_Chain object, replacing the old model of one fresh mpscpp.x
// process (and one scratch directory full of files) per DMRG task.
//
// This is the first slice of the migration -- ground state energy, generic
// VEV, and apply-operator, matching the "simplest/most-used, exercises the
// whole plumbing" template called out in the migration plan. Everything
// else (excited states, KPM correlators, time evolution, ...) follows the
// same pattern in later slices.
//
// Sites/Hamiltonian/wavefunction all live as ordinary members here instead
// of being written to sites.in/hamiltonian.in/psi_GS.mps and re-read by a
// freshly spawned process -- this is what "direct memory transfer" means in
// practice: Python calls a method on this object and gets a value or an
// opaque MPS handle back, with nothing touching disk in between.

// Result of the KPM dynamical-correlator computation (get_moments_dynamical_
// correlator in kpmmoments.h): the Chebyshev moments plus the energy-scale
// bookkeeping the Python side needs to unscale them (see kpmdmrg.py's
// get_dynamical_correlator, which today reads these back from KPM_MOMENTS.OUT/
// KPM_SCALE.OUT/KPM_NUM_POLYNOMIALS.OUT).
struct KPMResult
    {
    std::vector<std::complex<double>> moments;
    double emin, emax, scale;
    int num_polynomials;
    };

// Result of excited_states(): mirrors get_excited.h's EXCITED.OUT
// (energy + energy-fluctuation per state) plus the wavefunctions
// themselves (get_excited.h writes those to wavefunction_<i>.mps).
struct ExcitedResult
    {
    std::vector<double> energies;
    std::vector<double> fluctuations;
    std::vector<MPS> wavefunctions;
    };

// Result of quench()/evolve_and_measure(): the per-timestep overlap/
// expectation value (mirrors TIME_EVOLUTION.OUT) plus the final
// wavefunction (mirrors psi_time_evolution.mps/psi_evolve_and_measure.mps).
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
    num_sites() const { return sites_.N(); }

    // Local Hilbert-space dimension at a (1-based) site, e.g. for reshaping
    // reduced_dm()'s flat output without any float-precision guessing.
    int
    site_dim(int site) const { return sites_.si(site).m(); }

    // Mirrors get_random_mps.h exactly: a default-constructed MPS on this
    // Chain's sites (not a rigorously random state on its own -- see
    // randommps.py's random_mps(), which combines two of these with a
    // random complex phase and normalizes). This is what is_hermitian()/
    // gs_energy()'s Hermiticity check ultimately depends on, so without
    // it the old file-based backend stayed a hard dependency for every
    // single gs_energy() call even with the extension otherwise fully
    // wired up (one extra subprocess spawn per call, not just leftover
    // dead code).
    MPS
    random_mps() const { return MPS(sites_); }

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

    // Builds and caches the Hamiltonian MPO. Because this is a member
    // (rather than re-read from hamiltonian.in on every single call, as
    // get_hamiltonian.h does today), repeated gs_energy()/vev()/... calls
    // against the same Hamiltonian never rebuild it -- this is what fixes
    // the "get_hamiltonian.h always rereads from disk" inefficiency flagged
    // in the migration plan.
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
            wf0_ = MPS(sites_);
            have_wf0_ = true;
            }
        auto sweeps = make_sweeps();
        double energy = dmrg(wf0_,H_,sweeps);
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

    // Excited states via the state-targeting/Lagrange-multiplier technique,
    // mirroring get_excited.h's get_excited() exactly, with one deliberate
    // difference: the original computes a *second*, independent, unused
    // ground state (a fresh MPS(sites)+dmrg() whose result is immediately
    // discarded in favor of get_gs()'s already-cached one -- dead work, an
    // artifact of the file-based code having no memory of what an earlier
    // task already computed) before seeding wfs[0] from get_gs(); here
    // wfs[0] is seeded directly from gs_energy()'s cached ground state, one
    // fewer full DMRG run for exactly the same result.
    ExcitedResult
    excited_states(int n, double scale_lagrange=1.0, bool do_gram_schmidt=false)
        {
        if (!have_H_) Error("Chain::excited_states called before set_hamiltonian");
        if (!have_wf0_) gs_energy(); // ensure a ground state is available
        auto sweeps = make_sweeps();
        std::vector<MPS> wfs;
        auto psi0 = wf0_;
        normalize(psi0);
        wfs.push_back(psi0);
        double weight = bandwidth()*scale_lagrange;
        for (int i=1;i<n;i++)
            {
            MPS psi1(sites_);
            dmrg(psi1,H_,wfs,sweeps,{"Quiet",true,"Weight",weight});
            normalize(psi1);
            wfs.push_back(psi1);
            }
        if (do_gram_schmidt) wfs = gram_schmidt(wfs);
        ExcitedResult out;
        for (auto& wf : wfs)
            {
            out.fluctuations.push_back(energy_fluctuation(wf,H_));
            out.energies.push_back(overlap(wf,H_,wf));
            out.wavefunctions.push_back(wf);
            }
        return out;
        }

    // Generic vacuum expectation value <wf|A^npow|wf>, mirroring vev.h's
    // multi_vev logic exactly, just taking the operator/wavefunction
    // directly instead of reading vev_multioperator.in/wf_vev.mps.
    std::complex<double>
    vev(std::vector<MOTerm> const& terms, MPS const& wf, int npow=1)
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto psi = wf;
        psi /= sqrt(overlap(psi,psi)); // normalize
        Cplx c = 0;
        if (npow==1) { c = overlapC(psi,A,psi); }
        if (npow>1)
            {
            auto psi1 = psi;
            for (int i=0;i<npow-1;i++)
                psi1 = exactApplyMPO(psi1,A,{"Maxm",maxm_,"Cutoff",cutoff_});
            c = overlapC(psi,A,psi1);
            }
        return c;
        }

    // Applies an arbitrary operator to a wavefunction, mirroring
    // applyoperator.h's applyoperator() exactly.
    MPS
    apply_operator(std::vector<MOTerm> const& terms, MPS const& wf)
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        return exactApplyMPO(wf,A,args);
        }

    // Plain overlap <wf1|wf2>, mirroring compute_overlap.h's
    // compute_overlap() exactly (mpsalgebra.py::overlap_dmrg).
    std::complex<double>
    overlap_mps(MPS const& wf1, MPS const& wf2) const
        {
        return overlapC(wf1,wf2);
        }

    // Overlap <wf1|A|wf2> for an arbitrary operator, mirroring
    // applyoperator.h's overlap_aMb() exactly (mpsalgebra.py::
    // overlap_aMb_dmrg_MO).
    std::complex<double>
    overlap_aMb(MPS const& wf1, std::vector<MOTerm> const& terms, MPS const& wf2) const
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        return overlapC(wf1,A,wf2);
        }

    // Sum of two MPS, mirroring applyoperator.h's get_summps()/
    // mpsalgebra.h's sum_mps() exactly (mpsalgebra.py::summps_dmrg).
    MPS
    sum_mps(MPS const& wf1, MPS const& wf2) const
        {
        return sum(wf1,wf2,{"Maxm",maxm_,"Cutoff",cutoff_});
        }

    // Complex conjugate of an MPS, mirroring conjugatemps.h's
    // conjugate_mps() exactly (mpsalgebra.py::conjugate_mps). conjugatemps.h
    // itself delegates to cvm.h's conjMPS(), but cvm.h also pulls in
    // apply_inverse() (a plain, non-generic lambda that calls tasks.in-based
    // helpers like get_mpo_operator()/read_wf() directly in its body, so it
    // cannot be included without those); conjMPS's own logic is a few lines,
    // so it is inlined here instead of pulling in the rest of cvm.h.
    MPS
    conjugate(MPS const& wf) const
        {
        auto out = wf;
        for (int i=0;i<out.N();++i) out.Aref(i+1).conj();
        return out;
        }

    // Reduced density matrix, mirroring reduced_dm.h's reduced_dm() exactly
    // (including its psi /= overlap(psi,psi) normalization -- dividing by
    // the squared norm rather than its square root looks like a bug, but is
    // a no-op in practice since wf here is always already unit-normalized
    // coming out of dmrg()/gs_energy(), so it's preserved verbatim rather
    // than "fixed" as an unrelated change), just taking wf/site directly
    // and returning the flat dim*dim matrix instead of DM.OUT (which today
    // is written one tensor element at a time via a separate file
    // open/close per element). site is 1-based, matching every other
    // site-index convention at this extension boundary (MOTerm::factors,
    // build_ampo/build_mpo).
    std::vector<std::complex<double>>
    reduced_dm(MPS const& wf, int site) const
        {
        auto psi = wf;
        psi /= overlap(psi,psi); // normalize (see note above)
        psi.position(site);
        auto ir = commonIndex(psi.A(site),psi.A(site+1));
        auto rho = psi.A(site)*dag(prime(psi.A(site),Site,ir));
        for (int k=site+1;k<=psi.N();++k)
            {
            rho *= psi.A(k);
            rho *= dag(prime(psi.A(k),Link));
            }
        std::vector<std::complex<double>> out;
        auto collect = [&out](Cplx z) { out.push_back(z); };
        rho.visit(collect);
        return out;
        }

    // Applies exp(tau*H) to a wavefunction via nsteps repeated applications
    // of a 2nd-order Taylor-expanded exponential, mirroring time_evolution.h's
    // exponential_eMwf() task in its default (tevol_custom_exp=true) mode --
    // the only mode Python ever uses (mpsalgebra.py::exponential_dmrg always
    // sets tevol_custom_exp when self.tevol_custom_exp, which defaults to
    // true in manybodychain.py). H is passed in as a term list rather than
    // read from an implicit hamiltonian.in, since exponential_dmrg can be
    // called with an operator other than this chain's own Hamiltonian.
    MPS
    exponential_apply(std::vector<MOTerm> const& terms, MPS const& wf,
                       std::complex<double> tau, int nsteps) const
        {
        auto H = build_mpo(sites_,terms,mpomaxm_);
        auto taui = tau/double(nsteps);
        auto expH = custom_exp(H,taui);
        auto args = Args("Cutoff=",cutoff_,"Maxm=",maxm_);
        auto psi1 = wf;
        for (int it=1;it<=nsteps;it++) psi1 = exactApplyMPO(expH,psi1,args);
        return psi1;
        }

    // The methods below mirror multioperatortk/staticoperator.py's
    // StaticOperator: a *pre-compiled* operator handle (an ITensor MPO
    // held directly, no re-parsing of a term list on every use), matching
    // pureapplyoperator.h/multmpo.h's "pureoperator" family of tasks
    // (which serialize/deserialize a plain ITensor MPO to/from a .mpo
    // file, as distinct from get_ampo_operator.h's ampotk text format
    // used everywhere else in this class). build_operator() is the bridge
    // from a term list to such a handle, mirroring gen_pureoperator().

    MPO
    build_operator(std::vector<MOTerm> const& terms) const
        {
        return build_mpo(sites_,terms,mpomaxm_);
        }

    // Mirrors pureapplyoperator.h's pureapplyoperator() task.
    MPS
    apply_pure_operator(MPO const& A, MPS const& wf) const
        {
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        return exactApplyMPO(wf,A,args);
        }

    // Mirrors multmpo.h's multmpo_operator() task.
    MPO
    multiply_operators(MPO const& A, MPO const& B) const
        {
        return mult_mpo(A,B);
        }

    // Mirrors multmpo.h's trace_mpo_operator() task / operators.h's
    // trace_mpo() (Tr[A] = <Id|A>).
    std::complex<double>
    trace_operator(MPO const& A) const
        {
        auto ampo = AutoMPO(sites_);
        ampo += 1.0,"Id",1;
        auto ident = MPO(ampo);
        return overlapC(ident,A);
        }

    // Mirrors multmpo.h's hermitianmpo_operator() task / operators.h's
    // hermitian_mpo(): daggers each site tensor and swaps its Site-type
    // prime level 0<->1.
    MPO
    hermitian_operator(MPO const& A) const
        {
        auto out = A;
        for (auto j : range1(out.N())) out.Aref(j) = dag(swapPrime(out.A(j),0,1,Site));
        return out;
        }

    // Mirrors pureapplyoperator.h's overlap_aMb_static() task.
    std::complex<double>
    overlap_aMb_operator(MPS const& wf1, MPO const& A, MPS const& wf2) const
        {
        return overlapC(wf1,A,wf2);
        }

    // Von Neumann entanglement entropy across the bond between sites b and
    // b+1 (1-based), mirroring get_entropy.h's get_entropy() task /
    // entropy() helper exactly.
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

    // Real-time quench dynamical correlator, mirroring time_evolution.h's
    // quench() task (the "time_evolution" task): evolves two probe
    // wavefunctions (A1|GS>, A2|GS>) under exp(-i dt (H-EGS)) and records
    // their overlap at each step. H is passed as an explicit term list
    // (like exponential_apply/evolve_and_measure below) since
    // timedependent.py's evolution_dmrg_DC always takes an explicit h
    // parameter (defaulting to, but not necessarily, this chain's own
    // Hamiltonian). Only the tevol_custom_exp=true path (evoloperator(),
    // not toExpH<ITensor>) is implemented, matching every other
    // time-evolution method ported so far -- see exponential_apply's note.
    TimeEvolutionResult
    quench(std::vector<MOTerm> const& terms_h,
           std::vector<MOTerm> const& terms_i,
           std::vector<MOTerm> const& terms_j,
           int nt, double dt, bool fit_td=true)
        {
        if (!have_wf0_) gs_energy();
        auto H = build_mpo(sites_,terms_h,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        double EGS = overlap(wf0_,H,wf0_)/overlap(wf0_,wf0_);
        auto ampo = build_ampo(sites_,terms_h);
        ampo += -EGS,"Id",1;
        auto expH = evoloperator(MPO(ampo),dt);
        auto A1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto A2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto psi1 = exactApplyMPO(wf0_,A1,args);
        auto psi2 = exactApplyMPO(wf0_,A2,args);
        Cplx norm0 = std::sqrt(overlapC(psi1,psi1));
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (fit_td) fitApplyMPO(psi1,expH,psi1,args);
            else psi1 = exactApplyMPO(expH,psi1,args);
            normalize(psi1);
            psi1 *= norm0;
            out.correlator.push_back(overlapC(psi2,psi1));
            }
        out.final_wf = psi1;
        return out;
        }

    // Real-time evolution of an explicit wavefunction, measuring an
    // operator's expectation value at each step, mirroring
    // time_evolution.h's evolution_measure() task exactly (no energy
    // shift here, unlike quench() -- matches the original, which builds
    // expH straight from get_ampo(sites) with no "-EGS" term added).
    TimeEvolutionResult
    evolve_and_measure(std::vector<MOTerm> const& terms_h,
                        std::vector<MOTerm> const& terms_op,
                        MPS const& wf, int nt, double dt, bool fit_td=true)
        {
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        auto ampo = build_ampo(sites_,terms_h);
        auto expH = evoloperator(MPO(ampo),dt);
        auto A = build_mpo(sites_,terms_op,mpomaxm_);
        auto psi = wf;
        TimeEvolutionResult out;
        for (int it=0;it<nt;it++)
            {
            if (fit_td) fitApplyMPO(psi,expH,psi,args);
            else psi = exactApplyMPO(expH,psi,args);
            out.correlator.push_back(overlapC(psi,A,psi));
            }
        out.final_wf = psi;
        return out;
        }

    // Correction-vector-method (CVM) dynamical correlator at a single
    // frequency point, mirroring cvm_dynamical_correlator.h's task plus
    // cvm.h's spectral_function()/bicstab() exactly, but taking the two
    // operators/frequency/energy directly instead of reading
    // dc_multioperator_i.in/dc_multioperator_j.in/tasks.in, and using this
    // Chain's own cached ground state (wf0_) instead of an internal
    // get_gs() call (mirrors kpm_dynamical_correlator's same choice).
    // energy is taken as an explicit parameter (not read from wf0_energy_)
    // to match the original, which takes it from Python's self.e0 rather
    // than recomputing/reusing it internally.
    double
    cvm_dynamical_correlator(std::vector<MOTerm> const& terms_i,
                             std::vector<MOTerm> const& terms_j,
                             double omega, double eta, double energy,
                             double tol, int max_it) const
        {
        if (!have_wf0_) Error("Chain::cvm_dynamical_correlator called before gs_energy");
        auto S1 = build_mpo(sites_,terms_i,mpomaxm_);
        auto S2 = build_mpo(sites_,terms_j,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        const Cplx z(omega+energy,eta);
        auto ampo = AutoMPO(sites_);
        ampo += z,"Id",1;
        auto zId = MPO(ampo);
        auto A = sum(zId,(-1.0)*H_,args);
        auto b = exactApplyMPO(S2,wf0_,args);
        auto x = bicstab(A,b,tol,max_it,args);
        Cplx G = overlapC(wf0_,S1,x);
        return -G.imag()/M_PI;
        }

    // Solves A x = wf for x via BiCGSTAB, mirroring cvm.h's apply_inverse()
    // task (used by mpsalgebra.py's applyinverse_dmrg, e.g. for analytic-
    // continuation dynamical correlators).
    MPS
    apply_inverse(std::vector<MOTerm> const& terms, MPS const& wf,
                  double tol, int max_it) const
        {
        auto A = build_mpo(sites_,terms,mpomaxm_);
        auto args = Args("Cutoff",cutoff_,"Maxm",maxm_);
        return bicstab(A,wf,tol,max_it,args);
        }

    // Fermionic single-particle correlation matrix <Cdag_i C_j>, mirroring
    // correlationmatrix.h's get_correlation_matrix() exactly (same operator
    // names/site convention), just taking wf directly and returning the
    // flat, row-major N*N matrix instead of CORRELATION_MATRIX.OUT (today's
    // O(N^2) lines of text).
    std::vector<std::complex<double>>
    correlation_matrix(MPS const& wf) const
        {
        int N = sites_.N();
        std::vector<std::complex<double>> out(N*N);
        for (int i=0;i<N;i++)
            {
            for (int j=i;j<N;j++)
                {
                auto ampo = AutoMPO(sites_);
                ampo += 1.0,"Cdag",i+1,"C",j+1;
                auto op = MPO(ampo);
                auto c = overlapC(wf,op,wf);
                out[i*N+j] = c;
                out[j*N+i] = std::conj(c);
                }
            }
        return out;
        }

    // Four-fermion correlation tensor <Cdag_i C_j Cdag_k C_l>, mirroring
    // fourcorrelationtensor.h's get_four_correlation_tensor() exactly (both
    // the "accelerate" half-loop-plus-mirror mode and the brute-force
    // fallback), returning the flat, row-major N^4 tensor instead of
    // today's O(N^4) lines of text -- the single biggest I/O cost in the
    // whole file-based protocol.
    std::vector<std::complex<double>>
    four_correlation_tensor(MPS const& wf, bool accelerate=true) const
        {
        int N = sites_.N();
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
                    auto op = MPO(ampo);
                    auto c = overlapC(wf,op,wf);
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
                auto op = MPO(ampo);
                auto c = overlapC(wf,op,wf);
                out[idx(i,j,k,l)] = c;
                out[idx(l,k,j,i)] = std::conj(c);
                }
            }
        return out;
        }

    // KPM dynamical correlator: mirrors kpmmoments.h's
    // get_moments_dynamical_correlator(), i.e. the "dynamical_correlator"
    // task that kpmdmrg.py::get_moments_dynamical_correlator_dmrg always
    // drives (Python always supplies two MultiOperators -- the
    // "kpm_multioperator_i/j" branch is the only one ever taken from
    // Python, so the alternative single-named-operator branch is not
    // ported here, matching the plan's "drop rather than port unreachable
    // branches" guidance).
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
        auto psi1 = exactApplyMPO(wf0_,m1,{"Maxm",kpmmaxm,"Cutoff",kpm_cutoff});
        auto psi2 = exactApplyMPO(wf0_,m2,{"Maxm",kpmmaxm,"Cutoff",kpm_cutoff});
        KPMResult out;
        out.moments = kpm_moments(hs.scaled_H,psi1,psi2,n,kpmmaxm,kpm_cutoff,kpm_accelerate);
        out.emin = hs.emin; out.emax = hs.emax; out.scale = hs.scale;
        out.num_polynomials = n;
        return out;
        }

    // Generic KPM moments of an already-scaled operator between two
    // arbitrary wavefunctions, mirroring kpmmoments.h's general_kpm() (the
    // "general_kpm" task that kpmdmrg.py's general_kpm_moments/
    // kpm_moments_wfa_wfb drive). Unlike kpm_dynamical_correlator, the
    // scaling here is already done on the Python side (scale_operator() in
    // kpmdmrg.py), so this only builds the operator MPO and runs the same
    // Chebyshev recursion.
    std::vector<std::complex<double>>
    general_kpm(std::vector<MOTerm> const& terms_x, MPS const& wfa, MPS const& wfb,
                int kpmmaxm, bool kpm_accelerate, int num_polynomials, double kpm_cutoff)
        {
        auto m = build_mpo(sites_,terms_x,mpomaxm_);
        return kpm_moments(m,wfa,wfb,num_polynomials,kpmmaxm,kpm_cutoff,kpm_accelerate);
        }

    private:
    Sweeps
    make_sweeps() const
        {
        auto sweeps = Sweeps(nsweeps_);
        sweeps.maxm() = maxm_;
        sweeps.cutoff() = cutoff_;
        sweeps.noise() = noise_;
        // noise only in the first half, mirrors get_sweeps.h
        for (int i=nsweeps_/2;i<nsweeps_;i++) sweeps.setnoise(i,0.0);
        return sweeps;
        }

    // Ground/anti-ground energy of H_, used to set the KPM energy scale
    // (scaled_hamiltonian below). These replace bandwidth.h's
    // minimum_energy/maximum_energy: that file cached results in
    // process-lifetime function-local statics (safe only because each old
    // mpscpp.x task was a fresh process); here the cache is a member, so
    // two independent Chain instances (different Hamiltonians) never
    // observe each other's cached bandwidth, which a shared global would
    // have caused the moment this became a persistent, multi-session
    // process (the exact hazard flagged in Phase 0).
    //
    // minimum_energy() also deliberately does NOT reproduce
    // bandwidth.h's behavior of calling get_gs_energy(H) -> get_gs(), which
    // reruns a *full second* DMRG from a fresh random state even when a
    // ground state was already computed earlier in the same task -- that
    // was only ever "correct" because the file-based code had no memory of
    // anything computed earlier in the process. Here gs_energy() is already
    // cached on this Chain, so reusing it is strictly more correct, not a
    // behavior change worth preserving.
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
            auto psi = MPS(sites_);
            auto sweeps = make_sweeps();
            auto negH = (-1.0)*H_;
            bandwidth_emax_ = -dmrg(psi,negH,sweeps);
            have_bandwidth_max_ = true;
            }
        return bandwidth_emax_;
        }

    // emax-emin, mirroring bandwidth.h's bandwidth(), using this Chain's
    // own cached minimum_energy()/maximum_energy() rather than the
    // process-global bandwidth_cache in bandwidth.h (which is still used
    // by the old file-based executable and is fine there -- one process
    // per task -- but would leak state across independent chains here).
    double
    bandwidth() { return maximum_energy()-minimum_energy(); }

    // Mirrors get_excited.h's get_energy_fluctuation(): <H^2>-<H>^2 for a
    // (locally normalized) copy of wf -- psi1 is taken by value, matching
    // the original's local-copy-then-normalize semantics.
    double
    energy_fluctuation(MPS psi1, MPO const& H) const
        {
        normalize(psi1);
        auto psi2 = exactApplyMPO(psi1,H,{"Maxm",maxm_,"Cutoff",cutoff_});
        double de = overlap(psi1,psi2);
        de = overlap(psi2,psi2)-de*de;
        return de;
        }

    // Mirrors get_excited.h's gram_schmidt(): orthogonalizes a sequence of
    // MPS against each other in order, using this Chain's own sum_mps
    // instead of the free sum_mps() in mpsalgebra.h.
    std::vector<MPS>
    gram_schmidt(std::vector<MPS> wfs) const
        {
        for (size_t i=1;i<wfs.size();i++)
            {
            for (size_t j=0;j<i;j++)
                {
                auto proj = overlapC(wfs.at(j),wfs.at(i))*wfs.at(j);
                auto wf = sum_mps(wfs.at(i),(-1.0)*proj);
                normalize(wf);
                wfs.at(i) = wf;
                }
            }
        return wfs;
        }

    struct HamiltonianScale { MPO scaled_H; double emin, emax, scale; };

    // Mirrors scalehamiltonian.h's scale_hamiltonian(): shift+rescale H_ so
    // its spectrum fits in [-1,1], required for the Chebyshev/KPM
    // expansion. Builds the shift directly as a small identity MPO added to
    // the already-built H_ (an MPO sum, bond dimension 1 for the identity
    // term) instead of re-deriving the whole Hamiltonian a second time from
    // hamiltonian.in via get_ampo(), which is what the file-based version
    // does.
    HamiltonianScale
    scaled_hamiltonian(double kpm_scale)
        {
        double emin = minimum_energy();
        double emax = maximum_energy();
        double shift = -(emin+emax)/2.0;
        auto ampo = AutoMPO(sites_);
        ampo += shift,"Id",1;
        auto shift_mpo = MPO(ampo);
        auto m = sum(H_,shift_mpo,{"Maxm",mpomaxm_,"Cutoff",cutoff_});
        double scale = (emax-emin)*kpm_scale;
        scale = 1.0/scale;
        m = m*scale;
        HamiltonianScale out; out.scaled_H = m; out.emin = emin; out.emax = emax; out.scale = scale;
        return out;
        }

    // Whether two MPS are numerically the same state, mirroring
    // mpsalgebra.h's same_mps() (used to pick the accelerated single-vector
    // KPM recursion when vi==vj), just taking maxm/cutoff as explicit
    // parameters instead of reading tasks.in.
    bool
    same_mps(MPS const& vi, MPS const& vj, int maxm, double cutoff) const
        {
        auto d = sum(1.0*vi,-1.0*vj,{"Maxm",maxm,"Cutoff",cutoff});
        double dd = sqrt(overlap(d,d));
        return dd<1e-10;
        }

    // MPO sum/product, mirroring mpsalgebra.h's sum_mpo()/mult_mpo() exactly,
    // just taking maxm/mpomaxm/cutoff from this Chain's own members instead
    // of tasks.in. Used by custom_exp() below.
    MPO
    sum_mpo(MPO const& A1, MPO const& A2) const
        {
        return sum(A1,A2,{"Maxm",maxm_,"Cutoff",cutoff_});
        }

    MPO
    mult_mpo(MPO const& A1, MPO const& A2) const
        {
        auto ampo = AutoMPO(sites_);
        auto out = MPO(ampo);
        nmultMPO(A1,A2,out,{"Maxm",mpomaxm_,"Cutoff",cutoff_});
        return out;
        }

    // BiCGSTAB iterative solver for A x = b (MPO/MPS), mirroring cvm.h's
    // bicstab() exactly, using this Chain's own conjugate() instead of the
    // free conjMPS() (both do the same thing -- conjugate every tensor of
    // a copy of the MPS). Used by apply_inverse() and
    // cvm_dynamical_correlator() above.
    MPS
    bicstab(MPO const& A, MPS const& b, double tol, int max_it, Args const& args) const
        {
        MPS x = b;
        MPS r_old = sum(b,(-1.0)*exactApplyMPO(A,x,args),args);
        MPS r_ = r_old;
        MPS p = r_old;
        MPS s, Ap, As, r_new;
        Cplx alpha, beta, w;
        int k = 0;
        while (k<max_it)
            {
            Ap = exactApplyMPO(A,p,args);
            alpha = overlapC(conjugate(r_old),r_) / overlapC(conjugate(Ap),r_);
            s = sum(r_old,(-alpha)*Ap,args);
            As = exactApplyMPO(A,s,args);
            w = overlapC(conjugate(As),s) / overlapC(conjugate(As),As);
            x = sum(x,sum(alpha*p,w*s,args),args);
            r_new = sum(s,(-w)*As,args);
            double res = std::sqrt(std::abs(overlapC(conjugate(r_new),r_new).real()));
            if (res<=tol) break;
            beta = (alpha/w) * overlapC(conjugate(r_new),r_) / overlapC(conjugate(r_old),r_);
            p = sum(r_new,beta*sum(p,(-w)*Ap,args),args);
            r_old = r_new;
            k++;
            }
        return x;
        }

    // 2nd-order Taylor-expanded exp(z*H), mirroring time_evolution.h's
    // custom_exp() exactly (the variant exponential_eMwf() uses when
    // tevol_custom_exp is set, which is Python's default and only path --
    // see exponential_apply() above).
    MPO
    custom_exp(MPO const& H, Cplx z) const
        {
        auto ampo = AutoMPO(sites_);
        ampo += 1.0,"Id",1;
        auto Iden = MPO(ampo);
        auto out = sum_mpo(Iden,z*H);
        auto H2 = mult_mpo(H,H);
        out = sum_mpo(out,(0.5*z*z)*H2);
        return out;
        }

    // exp(-i*dt*H) Taylor-expanded to (nominally) 3rd order, mirroring
    // time_evolution.h's evoloperator() -- used by quench()/
    // evolve_and_measure() above -- *verbatim*, including a latent bug in
    // the original: H3 (H^3) is computed but the z^3/6 term multiplies H2
    // again instead of H3, so this is not actually a 3rd-order expansion
    // despite the name. Not fixed here: this migration changes how
    // Python and C++ exchange data, not the numerics any existing
    // calculation depends on, and silently correcting this would make
    // old-vs-new comparisons diverge for a reason unrelated to the
    // migration itself.
    MPO
    evoloperator(MPO const& H, double dt) const
        {
        auto ampo = AutoMPO(sites_);
        ampo += 1.0,"Id",1;
        Cplx z(0.0,-dt);
        auto Iden = MPO(ampo);
        auto out = sum_mpo(Iden,z*H);
        auto H2 = mult_mpo(H,H);
        auto H3 = mult_mpo(H,H2); // computed to match the original; unused below, see note
        (void)H3;
        out = sum_mpo(out,(0.5*z*z)*H2);
        out = sum_mpo(out,(z*z*z/6.0)*H2); // NOTE: original uses H2 here, not H3
        return out;
        }

    // Chebyshev moment recursion for two distinct vectors, mirroring
    // kpmcorrelator.h's moments_vi_vj_full() (minus the entanglement-entropy
    // side channel written to KPM_ENTROPY.OUT, which no Python call site
    // reads back).
    std::vector<std::complex<double>>
    kpm_moments_full(MPO const& m, MPS const& vi, MPS const& vj, int n,
                      int kpmmaxm, double kpmcutoff) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto v = 1.0*vi;
        auto am = 1.0*vi;
        auto a = exactApplyMPO(v,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
        auto ap = 1.0*a;
        out.push_back(overlapC(vj,v));
        out.push_back(overlapC(vj,a));
        for (int i=0;i<n;i++)
            {
            ap = exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
            out.push_back(overlapC(vj,ap));
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    // Chebyshev moment recursion when vi==vj, mirroring kpmcorrelator.h's
    // moments_vi_accelerated() (roughly 2x cheaper since it only ever
    // propagates one vector).
    std::vector<std::complex<double>>
    kpm_moments_accelerated(MPO const& m, MPS const& vi, int n,
                            int kpmmaxm, double kpmcutoff) const
        {
        std::vector<std::complex<double>> out;
        out.reserve(n+2);
        auto a = exactApplyMPO(vi,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
        auto am = 1.0*vi;
        auto ap = 1.0*a;
        Cplx mu0 = overlapC(vi,vi);
        Cplx mu1 = overlapC(vi,a);
        out.push_back(mu0);
        out.push_back(mu1);
        for (int i=0;i<n/2;i++)
            {
            ap = exactApplyMPO(a,m,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
            ap = sum(2.0*ap,-1.0*am,{"Maxm",kpmmaxm,"Cutoff",kpmcutoff});
            Cplx bk = 2.0*overlapC(a,a) - mu0;
            Cplx bk1 = 2.0*overlapC(a,ap) - mu1;
            out.push_back(bk);
            out.push_back(bk1);
            am = 1.0*a;
            a = 1.0*ap;
            }
        return out;
        }

    // Dispatcher, mirroring kpmmoments.h's moments_vi_vj(): picks the
    // accelerated recursion when the two vectors coincide and accelerate is
    // requested, otherwise falls back to the general two-vector recursion.
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
    };
