// Controlled bond expansion (CBE) for one-site DMRG and one-site TDVP on
// ITensor v3 -- a port of pyitensor/cbe.py. See that module's docstring
// for the algorithm, the references (Gleis, Li, von Delft, PRL 130, 246402
// (2023); Li, Gleis, von Delft, arXiv:2208.10972; FiniteMPS.jl, MIT) and
// the choice of thresholds; this file only says where the C++ differs.
//
// The kernels of cbe.py are written on dense arrays with fixed axis
// orders. Here they are written on ITensors, where an index is matched by
// identity rather than by position, so one function serves both sides of
// a bond: expand_isometry() is given the isometry to grow (A_l on a
// right-to-left sweep, B_{l+1} on a left-to-right one), the center next to
// it, the two environments and the two MPO tensors, and never asks which
// side it is on. cbe.py's expand_right_isometry() reflects the arrays to
// reach the same symmetry.
//
// Dense indices only. Keeping the D~ largest singular vectors is a column
// selection on the SVD index, which for a QN index would have to be done
// block by block; the session refuses CBE in a conserved sector instead
// (Chain::reject_sector), where cbe.py works on its charge grading.
//
// Environments follow ITensor's own convention,
//   L_i = L_{i-1} * psi(i) * H(i) * dag(prime(psi(i)))
// so an environment has the ket link, the MPO link and the primed bra link,
// and a null ITensor stands for the trivial environment at a chain end.
#pragma once

#include <cmath>
#include <vector>

namespace cbe {

// Singular values below this fraction of the scale of the object they
// come from are numerical zeros (cbe.py's _NUMERICAL_ZERO).
static const double NUMERICAL_ZERO = 1e-12;
// The thresholds of the CBE-TDVP paper, preselection and final selection.
static const double TDVP_EPS_PRESELECT = 1e-4;
static const double TDVP_EPS_FINAL = 1e-6;

// ------------------------------------------------------------------------
// small dense helpers
// ------------------------------------------------------------------------

inline std::vector<Index>
as_vector(IndexSet const& is)
    {
    std::vector<Index> v;
    for (auto const& i : is) v.push_back(i);
    return v;
    }

// The indices of T other than those in `drop`, as an IndexSet.
inline IndexSet
other_inds(ITensor const& T, std::vector<Index> const& drop)
    {
    std::vector<Index> keep;
    for (auto const& i : inds(T))
        {
        bool skip = false;
        for (auto const& d : drop) if (i==d) skip = true;
        if (!skip) keep.push_back(i);
        }
    return IndexSet(keep);
    }

// The weighted left factor U*S of T, U carrying `left`, from the
// eigendecomposition of the Gram matrix T T^dag instead of a full SVD,
// with the singular values s descending. Only for a factor used as
// weights: the Gram route resolves singular values down to about 1e-8 of
// the largest, below which the columns are noise of that size, harmless
// in a weight and not in an isometry (see cbe.py's svd()). Returns US with
// `left` and the new index u.
struct WeightedSplit
    {
    ITensor US;
    Index u;
    std::vector<double> s;
    };

inline WeightedSplit
gram_weighted_left(ITensor const& T, IndexSet const& left)
    {
    auto G = T*dag(prime(T,left)); // (left, left') = T T^dag
    ITensor U, D;
    // M == dag(U)*D*prime(U) (decomp.h), so dag(U) holds the left
    // singular vectors of T
    diagHermitian(G,U,D);
    auto ud = commonIndex(U,D);
    auto uo = uniqueIndex(D,U);
    WeightedSplit out;
    for (auto j : range1(dim(ud)))
        {
        double ev = elt(D,ud=j,uo=j);
        out.s.push_back(ev>0 ? std::sqrt(ev) : 0.0);
        }
    D.apply([](Real x){ return x>0 ? std::sqrt(x) : 0.0; });
    out.US = dag(U)*D;
    out.u = uo;
    return out;
    }

// A full, untruncated SVD of T with U carrying left_inds, and the
// singular values on the host, descending.
struct Split
    {
    ITensor U, S, V;
    Index u, v; // U's and V's bond index (S carries both)
    std::vector<double> s;
    };

inline Split
full_svd(ITensor const& T, IndexSet const& left_inds)
    {
    Split out;
    out.U = ITensor(left_inds);
    // no Cutoff/MaxDim: svd() then does not truncate (svd.cc)
    svd(T,out.U,out.S,out.V);
    out.u = commonIndex(out.U,out.S);
    out.v = commonIndex(out.S,out.V);
    for (auto j : range1(dim(out.u))) out.s.push_back(elt(out.S,out.u=j,out.v=j));
    return out;
    }

// How many of the descending singular values s to keep: those above
// max(eps,NUMERICAL_ZERO)*scale, at most cap (cap<0: no cap).
inline int
count_kept(std::vector<double> const& s, double scale, double eps, int cap)
    {
    if (scale<=0 || s.empty()) return 0;
    double floor = std::max(eps,NUMERICAL_ZERO)*scale;
    int keep = 0;
    for (auto x : s) if (x>floor) keep++;
    if (cap>=0 && keep>cap) keep = cap;
    return keep;
    }

// T with its index i replaced by a new index of dimension r, keeping the
// first r values of i. Returns the new index through `out`.
inline ITensor
keep_first(ITensor const& T, Index const& i, int r, Index& out)
    {
    out = Index(r,"Link,CBE");
    auto P = ITensor(i,out);
    for (auto j : range1(r)) P.set(i=j,out=j,1.0);
    return T*P;
    }

// U (outer x u) with the span of A's columns (outer x bond) removed, and
// made orthonormal again; numerically zero columns dropped (cbe.py's
// _orthonormal_complement_columns, the paper's step 14 done twice).
// `outer` are the indices U and A share. Returns a null ITensor when
// nothing is left.
inline ITensor
orthonormal_complement(ITensor U, ITensor const& A, IndexSet const& outer,
                       double scale=1.0)
    {
    for (int k=0;k<2;k++) U -= A*(dag(A)*U);
    auto sp = full_svd(U,outer);
    int r = count_kept(sp.s,scale,0.0,-1); // scale: U's largest singular value, or a bound on it
    if (r==0) return ITensor();
    Index q;
    return keep_first(sp.U,sp.u,r,q);
    }

// ------------------------------------------------------------------------
// environments and effective Hamiltonians
// ------------------------------------------------------------------------

inline ITensor
extend(ITensor const& E, ITensor const& M, ITensor const& W)
    {
    ITensor out = E ? E*M : M;
    out *= W;
    out *= dag(prime(M));
    return out;
    }

// One-site effective Hamiltonian L*W*R (davidson()/applyExp() BigMatrix).
struct OneSiteH
    {
    ITensor const& L;
    ITensor const& W;
    ITensor const& R;
    size_t n;
    OneSiteH(ITensor const& L_, ITensor const& W_, ITensor const& R_, ITensor const& phi)
        : L(L_), W(W_), R(R_), n(dim(inds(phi))) {}
    void
    product(ITensor const& phi, ITensor& out) const
        {
        out = L ? L*phi : phi;
        out *= W;
        if (R) out *= R;
        out.noPrime();
        }
    size_t size() const { return n; }
    };

// Zero-site (bond) effective Hamiltonian L*R.
struct ZeroSiteH
    {
    ITensor const& L;
    ITensor const& R;
    size_t n;
    ZeroSiteH(ITensor const& L_, ITensor const& R_, ITensor const& phi)
        : L(L_), R(R_), n(dim(inds(phi))) {}
    void
    product(ITensor const& phi, ITensor& out) const
        {
        out = L ? L*phi : phi;
        if (R) out *= R;
        out.noPrime();
        }
    size_t size() const { return n; }
    };

// ------------------------------------------------------------------------
// the expansion
// ------------------------------------------------------------------------

// New directions for the isometry `iso` of a bond whose center `cen` sits
// on its other side. Eo/Wo are the environment and MPO tensor on the
// isometry's side, Ei/Wi those on the center's side; `bond` is the link
// iso and cen share. The result carries iso's outer indices and one new
// index, orthonormal and orthogonal to iso's columns; a null ITensor when
// no direction is found. dtilde caps the number of directions, dpre the
// preselection rank ("shrewd" only).
inline ITensor
expand_isometry(ITensor const& Eo, ITensor const& Wo, ITensor const& Wi,
                ITensor const& Ei, ITensor const& iso, ITensor const& cen,
                Index const& bond, int dtilde, int dpre,
                double eps_pre, double eps_final, bool shrewd)
    {
    if (dtilde<=0) return ITensor();
    auto outer = other_inds(iso,{bond});
    // right half projected onto the center's own row space: B, an
    // orthonormal basis of the rows of cen (cbe.py's _row_basis)
    auto spc = full_svd(cen,IndexSet(bond));
    int rb = count_kept(spc.s,spc.s.empty()?0.0:spc.s[0],0.0,-1);
    Index kb;
    auto Bv = keep_first(spc.V,spc.v,rb,kb); // (kb, center outer indices)
    // X = Eo iso Wo: the isometry's half of H|psi>, open at the MPO bond
    // and at `bond`; Y = cen Wi Ei: the center's half
    auto X = Eo ? Eo*iso : iso;
    X *= Wo;
    X.noPrime();
    auto Y = cen*Wi;
    if (Ei) Y *= Ei;
    Y.noPrime();
    auto proj_center = [&](ITensor const& T) { return T - (T*dag(Bv))*Bv; };
    auto proj_iso = [&](ITensor const& T) { return T - iso*(dag(iso)*T); };

    if (!shrewd)
        {
        // H_2s psi, both legs projected onto the discarded spaces, then
        // the SVD across the bond (FiniteMPS.jl's NaiveCBE)
        auto T = X*Y;
        double scale = norm(T);
        T = proj_center(proj_iso(T));
        auto sp = full_svd(T,outer);
        int keep = count_kept(sp.s,scale,eps_final,dtilde);
        if (keep==0) return ITensor();
        Index nu;
        auto U = keep_first(sp.U,sp.u,keep,nu);
        return orthonormal_complement(U,iso,outer);
        }

    // (a) the center's half on (1 - B^dag B), split on `bond`. Only the
    // weights U1 S1 are used, so the Gram route replaces cbe.py's SVD
    // here and in (b): measured on the L=40 Heisenberg chain at maxm=100,
    // these two SVDs were 2.7 s of a 10 s solve, and the energy did not
    // move (1e-14).
    auto Yo = proj_center(Y);
    auto sp1 = gram_weighted_left(Yo,IndexSet(bond));
    int r1 = count_kept(sp1.s,norm(Y),0.0,-1);
    if (r1==0) return ITensor();
    Index n1;
    auto US = keep_first(sp1.US,sp1.u,r1,n1); // (bond, n1)
    // the isometry's half with those weights on the central MPS bond,
    // projected on (1 - A A^dag)
    auto Lt = X*US; // (outer, MPO bond, n1)
    double scale_pre = norm(Lt);
    auto Lo = proj_iso(Lt);
    // (b) truncate the central MPS bond r1 -> D': U2 S2 = Lo V2, with V2
    // the right singular vectors of Lo, the left ones of Lo^dag
    auto sp2 = gram_weighted_left(dag(Lo),IndexSet(n1));
    int r2 = count_kept(sp2.s,scale_pre,eps_pre,dpre);
    if (r2==0) return ITensor();
    Index n2;
    auto V2s = keep_first(sp2.US,sp2.u,r2,n2); // (n1, n2), columns V2 s
    auto inv = ITensor(n2,prime(n2));
    for (auto j : range1(r2)) inv.set(n2=j,prime(n2)=j,1.0/sp2.s[j-1]);
    auto us = (Lo*V2s*inv).noPrime(n2); // (outer, MPO bond, n2)
    // (c) regroup the MPO bond with the D' kept states: the preselected
    // complement A^pr, an orthonormal basis of the column space of us
    // (outer x (MPO bond, n2)) with iso's columns removed. cbe.py takes
    // the SVD of us first and then the complement; the complement's own
    // rank-revealing SVD spans the same space, so one SVD does both here,
    // with numerical zeros judged against the norm of us, which bounds its
    // largest singular value.
    auto Apr = orthonormal_complement(us,iso,outer,norm(us));
    if (!Apr) return ITensor();
    auto q = other_inds(Apr,as_vector(outer))(1);
    // (d) final selection: H|psi> on A^pr x (1 - B^dag B), MPO bond closed
    auto Ct = (dag(Apr)*X)*Y; // (q, center's outer indices)
    double scale_fin = norm(Ct);
    auto Co = proj_center(Ct);
    auto sp4 = full_svd(Co,IndexSet(q));
    int r4 = count_kept(sp4.s,scale_fin,eps_final,dtilde);
    if (r4==0) return ITensor();
    Index n4;
    return Apr*keep_first(sp4.U,sp4.u,r4,n4);
    }

// (dtilde, dpre) for bond b of an MPO H at current bond dimension k:
// expand to D_ex = ceil((1+delta)*maxdim), preselect D' = ceil(D_ex/w).
inline std::pair<int,int>
target_dims(MPO const& H, int b, int k, int maxdim, double delta)
    {
    int d_ex = (int)std::ceil((1.0+delta)*maxdim);
    int w = dim(linkIndex(H,b));
    int dpre = (int)std::ceil(d_ex/(double)w);
    return {std::max(0,d_ex-k),std::max(1,dpre)};
    }

// Grow bond (b,b+1) of psi by the directions expand_isometry() finds for
// the isometry at site `iso_site` (b+1 on a left-to-right sweep, b on a
// right-to-left one); the center tensor at the other site gets zero
// entries, so psi itself does not change. Eo/Ei as in expand_isometry().
// Returns the number of directions added.
inline int
expand_bond(MPS& psi, MPO const& H, int b, int iso_site,
            ITensor const& Eo, ITensor const& Ei, int maxdim, double delta,
            double eps_pre, double eps_final, bool shrewd)
    {
    int cen_site = (iso_site==b) ? b+1 : b;
    auto bond = linkIndex(psi,b);
    int k = dim(bond);
    auto dims = target_dims(H,b,k,maxdim,delta);
    if (dims.first==0) return 0;
    auto const& iso = psi(iso_site);
    auto const& cen = psi(cen_site);
    auto At = expand_isometry(Eo,H(iso_site),H(cen_site),Ei,iso,cen,bond,
                              dims.first,dims.second,eps_pre,eps_final,shrewd);
    if (!At) return 0;
    auto outer = other_inds(iso,{bond});
    auto nt_idx = other_inds(At,as_vector(outer))(1);
    int nt = dim(nt_idx);
    // embed bond and nt_idx into one new link (directSumITensors' dense
    // branch, itensor.cc)
    auto nl = Index(k+nt,tags(bond));
    auto E1 = itensor::delta(bond,nl);
    auto S = Matrix(nt,k+nt);
    for (auto j : range(nt)) S(j,k+j) = 1;
    auto E2 = matrixITensor(std::move(S),nt_idx,nl);
    auto iso_ex = iso*E1 + At*E2;
    auto cen_ex = cen*E1;
    psi.ref(iso_site) = iso_ex;
    psi.ref(cen_site) = cen_ex;
    return nt;
    }

// ------------------------------------------------------------------------
// CBE-DMRG
// ------------------------------------------------------------------------

// Lowest eigenpair of the one-site effective Hamiltonian, by davidson()
// as ITensor's own dmrg() runs it, from phi with a random part of
// norm sqrt(noise)*|phi| (cbe.py's _noisy_start). phi is overwritten.
inline double
local_ground_state(ITensor const& L, ITensor const& W, ITensor const& R,
                   ITensor& phi, int niter, double noise)
    {
    if (noise>0)
        {
        auto r = random(phi,{"Complex",isComplex(phi)}); // a complex draw would make a real solve complex, 4x the cost
        double nr = norm(r), np = norm(phi);
        if (nr>0) phi += (std::sqrt(noise)*np/nr)*r;
        }
    OneSiteH Heff(L,W,R,phi);
    return davidson(Heff,phi,{"MaxIter",niter,"DebugLevel",-1});
    }

struct SweepParams
    {
    int maxdim;
    double cutoff;
    double noise;
    int niter;
    };

// One left-to-right then right-to-left CBE sweep. psi enters with its
// center on site 1 and Renv valid for it, and leaves the same way.
// Lenv[i] is the environment through site i, Renv[i] the one from site
// i; Lenv[0] and Renv[n+1] are null. Returns the last local eigenvalue.
inline double
cbe_sweep(MPS& psi, MPO const& H, std::vector<ITensor>& Lenv,
          std::vector<ITensor>& Renv, SweepParams const& p, double delta,
          bool shrewd)
    {
    int n = length(psi);
    double energy = 0;
    auto args = Args("Cutoff",p.cutoff,"MaxDim",p.maxdim);
    for (int i=1;i<n;i++)
        {
        expand_bond(psi,H,i,i+1,Renv[i+2],Lenv[i-1],p.maxdim,delta,0.0,0.0,shrewd);
        auto R1 = extend(Renv[i+2],psi(i+1),H(i+1));
        auto phi = psi(i);
        energy = local_ground_state(Lenv[i-1],H(i),R1,phi,p.niter,p.noise);
        ITensor U = (i>1) ? ITensor(leftLinkIndex(psi,i),siteIndex(psi,i))
                          : ITensor(siteIndex(psi,i));
        ITensor S, V;
        svd(phi,U,S,V,args);
        psi.ref(i) = U;
        psi.ref(i+1) = (S*V)*psi(i+1);
        Lenv[i] = extend(Lenv[i-1],psi(i),H(i));
        }
    for (int i=n;i>1;i--)
        {
        expand_bond(psi,H,i-1,i-1,Lenv[i-2],Renv[i+1],p.maxdim,delta,0.0,0.0,shrewd);
        auto L1 = extend(Lenv[i-2],psi(i-1),H(i-1));
        auto phi = psi(i);
        energy = local_ground_state(L1,H(i),Renv[i+1],phi,p.niter,p.noise);
        ITensor U(leftLinkIndex(psi,i)), S, V;
        svd(phi,U,S,V,args);
        psi.ref(i) = V;
        psi.ref(i-1) = psi(i-1)*(U*S);
        Renv[i] = extend(Renv[i+1],psi(i),H(i));
        }
    return energy;
    }

inline void
right_environments(MPS const& psi, MPO const& H, std::vector<ITensor>& Renv)
    {
    int n = length(psi);
    Renv.assign(n+2,ITensor());
    for (int i=n;i>=2;i--) Renv[i] = extend(Renv[i+1],psi(i),H(i));
    }

inline void
left_environments(MPS const& psi, MPO const& H, std::vector<ITensor>& Lenv)
    {
    int n = length(psi);
    Lenv.assign(n+2,ITensor());
    for (int i=1;i<n;i++) Lenv[i] = extend(Lenv[i-1],psi(i),H(i));
    }

// Ground-state DMRG with one-site updates and controlled bond expansion,
// cbe.py's dmrg_cbe(). The noise goes into the first sweep only (it is
// passed in only for a fresh start). Mutates psi and returns the last
// local eigenvalue.
inline double
dmrg_cbe(MPS& psi, MPO const& H, std::vector<SweepParams> const& sweeps,
         double delta, bool shrewd, bool quiet)
    {
    int n = length(psi);
    psi.position(1);
    std::vector<ITensor> Lenv, Renv;
    right_environments(psi,H,Renv);
    Lenv.assign(n+2,ITensor());
    double energy = 0;
    for (size_t sw=0;sw<sweeps.size();sw++)
        {
        energy = cbe_sweep(psi,H,Lenv,Renv,sweeps[sw],delta,shrewd);
        if (!quiet) printfln("CBE sweep %d: energy = %.14f, maxdim = %d",
                             (int)sw+1,energy,maxLinkDim(psi));
        }
    psi.leftLim(0);
    psi.rightLim(2);
    return energy;
    }

// ------------------------------------------------------------------------
// CBE-TDVP
// ------------------------------------------------------------------------

inline void
forward(ITensor const& L, ITensor const& W, ITensor const& R, ITensor& phi,
        Cplx tau, int niter)
    {
    OneSiteH Heff(L,W,R,phi);
    applyExp(Heff,phi,Cplx(0,-1)*tau,{"MaxIter",niter});
    }

inline void
backward(ITensor const& L, ITensor const& R, ITensor& phi, Cplx tau, int niter)
    {
    ZeroSiteH Heff(L,R,phi);
    applyExp(Heff,phi,Cplx(0,1)*tau,{"MaxIter",niter});
    }

inline void
half_sweep_lr(MPS& psi, MPO const& H, Cplx tau, double cutoff, int maxdim,
              double delta, int niter, bool shrewd)
    {
    int n = length(psi);
    std::vector<ITensor> Lenv, Renv;
    right_environments(psi,H,Renv);
    Lenv.assign(n+2,ITensor());
    auto args = Args("Cutoff",cutoff,"MaxDim",maxdim);
    for (int i=1;i<=n;i++)
        {
        ITensor Rn;
        if (i<n)
            {
            expand_bond(psi,H,i,i+1,Renv[i+2],Lenv[i-1],maxdim,delta,
                        TDVP_EPS_PRESELECT,TDVP_EPS_FINAL,shrewd);
            Rn = extend(Renv[i+2],psi(i+1),H(i+1));
            }
        auto phi = psi(i);
        forward(Lenv[i-1],H(i),Rn,phi,tau,niter);
        if (i==n) { psi.ref(i) = phi; continue; }
        ITensor U = (i>1) ? ITensor(leftLinkIndex(psi,i),siteIndex(psi,i))
                          : ITensor(siteIndex(psi,i));
        ITensor S, V;
        svd(phi,U,S,V,args);
        auto C = S*V;
        psi.ref(i) = U;
        Lenv[i] = extend(Lenv[i-1],U,H(i));
        backward(Lenv[i],Rn,C,tau,niter);
        psi.ref(i+1) = C*psi(i+1);
        }
    }

inline void
half_sweep_rl(MPS& psi, MPO const& H, Cplx tau, double cutoff, int maxdim,
              double delta, int niter, bool shrewd)
    {
    int n = length(psi);
    std::vector<ITensor> Lenv, Renv;
    left_environments(psi,H,Lenv);
    Renv.assign(n+2,ITensor());
    auto args = Args("Cutoff",cutoff,"MaxDim",maxdim);
    for (int i=n;i>=1;i--)
        {
        ITensor Ln;
        if (i>1)
            {
            expand_bond(psi,H,i-1,i-1,Lenv[i-2],Renv[i+1],maxdim,delta,
                        TDVP_EPS_PRESELECT,TDVP_EPS_FINAL,shrewd);
            Ln = extend(Lenv[i-2],psi(i-1),H(i-1));
            }
        auto phi = psi(i);
        forward(Ln,H(i),Renv[i+1],phi,tau,niter);
        if (i==1) { psi.ref(i) = phi; continue; }
        ITensor U(leftLinkIndex(psi,i)), S, V;
        svd(phi,U,S,V,args);
        auto C = U*S;
        psi.ref(i) = V;
        Renv[i] = extend(Renv[i+1],V,H(i));
        backward(Ln,Renv[i],C,tau,niter);
        psi.ref(i-1) = psi(i-1)*C;
        }
    }

// One step exp(-i*dt*H) of one-site TDVP with controlled bond expansion
// before every local forward step, cbe.py's tdvp_cbe_step(): a
// left-to-right half sweep by dt/2, then a right-to-left one by dt/2. dt
// may be complex. The local exponentials are ITensor's Lanczos applyExp(),
// as in TDVP/tdvp.h, so H is taken as Hermitian (cbe.py switches to
// Arnoldi for a non-Hermitian MPO; v3's two-site TDVP does not either).
inline void
tdvp_cbe_step(MPS& psi, MPO const& H, Cplx dt, double cutoff, int maxdim,
              double delta, int niter, bool shrewd)
    {
    Cplx tau = dt/2.0;
    psi.position(1);
    half_sweep_lr(psi,H,tau,cutoff,maxdim,delta,niter,shrewd);
    half_sweep_rl(psi,H,tau,cutoff,maxdim,delta,niter,shrewd);
    psi.leftLim(0);
    psi.rightLim(2);
    }

// ------------------------------------------------------------------------
// Arnoldi exponential and two-site TDVP for a non-Hermitian MPO
// ------------------------------------------------------------------------

// exp(coeff*A)|phi> for a general BigMatrix A, pyitensor/tdvp.py's
// _arnoldi_expm_multiply(): the whole projection column is kept (two
// classical Gram-Schmidt passes), so the Hessenberg relation holds for a
// non-Hermitian A, where ITensor's Lanczos applyExp() does not converge
// on a Lindbladian. Stops on Saad's estimate relative to the result, and
// sub-steps when niter iterations do not reach the goal.
template<class BigM>
void
arnoldi_exp(BigM const& A, ITensor& phi, Cplx coeff, int niter,
            double errgoal=1e-10, double tol=1e-12, int max_substeps=1000)
    {
    double beta0 = norm(phi);
    if (beta0==0) return;
    int n = (int)A.size();
    int m = std::max(1,std::min(niter,n));
    auto inner = [](ITensor const& a, ITensor const& b)
        { return eltC(dag(a)*b); };
    auto first_column = [](CMatrix const& Hm, int k, Cplx z)
        {
        CMatrix Hk(k,k);
        for (int a=0;a<k;a++) for (int b=0;b<k;b++) Hk(a,b) = Hm(a,b);
        auto E = expMatrix(Hk,z);
        std::vector<Cplx> c(k);
        for (int a=0;a<k;a++) c[a] = E(a,0);
        return c;
        };
    auto estimate_ok = [](Cplx z, double hnext, std::vector<Cplx> const& c,
                          double goal)
        {
        double nc = 0;
        for (auto x : c) nc += std::norm(x);
        return std::abs(z)*hnext*std::abs(c.back()) < goal*std::sqrt(nc);
        };
    ITensor v = phi;
    double remaining = 1.0, frac = 1.0;
    for (int sub=0;sub<max_substeps;sub++)
        {
        std::vector<ITensor> Q;
        Q.push_back(v/beta0);
        CMatrix Hm(m+1,m);
        for (int a=0;a<m+1;a++) for (int b=0;b<m;b++) Hm(a,b) = 0;
        double scale = 0, hnext = 0;
        bool converged = false;
        std::vector<Cplx> c;
        int k = 0;
        for (int it=0;it<m;it++)
            {
            k = it+1;
            ITensor w;
            A.product(Q[it],w);
            std::vector<Cplx> h(k,0.);
            for (int pass=0;pass<2;pass++)
                for (int j=0;j<k;j++)
                    {
                    auto hj = inner(Q[j],w);
                    w -= hj*Q[j];
                    h[j] += hj;
                    }
            hnext = norm(w);
            for (int j=0;j<k;j++)
                {
                Hm(j,it) = h[j];
                scale = std::max(scale,std::abs(h[j]));
                }
            Hm(k,it) = hnext;
            scale = std::max(scale,hnext);
            c = first_column(Hm,k,coeff*frac);
            if (hnext<=tol*scale || k>=n || estimate_ok(coeff*frac,hnext,c,errgoal))
                { converged = true; break; }
            if (k==m) break;
            Q.push_back(w/hnext);
            }
        if (!converged)
            {
            double s = frac;
            bool found = false;
            for (int j=0;j<60;j++)
                {
                s /= 2.0;
                c = first_column(Hm,k,coeff*s);
                if (estimate_ok(coeff*s,hnext,c,errgoal)) { found = true; break; }
                }
            if (!found)
                throw std::runtime_error("arnoldi_exp: no sub-step of the "
                    "unconverged Krylov space meets the error goal");
            frac = s;
            }
        ITensor out = (beta0*c[0])*Q[0];
        for (int j=1;j<k;j++) out += (beta0*c[j])*Q[j];
        v = out;
        if (converged && frac>=remaining) { phi = v; return; }
        remaining -= frac;
        frac = std::min(frac,remaining);
        beta0 = norm(v);
        if (beta0==0) { phi = v; return; }
        }
    throw std::runtime_error("arnoldi_exp: the step needs more than the "
        "allowed number of sub-steps; reduce dt");
    }

// Two-site effective Hamiltonian L*W1*W2*R.
struct TwoSiteH
    {
    ITensor const& L;
    ITensor const& W1;
    ITensor const& W2;
    ITensor const& R;
    size_t n;
    TwoSiteH(ITensor const& L_, ITensor const& W1_, ITensor const& W2_,
             ITensor const& R_, ITensor const& phi)
        : L(L_), W1(W1_), W2(W2_), R(R_), n(dim(inds(phi))) {}
    void
    product(ITensor const& phi, ITensor& out) const
        {
        out = L ? L*phi : phi;
        out *= W1;
        out *= W2;
        if (R) out *= R;
        out.noPrime();
        }
    size_t size() const { return n; }
    };

// One step exp(-i*dt*H) of two-site TDVP whose local exponentials are
// arnoldi_exp(), pyitensor's tdvp_step() for a non-Hermitian MPO: a
// left-to-right half sweep by dt/2 (forward two-site step, SVD at
// maxdim/cutoff, backward one-site step) and its mirror. dt may be
// complex; the norm is not restored.
inline void
tdvp_arnoldi_step(MPS& psi, MPO const& H, Cplx dt, double cutoff,
                  int maxdim, int niter)
    {
    int n = length(psi);
    Cplx tau = dt/2.0;
    Cplx fwd = Cplx(0,-1)*tau, bwd = Cplx(0,1)*tau;
    auto args = Args("Cutoff",cutoff,"MaxDim",maxdim);
    psi.position(1);
    std::vector<ITensor> Lenv, Renv;
    right_environments(psi,H,Renv);
    Lenv.assign(n+2,ITensor());
    for (int b=1;b<n;b++)
        {
        auto sb = siteIndex(psi,b);
        ITensor U = (b>1) ? ITensor(leftLinkIndex(psi,b),sb) : ITensor(sb);
        auto phi = psi(b)*psi(b+1);
        arnoldi_exp(TwoSiteH(Lenv[b-1],H(b),H(b+1),Renv[b+2],phi),phi,fwd,niter);
        ITensor S, V;
        svd(phi,U,S,V,args);
        psi.ref(b) = U;
        Lenv[b] = extend(Lenv[b-1],U,H(b));
        auto C = S*V;
        if (b<n-1)
            arnoldi_exp(OneSiteH(Lenv[b],H(b+1),Renv[b+2],C),C,bwd,niter);
        psi.ref(b+1) = C;
        }
    Renv.assign(n+2,ITensor());
    for (int b=n-1;b>=1;b--)
        {
        auto sb = siteIndex(psi,b);
        ITensor U = (b>1) ? ITensor(leftLinkIndex(psi,b),sb) : ITensor(sb);
        auto phi = psi(b)*psi(b+1);
        arnoldi_exp(TwoSiteH(Lenv[b-1],H(b),H(b+1),Renv[b+2],phi),phi,fwd,niter);
        ITensor S, V;
        svd(phi,U,S,V,args);
        psi.ref(b+1) = V;
        Renv[b+1] = extend(Renv[b+2],V,H(b+1));
        auto C = U*S;
        if (b>1)
            arnoldi_exp(OneSiteH(Lenv[b-1],H(b),Renv[b+1],C),C,bwd,niter);
        psi.ref(b) = C;
        }
    psi.leftLim(0);
    psi.rightLim(2);
    }

} // namespace cbe
