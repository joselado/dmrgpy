// Real-time evolution via 2nd-order Trotter TEBD for strictly nearest-
// neighbor Hamiltonians -- a C++ port of pyitensor/tebd.py's TEBDEvolver
// onto ITensor v3's own BondGate (itensor/mps/bondgate.h, pulled in
// transitively via itensor/mps/tevol.h) and svd() primitives. See
// pyitensor/tebd.py's module docstring for the full algorithm/rationale
// (bond Hamiltonians, half-dt-odd/full-dt-even/half-dt-odd Trotter
// splitting, half-weight onsite-term splitting onto neighboring bonds,
// why fermionic terms need no special-casing beyond correct per-site
// operator resolution); this file only documents where the C++ port's
// *method* of getting each term's per-site operator diverges from the
// pure-Python reference.
//
// pyitensor/tebd.py resolves each MOTerm's per-site Jordan-Wigner-threaded
// operator by hand (its own HTerm.resolve(), mirroring autompo.cc's
// rewriteFermionic()). Here, that same algorithm is ported directly --
// same carried-parity/JW-threading logic, same per-site factor
// composition order -- rather than routed through ITensor's own
// AutoMPO/toMPO() automaton: a probe against real ITensor v3 confirmed a
// *single*, unsummed HTerm does NOT compile to a bond-dimension-1 MPO (a
// bare one-site "Sz" term already comes out with uniform bond dimension 2
// on every link, presumably automaton bookkeeping state rather than
// physical content), so there is no cheap, generic way to "slice out" a
// term's own local operator from its compiled MPO the way this file's
// other extraction tricks (e.g. Chain::build_product_state()'s setElt()
// idiom) work for genuinely bond-dimension-1 objects. Composing same-site
// factors uses ITensor's own multSiteOps(A,B) ("apply B then A", i.e.
// A@B) -- confirmed directly (multSiteOps(Cdag,C) == N on a Fermion
// site) -- instead of hand-rolled prime-level bridging.
bool inline
is_fermionic_opname(std::string const& name)
    {
    return !name.empty() && name[0]=='C';
    }

// The per-site standard-convention operator (out=primed,in=unprimed) for
// one MOTerm at one site, given the fermion parity carried in from sites
// <site -- mirrors pyitensor/autompo.py's HTerm.resolve() body exactly.
// Returns (operator, outgoing carry).
std::pair<ITensor,bool>
resolve_term_at_site(SiteSet const& sites, MOTerm const& term, int site, bool carry_in)
    {
    std::vector<std::string> names;
    for (auto const& f : term.factors)
        if (f.site==site) names.push_back(f.name);
    if (names.empty())
        return {sites.op(carry_in ? "F" : "Id", site), carry_in};
    bool is_site_fermionic = false;
    for (auto const& nm : names)
        if (is_fermionic_opname(nm)) is_site_fermionic = !is_site_fermionic;
    bool need_F = (carry_in != is_site_fermionic);
    ITensor combined = sites.op(names[0],site);
    for (size_t k=1; k<names.size(); ++k)
        combined = multSiteOps(combined,sites.op(names[k],site));
    if (need_F) combined = multSiteOps(combined,sites.op("F",site));
    return {combined, carry_in != is_site_fermionic};
    }

// The 2-site bond Hamiltonians h_1..h_{n-1} (bond b couples sites b,b+1),
// each an ITensor with indices s(b)',s(b),s(b+1)',s(b+1) (the standard
// output-primed/input-unprimed convention op(sites,name,site) itself
// uses, e.g. see ITensor's own tutorial/05_gates/gates.cc) -- ready to
// hand to BondGate's tReal/tImag constructor unchanged.
std::map<int,ITensor>
bond_hamiltonians(SiteSet const& sites, std::vector<MOTerm> const& terms)
    {
    int n = sites.length();
    if (n<2) throw ITError("bond_hamiltonians: TEBD needs at least 2 sites");
    std::map<int,ITensor> H;
    for (int b=1; b<n; ++b)
        H[b] = (sites.op("Id",b)*sites.op("Id",b+1)) * 0.0;

    for (auto const& term : terms)
        {
        std::vector<ITensor> mats(n+1); // 1-indexed, [0] unused
        bool carry = false;
        for (int i=1; i<=n; ++i)
            {
            auto resolved = resolve_term_at_site(sites,term,i,carry);
            mats[i] = resolved.first;
            carry = resolved.second;
            }

        // Touched span: where the resolved per-site operator differs
        // from plain Id -- mirrors tebd.py's _true_span() (determining
        // this from the *resolved* matrices, not the raw op-name list,
        // is what makes it immune to both of _true_span()'s documented
        // false-positive sources: MultiOperator.identity()'s placeholder
        // "Id" factor, and paired F-string factors that algebraically
        // cancel to Id). The threshold is matched to _true_span()'s own
        // np.allclose(m,eye) default (atol=1e-8,rtol=1e-5, i.e. an
        // effective tolerance of ~1e-8 to ~1e-5 depending on entry
        // magnitude) rather than a tighter one: this is a Frobenius norm
        // over the whole small (out,in) operator, not an element-wise
        // check, but staying below Python's own loosest bound keeps a
        // term that's borderline-touched from being classified
        // differently (nearest-neighbor vs "reject as long-range")
        // between the two backends.
        int lo=-1, hi=-1;
        for (int i=1; i<=n; ++i)
            {
            if (norm(mats[i]-sites.op("Id",i)) > 1E-8)
                {
                if (lo==-1) lo = i;
                hi = i;
                }
            }

        if (lo==-1)
            {
            // Untouched everywhere (a bare scalar*Identity term, e.g. a
            // ground-state-energy shift folded in as coef*Id): dump the
            // whole contribution on bond 1, exactly like tebd.py -- any
            // other bond would do equally well (sum over bonds is what
            // matters for Trotter correctness, not which bond carries a
            // global-phase constant), bond 1 is just tebd.py's own
            // convention, kept here for consistency.
            H[1] += term.coef * (sites.op("Id",1)*sites.op("Id",2));
            }
        else if (lo==hi)
            {
            // Onsite term: split half-and-half onto both neighboring
            // bonds (full weight at a chain boundary) -- TeNPy's
            // NearestNeighborModel convention, same as tebd.py.
            auto local = term.coef * mats[lo];
            bool has_left = (lo>1), has_right = (lo<n);
            double weight = (has_left && has_right) ? 0.5 : 1.0;
            if (has_left) H[lo-1] += weight * (sites.op("Id",lo-1) * local);
            if (has_right) H[lo] += weight * (local * sites.op("Id",lo+1));
            }
        else if (hi==lo+1)
            {
            H[lo] += term.coef * (mats[lo] * mats[hi]);
            }
        else
            {
            throw ITError("TEBD requires a nearest-neighbor Hamiltonian; "
                "found a term with non-trivial support spanning sites "
                +std::to_string(lo)+".."+std::to_string(hi));
            }
        }
    return H;
    }

// Precomputed odd/even-bond gates for one fixed dt -- built once from the
// (static, time-independent) bond Hamiltonians and reused unchanged for
// every subsequent tebd_step() call, mirroring TEBDEvolver.__init__: the
// odd bonds only ever need their half-dt gate, the even bonds only ever
// need their full-dt gate (see tebd_step()'s own comment), so only those
// are built. ITensor's BondGate computes exp(-i*tau*bondH) via its own
// Taylor series (bondgate.h, order 100) -- no expm() port needed here.
struct TebdGates
    {
    std::vector<int> odd, even;
    std::map<int,BondGate> half; // odd bonds, dt/2
    std::map<int,BondGate> full; // even bonds, dt
    };

TebdGates
build_tebd_gates(SiteSet const& sites, std::map<int,ITensor> const& h_bonds, double dt)
    {
    int n = sites.length();
    TebdGates g;
    for (int b=1; b<n; b+=2) g.odd.push_back(b);
    for (int b=2; b<n; b+=2) g.even.push_back(b);
    for (int b : g.odd)
        g.half.emplace(b,BondGate(sites,b,b+1,BondGate::tReal,dt/2.0,h_bonds.at(b)));
    for (int b : g.even)
        g.full.emplace(b,BondGate(sites,b,b+1,BondGate::tReal,dt,h_bonds.at(b)));
    return g;
    }

// Contracts `gate` onto bond (b,b+1) and SVD-truncates. Canonicalizes
// psi's orthogonality center to site b first (position(), the same
// primitive every other bond-update in this file -- e.g. bond_entropy()
// above -- uses), which is what makes the local SVD truncation-optimal
// regardless of which bond was last touched. U is pre-seeded from psi's
// own (pre-gate) tensor at site b, the same "seed the wanted left/right
// split via a template ITensor" idiom bond_entropy() already uses above,
// safe here because the gate's own noPrime() restores site b's physical
// index to that very same Index object. psi.set() (unlike position()'s
// own internal orthMPS()/svdBond() sweep helpers) does not update MPS's
// own l_orth_lim_/r_orth_lim_ orthogonality-center bookkeeping, so it is
// set explicitly here to (b,b+2) -- the (leftLim,rightLim) convention
// position() itself uses for a center at site b+1 (confirmed by reading
// MPS::position()'s own implementation) -- mirroring tebd.py's explicit
// `psi.center = i+1` after the same U/S*V split. Skipping this left the
// MPS's ortho-center metadata stuck wherever it was before this call,
// confirmed directly: MPS::normalize() (called after every tebd_step() in
// quench_tebd()) aborted with "MPS must have well-defined ortho center to
// compute norm" without it. Mutates psi in place.
void
apply_bond_gate(MPS& psi, int b, BondGate const& gate, double cutoff, int maxdim)
    {
    auto args = Args("Cutoff",cutoff,"MaxDim",maxdim);
    psi.position(b,args);
    ITensor U = psi.A(b);
    ITensor theta = noPrime(gate.gate() * (psi.A(b)*psi.A(b+1)));
    ITensor S,V;
    svd(theta,U,S,V,args);
    psi.set(b,U);
    psi.set(b+1,S*V);
    psi.leftLim(b);
    psi.rightLim(b+2);
    }

// One full step of size dt: exp(-i dt/2 H_odd) exp(-i dt H_even)
// exp(-i dt/2 H_odd) -- H_odd/H_even being the sums over odd-indexed and
// even-indexed bonds respectively, each internally commuting since every
// bond in one group acts on disjoint sites -- mirrors TEBDEvolver.step().
// Only ever applies the half-dt gate on odd bonds and the full-dt gate on
// even bonds, matching what build_tebd_gates() above actually builds.
// Mutates psi in place.
void
tebd_step(MPS& psi, TebdGates const& gates, double cutoff, int maxdim)
    {
    for (int b : gates.odd) apply_bond_gate(psi,b,gates.half.at(b),cutoff,maxdim);
    for (int b : gates.even) apply_bond_gate(psi,b,gates.full.at(b),cutoff,maxdim);
    for (int b : gates.odd) apply_bond_gate(psi,b,gates.half.at(b),cutoff,maxdim);
    }
