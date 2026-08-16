// Port of mpscpp2/mo_terms.h to the ITensor v3 API. The only actual change
// is toMPO(): v2 had a single-tensor-type-vs-QN-tensor-type template
// (toMPO<ITensor>(...)) because ITensor and IQTensor were distinct C++
// types; v3 merged them into one ITensor type (QN blocks live on the Index
// itself now), so toMPO() is no longer templated at all.
//
// In-memory replacement for the "ampotk" file format (get_ampo_operator.h +
// the generated extra/ampotk/ampotk.h switch). A MultiOperator on the
// Python side is already exactly a list of (coefficient, [(opname,site),...])
// terms (see multioperator.py's MultiOperator.op) -- this is the same shape,
// just handed over as in-memory objects instead of being serialized to text
// and reparsed. Builds directly on ITensor's own AutoMPO::add(HTerm)/
// HTerm::add(name,site) API, which the old ampotk switch was already calling
// indirectly through the "ampo += coef,name,site,..." accumulator syntax --
// no 99-way-switch code generation is needed at all.

#include <set>
#include <string>
#include <vector>

struct OpFactor
    {
    std::string name;
    int site; // 1-based, matching ITensor's AutoMPO/HTerm convention
    };

struct MOTerm
    {
    Cplx coef;
    std::vector<OpFactor> factors;
    };

// -- "Sy" elimination: build a *real*-valued MPO whenever the operator
// actually is real ------------------------------------------------------
//
// ITensor's toMPO() (autompo.cc's svdMPO) picks real-valued MPO tensors
// only when every AutoMPO *coefficient* is real -- but the per-site
// operator matrices go in unexamined, and "Sy" is purely imaginary
// (spinhalf.h et al. set +-0.5i). So the textbook dmrgpy Heisenberg
// Hamiltonian, Sx.Sx + Sy.Sy + Sz.Sz with coefficient 1, produces a
// *complex* MPO even though the operator it represents is real in the Sz
// basis. Everything downstream then inherits that: DMRG's ground state
// comes out complex, and every later applyMPO/sum/inner runs in complex
// arithmetic (zgemm, 4 real multiplies per element) for no reason at all.
// Measured directly on an L=30 S=1/2 Heisenberg chain at bond dimension
// 50 -- the KPM Chebyshev recursion's own working point -- against the
// identical Hamiltonian written with S+/S- instead:
//
//     applyMPO   0.316 s (complex)  vs  0.072 s (real)   4.4x
//     sum(MPS)   0.104 s            vs  0.035 s          3.0x
//     inner      5.5e-3 s           vs  1.7e-3 s         3.2x
//
// Sy always appears an even number of times in a real Hamiltonian, and
// Sy = (S+ - S-)/(2i) exactly (any spin S, not just 1/2), so such a term
// can be expanded into 2^k products of the *real* operators S+/S- with a
// real overall coefficient -- i^(even) is real. Sx = (S+ + S-)/2 is
// expanded in the same pass, not because Sx is complex (it isn't) but
// because doing both is what lets AutoMPO see the exact cancellation
// Sx.Sx + Sy.Sy = (S+.S- + S-.S+)/2: expanding Sy alone leaves S+S+/S-S-
// strings that don't cancel and *grow* the MPO bond dimension (5 -> 6),
// eating much of the win.
//
// This is a representation change only -- the operator, and hence every
// number dmrgpy reports, is mathematically identical. It is applied only
// when all of the following hold (otherwise the terms are used verbatim):
//
//   * some term actually contains an "Sy" (nothing to gain otherwise),
//   * every coefficient is already real (a genuinely complex coefficient
//     makes the MPO complex regardless),
//   * every term has an *even* number of "Sy" factors (an odd one is a
//     genuinely imaginary operator -- e.g. a bare <Sy> observable -- and
//     must stay complex),
//   * no term carries more than MAX_XY_FACTORS Sx/Sy factors (the
//     expansion is 2^k terms; this bounds it at 16),
//   * every site carrying an Sx/Sy factor actually defines "S+"/"S-".
//
// The last check is what keeps this away from site types where an
// "imaginary" operator is not a spin at all (Z3/Z4 clock operators, whose
// complexity is genuine and irremovable) -- those never define S+/S-.
// SiteSet::op() throws ITError for an unrecognized name, so the probe is a
// real try/catch, not a name whitelist that would silently rot as site
// types are added.

// Process-global opt-out, for tests that want to prove the rewritten and
// verbatim paths agree numerically (see bindings.cc's
// set_realify_spin_terms). On by default -- it is a pure speedup.
inline bool&
realify_spin_terms_enabled()
    {
    static bool enabled = true;
    return enabled;
    }

namespace motermsdetail {

constexpr int MAX_XY_FACTORS = 4; // 2^4 = 16 expanded terms, worst case

// One (coefficient, operator-name) alternative in the expansion of a
// single Sx/Sy factor.
struct XYPiece { Cplx coef; char const* name; };

inline bool
site_has_ladder_ops(SiteSet const& sites, int site)
    {
    try { sites.op("S+",site); sites.op("S-",site); }
    catch(...) { return false; }
    return true;
    }

// Whether the Sy-elimination rewrite is both worthwhile and safe for this
// particular term list -- see the comment block above for each condition.
inline bool
should_realify(SiteSet const& sites, std::vector<MOTerm> const& terms)
    {
    if (!realify_spin_terms_enabled()) return false;
    bool any_sy = false;
    std::set<int> xy_sites;
    for (auto const& term : terms)
        {
        if (term.coef.imag() != 0.0) return false;
        int n_sy = 0, n_xy = 0;
        for (auto const& f : term.factors)
            {
            if (f.name=="Sy") { ++n_sy; ++n_xy; xy_sites.insert(f.site); }
            else if (f.name=="Sx") { ++n_xy; xy_sites.insert(f.site); }
            }
        if (n_sy % 2 != 0) return false;
        if (n_xy > MAX_XY_FACTORS) return false;
        if (n_sy > 0) any_sy = true;
        }
    if (!any_sy) return false;
    for (auto site : xy_sites)
        if (!site_has_ladder_ops(sites,site)) return false;
    return true;
    }

// Sx = (S+ + S-)/2, Sy = (S+ - S-)/(2i) = -i/2 S+ + i/2 S-.
inline std::vector<XYPiece>
expansion_of(std::string const& name)
    {
    if (name=="Sx")
        return {{Cplx(0.5,0.0),"S+"},{Cplx(0.5,0.0),"S-"}};
    return {{Cplx(0.0,-0.5),"S+"},{Cplx(0.0,0.5),"S-"}};
    }

} // namespace motermsdetail

AutoMPO inline
build_ampo(SiteSet const& sites, std::vector<MOTerm> const& terms)
    {
    auto ampo = AutoMPO(sites);
    bool realify = motermsdetail::should_realify(sites,terms);
    for (auto const& term : terms)
        {
        if (!realify)
            {
            HTerm hterm;
            hterm.coef = term.coef;
            for (auto const& f : term.factors) hterm.add(f.name,f.site);
            ampo.add(hterm);
            continue;
            }
        // Expand this term over every Sx/Sy factor's two ladder-operator
        // alternatives, in place, so the surviving factor order (which
        // matters for non-commuting same-site products) is untouched.
        std::vector<int> xy_pos; // indices into term.factors
        for (size_t k=0;k<term.factors.size();++k)
            if (term.factors[k].name=="Sx" || term.factors[k].name=="Sy")
                xy_pos.push_back(int(k));
        int nxy = int(xy_pos.size());
        int ncomb = 1 << nxy;
        for (int mask=0;mask<ncomb;++mask)
            {
            HTerm hterm;
            Cplx coef = term.coef;
            std::vector<std::string> names;
            names.reserve(term.factors.size());
            for (auto const& f : term.factors) names.push_back(f.name);
            for (int b=0;b<nxy;++b)
                {
                auto pieces = motermsdetail::expansion_of(term.factors[xy_pos[b]].name);
                auto const& p = pieces.at((mask>>b)&1);
                coef *= p.coef;
                names[xy_pos[b]] = p.name;
                }
            // Even Sy count guarantees this; drop the residual rounding
            // dust so svdMPO's own is_real test (t.coef.imag() != 0.0,
            // an exact comparison) actually sees a real coefficient.
            coef = Cplx(coef.real(),0.0);
            if (coef == Cplx(0.0,0.0)) continue;
            hterm.coef = coef;
            for (size_t k=0;k<term.factors.size();++k)
                hterm.add(names[k],term.factors[k].site);
            ampo.add(hterm);
            }
        }
    return ampo;
    }

MPO inline
build_mpo(SiteSet const& sites, std::vector<MOTerm> const& terms, int mpomaxm)
    {
    auto ampo = build_ampo(sites,terms);
    if (mpomaxm<5) mpomaxm = 5000; // default, mirrors get_ampo_operator.h
    return toMPO(ampo,{"MaxDim",mpomaxm,"Exact",false});
    }
