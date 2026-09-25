// In-memory replacement for the "ampotk" file format (get_ampo_operator.h +
// the generated extra/ampotk/ampotk.h switch). A MultiOperator on the
// Python side is already exactly a list of (coefficient, [(opname,site),...])
// terms (see multioperator.py's MultiOperator.op) -- this is the same shape,
// just handed over as in-memory objects instead of being serialized to text
// and reparsed. Builds directly on ITensor's own AutoMPO::add(HTerm)/
// HTerm::add(name,site) API, which the old ampotk switch was already calling
// indirectly through the "ampo += coef,name,site,..." accumulator syntax --
// no 99-way-switch code generation is needed at all.

#include <algorithm>
#include <cmath>

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

AutoMPO inline
build_ampo(SiteSet const& sites, std::vector<MOTerm> const& terms)
    {
    auto ampo = AutoMPO(sites);
    for (auto const& term : terms)
        {
        HTerm hterm;
        hterm.coef = term.coef;
        for (auto const& f : term.factors) hterm.add(f.name,f.site);
        ampo.add(hterm);
        }
    return ampo;
    }

// -- unit scale for ITensor's absolute thresholds ------------------------
//
// The same helpers as mpscpp3/mo_terms.h's, see the comment there for the
// full account. In short: svdMPO's truncate() compares the discarded
// squared singular values against an absolute 1E-13 (and skips matrix
// elements below 1E-14), and davidson() randomizes a Krylov direction
// whose residual is below an absolute 1E-10, both calibrated for an
// operator of order one. On this backend a 6-site Heisenberg bond carries
// XX, YY and ZZ as three channels of weight J^2 each, so one goes below
// J=3.16e-7 and two below 2.24e-7, the third being kept by the minimum
// bond dimension of 1 (gs_energy()/J of -1.747 at 3e-7 and the Neel -1.25 at 1e-7, against
// -2.4936; 2026-09-25 audit, small-units). An operator whose largest
// coefficient is below 1 is therefore handed to ITensor multiplied by the
// power of two that brings it into [1,2), and divided back afterwards,
// exactly; at a largest coefficient of 1 or more nothing changes, byte for
// byte. One scale per operator does not reach a bond whose strongest
// crossing term is far below the operator's largest coefficient (an O(1)
// offset or field next to exchange below about 3e-7 of it, a weak link),
// since truncate() runs bond by bond; see the mpscpp3 comment.
inline double
unit_scale_up(double cmax)
    {
    if (!(cmax > 0.0) || cmax >= 1.0) return 1.0;
    int e = 0;
    std::frexp(cmax,&e); // cmax in [2^(e-1),2^e), and e <= 0
    return std::ldexp(1.0,std::min(1-e,1000)); // cmax*up in [1,2)
    }

// toMPO<ITensor>() with the operator at unit scale. AutoMPO keeps its terms
// ordered by operator string alone (LessNoCoef), so the scaled copy holds
// the same terms in the same order. The default Args are toMPO's own, so
// the MPO(ampo) conversions this replaces are unchanged at up==1.
MPO inline
to_mpo_unit(AutoMPO const& ampo, Args const& args = Args::global())
    {
    double cmax = 0.0;
    for (auto const& t : ampo.terms()) cmax = std::max(cmax,std::abs(t.coef));
    double up = unit_scale_up(cmax);
    if (up==1.0) return toMPO<ITensor>(ampo,args);
    auto scaled = AutoMPO(ampo.sites());
    for (auto t : ampo.terms()) { t.coef *= up; scaled.add(t); }
    auto W = toMPO<ITensor>(scaled,args);
    W *= 1.0/up; // exact: multiplies one tensor's elements by 2^-k
    return W;
    }

// The largest |coefficient| of a term list, which is what the session reads
// its Hamiltonian's unit scale from (Chain::hscale_up_).
inline double
max_abs_coef(std::vector<MOTerm> const& terms)
    {
    double cmax = 0.0;
    for (auto const& t : terms) cmax = std::max(cmax,std::abs(t.coef));
    return cmax;
    }

MPO inline
build_mpo(SiteSet const& sites, std::vector<MOTerm> const& terms, int mpomaxm)
    {
    auto ampo = build_ampo(sites,terms);
    if (mpomaxm<5) mpomaxm = 5000; // default, mirrors get_ampo_operator.h
    return to_mpo_unit(ampo,{"Maxm",mpomaxm,"Exact",false});
    }
