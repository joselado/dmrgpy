// In-memory replacement for the "ampotk" file format (get_ampo_operator.h +
// the generated extra/ampotk/ampotk.h switch). A MultiOperator on the
// Python side is already exactly a list of (coefficient, [(opname,site),...])
// terms (see multioperator.py's MultiOperator.op) -- this is the same shape,
// just handed over as in-memory objects instead of being serialized to text
// and reparsed. Builds directly on ITensor's own AutoMPO::add(HTerm)/
// HTerm::add(name,site) API, which the old ampotk switch was already calling
// indirectly through the "ampo += coef,name,site,..." accumulator syntax --
// no 99-way-switch code generation is needed at all.

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

MPO inline
build_mpo(SiteSet const& sites, std::vector<MOTerm> const& terms, int mpomaxm)
    {
    auto ampo = build_ampo(sites,terms);
    if (mpomaxm<5) mpomaxm = 5000; // default, mirrors get_ampo_operator.h
    return toMPO<ITensor>(ampo,{"Maxm",mpomaxm,"Exact",false});
    }
