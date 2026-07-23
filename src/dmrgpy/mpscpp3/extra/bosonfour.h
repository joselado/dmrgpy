
#ifndef __ITENSOR_BOSONFOUR_H
#define __ITENSOR_BOSONFOUR_H
#include "itensor/mps/siteset.h"
#include <algorithm>
#include <cctype>
#include <string>

// Port of mpscpp2/extra/bosonfour.h to the ITensor v3 API, see
// spinthreehalf.h in this directory for the general porting notes.
//
// Generalized (v3 only, not ported back to mpscpp2) from a fixed 4-level
// site to an arbitrary occupation cutoff, driven by the "MaxOcc" Args key
// (dimension = MaxOcc+1, defaulting to 3 i.e. the original 4-level site)
// -- see get_sites.h's boson type-code range (100+dim, 102..199) for how a
// dmrgpy-level Bosonic_Chain(maxnb=...) selects MaxOcc. Operator names N,
// A, Adag, and the occupation-projectors N0..N{MaxOcc} generalize the
// original hardcoded N0..N3 the same way pyboson/boson.py's ED BosonChain
// already did (that side was always dimension-generic, see its own
// per-site "maxnb" loop) -- this brings the DMRG side to the same
// generality instead of silently truncating to dim 4 regardless of the
// Python-level maxnb (a real, previously-silent DMRG/ED mismatch bug).

namespace itensor {

class BosonFour : public SiteSet
    {
    public:

    BosonFour() { }

    BosonFour(int N,
        Args const& args = Args::global());

    void
        read(std::istream& s);

    };


class BosonFourSite
	{
    Index s;
	public:

    BosonFourSite() { }

    BosonFourSite(Index I) : s(I) { }

    BosonFourSite(int n, Args const& args = Args::global())
		{
        auto maxOcc = args.getInt("MaxOcc",3); // dimension = maxOcc+1, default 4-level site
        auto dim = maxOcc+1;
        auto ts = TagSet("Site,Boson="+str(dim)+",n="+str(n));
        if(args.getBool("ConserveQNs",true))
            {
            auto qints = Index::qnstorage(dim);
            for(int occ : range(dim)) qints[occ] = QNInt(QN({"Nb",occ}),1);
            s = Index(std::move(qints),Out,ts);
            }
        else
            {
            s = Index(dim,ts);
            }
		}

    Index
    index() const { return s; }

    // Plain occupation-number strings ("0".."MaxOcc") only. The original
    // fixed-dim=4 site additionally accepted "Up"/"Upi"/"Dni"/"Dn" and
    // read its numeric strings ("3","1","-1","-3") as *Sz labels* from
    // its old fixed QN scheme (matching the +3,+1,-1,-3 Sz values that
    // scheme's QN block used, i.e. inverted relative to occupation
    // number: old "3" meant occupation 0, not 3) -- not carried forward
    // here since MaxOcc is now arbitrary and that convention doesn't
    // generalize. No caller in this codebase currently invokes
    // state() for a boson site (dmrg()'s default_mps() always starts
    // from randomMPS, never InitState), so this is a behavior change
    // with no live effect today, not a regression -- flagged here for
    // whoever adds the first real caller.
    IndexVal
    state(std::string const& state)
		{
        auto maxOcc = dim(s)-1;
        for(auto occ : range(1+maxOcc))
            {
            if(state == str(occ)) return s(1+occ);
            }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
		}

    ITensor
    op(std::string const& opname,
       Args const& args) const
		{
        auto sP = prime(s);
        auto maxOcc = dim(s)-1;

        auto Op = ITensor(dag(s),sP);

        if (opname == "N")
            {
            for(auto occ : range1(1+maxOcc)) Op.set(s=occ,sP=occ,occ-1);
            }
        else if (opname == "Adag")
            {
            for(auto occ : range1(maxOcc)) Op.set(s=occ,sP=occ+1,std::sqrt((Real)occ));
            }
        else if (opname == "A")
            {
            for(auto occ : range1(maxOcc)) Op.set(s=occ+1,sP=occ,std::sqrt((Real)occ));
            }
        else if (opname.size()>1 && opname[0]=='N' &&
                 std::all_of(opname.begin()+1,opname.end(),::isdigit))
            {
            // occupation-number projector N<k>, k=0..maxOcc (generalizes
            // the original hardcoded N0..N3)
            auto level = std::stoi(opname.substr(1));
            if(level<0 || level>maxOcc)
                throw ITError("Operator " + opname + " out of range for this site's MaxOcc");
            Op.set(s=level+1,sP=level+1,+1.0);
            }
        else
            {
            throw ITError("Operator " + opname + " name not recognized");
            }

        return Op;
		}
	}; //BosonFourSite

inline BosonFour::
BosonFour(int N,
        Args const& args)
	{
    auto sites = SiteStore(N);

    auto start = 1;

    for(int j = start; j < N; ++j)
        {
        sites.set(j,BosonFourSite(j,args));
        }

    SiteSet::init(std::move(sites));
	}

void inline BosonFour::
read(std::istream& s)
	{
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            auto I = Index{};
            I.read(s);
            store.set(j,BosonFourSite(I));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
