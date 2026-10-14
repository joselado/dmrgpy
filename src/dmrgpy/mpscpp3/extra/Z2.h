//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
// Z2 clock site (the Z=2 parafermion of parafermionchain.py, i.e. the
// transverse-field Ising model in clock language), the mod-2 member of the
// Z3.h/Z4.h family with the same operator names and the same matrix
// convention: Sig is the shift (Sig|1>=|0>, Sig|0>=|1>), Tau the clock
// diag(1,-1), N diag(0,1). Site-type code -1 (get_sites.h).
//
#ifndef __ITENSOR_Z2_H
#define __ITENSOR_Z2_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Z2Site;

using Z2 = BasicSiteSet<Z2Site>;

class Z2Site
    {
    Index s;
    public:

    Z2Site() { }

    Z2Site(Index I) : s(I) { }

    Z2Site(int n, Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Z2,n="+str(n));
        if(args.getBool("ConserveQNs",true))
            {
            s = Index(QN({"T",0,2}),1,
                      QN({"T",1,2}),1,Out,ts);
            }
        else
            {
            s = Index(2,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "0") { return s(1); }
        else
        if(state == "1") { return s(2); }
        else
            {
            throw ITError("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);

        auto Op = ITensor(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            }
        else
        if(opname == "Sig" || opname == "SigDag")
            {
            Op.set(Zer,OneP,1);
            Op.set(One,ZerP,1);
            }
        else
        if(opname == "Tau" || opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,-1);
            }
        else
        if(opname == "Proj0")
            {
            Op.set(Zer,ZerP,1);
            }
        else
        if(opname == "Proj1")
            {
            Op.set(One,OneP,1);
            }
        else
            {
            throw ITError("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

#endif
