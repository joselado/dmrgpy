//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
// Z2 clock site (the Z=2 parafermion of parafermionchain.py), the mod-2
// member of the Z3/Z4 family with the same operator names and matrix
// convention: Sig the shift, Tau the clock diag(1,-1), N diag(0,1).
// Site-type code -1 (get_sites.h).
//
#ifndef __ITENSOR_Z2_H
#define __ITENSOR_Z2_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Z2Site;

using Z2 = BasicSiteSet<Z2Site>;

class Z2Site
    {
    IQIndex s;
    public:

    Z2Site() { }

    Z2Site(IQIndex I) : s(I) { }

    Z2Site(int n, Args const& args = Args::global())
        {
        s = IQIndex{nameint("Z2 site=",n),
        Index(nameint("0|site",n),1,Site),QN({0,2}),
        Index(nameint("1|site",n),1,Site),QN({1,2})};
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0") { return s(1); }
        else
        if(state == "1") { return s(2); }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);

        auto Op = IQTensor(dag(s),sP);

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
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };

} //namespace itensor

#endif
