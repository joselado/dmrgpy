//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
// Port of mpscpp2/extra/Z4.h to the ITensor v3 API, following the same
// pattern as v3's own stock itensor/mps/sites/Z3.h (mod-3 -> mod-4).
//
#ifndef __ITENSOR_Z4_H
#define __ITENSOR_Z4_H
#include "itensor/mps/siteset.h"

namespace itensor {

class Z4Site;

using Z4 = BasicSiteSet<Z4Site>;

class Z4Site
    {
    Index s;
    public:

    Z4Site() { }

    Z4Site(Index I) : s(I) { }

    Z4Site(int n, Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Z4,n="+str(n));
        if(args.getBool("ConserveQNs",true))
            {
            s = Index(QN({"T",0,4}),1,
                      QN({"T",1,4}),1,
                      QN({"T",2,4}),1,
                      QN({"T",3,4}),1,Out,ts);
            }
        else
            {
            s = Index(4,ts);
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
        if(state == "2") { return s(3); }
        else
        if(state == "3") { return s(4); }
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
        auto Two = s(3);
        auto TwoP = sP(3);
        auto Three = s(4);
        auto ThreeP = sP(4);

        auto Op = ITensor(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,2);
            Op.set(Three,ThreeP,3);
            }
        else
        if(opname == "Sig")
            {
            Op.set(Zer,ThreeP,1);
            Op.set(One,ZerP,1);
            Op.set(Two,OneP,1);
            Op.set(Three,TwoP,1);
            }
        else
        if(opname == "SigDag")
            {
            Op.set(Three,ZerP,1);
            Op.set(Zer,OneP,1);
            Op.set(One,TwoP,1);
            Op.set(Two,ThreeP,1);
            }
        else
        if(opname == "Tau")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/4.)+sin(2.*Pi/4.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/4.)+sin(4.*Pi/4.)*1_i);
            Op.set(Three,ThreeP,cos(6.*Pi/4.)+sin(6.*Pi/4.)*1_i);
            }
        else
        if(opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/4.)-sin(2.*Pi/4.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/4.)-sin(4.*Pi/4.)*1_i);
            Op.set(Three,ThreeP,cos(6.*Pi/4.)-sin(6.*Pi/4.)*1_i);
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
        if(opname == "Proj2")
            {
            Op.set(Two,TwoP,1);
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
