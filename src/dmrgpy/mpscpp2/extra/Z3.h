//
// Distributed under the ITensor Library License, Version 1.2
// (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPINTWO_H
#define __ITENSOR_SPINTWO_H
#include "itensor/mps/siteset.h"
//#include "itensor/mps/sites/spinhalf.h"


namespace itensor {

class Z3 : public SiteSet
    {
    public:

    Z3() { }

    Z3(int N, 
        Args const& args = Args::global());

    void
        read(std::istream& s);

    };


class Z3Site
	{
    IQIndex s;
	public:

    Z3Site() { }

    Z3Site(IQIndex I) : s(I) { }

    Z3Site(int n, Args const& args = Args::global())
		{
        s = IQIndex{nameint("S=2 site=",n),
            Index(nameint("Up:site",n),1,Site),QN("Sz=",+4),
            Index(nameint("Upi:site",n),1,Site),QN("Sz=",+2),
            Index(nameint("Z0:site",n),1,Site),QN("Sz=",0),
            Index(nameint("Dni:site",n),1,Site),QN("Sz=",-2),
            Index(nameint("Dn:site",n),1,Site),QN("Sz=",-4)};
		}

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
		{
        if (state == "Up" || state == "4") 
            {
            return s(1);
            }
        else if (state == "Upi" || state == "2")
            {
            return s(2);
            }
        else if (state == "Z0" || state == "0")
            {
            return s(3);
            }
        else if (state == "Dni" || state == "-2")
            {
            return s(4);
            }
        else if (state == "Dn" || state == "-4")
            {
            return s(5);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
		}

    IQTensor
    op(std::string const& opname,
       Args const& args) const
        auto sP = prime(s);
        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);
        auto Two = s(3);
        auto TwoP = sP(3);
        auto Op = ITensor(dag(s),sP);
        if(opname == "N")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,2);
            }
        else
        if(opname == "Sig")
            {
            Op.set(Zer,TwoP,1);
            Op.set(One,ZerP,1);
            Op.set(Two,OneP,1);
            }
        else
        if(opname == "SigDag")
            {
            Op.set(Two,ZerP,1);
            Op.set(Zer,OneP,1);
            Op.set(One,TwoP,1);
            }
        else
        if(opname == "Tau")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/3.)+sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)+sin(4.*Pi/3.)*1_i);
            }
        else
        if(opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/3.)-sin(2.*Pi/3.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/3.)-sin(4.*Pi/3.)*1_i);
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
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        else
            {
            Error("Operator " + opname + " name not recognized");
            }

        return Op;
		}
	}; //Z3Site

inline Z3::
Z3(int N, 
        Args const& args)
	{
    auto sites = SiteStore(N);
    auto start = 1;

    SiteSet::init(std::move(sites));
	}

void inline Z3::
read(std::istream& s)
	{
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        for(int j = 1; j <= N; ++j) 
            {
            auto I = IQIndex{};
            I.read(s);
            store.set(j,Z3Site(I));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
