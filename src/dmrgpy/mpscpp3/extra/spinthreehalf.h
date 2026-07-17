
//
// Distributed under the ITensor Library License, Version 1.2
// (See accompanying LICENSE file.)
//
// Port of mpscpp2/extra/spinthreehalf.h to the ITensor v3 API: IQIndex ->
// Index (QN blocks carried directly on the Index), IQTensor -> ITensor,
// IndexType Site -> TagSet("Site"), mixedIQTensor(s,sP) -> plain
// ITensor(dag(s),sP) (v3 does not need the "mixed" workaround for
// non-flux-conserving operators like Sx/Sy). Same operator names/values as
// the v2 version.
//

#ifndef __ITENSOR_SPINTHREEHALF_H
#define __ITENSOR_SPINTHREEHALF_H
#include "itensor/mps/siteset.h"

namespace itensor {

class SpinThreeHalf : public SiteSet
    {
    public:

    SpinThreeHalf() { }

    SpinThreeHalf(int N,
        Args const& args = Args::global());

    void
        read(std::istream& s);

    };


class SpinThreeHalfSite
	{
    Index s;
	public:

    SpinThreeHalfSite() { }

    SpinThreeHalfSite(Index I) : s(I) { }

    SpinThreeHalfSite(int n, Args const& args = Args::global())
		{
        auto ts = TagSet("Site,S=3/2,n="+str(n));
        if(args.getBool("ConserveQNs",true))
            {
            s = Index(QN({"Sz",+3}),1,
                      QN({"Sz",+1}),1,
                      QN({"Sz",-1}),1,
                      QN({"Sz",-3}),1,Out,ts);
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
        if (state == "Up" || state == "3")
            {
            return s(1);
            }
        else if (state == "Upi" || state == "1")
            {
            return s(2);
            }
        else if (state == "Dni" || state == "-1")
            {
            return s(3);
            }
        else if (state == "Dn" || state == "-3")
            {
            return s(4);
            }
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
        const Real val1 = std::sqrt(3.0)/2.0;

        auto sP = prime(s);

        auto Up  = s(1);
        auto UpP = sP(1);
        auto Upi = s(2);
        auto UpiP = sP(2);
        auto Dni = s(3);
        auto DniP = sP(3);
        auto Dn  = s(4);
        auto DnP = sP(4);

        auto Op = ITensor(dag(s),sP);

        if (opname == "Sz")
            {
            Op.set(Up,UpP,+1.5);
            Op.set(Upi,UpiP,+0.5);
            Op.set(Dni,DniP,-0.5);
            Op.set(Dn,DnP,-1.5);
            }
        else if (opname == "Sx")
            {
            Op.set(Up,UpiP,val1);
            Op.set(Upi,UpP,val1);
            Op.set(Upi,DniP,1.0);
            Op.set(Dni,UpiP,1.0);
            Op.set(Dni,DnP,val1);
            Op.set(Dn,DniP,val1);
            }
        else if (opname == "Sy")
            {
            Op.set(Up,UpiP,val1*Cplx_i);
            Op.set(Upi,UpP,-val1*Cplx_i);
            Op.set(Upi,DniP,1.0*Cplx_i);
            Op.set(Dni,UpiP,-1.0*Cplx_i);
            Op.set(Dni,DnP,val1*Cplx_i);
            Op.set(Dn,DniP,-val1*Cplx_i);
            }
        else
            {
            throw ITError("Operator " + opname + " name not recognized");
            }

        return Op;
		}
	}; //SpinThreeHalfSite

inline SpinThreeHalf::
SpinThreeHalf(int N,
        Args const& args)
	{
    auto sites = SiteStore(N);

    auto start = 1;

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinThreeHalfSite(j));
        }

    SiteSet::init(std::move(sites));
	}

void inline SpinThreeHalf::
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
            if(dim(I) == 4) store.set(j,SpinThreeHalfSite(I));
            else if(dim(I) == 2) store.set(j,SpinHalfSite(I));
            else throw ITError(tinyformat::format("SpinThreeHalf cannot read index of size %d",dim(I)));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
