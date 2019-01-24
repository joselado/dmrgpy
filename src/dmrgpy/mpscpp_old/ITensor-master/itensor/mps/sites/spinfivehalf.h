
//
// Distributed under the ITensor Library License, Version 1.2
// (See accompanying LICENSE file.)
//


#ifndef __ITENSOR_SPINFIVEHALF_H
#define __ITENSOR_SPINFIVEHALF_H
#include "itensor/mps/siteset.h"


namespace itensor {

class SpinFiveHalf : public SiteSet
    {
    public:

    SpinFiveHalf() { }

    SpinFiveHalf(int N, 
        Args const& args = Args::global());

    void
        read(std::istream& s);

    };


class SpinFiveHalfSite
	{
    IQIndex s;
	public:

    SpinFiveHalfSite() { }

    SpinFiveHalfSite(IQIndex I) : s(I) { }

    SpinFiveHalfSite(int n, Args const& args = Args::global())
		{
        s = IQIndex{nameint("S=2 site=",n),
            Index(nameint("Up:site",n),1,Site),QN("Sz=",+5),
            Index(nameint("Upi:site",n),1,Site),QN("Sz=",+3),
            Index(nameint("Upii:site",n),1,Site),QN("Sz=",+1),
            Index(nameint("Dnii:site",n),1,Site),QN("Sz=",-1),
            Index(nameint("Dni:site",n),1,Site),QN("Sz=",-3),
            Index(nameint("Dn:site",n),1,Site),QN("Sz=",-5)};
		}

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
		{
        if (state == "Up" || state == "5") 
            {
            return s(1);
            }
        else if (state == "Upi" || state == "3")
            {
            return s(2);
            }
        else if (state == "Upii" || state == "1")
            {
            return s(3);
            }
        else if (state == "Dnii" || state == "-1")
            {
            return s(4);
            }
        else if (state == "Dni" || state == "-3")
            {
            return s(5);
            }
        else if (state == "Dn" || state == "-5")
            {
            return s(6);
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
		{
        const Real val1 = std::sqrt(5.0)/2.0;
        const Real val2 = std::sqrt(8.0)/2.0;
        const Real val3 = std::sqrt(9.0)/2.0;

        auto sP = prime(s);

        auto Up  = s(1);
        auto UpP = sP(1);
        auto Upi = s(2);
        auto UpiP = sP(2);
        auto Upii = s(3);
        auto UpiiP = sP(3);
        auto Dnii = s(4);
        auto DniiP = sP(4);
        auto Dni = s(5);
        auto DniP = sP(5);
        auto Dn  = s(6);
        auto DnP = sP(6);

        auto Op = IQTensor(dag(s),sP);

        if (opname == "Sz")
            {
            Op.set(Up,UpP,+2.5);
            Op.set(Upi,UpiP,+1.5);
            Op.set(Upii,UpiiP,+0.5);
            Op.set(Dnii,DniiP,-0.5);
            Op.set(Dni,DniP,-1.5);
            Op.set(Dn,DnP,-2.5);
            }
        else if (opname == "Sx")
            {
            //mixedIQTensor call needed here
            //because as an IQTensor, Op would
            //not have a well defined QN flux
            Op = mixedIQTensor(s,sP);
            Op.set(Up,UpiP,val1);
            Op.set(Upi,UpP,val1);
            Op.set(Upi,UpiiP,val2);
            Op.set(Upii,UpiP,val2);
            Op.set(Upii,DniiP,val3); 
            Op.set(Dnii,UpiiP,val3);
            Op.set(Dnii,DniP,val2);
            Op.set(Dni,DniiP,val2);
            Op.set(Dni,DnP,val1);
            Op.set(Dn,DniP,val1);
            }
        else if (opname == "Sy")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(Up,UpiP,val1*Cplx_i);
            Op.set(Upi,UpP,-val1*Cplx_i);
            Op.set(Upi,UpiiP,val2*Cplx_i);
            Op.set(Upii,UpiP,-val2*Cplx_i);
            Op.set(Upii,DniiP,val3*Cplx_i);
            Op.set(Dnii,UpiiP,-val3*Cplx_i);
            Op.set(Dnii,DniP,val2*Cplx_i);
            Op.set(Dni,DniiP,-val2*Cplx_i);
            Op.set(Dni,DnP,val1*Cplx_i);
            Op.set(Dn,DniP,-val1*Cplx_i);
            }
        else
            {
            Error("Operator " + opname + " name not recognized");
            }

        return Op;
		}
	}; //SpinFiveHalfSite

inline SpinFiveHalf::
SpinFiveHalf(int N, 
        Args const& args)
	{
    auto sites = SiteStore(N);

    auto start = 1;

    for(int j = start; j < N; ++j)
        {
        sites.set(j,SpinFiveHalfSite(j));
        }

    SiteSet::init(std::move(sites));
	}

void inline SpinFiveHalf::
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
            if(I.m() == 3) store.set(j,SpinFiveHalfSite(I));
            else if(I.m() == 2) store.set(j,SpinHalfSite(I));
            else Error(format("SpinFiveHalf cannot read index of size %d",I.m()));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
