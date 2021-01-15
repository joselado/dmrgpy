
#ifndef __ITENSOR_BOSONFOUR_H
#define __ITENSOR_BOSONFOUR_H
#include "itensor/mps/siteset.h"


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
    IQIndex s;
	public:

    BosonFourSite() { }

    BosonFourSite(IQIndex I) : s(I) { }

    BosonFourSite(int n, Args const& args = Args::global())
		{
        s = IQIndex{nameint("Boson=4 site=",n),
            Index(nameint("Up:site",n),1,Site),QN("Sz=",+3),
            Index(nameint("Upi:site",n),1,Site),QN("Sz=",+1),
            Index(nameint("Dni:site",n),1,Site),QN("Sz=",-1),
            Index(nameint("Dn:site",n),1,Site),QN("Sz=",-3)};
		}

    IQIndex
    index() const { return s; }

    IQIndexVal
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
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
		}

    IQTensor
    op(std::string const& opname,
       Args const& args) const
		{

	const Real val1 = std::sqrt(1.0);
	const Real val2 = std::sqrt(2.0);
	const Real val3 = std::sqrt(3.0);
	const Real val4 = std::sqrt(4.0);
        auto sP = prime(s);

        auto A0  = s(1);
        auto A0P = sP(1);
        auto A1  = s(2);
        auto A1P = sP(2);
        auto A2  = s(3);
        auto A2P = sP(3);
        auto A3  = s(4);
        auto A3P = sP(4);
        auto Op = IQTensor(dag(s),sP);

        if (opname == "N")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A0,A0P,+0.0);
            Op.set(A1,A1P,+1.0);
            Op.set(A2,A2P,+2.0);
            Op.set(A3,A3P,+3.0);
            }
        else if (opname == "N0")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A0,A0P,+1.0);
            }
        else if (opname == "N1")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A1,A1P,+1.0);
            }
        else if (opname == "N2")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A2,A2P,+1.0);
            }
        else if (opname == "N3")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A3,A3P,+1.0);
            }
        else if (opname == "Adag")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A0,A1P,val1);
            Op.set(A1,A2P,val2);
            Op.set(A2,A3P,val3);
            }
        else if (opname == "A")
            {
            Op = mixedIQTensor(s,sP);
            Op.set(A1,A0P,val1);
            Op.set(A2,A1P,val2);
            Op.set(A3,A2P,val3);
            }
        else
            {
            Error("Operator " + opname + " name not recognized");
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
        sites.set(j,BosonFourSite(j));
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
            auto I = IQIndex{};
            I.read(s);
            store.set(j,BosonFourSite(I));
            }
        init(std::move(store));
        }
	}

} //namespace itensor

#endif
