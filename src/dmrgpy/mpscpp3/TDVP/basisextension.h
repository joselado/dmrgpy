#include "itensor/decomp.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/mpsalgs.cc"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/cputime.h"
#include "itensor/tensor/slicemat.h"

namespace itensor{

void
denmatSumDecomp(std::vector<MPS> const& psis,
                MPS & res,
                std::vector<ITensor> & Bs,
                int b,
                Direction dir,
                Args args = Args::global())
	{
	// NumCenter can only be 1 if not want to treat res exactly
	const int numCenter = args.getInt("NumCenter",1);
	const bool quiet = args.getBool("Quiet",false);

	// SVD site tensor of res without truncaion
	auto [V1,S1,U1] = svd(Bs.front(), dir == Fromleft? rightLinkIndex(res,b): leftLinkIndex(res,b));

	// Find the indices to be left to the density matrix
	auto& to_orth = Bs.front();
	auto& newoc   = (dir==Fromleft? res(b+1) : res(b-1));
	auto& activeInds = to_orth.inds();
	auto cinds = stdx::reserve_vector<Index>(activeInds.r());
	for(auto& I : activeInds)
   		{
   		if(!hasIndex(newoc,I)) cinds.push_back(I);
   		}

	auto [cmb,mid] = combiner(std::move(cinds));

	if(dim(mid) <= dim(commonIndex(U1,S1)))
		{
		res.ref(b) = U1;
		if(!quiet)
			printfln("warning: at bond %d, already reach maximum bond dimension.", b);
		}
	else
		{
		// Density matrix summation to be truncated
		ITensor rho2;
		if(numCenter == 1)
			{
			auto psi = psis.begin();
			for(auto B = Bs.begin()+1; B != Bs.end(); ++B, ++psi)
				{
				rho2 += prime(*B)*dag(prime(*B,dir == Fromleft? rightLinkIndex(*psi,b): leftLinkIndex(*psi,b)));
				}
			rho2.swapPrime(0,1);
			}
		else
			{
			Error("numCenter can only be one");
			}

		// Form the density matrix with only 2 indices
		U1 *= cmb;
		rho2 *= cmb;
		cmb.prime();
		cmb.dag();
		rho2 *= cmb;
		cmb.noPrime();

		// Project rho2c to the orthogonal complement of U1
		auto proj2 = toDense(delta(dag(mid),prime(mid)))-dag(U1)*prime(U1,mid); 
		auto normrho2 = norm(rho2);
		rho2 *= mapPrime(proj2,0,2);
		rho2 *= proj2;
		rho2.mapPrime(2,0);
		rho2.swapPrime(0,1);
		auto normPrho2P = norm(rho2);
		if(normPrho2P/normrho2 < 1E-12)// TODO: changed to calculate the trace will have less complexity and have the same effect!
			{
			res.ref(b) = cmb * U1;
			if(!quiet)
				printfln("warning: at bond %d, not adding any new basis.", b);
			}
		else
			{
			// Diagonalize rho2c to obtain U2
			ITensor U2, D2;
			args.add("Truncate",true);
			diag_hermitian(rho2,U2,D2,args);// T==prime(U)*D*dag(U)
			U2.dag();

			// Direct sum of U1 and U2
			auto i1 = commonIndex(U1,S1);
			auto i2 = commonIndex(U2,D2);
			auto sumind = Index(dim(i1)+dim(i2),"Link");
			sumind.setDir(i1.dir());
			ITensor expand1,expand2;
			plussers(i1,i2,sumind,expand1,expand2);
			auto U = U1*expand1 + U2*expand2;

			res.ref(b) = cmb * U;
			}

		}

	// Obtain the new Bs for the operation of the next site
        Bs.front() *= dag(res(b));
        Bs.front() *= (dir == Fromleft? res(b+1): res(b-1));
        auto psi = psis.begin();
	for(auto B = Bs.begin()+1; B != Bs.end(); ++B, ++psi)
		{
		(*B) *= dag(res(b));
		(*B) *= (dir == Fromleft? (*psi).A(b+1): (*psi).A(b-1));
		}
	
	}

void
addBasisWorker(std::vector<MPS> const& psis,
               MPS & res,
               Direction dir,
               const Args & args = Args::global())
	{
	int N = length(res);
	int nt = psis.size()+1;

	if(dir == Fromleft)
		{
		if(orthoCenter(res) != 1)
			Error("OC need set to be 1");
		for(auto& psi : psis)
			{
			if(orthoCenter(psi) != 1)
				Error("OC need set to be 1");
			}

		auto Bs = std::vector<ITensor>(nt);
		Bs.front() = res(1);
		auto psi = psis.begin();
		for(auto B = Bs.begin()+1; B != Bs.end(); ++B, ++psi)
			{
			(*B) = (*psi).A(1);
			}

		for(int b = 1; b < N ; ++b)
			{
			denmatSumDecomp(psis,res,Bs,b,Fromleft,args);
			}

		res.Aref(N) = Bs.front();
		}
	else
		{
		if(orthoCenter(res) != N)
			Error("OC need set to be N");
		for(auto& psi : psis)
			{
			if(orthoCenter(psi) != N)
				Error("OC need set to be N");
			}

		auto Bs = std::vector<ITensor>(nt);
		Bs.front() = res(N);
		auto psi = psis.begin();
		for(auto B = Bs.begin()+1; B != Bs.end(); ++B, ++psi)
			{
			(*B) = (*psi).A(N);
			}

		for(int b = N; b > 1 ; --b)
			{
			denmatSumDecomp(psis,res,Bs,b,Fromright,args);
			}

		res.Aref(1) = Bs.front();
		}
	}

void addBasis(MPS& phi,
              const MPO& H,
              std::vector<Real> const& truncK,
              const Args& args0 = Args::global())
	{
	auto quiet = args0.getBool("Quiet",false);
	auto dk = args0.getInt("KrylovOrd",2);
	auto method = args0.getString("Method","DensityMatrix");
	auto nsw = args0.getInt("Nsweep",2);
	auto donormalize = args0.getBool("DoNormalize",false);//TODO: add a function to construct general 1-tauH
	
	auto psis = std::vector<MPS>(dk-1);	

	cpu_time expand_time;
	for(int i = 0; i < dk-1; ++i)
		{
		auto args1 = Args("Method=",method,"Cutoff=",truncK.at(i),"Nsweep=",nsw);
		if(args0.defined("WriteDim"))
			{
			args1.add("WriteDim",args0.getInt("WriteDim"));
			if(args0.defined("WriteDir"))  args1.add("WriteDir",args0.getString("WriteDir"));
			}
		
		if(i==0)
			psis.at(i) = applyMPO(H,phi,args1);
		else
			psis.at(i) = applyMPO(H,psis.at(i-1),args1);
		
		psis.at(i).noPrime();
		if(donormalize)
			psis.at(i).normalize();
		
		if(!quiet)
			{
			printfln("norm(psi%d)=%.20f",i+1,norm(psis.at(i)));
			printfln("maxLinkDim(psi%d) = %d",i+1,maxLinkDim(psis.at(i)));
			}
		}

	int N = length(phi);
	for(int i = 0; i < dk-1; ++i)
		{
		psis.at(i).position(N);
		}
	phi.position(N);

	//TODO: adjustable weight for each psi
	addBasisWorker(psis,phi,Fromright,args0);

	auto sm = expand_time.sincemark();
	printfln("\nmaxLinkDim after global subspace expansion = %d",maxLinkDim(phi));
	printfln("Global subspace expansion: cputime = %s, walltime = %s",showtime(sm.time),showtime(sm.wall));
	}

void addBasis(MPS& phi,
              const MPO& H,
              std::vector<int> const& maxdimK,
              const Args& args0 = Args::global())
	{
        auto dk = args0.getInt("KrylovOrd",2);
        auto method = args0.getString("Method","DensityMatrix");
        auto nsw = args0.getInt("Nsweep",2);
        auto donormalize = args0.getBool("DoNormalize",false);

	auto psis = std::vector<MPS>(dk-1);	
        
        cpu_time expand_time;
	for(int i = 0; i < dk-1; ++i)
		{
		auto args1 = Args("Method=",method,"MaxDim=",maxdimK.at(i),"Nsweep=",nsw);
		if(args0.defined("WriteDim"))
			{
			args1.add("WriteDim",args0.getInt("WriteDim"));
			if(args0.defined("WriteDir"))  args1.add("WriteDir",args0.getString("WriteDir"));
			}
		
		if(i==0)
			psis.at(i) = applyMPO(H,phi,args1);
		else
			psis.at(i) = applyMPO(H,psis.at(i-1),args1);
		
		psis.at(i).noPrime();
		if(donormalize)
			psis.at(i).normalize();

                printfln("norm(psi%d)=%.20f",i+1,norm(psis.at(i)));
                printfln("maxLinkDim(psi%d) = %d",i+1,maxLinkDim(psis.at(i)));
		}

	int N = length(phi);
	for(int i = 0; i < dk-1; ++i)
		{
		psis.at(i).position(N);
		}
	phi.position(N);

	addBasisWorker(psis,phi,Fromright,args0);
        
        auto sm = expand_time.sincemark();
        printfln("\nmaxLinkDim after global subspace expansion = %d",maxLinkDim(phi));
        printfln("Global subspace expansion: cputime = %s, walltime = %s",showtime(sm.time),showtime(sm.wall));
	}

}// namespace itensor
