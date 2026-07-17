#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"


namespace itensor {

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
           Cplx t,
           const Sweeps& sweeps,
           const Args& args = Args::global());

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
           Cplx t,
           const Sweeps& sweeps,
           DMRGObserver & obs,
           Args args = Args::global());

//
// Available TDVP methods:
// second order integrator: sweep left-to-right and right-to-left
//

//
//TDVP with an MPO
//
Real inline
tdvp(MPS & psi, 
     MPO const& H,
     Cplx t, 
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with an MPO and custom DMRGObserver
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
     Cplx t,
     const Sweeps& sweeps, 
     DMRGObserver & obs,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

//
//TDVP with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
     Cplx t,
     ITensor const& LH, 
     ITensor const& RH,
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with an MPO and boundary tensors LH, RH
//and a custom observer
//
Real inline
tdvp(MPS & psi, 
     MPO const& H, 
     Cplx t,
     ITensor const& LH, 
     ITensor const& RH,
     const Sweeps& sweeps, 
     DMRGObserver& obs,
     const Args& args = Args::global())
    {
    LocalMPO PH(H,LH,RH,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

//
//TDVP with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
Real inline
tdvp(MPS& psi, 
     std::vector<MPO> const& Hset,
     Cplx t, 
     const Sweeps& sweeps,
     const Args& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,args);
    return energy;
    }

//
//TDVP with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
Real inline
tdvp(MPS & psi, 
     std::vector<MPO> const& Hset, 
     Cplx t,
     const Sweeps& sweeps, 
     DMRGObserver& obs,
     const Args& args = Args::global())
    {
    LocalMPOSet PH(Hset,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }


//
// TDVPWorker
//

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& PH,
           Cplx t,
           Sweeps const& sweeps,
           Args const& args)
    {
    DMRGObserver obs(psi,args);
    Real energy = TDVPWorker(psi,PH,t,sweeps,obs,args);
    return energy;
    }

template <class LocalOpT>
Real
TDVPWorker(MPS & psi,
           LocalOpT& H,
           Cplx t,
           Sweeps const& sweeps,
           DMRGObserver& obs,
           Args args)
    { 
    // Truncate blocks of degenerate singular values (or not)
    args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

    const bool silent = args.getBool("Silent",false);
    if(silent)
        {
        args.add("Quiet",true);
        args.add("PrintEigs",false);
        args.add("NoMeasure",true);
        args.add("DebugLevel",-1);
        }
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? -1 : 0));
    const int numCenter = args.getInt("NumCenter",2);
    if(numCenter != 1)
        args.add("Truncate",args.getBool("Truncate",true));
    else
        args.add("Truncate",args.getBool("Truncate",false));

    const int N = length(psi);
    Real energy = NAN;

    if((!isOrtho(psi)) || (psi.leftLim() != 0))
        psi.position(1);

    args.add("DebugLevel",debug_level);

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        cpu_time sw_time;
        args.add("Sweep",sw);
        args.add("NSweep",sweeps.nsweep());
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("MaxIter",sweeps.niter(sw));
 
        if(!H.doWrite()
           && args.defined("WriteDim")
           && maxLinkDim(psi) >= args.getInt("WriteDim"))
            {
            if(!quiet)
                {
                println("\nTurning on write to disk, write_dir = ",
                        args.getString("WriteDir","./"));
                }
  
            //psi.doWrite(true);
            H.doWrite(true,args);
            }

        // 0, 1 and 2-site wavefunctions
        ITensor phi0,phi1;
        Spectrum spec;
        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N,{"NumCenter=",numCenter}))
            {
            if(!quiet)
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
  
            H.numCenter(numCenter);
            H.position(b,psi);

            if(numCenter == 2)
                phi1 = psi(b)*psi(b+1);
            else if(numCenter == 1)
                phi1 = psi(b);

            applyExp(H,phi1,t/2,args);

            if(args.getBool("DoNormalize",true))
                phi1 /= norm(phi1);
   
            if(numCenter == 2)
                spec = psi.svdBond(b,phi1,(ha==1 ? Fromleft : Fromright),H,args);
            else if(numCenter == 1)
                psi.ref(b) = phi1;
 
            if((ha == 1 && b+numCenter-1 != N) || (ha == 2 && b != 1))
                {
                auto b1 = (ha == 1 ? b+1 : b);
 
                if(numCenter == 2)
                    {
                    phi0 = psi(b1);
                    }
                else if(numCenter == 1)
                    {
                    Index l;
                    if(ha == 1) l = commonIndex(psi(b),psi(b+1));
                    else        l = commonIndex(psi(b-1),psi(b));
                    ITensor U,S,V(l);
                    spec = svd(phi1,U,S,V,args);
                    psi.ref(b) = U;
                    phi0 = S*V;
                    }
 
                H.numCenter(numCenter-1);
                H.position(b1,psi);
 
                applyExp(H,phi0,-t/2,args);
 
                if(args.getBool("DoNormalize",true))
                    phi0 /= norm(phi0);
                
                if(numCenter == 2)
                    {
                    psi.ref(b1) = phi0;
                    }
                if(numCenter == 1)
                    {
                    if(ha == 1) psi.ref(b+1) *= phi0;
                    else        psi.ref(b-1) *= phi0;
                    }
 
                // Calculate energy
                ITensor H_phi0;
                H.product(phi0,H_phi0);
                energy = real(eltC(dag(phi0)*H_phi0));
                }
            else
                {
                // Calculate energy
                ITensor H_phi1;
                H.product(phi1,H_phi1);
                energy = real(eltC(dag(phi1)*H_phi1));
                }
 
             if(!quiet)
                 { 
                 printfln("    Truncated to Cutoff=%.1E, Min_dim=%d, Max_dim=%d",
                          sweeps.cutoff(sw),
                          sweeps.mindim(sw), 
                          sweeps.maxdim(sw) );
                 printfln("    Trunc. err=%.1E, States kept: %s",
                          spec.truncerr(),
                          showDim(linkIndex(psi,b)) );
                 }

            obs.lastSpectrum(spec);

            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Energy",energy); 
            args.add("Truncerr",spec.truncerr()); 

            obs.measure(args);

            } //for loop over b

        if(!silent)
            {
            auto sm = sw_time.sincemark();
            printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                     sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));
            }
        
        if(obs.checkDone(args)) break;
  
        } //for loop over sw
  
    psi.rightLim(psi.leftLim()+2);
    if(args.getBool("DoNormalize",true))
        psi.normalize();

    return energy;
    }

} //namespace itensor

#endif
