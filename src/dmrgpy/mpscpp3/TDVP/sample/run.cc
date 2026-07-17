#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

int main()
    {
    // The system will be a 2x20 ladder
    int Nx = 20;
    int Ny = 2;
    int N = Nx*Ny;
    auto t = 0.1;
    auto tend = 0.5;
    int nsw = tend/t;
    auto t0 = 0.01;

    // Make N spin 1/2's
    auto sites = SpinHalf(N);

    // Make the Hamiltonian for rung-decoupled Heisenberg ladder
    auto ampo = AutoMPO(sites);    
    for(int i = 1; i <= N-2; ++ i)
        {
        ampo += 0.5,"S+",i,"S-",i+2;
        ampo += 0.5,"S-",i,"S+",i+2;
        ampo +=     "Sz",i,"Sz",i+2;
        }
    auto H = toMPO(ampo);
    printfln("Maximum bond dimension of H is %d",maxLinkDim(H));

    // Set the initial state to be Neel state
    auto state = InitState(sites);
    state.set(1,"Up");
    for(int i = 2; i < N; i=i+2)
        {
        if((i/2)%2==1)
            {
            state.set(i,"Dn");
            state.set(i+1,"Dn");
            }
        else
            {
            state.set(i,"Up");
            state.set(i+1,"Up");
            }
        }
    state.set(N,"Up");
    
    auto psi1 = MPS(state);
    auto psi2 = psi1;

    // start TDVP, either one site or two site algorithm can be used by adjusting the "NumCenter" argument
    println("----------------------------------------GSE-TDVP---------------------------------------");

    auto energy = real(innerC(psi1,H,psi1));
    printfln("Initial energy = %.5f", energy);

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 10;

    for(int n = 1; n <= nsw; ++n)
        {
        if(n < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi1,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",true,
                                      "Quiet",true});
            }
        
        // TDVP sweep
        energy = tdvp(psi1,H,-t,sweeps,{"Truncate",true,
                                        "DoNormalize",true,
                                        "Quiet",true,
                                        "NumCenter",1,
                                        "ErrGoal",1E-7});
        }

    printfln("\nEnergy after imaginary time evolution = %.10f",energy);
    printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );

    println("-------------------------------------MPO W^I 2nd order---------------------------------------");

    auto expH1 = toExpH(ampo,(1-1_i)/2*t0);
    auto expH2 = toExpH(ampo,(1+1_i)/2*t0);
    printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","DensityMatrix","Cutoff=",1E-12,"MaxDim=",2000);
    for(int n = 1; n <= nsw*std::real(t/t0); ++n)
        {
        psi2 = applyMPO(expH1,psi2,args);
        psi2.noPrime();
        psi2 = applyMPO(expH2,psi2,args);
        psi2.noPrime().normalize();
        if(n%int(std::real(t/t0)) == 0)
            {
            printfln("\nMaximum bond dimension at time %.1f is %d ", n*t0, maxLinkDim(psi2));
            printfln("Energy using overlap at time %.1f is %.10f", n*t0, real(innerC(psi2,H,psi2)) );
            }
        }

    return 0;
    }
