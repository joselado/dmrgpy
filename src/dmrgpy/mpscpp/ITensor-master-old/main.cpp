#include "itensor/all.h"

using namespace itensor;

int 
main()
    {
    int N = 100;

    auto sites = SpinOne(N);

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);

    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;

    auto psi = MPS(sites);

    auto energy = dmrg(psi,H,sweeps);

    println("Ground State Energy = ",energy);

    return 0;
    }
