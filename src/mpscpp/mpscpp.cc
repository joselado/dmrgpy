#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace itensor;
using namespace std;



#include"check_task.h" // read the different tasks
#include"get_sweeps.h" // get the sweep info
#include"get_gap.h" // compute the gap
#include"get_sites.h" // get the sites from a file
#include"spinoperators.h" // read the different tasks
#include"get_hamiltonian.h" // get the hoppings (in case there are)
#include"read_wf.h" // this does not work yet
#include"get_gs.h" // compute ground state energy and wavefunction
#include"get_dos.h" // compute the DOS
#include"get_field.h" // get local magnetic fields
#include"get_correlator.h" // compute correaltors betwee sites
#include"measure.h" // compute expectation values
#include"get_entropy.h" // compute entanglement entropy
#include"get_excited.h" // compute excited states
#include"kpm.h" // KPM routines
#include"compute_overlap.h" // Compute overlap


int 
main()
    {
    system("touch ERROR") ; // create error file


    // read the number of sites
    ifstream sfile; // file to read
    auto sites = get_sites(); // Get the different sites
    auto H = get_hamiltonian(sites) ; // get the Hamiltonian
//    test_hopping(H,sites); // test the hoppings
    auto sweeps = get_sweeps(); // get the DMRG sweeps

    auto psi = MPS(sites);   // create random wavefunction
//    auto psi = MPS(sites);   // create random wavefunction
    if (check_task("GS")) {
      psi = get_gs(sites,H) ; // get ground state wavefunction
      measure(psi,sites) ; // and compute expectation values
    } ;
//    if (not(check_task("GS"))) psi = read_wf() ; // read the wavefunction from a file
    if (check_task("correlator")) get_correlator(psi,sites) ; // write correlators 
    if (check_task("gap")) get_gap(H,sites,sweeps); // calculate the gap 
    if (check_task("entropy")) get_entropy(psi,2); 
    if (check_task("excited")) {
      int nexcited = get_int_value("nexcited") ; // number of excited states
      auto wfs = get_excited(H,sites,sweeps,nexcited); // compute states 
    } ;
    if (check_task("dos")) get_dos(H,sites) ; // get the DOS
    if (check_task("spismj"))  { // dynamical correlation function
       if (check_task("smart_kpm_window"))
         get_moments_spismj(sites,H,get_int_value("nkpm"),
           get_int_value("site_i_kpm")+1,get_int_value("site_j_kpm")+1) ;
       else 
         get_moments_spismj_brute(sites,H,get_int_value("nkpm"),
           get_int_value("site_i_kpm")+1,get_int_value("site_j_kpm")+1) ; 
    }
    if (check_task("dynamical_correlator"))  { // dynamical correlation
       get_moments_dynamical_correlator(sites,H,get_int_value("nkpm"),
         get_int_value("site_i_kpm"),get_int_value("site_j_kpm"),
         get_str("kpm_operator_i"),get_str("kpm_operator_j")) ; 
    } ;
// overlap task
    if (check_task("overlap"))  compute_overlap() ; // compute overlap
    system("rm -f ERROR") ; // create error file
    return 0;
    }
