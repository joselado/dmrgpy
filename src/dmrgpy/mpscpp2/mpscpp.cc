#include "itensor/all.h"
#include "extra/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace itensor;
using namespace std;



#include"check_task.h" // read the different tasks
#include"mpsalgebra.h" // functions to deal with MPS
#include"get_sweeps.h" // get the sweep info
#include"bandwidth.h"  // return the bandwidth of the hamiltonian
#include"get_gap.h" // compute the gap
#include"get_sites.h" // get the sites from a file
#include"get_ampo_operator.h" // get an arbitrary AMPO operator
#include"operators.h" // read the different tasks
#include"get_hamiltonian.h" // get the hoppings (in case there are)
#include"read_wf.h" // this does not work yet
#include"get_gs.h" // compute ground state energy and wavefunction
#include"get_excited.h" // compute excited states
//#include"get_dos.h" // compute the DOS
//#include"get_correlator.h" // compute correaltors betwee sites
#include"measure.h" // compute expectation values
#include"get_entropy.h" // compute entanglement entropy
#include"kpm.h" // KPM routines
#include"compute_overlap.h" // Compute overlap
#include"cvm_dynamical_correlator.h" // CVM dynamical correlator
#include"time_evolution.h" // Time evolution
#include"reduced_dm.h" // Reduced density matrix
#include"dynamical_correlator_excited.h" // dynamical correlator with exited
#include"vev.h" // Reduced density matrix


int 
main()
    {
    system("touch ERROR") ; // create error file


    // read the number of sites
    ifstream sfile; // file to read
//    test_hopping(H,sites); // test the hoppings
    if (check_task("GS")) {
      get_gs() ; // get ground state wavefunction
//      measure(psi,sites) ; // and compute expectation values
    } ;
//    if (not(check_task("GS"))) psi = read_wf() ; // read the wavefunction from a file
//    if (check_task("correlator")) get_correlator() ; // write correlators 
//    if (check_task("gap")) get_gap(H,sites,sweeps); // calculate the gap 
//    if (check_task("entropy")) get_entropy(psi,2); 
    if (check_task("excited")) {
      get_excited(); // compute states 
    } ;
//    if (check_task("dos")) get_dos(H,sites) ; // get the DOS
    if (check_task("dynamical_correlator"))  { // dynamical correlation
       get_moments_dynamical_correlator();
    } ;
//    if (check_task("dos"))  { // global DOS
//       get_moments_dos(sites,H);
//    } ;
    if (check_task("cvm"))  cvm_dynamical_correlator() ; // CVM
// overlap task
    if (check_task("overlap"))  compute_overlap() ; // compute overlap
    if (check_task("time_evolution"))  quench() ; // time evolution
    if (check_task("density_matrix"))  reduced_dm() ; // DM
    if (check_task("vev"))  vev() ; // Vacuum expectation value
//    if (check_task("excited_vev"))  excited_vev() ; // VEV excited
    if (check_task("dynamical_correlator_excited"))  
	    dynamical_correlator_excited(); // DM
    system("rm -f ERROR") ; // remove error file
    return 0;
    }
