#include "itensor/all.h"
#include "extra/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace itensor;
using namespace std;



#include"check_task.h" // read the different tasks
#include"get_sites.h" // get the sites from a file
#include"mpsalgebra.h" // functions to deal with MPS
#include"get_sweeps.h" // get the sweep info
#include"get_gap.h" // compute the gap
#include"get_ampo_operator.h" // get an arbitrary AMPO operator
#include"operators.h" // read the different tasks
#include"get_hamiltonian.h" // get the hoppings (in case there are)
#include"read_wf.h" // this does not work yet
#include"get_gs.h" // compute ground state energy and wavefunction
#include"bandwidth.h"  // return the bandwidth of the hamiltonian
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
#include"dynamical_correlator_excited.h" // dynamical correlator with excited
#include"vev.h" // VEV
#include"applyoperator.h" // apply operator to a vector
#include"pureapplyoperator.h" // apply a pure operator to a vector
#include"multmpo.h" // apply a pure operator to a vector
#include"get_random_mps.h" // get a random MPS
#include"conjugatemps.h" // conjugate the MPS


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
    if (check_task("entropy")) get_entropy(); 
    if (check_task("write_sites")) write_sites(); 
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
    if (check_task("exponential_eMwf"))  exponential_eMwf() ; 
    if (check_task("evolution_AeiHtB"))  evolution_AeiHtB() ; // time evolution
    if (check_task("evolution_measure"))  evolution_measure() ; // time evol
    if (check_task("density_matrix"))  reduced_dm() ; // DM
    if (check_task("vev"))  vev() ; // Vacuum expectation value
    if (check_task("applyoperator"))  applyoperator() ; 
    if (check_task("pureapplyoperator"))  pureapplyoperator() ; 
    if (check_task("gen_pureoperator"))  gen_pureoperator() ; 
    if (check_task("multmpo_operator"))  multmpo_operator() ; 
    if (check_task("trace_mpo_operator"))  trace_mpo_operator() ; 
    if (check_task("hermitian_mpo_operator"))  hermitianmpo_operator() ; 
    if (check_task("overlap_aMb"))  overlap_aMb() ; 
    if (check_task("overlap_aMb_static"))  overlap_aMb_static() ; 
    if (check_task("summps"))  get_summps() ; 
    if (check_task("random_mps"))  get_random_mps() ; 
    if (check_task("distribution"))  get_moments_distribution() ; 
    if (check_task("general_kpm"))  general_kpm() ; 
    if (check_task("apply_inverse"))  apply_inverse() ; 
    if (check_task("conjugate_mps"))  conjugate_mps() ; 
//    if (check_task("excited_vev"))  excited_vev() ; // VEV excited
    if (check_task("dynamical_correlator_excited"))  
	    dynamical_correlator_excited(); // DM
    system("rm -f ERROR") ; // remove error file
    return 0;
    }
