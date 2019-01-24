// compute dynamical correlators using the correction vector algorithm

#include"cg.h" // conjugate gradient function

static auto dynamical_correlator=[]() {
  auto sites = get_sites(); // get the sites
  auto H = get_hamiltonian(sites) ; // get the Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  // get the two operator
  auto opi = get_operator(sites,i,get_str("dynamical_correlator_operator_i")) ;
  auto opj = get_operator(sites,j,get_str("dynamical_correlator_operator_j")) ;
  auto delta = get_float_value("delta_dynamical_correlator") ; // smearing
  auto pp = {"Maxm",maxm,"Cutoff",cutoff} ; // DMRG parameters
  auto wfi =  exactApplyMPO(psi,opj,pp) ; // Apply operator
  auto v = solveHwdb(H,energy,delta,wfi) ; // get the correction vector
  // this is not finished

};

