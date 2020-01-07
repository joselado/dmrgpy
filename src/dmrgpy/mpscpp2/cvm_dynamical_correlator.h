
// CVM library
#include"cvm.h" // Time evolution



// compute the dynamical correlator with CVM
static auto cvm_dynamical_correlator=[]() {
  auto sites = get_sites();
  // now get the operators
  auto A1 = get_mpo_operator("dc_multioperator_i.in");
  auto A2 = get_mpo_operator("dc_multioperator_j.in");
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = get_gs() ; // get the ground state
  int maxm = get_int_value("maxm") ; // bond dimension
  int max_it = get_int_value("cvm_nit") ; // number of iterations
  auto cutoff = get_float_value("cutoff") ; // cutoff of DMRG
  auto delta = get_float_value("cvm_delta") ; // delta
  auto e0 = get_float_value("cvm_e0") ; // GS energy
  auto tol = get_float_value("cvm_tol") ; // GS energy
  auto omega = get_float_value("cvm_energy") ; // frequency
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm); // MPS arguments
  // Call the CVM function
  auto z = spectral_function(psi,H,A1,A2,omega,delta,
		  e0,tol,max_it,maxm,cutoff,sites); // return the correlator
  // end of the call
  ofstream filecvm;
  filecvm.open("CVM.OUT"); // open file
  filecvm << std::setprecision(8) << real(z) << "  "
                       << std::setprecision(8)<< imag(z) << endl;
  filecvm.close(); // close file
} ;
