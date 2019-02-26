
// STEFAN: uncomment this to put your header
#include"cvm.h" // Time evolution

// STEFAN: comment this dummy function to put yours
//static auto spectral_function=[](auto a1,auto a2,auto a3,auto a4,
//		auto a5, auto a6, auto a7, auto a8) {
//	return 1.0 ; } ;


// compute the dynamical correlator with CVM
static auto cvm_dynamical_correlator=[](auto sites) {
  auto name1 = get_str("cvm_operator_i");  // first operator
  auto name2 = get_str("cvm_operator_j");  // second operator
  auto i1 = get_int_value("cvm_site_i");  // second operator
  auto i2 = get_int_value("cvm_site_j");  // second operator
  // now get the operators
  auto A1 = get_operator(sites,i1,name1); // first operator
  auto A2 = get_operator(sites,i2,name2); // second operator
  auto H = get_hamiltonian(sites) ; // get the ampo for the Hamiltonian
  auto psi = get_gs(sites,H) ; // get the ground state
  int maxm = get_int_value("maxm") ; // bond dimension
  int max_it = get_int_value("cvm_nit") ; // number of iterations
  auto cutoff = get_float_value("cutoff") ; // cutoff of DMRG
  auto delta = get_float_value("cvm_delta") ; // delta
  auto e0 = get_float_value("cvm_e0") ; // GS energy
  auto tol = get_float_value("cvm_tol") ; // GS energy
  auto omega = get_float_value("cvm_energy") ; // frequency
  auto args = Args("Cutoff=",cutoff,"Maxm=",maxm); // MPS arguments
  // STEFAN: Here I am calling your function
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
