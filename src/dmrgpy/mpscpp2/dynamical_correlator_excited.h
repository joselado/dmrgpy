// compute a dynamical correlator using excited states
//
//

int dynamical_correlator_excited() {
	auto sites = get_sites(); // Get the different sites
	auto H = get_hamiltonian(sites) ; // get the Hamiltonian
	auto sweeps = get_sweeps(); // get sweeps
	auto nexcited = get_int_value("nexcited") ; // get excited states
	auto wfs = get_excited(); // get excited states
	nexcited = wfs.size() ; // get excited states
        // now compute the different matrix elements
	auto A1 = get_mpo_operator("dc_multioperator_i.in");
	auto A2 = get_mpo_operator("dc_multioperator_j.in");
	auto psi0 = get_gs(); // ground state
	// open files
	ofstream fileoverlap;
        fileoverlap.open("EXCITED_OVERLAPS.OUT"); // open file
	for(int i=1;i<nexcited;i++) {
		auto c1 = overlapC(wfs.at(i),A1,psi0); // compute overlap
		auto c2 = overlapC(wfs.at(i),A2,psi0); // compute overlap
		fileoverlap << std::setprecision(20) << real(c1) << "  "; 
		fileoverlap << std::setprecision(20) << imag(c1) << "  "; 
		fileoverlap << std::setprecision(20) << real(c2) << "  "; 
		fileoverlap << std::setprecision(20) << imag(c2) << "  "; 
		fileoverlap << endl ; // next line 
	}
	fileoverlap.close(); // close file
	// now that we have the matrix elements, compute the correlator
	return 0; // return
}
