// compute a dynamical correlator using excited states
//
//

int dynamical_correlator_excited() {
	auto sites = get_sites(); // Get the different sites
	auto H = get_hamiltonian(sites) ; // get the Hamiltonian
	auto sweeps = get_sweeps(); // get sweeps
	auto nexcited = get_int_value("nexcited") ; // get excited states
	auto wfs = get_excited(H,sites,sweeps,nexcited); // get excited states
        // now compute the different matrix elements
	auto namei = get_str("operator_i");
        auto namej = get_str("operator_j");
	auto ii = get_int_value("site_i");
        auto jj = get_int_value("site_j");
	auto A1 = get_operator(sites,ii,namei); // first operator
	auto A2 = get_operator(sites,jj,namej); // second operator
	auto psi0 = wfs.at(0); // ground state
	// open files
	ofstream fileoverlap;
        fileoverlap.open("EXCITED_OVERLAPS.OUT"); // open file
	for(int i=1;i<nexcited;i++) {
		auto c1 = overlapC(wfs.at(i),A1,psi0); // compute overlap
		auto c2 = overlapC(wfs.at(i),A2,psi0); // compute overlap
		fileoverlap << std::setprecision(8) << real(c1) << "  "; 
		fileoverlap << std::setprecision(8) << imag(c1) << "  "; 
		fileoverlap << std::setprecision(8) << real(c2) << "  "; 
		fileoverlap << std::setprecision(8) << imag(c1) << "  "; 
		fileoverlap << endl ; // next line 
	}
	fileoverlap.close(); // close file
	// now that we have the matrix elements, compute the correlator
	
}
