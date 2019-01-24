static auto get_dos=[](auto H, auto sites) {
	int nexcited = get_int_value("nexcited") ; // number of excited states
	int isite = get_int_value("dos_site") ; // number of excited states
	auto sweeps = get_sweeps(); // get the sweeps
	auto wfs = get_excited(H,sites,sweeps,nexcited) ; // get states
        auto op = get_up_creation_operator(sites,isite) ; // get operator 
	ofstream myfile;
        myfile.open("DOS_ELEMENT.OUT"); // open file
	// compute the neccesary matrix elements
	for (int i=0;i<nexcited;i++) {
	   for (int j=0;j<nexcited;j++) {
		auto c = overlapC(wfs.at(i),op,wfs.at(j)) ;
		myfile << real(c) << "  " << imag(c) << endl ;
	   };
	}
	myfile.close(); // close file
	return 0;
};
