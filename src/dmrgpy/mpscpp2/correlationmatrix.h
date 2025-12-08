

static auto get_correlation_matrix=[]() {
	auto psi = read_wf("wavefunction.mps"); // read the wavefunction
	auto sites = get_sites(); // get the sites
	auto N = sites.N() ; // number of sites
	// create storage matrix
	std::vector<std::vector<std::complex<double>>> matrix(
        N, std::vector<std::complex<double>>(N)
    );
	for (int i=0;i<N;i++) {
	   for (int j=i;j<N;j++) {
		auto ampo = AutoMPO(sites);  // create MPO
		cout << i << "  " << j << endl ;
		ampo += 1.0,"Cdag",i+1,"C",j+1 ; // add contribution
		auto op = MPO(ampo) ; // convert to MPO
		auto c = overlapC(psi,op,psi) ;
		matrix[i][j] = c ; // store
		matrix[j][i] = std::conj(c) ; // store complex conjugate
	   };
	}
	// write in a file
	ofstream myfile; // open file
        myfile.open("CORRELATION_MATRIX.OUT"); // open file
	for (int i=0;i<N;i++) {
	   for (int j=0;j<N;j++) {
                auto c = matrix[i][j] ; // get element
		myfile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
	   }
	}
	myfile.close(); // close file
	return 0;
};
