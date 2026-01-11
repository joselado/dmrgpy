

static auto get_four_correlation_tensor=[]() {
	auto psi = read_wf("wavefunction.mps"); // read the wavefunction
	auto sites = get_sites(); // get the sites
	auto N = sites.N() ; // number of sites
	// create storage mfour index tensor
	std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> tensor(
    N, std::vector<std::vector<std::vector<std::complex<double>>>>(
        N, std::vector<std::vector<std::complex<double>>>(
            N, std::vector<std::complex<double>>(N)
        )
    )
);
	for (int i=0;i<N;i++) {
	   for (int j=0;j<N;j++) {
	     for (int k=0;k<N;k++) {
	       for (int l=0;l<N;l++) {
        		auto ampo = AutoMPO(sites);  // create MPO
			// add contribution
        		ampo += 1.0,"Cdag",i+1,"C",j+1,"Cdag",k+1,"C",l+1 ; 
        		auto op = MPO(ampo) ; // convert to MPO
        		auto c = overlapC(psi,op,psi) ;
        		tensor[i][j][k][l] = c ; // store
	   };
	   };
	   };
	}
	// write in a file
	ofstream myfile; // open file
        myfile.open("FOUR_CORRELATION_TENSOR.OUT"); // open file
	for (int i=0;i<N;i++) {
	   for (int j=0;j<N;j++) {
	     for (int k=0;k<N;k++) {
	       for (int l=0;l<N;l++) {
                auto c = tensor[i][j][k][l] ; // get element
		myfile << std::setprecision(20) << real(c) << "  " << imag(c) << endl ;
	   };
	   };
	   };
	}
	myfile.close(); // close file
	return 0;
};
