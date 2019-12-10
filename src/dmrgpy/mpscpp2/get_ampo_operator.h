// read an auto Hamiltonian from file

static auto get_ampo_operator=[](auto ampo,std::string filename) {
    ifstream hfile; // file to read
    hfile.open(filename); // file with the operator
    int numterms ; // number of terms in the sum
    auto cr=0.0;
    auto ci=0.0*1i;
    int numprod; // number of terms in the product
    hfile >> numterms; // read the number of terms
    string op0,op1,op2,op3; // names of the operators
    int i0,i1,i2,i3; // index of the operators
    for (int ite=0;ite<numterms;ite++) { // loop over terms
	    hfile >> numprod; // read number of factors and coup
	    if (numprod==1) {  // two terms
		    hfile >> cr >> ci >> op0 >> i0; // read
	            auto cz = cr + ci*1i; // redefine
		    ampo += cz,op0,i0; // add
	    }
	    else if (numprod==2) {  // two terms
		    hfile >> cr >> ci >> op0 >> i0 >> op1 >> i1; // read
	            auto cz = cr + ci*1i; // redefine
		    ampo += cz,op0,i0,op1,i1; // add
	    }
	    else if (numprod==3) {  // two terms
		    hfile >> cr>> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2; 
	            auto cz = cr + ci*1i; // redefine
		    ampo += cz,op0,i0,op1,i1,op2,i2; // add
	    }
	    else if (numprod==4) {  // two terms
		    hfile >> cr>> ci >> op0 >> i0 >> op1 >> i1 >> op2 >> i2 >> op3 >>i3; 
	            auto cz = cr + ci*1i; // redefine
		    ampo += cz,op0,i0,op1,i1,op2,i2,op3,i3; // add
	    }
	    else exit(EXIT_FAILURE) ;
    };
    hfile.close();
    return ampo ;
};



static auto get_mpo_operator=[](std::string filename) {
	auto sites = get_sites(); // read sites
	auto ampo = AutoMPO(sites); // generate ampo
	ampo = get_ampo_operator(ampo,filename) ; // generate ampo
	return MPO(ampo); // return mpo
};



