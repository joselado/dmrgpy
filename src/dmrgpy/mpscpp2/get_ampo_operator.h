// read an auto Hamiltonian from file

static auto get_ampo_operator=[](auto ampo,std::string filename) {
    ifstream hfile; // file to read
    hfile.open(filename); // file with the operator
    int numterms ; // number of terms in the sum
    auto cr=0.0;
    auto ci=0.0*1i;
    int numprod; // number of terms in the product
    hfile >> numterms; // read the number of terms
    for (int ite=0;ite<numterms;ite++) { // loop over terms
	    hfile >> numprod; // read number of factors and coup
            #include"extra/ampotk/ampotk.h"
    };
    hfile.close();
    return ampo ;
};



static auto get_mpo_operator=[](std::string filename) {
	auto sites = get_sites(); // read sites
	auto ampo = AutoMPO(sites); // generate ampo
	ampo = get_ampo_operator(ampo,filename) ; // generate ampo
        int mpomaxm = get_int_value("mpomaxm") ; // get bond dimension
	if (mpomaxm<5) { mpomaxm = 5000 ;}; // default value
	return toMPO<ITensor>(ampo,{"Maxm",mpomaxm,"Exact",false}); // return mpo
};

