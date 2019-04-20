
// return a multioperator of a certain name

static auto get_multioperator=[](std::string name) {
	auto sites = get_sites(); //
	auto A = Iden(sites); // get the identity
        auto n = get_int_value(name+"_n") ; // number of operators
        int maxm = get_int_value("maxm") ; // bond dimension
        auto cutoff = get_float_value("cutoff") ;
	for (int i=0;i<n;i++) { // loop
		auto nameop = get_str(name+"_operator_"+
				std::to_string(i)+"_name");
		int iop = get_int_value(name+"_operator_"+
                                std::to_string(i)+"_site");
		// now read the operator
		auto Ai = get_operator(sites,iop,nameop); // read it
		auto ampot = AutoMPO(sites); // temporal one
		auto mt = MPO(ampot) ; // create MPO
		nmultMPO(Ai,A,mt,{"Maxm",maxm,"Cutoff",cutoff}); // multiply MPO
		A = mt; // reassign
	}; // end loop
	return A;
};
