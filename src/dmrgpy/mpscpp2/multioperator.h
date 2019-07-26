
// return a multioperator of a certain name

static auto get_multioperator=[](std::string name, auto sites) {
//	auto sites = get_sites(); //
	auto out = Iden(sites); // get the identity
	out = 0.*out; // initialize operator
        int maxm = get_int_value("maxm") ; // bond dimension
        auto cutoff = get_float_value("cutoff") ;
	auto args = Args({"Maxm", maxm, "Cutoff", cutoff});
        auto nop = get_int_value(name+"_n") ; // number of operators
	// loop over the different terms in the sum
	for (int j=0;j<nop;j++) { // loop over operators
      	  auto namej = name+"_operator_"+std::to_string(j); // get this one
	  int n = get_int_value(namej+"_nterms"); // number of terms
	  auto A = Iden(sites); // get the identity as starting point
	  auto cr = get_float_value(namej+"_coefficient_real") ; // real
	  auto ci = get_float_value(namej+"_coefficient_imag") ; // imaginary
	  auto cz = cr + 1i*ci; // complex coefficient
	  A = cz*A; // multiply by the coefficient
  	  for (int i=0;i<n;i++) { // loop over the term of this coefficient
		auto nameop = get_str(namej+"_term_"+
				std::to_string(i)+"_name");
		int iop = get_int_value(namej+"_term_"+
                                std::to_string(i)+"_site");
		// now read the operator
		auto Ai = get_operator(sites,iop,nameop); // read it
		auto ampot = AutoMPO(sites); // temporal one
		auto mt = MPO(ampot) ; // create MPO
		nmultMPO(Ai,A,mt,args); // multiply MPO
		A = mt; // reassign
	  };
	  // now add the contribution to the output
	auto tmp = sum(out,A, args); // add contribution
	out = tmp; // reassign total result
	}; // end loop
	return out; // return output
};
