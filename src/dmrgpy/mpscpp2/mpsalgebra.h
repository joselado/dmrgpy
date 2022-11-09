// chack that two MPS are the same
bool same_mps(auto vi, auto vj) {
        // check if the default should be used
	auto wi = 1.0*vi;
//	wi = wi*1/sqrt(overlap(vi,vi)); // normalize
	auto wj = 1.0*vj;
//	wj = wj*1./sqrt(overlap(vj,vj)); // normalize
	int maxm = get_int_value("maxm") ; // bond dimension for KPM
        auto cutoff = get_float_value("cutoff") ; // bond dimension for KPM
        auto d = sum(wi,-1.0*wj,{"Maxm",maxm,"Cutoff",cutoff}) ; // diff
	auto dd = 1.0; // initialize
	dd = sqrt(overlap(d,d)) ; // compute the norm
	cout << "Vector difference " << dd << endl;
	if (real(dd)<1e-10) return true;
	return false;
};



static auto sum_mpo=[](auto A1,auto A2) {
        int maxm = get_int_value("maxm") ; // bond dimension
        auto cutoff = get_float_value("cutoff") ;
	auto args = Args({"Maxm", maxm, "Cutoff", cutoff});
	auto out = sum(A1,A2, args); // add contribution
	return out;
}
;


static auto mult_mpo=[](auto A1,auto A2) {
        int maxm = get_int_value("mpomaxm") ; // bond dimension
        auto cutoff = get_float_value("cutoff") ;
	auto args = Args({"Maxm", maxm, "Cutoff", cutoff});
	auto sites = get_sites();
	auto ampo = AutoMPO(sites); // temporal one
        auto out = MPO(ampo) ; // create MPO
	nmultMPO(A1,A2,out,args); // multiply MPO
	return out;
}
;



static auto sum_mps=[](auto A1,auto A2) {
	return sum_mpo(A1,A2);
//	auto args = Args({"Maxm", maxm, "Cutoff", cutoff});
//	auto out = sum(A1,A2, args); // add contribution
//	return out;
}
;



