auto get_spin_operator(auto sites, int i, auto name) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)!=1)  // spin site
		ampo += 1.0,name,i+1 ;
	else if (site_type(i)==1) { // fermionic site
        	if (name=="Sx") {
        		ampo += 0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += 0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
		if (name=="Sy") {
        		ampo += -0.5*1i,"Cdagdn",i+1,"Cup",i+1;
        		ampo += 0.5*1i,"Cdagup",i+1,"Cdn",i+1;
        	} ;
		if (name=="Sz") {
        		ampo += 0.5,"Cdagup",i+1,"Cup",i+1;
        		ampo += -0.5,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
		if (name=="Cdag") {
        		ampo += 1.0,"Cdagup",i+1;
        		ampo += 1.0,"Cdagdn",i+1;
        	} ;
		if (name=="C") {
        		ampo += 1.0,"Cup",i+1;
        		ampo += 1.0,"Cdn",i+1;
        	} ;
		if (name=="Cup") {
        		ampo += 1.0,"Cup",i+1;
        	} ;
		if (name=="Cdn") {
        		ampo += 1.0,"Cdn",i+1;
        	} ;
		if (name=="Cdagup") {
        		ampo += 1.0,"Cdagup",i+1;
        	} ;
		if (name=="Cdagdn") {
        		ampo += 1.0,"Cdagdn",i+1;
        	} ;
		// superconducting terms
		if (name=="delta") {
        		ampo += 1.0,"Cdn",i+1,"Cup",i+1;
        	} ;
		if (name=="deltad") {
        		ampo += 1.0,"Cdagdn",i+1,"Cdagup",i+1;
        	} ;
		// density terms
		if (name=="density") {
        		ampo += 1.0,"Cdagdn",i+1,"Cdn",i+1;
        		ampo += 1.0,"Cdagup",i+1,"Cup",i+1;
        	} ;
		if (name=="density_up") {
        		ampo += 1.0,"Cdagup",i+1,"Cup",i+1;
        	} ;
		if (name=="density_dn") {
        		ampo += 1.0,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
        }
        auto m = MPO(ampo) ;	
return m ;
}



auto get_occupation_operator(auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Ntot",i+1;
        auto m = MPO(ampo) ;	
return m ;
}




auto get_occupation2_operator(auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Ntot",i+1,"Ntot",i+1;
        auto m = MPO(ampo) ;	
return m ;
}





auto get_delta_operator(auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Cdn",i+1,"Cup",i+1;
        	ampo += 1.0,"Cdagup",i+1,"Cdagdn",i+1;
        auto m = MPO(ampo) ;	
return m ;
}


auto get_up_creation_operator(auto sites, int i) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Cdagup",i+1;
        auto m = MPO(ampo) ;	
return m ;
}



auto get_hopping_operator(auto sites, int i, int j) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)==1)  // fermionic site
        	ampo += 1.0,"Cdagup",i+1,"Cup",j+1;
        	ampo += 1.0,"Cdagdn",i+1,"Cdn",j+1;
        auto m = MPO(ampo) ;	
return m ;
}



auto get_sidotsj_operator(auto sites, int i, int j) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)!=1)  // spin site
        	ampo += 1.0,"Sx",i+1,"Sx",j+1;
        	ampo += 1.0,"Sy",i+1,"Sy",j+1;
        	ampo += 1.0,"Sz",i+1,"Sz",j+1;
        auto m = MPO(ampo) ;	
return m ;
}




auto add_spin_operator(auto ampo, auto sites, float v, int i, auto name) {
	if (site_type(i)!=1)  // spin site
		ampo += v,name,i+1 ;
	else if (site_type(i)==1) { // fermionic site
        	if (name=="Sx") {
        		ampo += v*0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += v*0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
        	if (name=="Sy") {
        		ampo += -v*0.5*1i,"Cdagdn",i+1,"Cup",i+1;
        		ampo += v*0.5*1i,"Cdagup",i+1,"Cdn",i+1;
        	} ;
        	if (name=="Sz") {
        		ampo += v*0.5,"Nup",i+1;
        		ampo += -v*0.5,"Ndn",i+1;
        	} ;
        }
	return ampo ;
}




