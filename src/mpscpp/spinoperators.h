auto get_spin_operator(auto sites, int i, auto name) {
	auto ampo = AutoMPO(sites);
	if (site_type(i)!=1)  // spin site
		ampo += 1.0,name,i+1 ;
	else if (site_type(i)==1) { // fermionic site
        	if (name=="Sx") {
        		ampo += 0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += 0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
		else if (name=="Sy") {
        		ampo += -0.5,"Cdagdn",i+1,"CupI",i+1;
        		ampo += 0.5,"Cdagup",i+1,"CdnI",i+1;
        	} ;
		else if (name=="Sz") {
        		ampo += 0.5,"Cdagup",i+1,"Cup",i+1;
        		ampo += -0.5,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
		else if (name=="Cdag") {
        		ampo += 1.0,"Cdagup",i+1;
        		ampo += 1.0,"Cdagdn",i+1;
        	} ;
		else if (name=="C") {
        		ampo += 1.0,"Cup",i+1;
        		ampo += 1.0,"Cdn",i+1;
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








auto add_spin_operator(auto ampo, auto sites, float v, int i, auto name) {
	if (site_type(i)!=1)  // spin site
		ampo += v,name,i+1 ;
	else if (site_type(i)==1) { // fermionic site
        	if (name=="Sx") {
        		ampo += v*0.5,"Cdagdn",i+1,"Cup",i+1;
        		ampo += v*0.5,"Cdagup",i+1,"Cdn",i+1;
        	} ;
        	if (name=="Sy") {
        		ampo += -v*0.5,"Cdagdn",i+1,"CupI",i+1;
        		ampo += v*0.5,"Cdagup",i+1,"CdnI",i+1;
        	} ;
        	if (name=="Sz") {
        		ampo += v*0.5,"Cdagup",i+1,"Cup",i+1;
        		ampo += -v*0.5,"Cdagdn",i+1,"Cdn",i+1;
        	} ;
        }
	return ampo ;
}




