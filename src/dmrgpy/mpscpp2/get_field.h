
static auto get_field = [](auto sites, auto ampo) {
    ifstream bfile; // file to read
    bfile.open("fields.in"); // file with the fields
    int nb=0;
    bfile >> nb; // read the number of fields
    int index; // declare index
    float bx=0.0,by=0.0,bz=0.0; // declare fields
    for (int i=0;i<nb;++i) {
      bfile >> index >> bx >> by >> bz; // read this coupling
//      if ((site_type(index-1)==1)  ;
      ampo = add_spin_operator(ampo,sites,bx,index,"Sx"); // 
      ampo = add_spin_operator(ampo,sites,by,index,"Sy"); // 
      ampo = add_spin_operator(ampo,sites,bz,index,"Sz"); // 
//      ampo += bx,"Sx",index+1; // add magnetic field
//      ampo += by,"Sy",index+1; // add magnetic field
//      ampo += bz,"Sz",index+1; // add magnetic field
    } ;
    bfile.close() ; // close file
    return ampo ;  // return the Hamiltonian with exchange added
}
;