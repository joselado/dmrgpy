class SpinX : public SiteSet
    {
    public:


    SpinX(Args const& args = Args::global());


    void
    read(std::istream& s);

    };


/// this does not work properly
void inline SpinX::
read(std::istream& s)
    {
    int N = itensor::read<int>(s);
    if(N > 0)
        {
        auto store = SiteStore(N);
        ifstream sfile; // file to read
        sfile.open("sites.in"); // file with the coupling
	int N2 = 0;
        sfile >> N2; // read the number of sites from file
        auto sites = SiteStore(N); // get an empty list of sites
	int nm = 0; // initialize
        for(int j = 1; j <= N; ++j)
            {
            auto I = IQIndex{};
            sfile >> nm ; // read this spin
            I.read(s);
	    // S=3/2 is being confusing by spinful fermion!!!!
            if(nm == 3) store.set(j,SpinOneSite(I));
            else if(nm == 1) store.set(j,HubbardSite(I));
            else if(nm == 2) store.set(j,SpinHalfSite(I));
            else if(nm == 4) store.set(j,SpinThreeHalfSite(I));
            else if(nm == 5) store.set(j,SpinTwoSite(I));
            else if(nm == 6) store.set(j,SpinFiveHalfSite(I));
            else Error(format("SpinX cannot read index of size %d",nm));
            }
	sfile.close() ;
        init(std::move(store));
        }
    }

//} //namespace itensor







inline SpinX::
SpinX(Args const& args)
    {
    ifstream sfile; // file to read
    sfile.open("sites.in"); // file with the coupling
    int N, nm;
    sfile >> N; // read the number of sites and number of projections
    auto sites = SiteStore(N); // get an empty list of sites
    for (int i=1;i<=N;i++)  {
      sfile >> nm ; // read this spin
      cout << "Reading  " << nm << endl;
      if (nm==2) sites.set(i,SpinHalfSite(i)); // use spin=1/2
      else if (nm==1) sites.set(i,HubbardSite(i)); // use fermions
      else if (nm==3) sites.set(i,SpinOneSite(i)); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i)); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite(i)); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i)); // use spin=5/2
      else Error(format("SpinX cannot read index of size "));
    } ;
    sfile.close(); // close file

    SiteSet::init(std::move(sites));
    }




auto generate_sites() { // function to generate the sites
    auto sites = SpinX() ; // read from file
    return sites ;
}



auto get_sites() { // function to get the sites
    auto sites = generate_sites() ;  // generate the sites
    if (check_task("gs_from_file")) readFromFile("sites.sites",sites);
    cout << "Number of sites " << sites.N() << endl ;
    return sites ;
}



int site_type(int index) {
    ifstream sfile; // file to read
    sfile.open("sites.in"); // file with the sites
    int N, nm, out=-1;
    sfile >> N; // read the number of sites and number of projections
    for (int i=1;i<=N;i++)  {
      sfile >> nm ; // read this spin
      if (i-1==index) out = nm ; }
    sfile.close() ;
    cout << index << "  " << out << endl ;
    return out ;
}





