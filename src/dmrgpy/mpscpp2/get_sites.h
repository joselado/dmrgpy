class SpinX : public SiteSet
    {
    public:


    SpinX(Args const& args = Args::global());

    // Build directly from in-memory per-site type codes, with no file I/O
    // at all (same type-code convention as the sites.in-based constructor
    // below). This is what the in-process Chain session (chain_session.h)
    // uses instead of writing/reading sites.in + sites.sites.
    SpinX(std::vector<int> const& site_types);


    void
    read(std::istream& s);

    };


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
            auto I = IQIndex{}; // generate an identity for this index
            sfile >> nm ; // read the label of the site from sites.in
            I.read(s); // read the site from sites.sites
	    cout << "Reading site from sites.sites "<< nm << endl ;
	    cout << "Read from sites.sites "<< I << endl ;
            if(nm == 3) store.set(j,SpinOneSite(I));
            else if (nm==0) store.set(j,SpinlessSite(I)); // use fermions
            else if(nm == 1) store.set(j,HubbardSite(I));
            else if(nm == 2) store.set(j,SpinHalfSite(I));
            else if(nm == 4) store.set(j,SpinThreeHalfSite(I));
            else if(nm == 5) store.set(j,SpinTwoSite(I));
            else if(nm == 6) store.set(j,SpinFiveHalfSite(I));
            else if(nm == -2) store.set(j,Z3Site(I));
            else if(nm == -3) store.set(j,Z4Site(I));
            else if(nm == 104) store.set(j,BosonFourSite(I));
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
      sfile >> nm ; // read the name of that site
      cout << "Reading  " << nm << endl; // record this
      if (nm==2) sites.set(i,SpinHalfSite(i)); // use spin=1/2
      else if (nm==0) sites.set(i,SpinlessSite(i)); // use fermions
      else if (nm==1) sites.set(i,HubbardSite(i)); // use spinful fermions
      else if (nm==3) sites.set(i,SpinOneSite(i)); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i)); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite(i)); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i)); // use spin=5/2
      else if (nm==104) sites.set(i,BosonFourSite(i)); // use Boson
      else if (nm==(-2)) sites.set(i,Z3Site(i)); // use Z3
      else if (nm==(-3)) sites.set(i,Z4Site(i)); // use Z4
      else Error(format("SpinX cannot read index of size "));
    } ;
    sfile.close(); // close file

    SiteSet::init(std::move(sites));
    }


inline SpinX::
SpinX(std::vector<int> const& site_types)
    {
    int N = site_types.size();
    auto sites = SiteStore(N); // get an empty list of sites
    for (int i=1;i<=N;i++)  {
      int nm = site_types.at(i-1); // read the name of that site
      if (nm==2) sites.set(i,SpinHalfSite(i)); // use spin=1/2
      else if (nm==0) sites.set(i,SpinlessSite(i)); // use fermions
      else if (nm==1) sites.set(i,HubbardSite(i)); // use spinful fermions
      else if (nm==3) sites.set(i,SpinOneSite(i)); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i)); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite(i)); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i)); // use spin=5/2
      else if (nm==104) sites.set(i,BosonFourSite(i)); // use Boson
      else if (nm==(-2)) sites.set(i,Z3Site(i)); // use Z3
      else if (nm==(-3)) sites.set(i,Z4Site(i)); // use Z4
      else Error(format("SpinX cannot read index of size "));
    } ;

    SiteSet::init(std::move(sites));
    }




auto generate_sites() { // function to generate the sites
    auto sites = SpinX() ; // read from file
    return sites ;
}



auto get_sites() { // function to get the sites
    auto sites = generate_sites() ;  // generate the sites
    // overwrite the sites
    if (check_task("gs_from_file")) readFromFile("sites.sites",sites);
    if (check_task("sites_from_file")) readFromFile("sites.sites",sites);
    cout << "Number of sites " << sites.N() << endl ;
    return sites ;
}


auto write_sites() { // function to get the sites
    auto sites = get_sites(); // Get the different sites
    writeToFile("sites.sites",sites); // write the sites
}



// Explicit cache replacing the previous hidden per-function statics.
// Same lifetime/semantics as before (populated on first use, reused for the
// rest of one mpscpp.x invocation), but named and resettable so a future
// persistent-process backend can give each independent chain/session its
// own instance instead of silently sharing this one across chains of
// different site-type composition.
struct SiteTypeCache
    {
    bool loaded = false; // whether sites.in has been read
    int N = 0; // number of sites
    std::vector<int> stypes; // site type code per site
    };

static SiteTypeCache site_type_cache;

void inline
reset_site_type_cache()
    {
    site_type_cache = SiteTypeCache(); // drop any cached site types
    }

int site_type(int index) {
    if (!site_type_cache.loaded) {
      ifstream sfile; // file to read
      int nm;
      sfile.open("sites.in"); // file with the sites
      sfile >> site_type_cache.N; // read the number of sites and number of projections
      site_type_cache.stypes.resize(site_type_cache.N); // resize the array
      for (int i=1;i<=site_type_cache.N;i++)  {
        sfile >> nm ; // read this spin
	site_type_cache.stypes.at(i-1) = nm ;// store
        }
      sfile.close() ;
      site_type_cache.loaded = true; // next time do not read
      };
    int out = site_type_cache.stypes.at(index); // get the value
    cout << index << " site is of type  " << out << endl ;
    return out ;
}





