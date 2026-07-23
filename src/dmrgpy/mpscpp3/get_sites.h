// Port of mpscpp2/get_sites.h to the ITensor v3 API.
//
// Site-type code convention is unchanged from v2 (see chain_session.h's
// Chain constructor / manybodychain.py's callers): 2=spin-1/2, 0=spinless
// fermion, 1=spinful fermion (Hubbard), 3=spin-1, 4=spin-3/2, 5=spin-2,
// 6=spin-5/2, -2=Z3, -3=Z4. Boson is a *range* of codes, 102..199,
// v3-only (mpscpp2/pyitensor still only understand the single 104 code):
// code = 100+dim, dim being the local Hilbert space dimension (dmrgpy's
// own "maxnb" convention, see bosonchain.py's Bosonic_Chain) so 104 is
// the original fixed 4-level site and e.g. 108 is an 8-level site. This
// generalizes BosonFourSite's fixed dim=4 to an arbitrary MaxOcc=dim-1
// (see extra/bosonfour.h) so DMRG can actually honor a non-default
// maxnb instead of silently truncating to dim 4 regardless of what the
// Python layer requested.
//
// Class-name changes vs v2: SpinlessSite/HubbardSite are v3's own
// deprecated-but-supported aliases for FermionSite/ElectronSite (same
// operator names, see mpscpp3/ITensor/itensor/mps/sites/{fermion,electron}.h);
// SpinTwoSite is v3's own stock class (itensor/mps/sites/spintwo.h) with the
// same operator names/values dmrgpy's v2-only extra/spintwo.h had; the
// IQIndex-based read() path uses Index/dim() instead of IQIndex/.m().
//
// Every site is built with ConserveQNs=false. This isn't a v2/v3 API
// difference to paper over -- it reproduces a real v2 *behavior*: v2's
// Chain::gs_energy() seeds dmrg() from a plain default-constructed
// "MPS(sites_)" with no InitState, which in v2's IQIndex-based MPS ends up
// with ordinary (non-QN) auxiliary link indices, so despite every site
// being built from a QN-carrying IQIndex, the actual DMRG run v2 performs
// is an *unconstrained* search of the full Hilbert space, not one pinned to
// whatever total-Sz/particle-number sector a specific starting product
// state would fix. Verified directly against the compiled v2 extension: a
// Heisenberg chain plus a strong field favors a fully-polarized ground
// state (total Sz = N/2, not 0), and v2's gs_energy() finds exactly that
// energy -- something only possible without a fixed QN sector. v3's own
// "MPS(sites_)" with no InitState instead leaves every site tensor
// completely unallocated (a hard error the moment anything tries to
// contract it, not a silently-wrong-sector result like v2's), so simply
// not passing an InitState isn't an option at all in v3; building sites
// with ConserveQNs=false and starting from MPS(InitState(sites_)) (see
// chain_session.h's default_mps()) is what actually reproduces v2's
// unconstrained-search behavior, at the cost of the QN block-sparsity
// speedup -- a real perf/memory tradeoff of this backend, not a bug.

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
            auto I = Index{}; // generate an identity for this index
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
            else if(nm > 100 && nm < 200) store.set(j,BosonFourSite(I)); // dim baked into serialized Index I
            else Error(tinyformat::format("SpinX cannot read index of size %d",nm));
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
      if (nm==2) sites.set(i,SpinHalfSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=1/2
      else if (nm==0) sites.set(i,SpinlessSite({"SiteNumber=",i,"ConserveQNs=",false})); // use fermions
      else if (nm==1) sites.set(i,HubbardSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spinful fermions
      else if (nm==3) sites.set(i,SpinOneSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i,{"ConserveQNs",false})); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i,{"ConserveQNs",false})); // use spin=5/2
      else if (nm>100 && nm<200) sites.set(i,BosonFourSite(i,{"ConserveQNs",false,"MaxOcc",nm-101})); // Boson, dim=nm-100
      else if (nm==(-2)) sites.set(i,Z3Site({"SiteNumber=",i,"ConserveQNs=",false})); // use Z3
      else if (nm==(-3)) sites.set(i,Z4Site(i,{"ConserveQNs",false})); // use Z4
      else Error(tinyformat::format("SpinX cannot read index of size "));
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
      if (nm==2) sites.set(i,SpinHalfSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=1/2
      else if (nm==0) sites.set(i,SpinlessSite({"SiteNumber=",i,"ConserveQNs=",false})); // use fermions
      else if (nm==1) sites.set(i,HubbardSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spinful fermions
      else if (nm==3) sites.set(i,SpinOneSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i,{"ConserveQNs",false})); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite({"SiteNumber=",i,"ConserveQNs=",false})); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i,{"ConserveQNs",false})); // use spin=5/2
      else if (nm>100 && nm<200) sites.set(i,BosonFourSite(i,{"ConserveQNs",false,"MaxOcc",nm-101})); // Boson, dim=nm-100
      else if (nm==(-2)) sites.set(i,Z3Site({"SiteNumber=",i,"ConserveQNs=",false})); // use Z3
      else if (nm==(-3)) sites.set(i,Z4Site(i,{"ConserveQNs",false})); // use Z4
      else Error(tinyformat::format("SpinX cannot read index of size "));
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
