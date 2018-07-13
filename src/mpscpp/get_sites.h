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
        for(int j = 1; j <= N; ++j)
            {
            auto I = IQIndex{};
            I.read(s);
            if(I.m() == 3) store.set(j,SpinOneSite(I));
            else if(I.m() == 2) store.set(j,SpinHalfSite(I));
            else if(I.m() == 4) store.set(j,SpinThreeHalfSite(I));
            else if(I.m() == 5) store.set(j,SpinTwoSite(I));
            else if(I.m() == 6) store.set(j,SpinFiveHalfSite(I));
            else Error(format("SpinX cannot read index of size %d",I.m()));
            }
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
      if (nm==2) sites.set(i,SpinHalfSite(i)); // use spin=1/2
      else if (nm==3) sites.set(i,SpinOneSite(i)); // use spin=1
      else if (nm==4) sites.set(i,SpinThreeHalfSite(i)); // use spin=3/2
      else if (nm==5) sites.set(i,SpinTwoSite(i)); // use spin=2
      else if (nm==6) sites.set(i,SpinFiveHalfSite(i)); // use spin=5/2
    } ;
    sfile.close(); // close file

    SiteSet::init(std::move(sites));
    }







auto get_sites() { // function to get the sites
    auto sites = SpinX() ; // read from file
    return sites ;
}
