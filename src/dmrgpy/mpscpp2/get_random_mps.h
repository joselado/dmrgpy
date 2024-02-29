static auto get_random_mps=[]() {
    // read the GS from a file
    auto sites = get_sites(); // Get the different sites
    auto psi = MPS(sites);
 //   psi.orthogonalize(); //
//    normalize(psi); // to orthogonalize, but this throws the Zero norm error
		    // for large systems
    writeToFile("random.mps",psi); // write to file
};
