auto get_sweeps() {
//    auto input = InputGroup("sweeps.in","sweeps");
//    auto N = input.getInt("nsweeps",3); // number of sweeps
//    auto maxm = input.getInt("maxm",30); // bond dimension
//    auto cutoff = input.getReal("cutoff",1E-07); // cutoff
    auto maxm = get_int_value("maxm"); // bond dimension
    auto N = get_int_value("nsweeps"); // bond dimension
    auto cutoff = get_float_value("cutoff"); // bond dimension
    auto moise = get_float_value("moise"); // noise
    auto sweeps = Sweeps(N); //number of sweeps
    sweeps.maxm() = maxm;
    sweeps.cutoff() = cutoff;
    sweeps.noise() = moise;
    // noise only in the first half
    for (int i=N/2;i<N;i++) sweeps.setnoise(i,0.0) ;
    return sweeps ;
}

