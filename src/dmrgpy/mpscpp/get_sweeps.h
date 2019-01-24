auto get_sweeps() {
//    auto input = InputGroup("sweeps.in","sweeps");
//    auto N = input.getInt("nsweeps",3); // number of sweeps
//    auto maxm = input.getInt("maxm",30); // bond dimension
//    auto cutoff = input.getReal("cutoff",1E-07); // cutoff
    auto maxm = get_int_value("maxm"); // bond dimension
    auto N = get_int_value("nsweeps"); // bond dimension
    auto cutoff = get_float_value("cutoff"); // bond dimension
    auto sweeps = Sweeps(N); //number of sweeps
    sweeps.maxm() = maxm;
    sweeps.cutoff() = cutoff;
    return sweeps ;
}

