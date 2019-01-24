auto get_sweeps() {
    auto input = InputGroup("sweeps.in","sweeps");
    auto N = input.getInt("nsweeps",3); // number of sweeps
    auto maxm = input.getInt("maxm",30); // number of sweeps
    auto cutoff = input.getReal("cutoff",1E-07); // cutoff
    auto sweeps = Sweeps(N); //number of sweeps
    sweeps.maxm() = maxm;
    sweeps.cutoff() = cutoff;
    return sweeps ;
}

