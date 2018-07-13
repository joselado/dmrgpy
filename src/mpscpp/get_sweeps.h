auto get_sweeps() {
    auto input = InputGroup("sweeps.in","sweeps");
    auto N = input.getInt("nsweeps",3); // number of sweeps
    auto maxm = input.getInt("maxm",30); // number of sweeps
    auto sweeps = Sweeps(N); //number of sweeps
    sweeps.maxm() = maxm;
    sweeps.cutoff() = 1E-10;
    return sweeps ;
}

