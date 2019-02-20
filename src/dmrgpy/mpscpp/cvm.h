MPO I() {

    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N; ++i){

        ampo += "Id", i;

    }

    return (1./N) * MPO(ampo);

}

MPS conjugate_(MPS psi){

    auto tmp = psi;

    for(int i = 1; i <= N; ++i){

        tmp.Aref(i).conj();

    }

    return tmp;
}

MPS conjugate_gradient_squared(MPO A, MPS b, int iter) {

    auto x = b;
    MPS r_old = sum(b, -1 * applyMPO(A, x, args));
    MPS r_ = r_old;
    MPS p = r_old;
    MPS u = r_old;
    MPS q;
    MPS r_new;
    MPS Ap;
    std::complex<double> alpha;
    std::complex<double> beta;
    MPS u_q;

    // Print(overlapC(r_old, r_)); // This must not be zero

    for(int i = 0; i < iter; ++i){

        Ap = applyMPO(A, p, args);
        alpha = overlapC(conjugate_(r_old), r_) / overlapC(conjugate_(Ap), r_);
        q = sum(u, -alpha * Ap, args);
        u_q = sum(u, q, args);
        x = sum(x, alpha * u_q, args);
        r_new = sum(r_old, -alpha * applyMPO(A, u_q, args), args);
        beta = overlapC(conjugate_(r_new), r_) / overlapC(conjugate_(r_old), r_);
        u = sum(r_new, beta * q, args);
        p = sum(u, beta * sum(q, beta * p, args), args);
        r_old = r_new;

    }

    std::cout << "\nResidue = " << norm(r_old) << std::endl; // This checks convergence

    return x;
}

double spectral_function(MPS psi, MPO H, MPO S1, MPO S2, double omega, double eta, double energy, int iter, int i, int j) {

    const std::complex<double> z(omega + energy, eta);
    auto A = sum(z * I(), -1. * H, args);
    auto b =  applyMPO(S2j, psi, args);
    auto x = conjugate_gradient_squared(A, b, iter);

    std::complex<double> G = overlapC(psi, S1, x);

    return -G.imag() / M_PI;
};
