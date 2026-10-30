import sys, time; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.krylovspectrum import krylov_spectrum_dmrg
from dmrgpy.pyitensor.clik import matrix_elements
L, N = int(sys.argv[1]), int(sys.argv[2])
bc, lc = bh(L, N)
dc = lc.doubled; S = dc._session
t = time.time()
ks = krylov_spectrum_dmrg(dc, lc.L, nk=20, dt=0.5, alpha=0.2, starts=lc._seeds(*bh_start(L, N)), traceless=True, sector_charge=N)
print("time %.0f s substeps %d" % (time.time() - t, ks.substeps))
nk_ = sum(dc.N[2*i] for i in range(L)); nb_ = sum(dc.N[2*i+1] for i in range(L))
Q = S._mpo((((nk_ - N)*(nk_ - N)) + ((nb_ - N)*(nb_ - N))).to_terms())
Qm = matrix_elements(ks.states, Q)
bc.set_conserved_sector(Nb=N)
lamE = lc.get_spectrum(mode="ED")[0]
for k in range(16):
    c = ks.coefficients[:, k]
    dev = (np.vdot(c, Qm @ c) / np.vdot(c, ks.M @ c)).real
    print("%2d %22s  dev %.1e  dist to ED %.1e" % (k, np.round(ks.eigenvalues[k], 5), dev, np.min(np.abs(lamE - ks.eigenvalues[k]))))
