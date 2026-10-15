# Self-consistent mean-field decoupling of the spin-1/2 Heisenberg chain, as
# the fraction p of the exchange kept exactly goes from 0 (pure mean field)
# to 1 (the exact chain), started from a Neel state.
import numpy as np
from dmrgpy import spinchain, meanfield

n = 24
ps = [0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]

sc = spinchain.Spin_Chain(["S=1/2"]*n)
sc.maxm, sc.nsweeps = 60, 12
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
e_exact = sc.gs_energy()

m0 = [[0, 0, 0.5*(-1)**i] for i in range(n)]  # staggered initial guess
mz, e_var, e_mf = [], [], []
for p in ps:
    mf = meanfield.spinchain_meanfield(sc, p=p, m0=m0, mix=0.5,
                                       maxerror=1e-6, maxite=500)
    mz.append(np.real(mf.get_magnetization())[2])  # <Sz_i>
    e_var.append(mf.vev(h).real)  # <H> of the original chain
    e_mf.append(mf.gs_energy())   # energy of the mean-field chain itself

np.savez("data.npz", n=n, p=ps, mz=mz, e_var=e_var, e_mf=e_mf,
         e_exact=e_exact)
