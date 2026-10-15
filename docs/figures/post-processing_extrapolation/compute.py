# Energy per site of open spin-1/2 Heisenberg chains against 1/L, extrapolated
# to L -> infinity with extrapolate.size_extrapolation.
import numpy as np
from dmrgpy import spinchain, extrapolate

Ls = np.array([16, 24, 32, 48, 64, 96, 128])

energy = []
for L in Ls:
    sc = spinchain.Spin_Chain(["S=1/2"]*L)
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                           + sc.Sz[i]*sc.Sz[i+1] for i in range(L-1)))
    sc.maxm, sc.nsweeps = 100, 16
    energy.append(sc.gs_energy())
energy = np.array(energy)

# E/L = e_inf + a/L (n=1), plus b/L^2 (n=2)
e_inf_n1 = extrapolate.size_extrapolation(Ls, energy/Ls, n=1)
e_inf_n2 = extrapolate.size_extrapolation(Ls, energy/Ls, n=2)
e_exact = 0.25 - np.log(2)  # Bethe ansatz

np.savez("data.npz", L=Ls, energy=energy, e_inf_n1=e_inf_n1,
         e_inf_n2=e_inf_n2, e_exact=e_exact)
