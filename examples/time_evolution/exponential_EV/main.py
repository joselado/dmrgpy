# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Exponential of an operator acting on a wavefunction, exp(z Mx)|GS>.
# The ground state of H = Mz is fully polarized along z, so that
# <GS|exp(z Mx)|GS> = cosh(z/2)^n for n spins 1/2. We compute it with
# DMRG and with exact diagonalization (ED), and plot its logarithm.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
Mz = sum(sc.Sz) # total magnetization along z
Mx = sum(sc.Sx) # total magnetization along x
sc.set_hamiltonian(Mz) # only a magnetic field
sc.maxm = 20 # bond dimension

def overlap(z, mode):
    wf = sc.get_gs(mode=mode) # polarized ground state
    wf1 = sc.exponential(z*Mx, wf, mode=mode) # exp(z Mx)|GS>
    return wf.dot(wf1).real

zs = np.linspace(-1.0, 1.0, 10)
u_dmrg = np.log([overlap(z, "DMRG") for z in zs])
u_ed = np.log([overlap(z, "ED") for z in zs])
print("Max difference DMRG-ED =", np.max(np.abs(u_dmrg-u_ed)))

import matplotlib.pyplot as plt
plt.plot(zs, u_ed, linewidth=4, label="ED")
plt.plot(zs, u_dmrg, "--", linewidth=2, label="DMRG")
plt.plot(zs, n*np.log(np.cosh(zs/2)), ":", color="black", label="exact")
plt.xlabel("z")
plt.ylabel("log <GS|exp(z Mx)|GS>")
plt.legend()
plt.show()
