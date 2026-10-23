# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Infinite chain of spinless fermions with two sites A,B per cell,
# H = t1 (A+_n B_n + h.c.) + t2 (B+_n A_{n+1} + h.c.) + t3 (A+_n A_{n+1} + h.c.).
# The t3 term jumps over the site B, so it needs a Jordan-Wigner string.
# The model is quadratic, so the exact energy per site is the sum of the
# negative bands. We compare it with iDMRG.
import numpy as np
from dmrgpy import infinitechain

t1, t2 = 1.0, 0.4

def exact(t3):
    k = np.linspace(-np.pi, np.pi, 2000, endpoint=False)
    e = 0.
    for ki in k: # Bloch Hamiltonian at each k
        hk = [[2*t3*np.cos(ki), t1 + t2*np.exp(-1j*ki)], [0, 0]]
        hk[1][0] = np.conj(hk[0][1])
        w = np.linalg.eigvalsh(np.array(hk))
        e = e + np.sum(w[w < 0])
    return e/len(k)/2 # per site

def idmrg(t3):
    ic = infinitechain.Infinite_Many_Body_Chain([0, 0], itensor_version="python")
    ic.gs_method = "idmrg"
    ic.maxm = 12 # bond dimension
    ic.maxiter = 40 # number of iterations
    C = [ic.get_operator("C", i, "C") for i in range(2)] # in the cell
    Cd = [ic.get_operator("Cdag", i, "C") for i in range(2)]
    CR = [ic.get_operator("C", i, "R") for i in range(2)] # in the next cell
    CdR = [ic.get_operator("Cdag", i, "R") for i in range(2)]
    h = t1*(Cd[0]*C[1] + Cd[1]*C[0])
    h = h + t2*(Cd[1]*CR[0] + CdR[0]*C[1])
    h = h + t3*(Cd[0]*CR[0] + CdR[0]*C[0])
    ic.set_hamiltonian(h)
    return ic.gs_energy()

t3s = np.linspace(0.0, 0.3, 7)
e0 = np.array([exact(t3) for t3 in t3s])
e1 = np.array([idmrg(t3) for t3 in t3s])
print("Max difference =", np.max(np.abs(e0 - e1)))

import matplotlib.pyplot as plt
plt.plot(t3s, e0, linewidth=4, label="exact")
plt.plot(t3s, e1, "--", linewidth=2, label="iDMRG")
plt.xlabel("t3")
plt.ylabel("Energy per site")
plt.legend()
plt.show()
