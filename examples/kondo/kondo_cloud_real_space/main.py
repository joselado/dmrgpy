# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo screening cloud of an Anderson impurity at the end of a
# tight-binding lead (Holzner et al., PRB 80, 205114). The integrated spin
# correlation Sigma(x) = 1 + sum_{i<=x} <S_d.S_i>/<S_d.S_d> says how much of
# the impurity spin is screened by the first x lead sites. In a singlet it
# goes from 1 to exactly 0 at the end of the chain. The screening length
# xi_0.9 is where Sigma drops below 0.1, and it saturates for long chains.
import numpy as np
from dmrgpy import kondocloud

U, Gamma = 1.0, 0.5 # interaction and hybridization, in units of the hopping
TK = kondocloud.kondo_temperature(U, Gamma) # Bethe ansatz estimate
print("T_K =", TK)

import matplotlib.pyplot as plt
for L in [16, 24, 40]: # total number of sites, even
    fc = kondocloud.anderson_lead(L, ed=-U/2, U=U, Gamma=Gamma)
    fc.maxm = 80 # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    fc.cutoff = 1e-10
    r = fc.get_kondo_cloud(a=0.9)
    print("L =", L, "xi_0.9 =", r["xi"], "xi_0.9*T_K =", r["xi"]*TK,
          "Sigma(L-1) =", r["sum_rule"])
    plt.plot(r["x"], r["Sigma"], label="L = %d" % L)

plt.axhline(0.1, color="gray", linestyle="--", label="1 - a, a = 0.9")
plt.yscale("log")
plt.ylim(1e-3, 1.1)
plt.xlabel("distance x from the impurity")
plt.ylabel("Sigma(x)")
plt.legend()
plt.show()
