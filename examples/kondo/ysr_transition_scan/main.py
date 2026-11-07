# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A spin 1/2 on a superconductor with gap Delta, coupled by J S.s to the
# Wilson chain of the substrate. E_ysr = E_odd - E_even is the energy of
# the Yu-Shiba-Rusinov peak in dI/dV: positive for a screened singlet,
# negative for a free doublet. It changes sign at Delta_c, a few T_K.
# The scan solves each gap starting from the states of the previous one,
# and the bisection finds Delta_c the same way.
import numpy as np
from dmrgpy.ysrtk import subgapscan
from dmrgpy import ysr

J, Lambda, n = 0.4, 4.0, 14 # exchange (continuum), discretization, sites
Jchain = J*ysr.lambda_factor(Lambda) # the same exchange on the chain
Deltas = np.logspace(-3, -1, 9) # superconducting gaps
r = subgapscan.subgap_scan(n, Deltas=Deltas, J=Jchain, Lambda=Lambda,
                           maxm=64, nsweeps=20, # bond dimension, sweeps
                           nsweeps_warm=6) # sweeps from a warm start
es = r["E_ysr"]/Deltas
for Delta, e, gs in zip(Deltas, es, r["ground_state"]):
    print("Delta =", Delta, "(E_odd-E_even)/Delta =", e, gs)

# bisection of the transition between the two gaps where E_ysr changes sign
i = np.argmax(es < 0)
t = subgapscan.find_transition(n, Jchain, Deltas[i-1], Deltas[i],
                               Lambda=Lambda, maxm=64, nsweeps=20,
                               nsweeps_warm=6)
print("Delta_c =", t["Delta_c"])

import matplotlib.pyplot as plt
plt.semilogx(Deltas, es, "o-", label="warm-started scan")
plt.semilogx(t["Delta"], t["x"], "x", label="bisection points")
plt.axvline(t["Delta_c"], color="red", linestyle="--", label="Delta_c")
plt.axhline(0, color="gray", linewidth=0.5)
plt.xlabel("superconducting gap Delta")
plt.ylabel("(E_odd - E_even)/Delta")
plt.legend()
plt.show()
