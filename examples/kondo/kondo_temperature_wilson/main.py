# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo temperature of a spin 1/2 on a metal, in Wilson's definition.
# T chi_imp(T) is the impurity part of the susceptibility to a uniform
# field. It goes from 1/4 (free moment) to w T/(4 T_K) (screened moment),
# w = 0.4128, and this slope defines T_K. Two couplings give the same curve
# against T/T_K, and it crosses 0.0701 near T = T_K.
import numpy as np
from dmrgpy import nrgthermo, ysr

Lambda = 4.0 # discretization of the band
A = ysr.lambda_factor(Lambda) # continuum J is J*A on this Wilson chain
curves = []
for rhoJ in [0.1, 0.15]:
    r = nrgthermo.kondo_nrg(40, J=2*rhoJ*A, Lambda=Lambda, keep=400)
    k = nrgthermo.kondo_temperature(r) # T_K from the slope at low T
    T, Tchi = nrgthermo.even_odd_average(r.T, r.Tchi_imp)
    curves.append((rhoJ, T/k["T_K"], Tchi))
    print("rho J =", rhoJ, " T_K =", "%.4e" % k["T_K"],
          " T(0.0701)/T_K =", round(k["T_cross"]/k["T_K"], 3))
    B = 0.01*k["T_K"] # a weak field on the impurity only
    rb = nrgthermo.kondo_nrg(40, J=2*rhoJ*A, Lambda=Lambda, keep=400, B=B)
    kl = nrgthermo.kondo_temperature(rb, definition="local")
    print("   T_K from the local susceptibility =", "%.4e" % kl["T_K"])

import matplotlib.pyplot as plt
for rhoJ, x, y in curves:
    plt.plot(x, y, "o-", label="rho J = %g" % rhoJ)
x = np.logspace(-3, 0, 20)
plt.plot(x, 0.4128*x/4, "--", color="gray", label="w T/(4 T_K)")
plt.axhline(0.0701, color="gray", linewidth=0.5)
plt.axvline(1.0, color="gray", linewidth=0.5)
plt.xscale("log")
plt.xlabel("T/T_K")
plt.ylabel("T chi_imp")
plt.legend()
plt.show()
