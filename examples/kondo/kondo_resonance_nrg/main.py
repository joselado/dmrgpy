# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The Kondo resonance of one spin 1/2 on a metal, by the numerical
# renormalization group (NRG) with the full density matrix.
# STM measures A(w), and t(w) = pi^2 rho J^2 A(w)/2 is the T-matrix:
# at T = 0 it reaches the unitary value 1 at w = 0, and at T above
# the Kondo temperature T_K the resonance melts.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import nrgspectral, nrgthermo, ysr

Lambda = 4.0 # discretization of the Wilson chain
rhoJ = 0.2 # Kondo coupling in the continuum, rho = 1/(2D), D = 1
J = 2*rhoJ*ysr.lambda_factor(Lambda) # the same coupling on the chain
TK = nrgthermo.kondo_temperature(nrgthermo.kondo_nrg(40, J=J, Lambda=Lambda))["T_K"]
print("Wilson T_K =", TK)

es = np.logspace(-8, 0, 60) # positive energies, the spectrum is symmetric
for T in [0.0, TK, 10*TK]:
    (e, A, info) = nrgspectral.stm_spectrum(es, ns=1, n=40, JK=J, Lambda=Lambda,
                                            nz=4, b=0.4, keep=400, T=T)
    rho, Jc = 0.5, 2*rhoJ # density of states and continuum coupling
    t = np.pi**2*rho*Jc**2*A/2
    if T == 0.0:
        print("t(w) at w = 1e-6:", t[np.argmin(abs(es-1e-6))])
        assert abs(t[np.argmin(abs(es-1e-6))] - 1.0) < 0.02 # unitary limit
    plt.semilogx(es/TK, t, label="T = %.3g T_K" % (T/TK))
plt.xlabel("eV/T_K")
plt.ylabel("t(w)")
plt.legend()
plt.show()
