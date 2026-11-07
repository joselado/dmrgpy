# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The density of states that an impurity sees on a Wilson chain of a
# superconductor, averaged over 16 twisted meshes z. The usual mesh puts
# each level at the mean energy of its interval: the bulk is rho/A_Lambda
# and the levels pile up near the band edge D. The mesh of Zitko and
# Pruschke (mesh="zp") gives rho = 1/(2D) with no factor up to the band
# edge, where the width b of the broadening smears the step.
import numpy as np
from dmrgpy import wilsonchain, ysr
from dmrgpy.ysrtk import zaverage

n, Lambda, Delta, nz, b = 20, 4.0, 1e-3, 16, 0.15 # sites, mesh, gap, twists, width
xi = np.logspace(-4, 0.5, 400) # normal-state energy, w = sqrt(xi^2 + Delta^2)
w = np.sqrt(xi**2 + Delta**2)
dos = {}
for mesh in ["mean", "zp"]:
    A = 0*xi
    for z in zaverage.twists(nz):
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, z=z, mesh=mesh)
        e, weight = zaverage.bdg_poles(eps, t, Delta) # levels of the BdG chain
        w0 = zaverage._default_omega0(e, Delta, t[-1])
        A = A + zaverage.broaden_poles(w, e, weight, b=b, omega0=w0, gap=Delta)[0]
    if mesh == "mean":
        A = A*ysr.lambda_factor(Lambda) # the usual A_Lambda correction
    dos[mesh] = A/nz*xi/w # density per spin in the variable xi
exact = 0.5*(xi < 1.0) # flat band of half-width D=1

for x in [0.01, 0.1, 0.5, 0.7, 0.9, 2.0]:
    i = np.argmin(np.abs(xi - x))
    print("xi = %.3g: mean mesh %.4f, zp mesh %.4f, exact %.1f"
          % (xi[i], dos["mean"][i], dos["zp"][i], exact[i]))

import matplotlib.pyplot as plt
plt.semilogx(xi, exact, "k", linewidth=3, label="continuum")
plt.semilogx(xi, dos["mean"], "--", label="mean mesh, times A_Lambda")
plt.semilogx(xi, dos["zp"], "--", label="zp mesh, no factor")
plt.xlabel("xi/D")
plt.ylabel("density of states per spin")
plt.legend()
plt.show()
