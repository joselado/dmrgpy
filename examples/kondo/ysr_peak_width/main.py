# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Height of the YSR peaks measured with a superconducting tip of gap Delta_t.
# The tip fills the YSR state with an electron (weight u2) or a hole (v2),
# and the state relaxes to the continuum with the rate Gamma1. At weak
# tunneling the peaks at eV = +-(Delta_t + eps) are single-electron peaks
# with heights u2 : v2. At strong tunneling Andreev reflection dominates
# and the asymmetry inverts (Ruby et al., arXiv:1502.05048).
import numpy as np
from dmrgpy import ysr
from dmrgpy.ysrtk import zaverage, ysrcurrent

Lambda, Delta = 4.0, 0.1 # Wilson mesh and gap of the substrate
fc = ysr.kondo_chain(8, J=0.8*ysr.lambda_factor(Lambda), Delta=Delta,
                     Lambda=Lambda) # spin 1/2 + 8 sites, this reaches the gap
fc.maxm = 60 # bond dimension
fc.nsweeps = 15 # number of DMRG sweeps
poles = zaverage.stm_poles(fc, U=0.3, nkry=100) # poles of the tip vertex
p = ysrcurrent.rate_parameters(poles, Delta, Lambda=Lambda)
eps, u2, v2 = p["eps"], p["u2"], p["v2"] # YSR energy and weights
print("eps/Delta =", eps/Delta, " u2 =", u2, " v2 =", v2)

Dt = Delta # gap of the tip
G1, G2 = 1e-3*Delta, 2.5e-4*Delta # rates that empty and fill the state
GNs = np.logspace(-8, -1, 15) # normal-state conductance, units of 2e^2/h
names = ["alpha+", "alpha-", "beta+", "beta-"]
heights = []
for GN in GNs:
    g = ysrcurrent.gamma_tip_from_GN(GN) # tip coupling 2 pi t^2 nu_t
    h = ysrcurrent.peak_heights(eps, u2, v2, Dt, g, G1, G2)
    heights.append([h[k]["total"] for k in names])
heights = np.array(heights)
r = heights[:, 0]/heights[:, 1]
print("alpha+/alpha- at the lowest and highest GN:", r[0], r[-1])

import matplotlib.pyplot as plt
for k in range(4):
    plt.loglog(GNs, heights[:, k], "o-", label=names[k])
plt.xlabel("normal-state conductance (2e^2/h)")
plt.ylabel("peak height of dI/dV (2e^2/h)")
plt.legend()
plt.show()
