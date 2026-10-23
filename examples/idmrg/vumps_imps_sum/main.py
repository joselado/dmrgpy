# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Direct sum of two infinite MPS, vumps.imps_sum. We sum the ground state
# of the transverse field Ising chain at g=1.5 with a copy of the g=2.0
# ground state whose norm is 0.9 times smaller. The state with the larger
# norm dominates, so the sum must give back the g=1.5 state. We check it
# with the correlator <Sz_0 Sz_r>, which an STM or neutrons would probe.
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps

def tfim(g): # H = -sigma^x_i sigma^x_i+1 - g sigma^z_i
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 2 # bond dimension
    ic.etol = 1e-11
    ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.gs_energy()
    return ic._vumps_result

a, b = tfim(1.5), tfim(2.0)
f = 0.9 # rescale the norm of the second state
b = vumps.UniformMPS(b.sites_uc, b.n_uc, b.D, b.d_g, b.AL*f, b.AR*f,
        b.C*f, b.AC*f**2, eta=f**2)
s = vumps.imps_sum(a, b, cutoff=1e-12, maxdim=None)

rs = range(1, 5) # distances
ca = [vumps.two_point_correlator(a, "Sz", 0, "Sz", r).real for r in rs]
cs = [vumps.two_point_correlator(s, "Sz", 0, "Sz", r).real for r in rs]
print("Max difference =", np.max(np.abs(np.array(ca)-cs)))

import matplotlib.pyplot as plt
plt.plot(rs, ca, "o-", linewidth=4, label="g=1.5 state")
plt.plot(rs, cs, "x--", linewidth=2, label="sum")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.legend()
plt.show()
