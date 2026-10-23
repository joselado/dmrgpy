# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Applying an operator, written as an MPO, to an infinite MPS found with VUMPS.
# We take the ground state of the transverse field Ising chain,
# H = -4 sum_i Sx_i Sx_{i+1} - 2 g sum_i Sz_i, and apply the Pauli matrix X
# on every site. This flips <Sz> but leaves <Sz_0 Sz_r> unchanged.
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
ic.gs_method = "vumps"
ic.maxm = 2 # bond dimension
ic.etol = 1e-11 # energy tolerance
g = 1.5 # transverse field
ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
ic.gs_energy()
psi = ic._vumps_result # the infinite MPS

# the Pauli matrix X as an MPO with bond dimension 1
sites = psi.sites_uc
d = sites.dim(1)
s = sites.si(1)
X = 2*sites.site_type(1).matrix("Sx")
links = (Index(1, tags="Link"), s, s.prime(1), Index(1, tags="Link"))
W = [ITensor(links, X.reshape(1, d, d, 1))]
flipped = vumps.apply_mpo(psi, W, cutoff=1e-12, maxdim=None)

rs = range(1, 6) # distances
c0 = [vumps.two_point_correlator(psi, "Sz", 0, "Sz", r).real for r in rs]
c1 = [vumps.two_point_correlator(flipped, "Sz", 0, "Sz", r).real for r in rs]
sz0 = vumps.onsite_expectation(psi, "Sz", 0)
sz1 = vumps.onsite_expectation(flipped, "Sz", 0)
print("<Sz> before and after =", sz0, sz1)
assert abs(sz1 + sz0) < 1e-6

import matplotlib.pyplot as plt
plt.plot(rs, c0, "o-", linewidth=4, label="before X")
plt.plot(rs, c1, "x--", linewidth=2, label="after X")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.legend()
plt.show()
