# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Apply the Pauli matrix X on every site of the VUMPS ground state of the
# infinite transverse field Ising chain, H = -4 Sx_i Sx_{i+1} - 2 g Sz_i.
# X is unitary and flips Sz, so <Sz_0 Sz_r> must not change.
# We apply it with ITensor v3 and with the python backend and compare.
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

def chain(version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.gs_method = "vumps"
    ic.maxm = 2 # bond dimension
    ic.etol = 1e-11
    ic.set_hamiltonian(-4*ic.SxC[0]*ic.SxR[0] - 3*ic.SzC[0])
    ic.gs_energy()
    return ic

X = np.array([[0, 1], [1, 0]], dtype=complex)
rs = list(range(1, 6)) # distances
icp, ic3 = chain("python"), chain(3)

# python backend, the operator is an MPO with one tensor
s = icp._vumps_result.sites_uc.si(1)
W = ITensor((Index(1, tags="Link"), s, s.prime(1), Index(1, tags="Link")), X.reshape(1, 2, 2, 1))
flipped = vumps.apply_mpo(icp._vumps_result, [W], cutoff=1e-12, maxdim=None)
cp = [vumps.two_point_correlator(flipped, "Sz", 0, "Sz", r).real for r in rs]

# ITensor v3, apply the same operator and load the new state
D, d, AL, AR, C, AC, eta = ic3._session3.vumps_apply_mpo([X.flatten().tolist()], [1], [1], 1e-12, 0)
ic3._session3.vumps_load_uniform_state(D, d, AL.flatten().tolist(), AR.flatten().tolist(), C.flatten().tolist())
c3 = [ic3.correlator("Sz", 0, "Sz", r).real for r in rs]
print("Max difference v3-python =", np.max(np.abs(np.array(c3)-cp)))

import matplotlib.pyplot as plt
plt.plot(rs, c3, "o-", linewidth=4, label="ITensor v3")
plt.plot(rs, cp, "x--", linewidth=2, label="python")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r> after X")
plt.legend()
plt.show()
