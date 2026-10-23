# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlator <Sz_0 Sz_r> in the ground state of the infinite transverse
# field Ising chain, H = -J sum sx_i sx_{i+1} - h sum sz_i, with h>J so
# that the chain is gapped and the correlator decays exponentially.
# We compute it with iDMRG on the python and the ITensor v3 backends.
import numpy as np
from dmrgpy import infinitechain

J, hz = 1.0, 2.0 # exchange and transverse field
rs = range(1,8) # distances
cs = {}
for version in ["python", 3]:
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(-4*J*(ic.SxC[0]*ic.SxR[0]) - 2*hz*ic.SzC[0])
    ic.maxm = 16 # bond dimension
    ic.maxiter = 120 # maximum number of iterations
    ic.etol = 1e-12 # energy tolerance
    ic.gs_energy()
    cs[version] = np.array([ic.correlator("Sz",0,"Sz",r).real for r in rs])
diff = np.max(np.abs(cs["python"]-cs[3]))
print("Max difference python-v3 =", diff)
assert diff < 1e-5

import matplotlib.pyplot as plt
plt.plot(rs, np.abs(cs["python"]), linewidth=4, label="python")
plt.plot(rs, np.abs(cs[3]), "--", linewidth=2, label="v3")
plt.yscale("log")
plt.xlabel("distance r")
plt.ylabel("|<Sz_0 Sz_r>|")
plt.legend()
plt.show()
