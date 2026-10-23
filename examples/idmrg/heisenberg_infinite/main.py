# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Infinite S=1/2 Heisenberg chain with infinite DMRG (iDMRG).
# The energy per site is exactly 1/4 - log(2) (Bethe ansatz). We also
# compute the correlator <Sz_0 Sz_r>, which alternates in sign and
# decays slowly with the distance r.
import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"]) # one site per unit cell
ic.gs_method = "idmrg"

# H = sum_i S_i . S_{i+1}, C is the central cell and R the next one
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm = 30 # bond dimension
ic.maxiter = 40 # number of iDMRG steps

e = ic.gs_energy()
print("Energy per site =", e, " exact =", 0.25 - np.log(2))

rs = range(1, 8) # distances
cs = [ic.correlator("Sz", 0, "Sz", r) for r in rs]

import matplotlib.pyplot as plt
plt.plot(rs, cs, "o-")
plt.axhline(0, color="gray", linewidth=0.5)
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.show()
