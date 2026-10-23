# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Apply an operator to every unit cell of an infinite S=1/2 Heisenberg chain.
# We apply the Pauli matrix X = 2*Sx on each site. This operator is
# unitary and flips Sz, so <Sz> changes sign, the correlator
# <Sz_0 Sz_r> stays the same, and the norm eta stays equal to 1.
import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg"
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9 # bond dimension, iterations, tolerance
ic.gs_energy()

flipped = ic.apply_operator(2*ic.SxC[0]) # new chain, X applied to every site
print("Norm eta =", flipped.eta)
print("<Sz> before and after =", ic.vev("Sz", 0), flipped.vev("Sz", 0))

rs = range(1, 6) # distances
c0 = [ic.correlator("Sz", 0, "Sz", r) for r in rs]
c1 = [flipped.correlator("Sz", 0, "Sz", r) for r in rs]
print("Max difference =", np.max(np.abs(np.array(c0)-c1)))
assert np.max(np.abs(np.array(c0)-c1)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(rs, c0, "o-", linewidth=4, label="before")
plt.plot(rs, c1, "x--", linewidth=2, label="after X")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.legend()
plt.show()
