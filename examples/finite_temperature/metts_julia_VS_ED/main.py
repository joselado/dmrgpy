# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Thermal average <Sz_0> of the Heisenberg chain in a field, at
# temperature T, <A> = Tr(A exp(-H/T))/Z. METTS samples it with
# imaginary-time evolution of product states (Julia backend). We
# compare with exact diagonalization (ED) at each temperature.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_julia() # live Julia backend
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n): h = h + 0.3*sc.Sz[i] # magnetic field
sc.set_hamiltonian(h)

Ts = [0.5, 1.0, 2.0, 3.0] # temperatures
ed = [sc.vev(sc.Sz[0], mode="ED", T=T).real for T in Ts]
metts = []
for T in Ts:
    (m, err) = sc.metts_vev([sc.Sz[0]], T, nsamples=300, nwarmup=30,
            dbeta_half_step=0.08, seed=42)[0]
    metts.append(m.real)
print("Max difference METTS-ED =", np.max(np.abs(np.array(metts)-ed)))

import matplotlib.pyplot as plt
plt.plot(Ts, ed, linewidth=4, label="ED")
plt.plot(Ts, metts, "--o", linewidth=2, label="METTS")
plt.xlabel("temperature T")
plt.ylabel("<Sz_0>")
plt.legend()
plt.show()
