# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Levels of smallest real part of a non-Hermitian chain: the XX chain with
# a staggered imaginary field i g (-1)^i Sz_i, which is PT symmetric. Four
# levels share the smallest real part, a real one twice and a
# complex-conjugate pair. get_nh_levels() finds them by imaginary-time
# evolution with H itself. We compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
g = 0.3 # strength of the imaginary field
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
for i in range(n):
    h = h + 1j*g*(-1)**i*sc.Sz[i]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.cutoff = 1e-10 # truncation of the states

es, res = sc.get_nh_levels() # levels and their relative residuals
eed = np.linalg.eigvals(np.array(sc.get_ED_obj().get_hamiltonian().todense()))
low = es[es.real < eed.real.min() + 1e-3] # the levels of smallest real part
print("Levels of smallest real part:", np.round(low, 6))
diff = max(np.min(np.abs(eed - e)) for e in low)
print("Max distance to the exact spectrum =", diff)
assert diff < 1e-5

import matplotlib.pyplot as plt
plt.plot(eed.real, eed.imag, "o", markersize=6, label="ED")
plt.plot(es.real, es.imag, "x", markersize=14, label="get_nh_levels")
plt.xlabel("Re E")
plt.ylabel("Im E")
plt.legend()
plt.show()
