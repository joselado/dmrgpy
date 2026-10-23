# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time correlator S(x,t) = <C^dagger_x(t) C_0> of an infinite dimerized
# chain of spinless fermions (hoppings t1, t2, t3). The fermionic operators
# carry a Jordan-Wigner string. We use the infinite-boundary window method
# and compare with the exact free-fermion result.
import numpy as np, scipy.linalg
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_window
t1, t2, t3 = 1.0, 0.4, 0.1
ic = infinitechain.Infinite_Many_Body_Chain([0]*2, itensor_version="python")
ic.gs_method = "idmrg"
ic.maxm, ic.maxiter = 20, 40 # bond dimension, iterations
C, Cd = [[ic.get_operator(o, i, "C") for i in range(2)] for o in ("C", "Cdag")]
CR, CdR = [[ic.get_operator(o, i, "R") for i in range(2)] for o in ("C", "Cdag")]
ic.set_hamiltonian(t1*(Cd[0]*C[1] + Cd[1]*C[0]) + t2*(Cd[1]*CR[0] + CdR[0]*C[1])
                   + t3*(Cd[0]*CR[0] + CdR[0]*C[0]))
ic.gs_energy()
dt, nt = 0.1, 16 # time step and number of steps
xs = list(range(-5, 6)) # distances
ts, xarr, S = idmrg_window.dynamical_correlator_td(ic._result, 8, "Cdag", "C",
        dt=dt, nt=nt, cutoff=1e-10, maxdim=60, x_values=xs, connected=False)
# exact answer, a ring of L cells with the negative levels filled
L = 100
h = np.zeros((2*L, 2*L))
for n in range(L):
    a, b, a2 = 2*n, 2*n+1, (2*n+2) % (2*L)
    for (p, q, t) in ((a, b, t1), (b, a2, t2), (a, a2, t3)):
        h[p,q] += t ; h[q,p] += t
w, v = np.linalg.eigh(h)
v = v[:, w < 0] # filled levels
P = (v.conj() @ v.T).real # <c^dagger_l c_m>
i0 = L # a site of the A sublattice, in the bulk
Se = np.zeros_like(S)
for it, t in enumerate(ts):
    U = scipy.linalg.expm(1j*h*t)
    Se[it] = [U[i0+int(x), :] @ P[:, i0] for x in xarr]
print("Max difference =", np.max(np.abs(S-Se)))

import matplotlib.pyplot as plt
for x in (0, 1, 3):
    ix = list(xarr).index(x)
    plt.plot(ts, Se[:, ix].real, linewidth=4, label="exact, x=%d" % x)
    plt.plot(ts, S[:, ix].real, "k--")
plt.xlabel("time t")
plt.ylabel("Re S(x,t)")
plt.legend()
plt.show()
