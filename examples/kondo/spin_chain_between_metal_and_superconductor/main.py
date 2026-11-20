# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A chain of 4 spins 1/2 between two semi-infinite leads, each a Wilson
# chain coupled by Kondo exchange JK to one end spin. The left lead is a
# normal metal, the right lead a superconductor with gap Delta. We compute
# <S_end . s_0>, the binding of each end spin to the first site of its
# lead. The gap removes the states that screen the spin, so the binding
# of the right spin goes down with Delta, while the left spin does not move.
import numpy as np
from dmrgpy import ysr

Deltas = np.linspace(0.0, 0.5, 6) # gap of the right lead
left, right = [], []
for Delta in Deltas:
    fc = ysr.spin_chain_on_metal(4, 8, site=[0, 3], Js=0.3, JK=0.5,
                                 Delta=[0.0, Delta], Lambda=4.0)
    fc.maxm = 64 # bond dimension
    fc.nsweeps = 16 # number of DMRG sweeps
    for k, out in ((0, left), (1, right)): # the left and the right lead
        p, w = fc.impurities[k], fc.bath_sites[k][0] # end spin, first lead site
        o = fc.Sx[p]*fc.Sx[w] + fc.Sy[p]*fc.Sy[w] + fc.Sz[p]*fc.Sz[w]
        out.append(fc.vev(o).real) # <S_end . s_0>
    print("Delta =", round(Delta, 2), " left =", round(left[-1], 4),
          " right =", round(right[-1], 4))

import matplotlib.pyplot as plt
plt.plot(Deltas, left, "o-", linewidth=3, label="left spin, normal metal")
plt.plot(Deltas, right, "s--", linewidth=2, label="right spin, superconductor")
plt.xlabel("gap of the right lead, Delta")
plt.ylabel("<S_end . s_0>")
plt.legend()
plt.show()
