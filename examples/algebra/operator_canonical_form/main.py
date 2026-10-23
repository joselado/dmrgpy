# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Is an operator Hermitian? dmrgpy can settle it in two ways.
# Symbolic: write H - H^dagger in canonical form (terms sorted by site,
# with the fermionic signs). If everything cancels, H is Hermitian.
# Numerical: apply H - H^dagger to a random state and look at the norm.
# We compare the time of the two for a Heisenberg chain.
import time
import numpy as np
from dmrgpy import spinchain, mpsalgebra
from dmrgpy.multioperator import MultiOperator

ns = [4, 8, 16, 24, 32] # chain lengths
proof = MultiOperator.is_hermitian # the symbolic proof
tsym, tnum = [], []
for n in ns:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    t0 = time.time()
    ra = h.is_hermitian() # symbolic proof
    tsym.append(time.time() - t0)
    MultiOperator.is_hermitian = lambda self: False # switch the proof off
    t0 = time.time()
    rb = mpsalgebra.is_hermitian(sc, h) # numerical test
    tnum.append(time.time() - t0)
    print("n =", n, " symbolic:", ra, " numerical:", rb)
    assert ra and rb
    MultiOperator.is_hermitian = proof # switch the proof on again

import matplotlib.pyplot as plt
plt.plot(ns, np.array(tnum)*1e3, "o-", label="numerical")
plt.plot(ns, np.array(tsym)*1e3, "s-", label="canonical form")
plt.yscale("log")
plt.xlabel("number of sites")
plt.ylabel("time (ms)")
plt.legend()
plt.show()
