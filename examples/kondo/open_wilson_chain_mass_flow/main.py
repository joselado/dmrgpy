# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Mass flow on a Wilson chain (Bruognolo et al., arXiv:1611.05291).
# A harmonic oscillator is coupled to a sub-Ohmic bath, Gamma ~ omega^s.
# Its susceptibility chi = -d<x>/d eps is 1/(Omega - 2 alpha/s), for any
# temperature. A chain of N sites mimics T ~ Lambda^-N. The standard Wilson
# chain gives a chi that drifts with N, the renormalized chain C2 does not.
import numpy as np
from dmrgpy import wilsonchain

alpha, s, Omega = 0.1, 0.4, 1.0
exact = 1.0/(Omega - 2*alpha/s)
ns = range(1, 21) # chain lengths
chi = {}
for kind in ["standard", "C2"]: # the model is quadratic, chi from the chain
    chi[kind] = []
    for n in ns:
        timp, e, t = wilsonchain.boson_bath_chain(n, alpha=alpha, s=s, kind=kind)
        chi[kind].append(1/(Omega + timp**2*wilsonchain.chain_green(e, t).real))

# the same on the C2 chain with DMRG, chi = -<x>/eps
nd = [2, 4, 6]
chid = []
for n in nd:
    bc = wilsonchain.oscillator_chain(n, alpha=alpha, s=s, Omega=Omega,
                                      eps=0.05, maxnb=8)
    bc.maxm = 40 # bond dimension
    bc.nsweeps = 20 # number of DMRG sweeps
    chid.append(-bc.vev(bc.x_imp).real/0.05)
print("Max difference DMRG-exact =", np.max(np.abs(np.array(chid)-exact)))
assert np.max(np.abs(np.array(chid)-exact)) < 1e-4

import matplotlib.pyplot as plt
plt.axhline(exact, color="black", linewidth=3, label="continuum")
plt.plot(ns, chi["standard"], "o-", label="standard Wilson chain")
plt.plot(ns, chi["C2"], "s--", label="renormalized chain C2")
plt.plot(nd, chid, "x", markersize=12, label="DMRG on C2")
plt.xlabel("chain length N")
plt.ylabel("chi")
plt.legend()
plt.show()
