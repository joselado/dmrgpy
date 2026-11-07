# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Two magnetic atoms on a superconductor with gap Delta. Each spin 1/2 is
# coupled by J S.s to its own BCS bath of L=2 levels (the surrogate model,
# fitted to the hybridization of a flat BCS band). A hopping t_b between
# the two baths hybridizes the two Yu-Shiba-Rusinov (YSR) states, and the
# subgap peak in dI/dV splits in two. A potential scattering V is needed:
# at particle-hole symmetry the two YSR states do not split.
import numpy as np
from dmrgpy.ysrtk import surrogate

Delta = 0.1 # superconducting gap
es = np.linspace(-1.2*Delta, 1.2*Delta, 241) # bias eV
import matplotlib.pyplot as plt
for tb in [0.0, 0.05, 0.1]: # hopping between the two baths
    fc = surrogate.surrogate_kondo_chain(2, J=2.0, V=1.0, Delta=Delta, L=2,
                                         t_b=tb)
    fc.maxm = 64 # bond dimension
    fc.nsweeps = 20 # number of DMRG sweeps
    # exact energies of the two lowest odd states (2Sz=1), for reference
    r = surrogate.sector_energies(fc, mode="ED", sectors=(0, 1), nstates=2)
    eysr = (r["sectors"][1] - r["E0"])/Delta
    print("t_b/Delta =", tb/Delta, " YSR states at E/Delta =", eysr)
    # dI/dV with the tip over the first spin, from DMRG
    es, didv = surrogate.get_stm_spectrum(fc, es, spin=0, delta=0.01*Delta,
                                          U=0.3)
    p = plt.plot(es/Delta, didv, label="t_b/Delta = "+str(tb/Delta))
    for e in eysr:
        plt.axvline(e, color=p[0].get_color(), linestyle=":")
        plt.axvline(-e, color=p[0].get_color(), linestyle=":")

plt.xlabel("bias eV/Delta")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
