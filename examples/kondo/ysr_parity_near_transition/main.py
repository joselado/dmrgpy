# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# STM spectrum of a spin 1/2 on a superconductor near the singlet-doublet
# transition, where the two fermion parities are almost degenerate. A
# plain DMRG solve can stop at a mixture of the two parities, which gives
# a false peak at zero bias and moves the YSR peak. Solving inside the
# parity sector (parity_ground_state) gives one YSR peak at +-E_YSR.
import numpy as np
from dmrgpy import ysr
from dmrgpy.ysrtk import paritystate

Lambda, Delta, U = 4.0, 0.1, 0.3 # Wilson mesh, gap, direct tunneling
fc = ysr.kondo_chain(8, J=0.6*ysr.lambda_factor(Lambda), Delta=Delta,
                     Lambda=Lambda) # spin 1/2 + 8 Wilson chain sites
fc.maxm = 60 # bond dimension
fc.nsweeps = 15 # number of DMRG sweeps

E_ysr = ysr.subgap_energies(fc, mode="ED")["E_ysr"] # exact YSR energy
print("E_YSR/Delta from ED =", E_ysr/Delta)
es = np.linspace(-0.004, 0.004, 401) # bias eV
(e1, d1) = ysr.get_stm_spectrum(fc, es, delta=1e-4, U=U) # plain solve
gs, info = paritystate.parity_ground_state(fc) # solve in the parity sector
print("ground state 2Sz =", info["Sz_ground"], "parity", info["parity"])
(e2, d2) = ysr.get_stm_spectrum(gs, es, delta=1e-4, U=U)
for name, d in (("plain", d1), ("parity", d2)):
    print(name, ": peak at eV/Delta =", es[es > 0][np.argmax(d[es > 0])]/Delta,
          " dI/dV at zero bias =", d[len(es)//2])
assert abs(es[es > 0][np.argmax(d2[es > 0])] - abs(E_ysr)) < 2e-5

import matplotlib.pyplot as plt
plt.plot(e2/Delta, d2, linewidth=4, label="parity sector")
plt.plot(e1/Delta, d1, "--", linewidth=2, label="plain DMRG")
plt.axvline(abs(E_ysr)/Delta, color="gray", linestyle=":")
plt.axvline(-abs(E_ysr)/Delta, color="gray", linestyle=":")
plt.xlabel("bias eV/Delta")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
