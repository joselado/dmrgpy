# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# dI/dV of a spin 1/2 on a superconductor. One Wilson mesh gives a comb of
# levels above the gap, so we average 4 twisted meshes z and broaden each
# level with a log-Gaussian of width b, as in NRG. The YSR poles in the gap
# are kept apart. Dashed: J=0, the BCS coherence peaks; with J the YSR state
# takes weight from the gap edge. 10 sites give the YSR energy, but not the
# continuum below about 1.2 Delta (info["reaches_gap"]; 14 sites reach it).
import numpy as np
from dmrgpy import ysr
from dmrgpy.ysrtk import zaverage

J, Lambda, n, Delta, U = 0.4, 4.0, 10, 0.01, 0.3 # continuum J, mesh, sites, gap, direct tunneling
es = Delta*np.linspace(-4, 4, 801) # bias eV
es, dIdV, info = zaverage.zaverage_stm_spectrum(n, es,
            J=J*ysr.lambda_factor(Lambda), Delta=Delta, Lambda=Lambda,
            nz=4, b=0.4, U=U, delta_sub=0.03*Delta, # 4 meshes, width b
            maxm=48, nsweeps=30, nkry=100) # bond dimension, sweeps, Krylov size
for z, sub, ok in zip(info["z"], info["subgap"], info["reaches_gap"]):
    big = np.abs(sub["weight"]) > 1e-6
    print("z =", z, "YSR poles at eV/Delta =",
          np.unique(np.round(sub["omega"][big]/Delta, 4)), "reaches gap:", ok)

# reference: the spin decoupled (J=0), same meshes, from the BdG levels
free = (U**2/2 + 3/8)*zaverage.bdg_chain_dos(es, n, Delta, Lambda=Lambda,
                                           nz=4, b=0.4)

import matplotlib.pyplot as plt
plt.plot(es/Delta, info["A_z"][0], label="one mesh, z = 0")
plt.plot(es/Delta, dIdV, linewidth=3, label="average of 4 meshes")
plt.plot(es/Delta, free, "k--", label="J = 0, coherence peaks")
plt.ylim(0, 3*dIdV[-1])
plt.xlabel("bias eV/Delta")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
