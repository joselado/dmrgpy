# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Tangent-space Krylov (TaSK, arXiv:2510.07279) dynamical correlator of the
# spin-1/2 Heisenberg chain, against the exact Lehmann sum. With maxm large
# enough to hold the whole Hilbert space the tangent space is the full
# complement of the ground state and TaSK is exact (asserted below); with a
# small maxm the poles move, and the two-site variance of each pole, drawn
# as an error bar on its energy, says by how much.

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

n = 10
es = np.linspace(0., 3.5, 400)
delta = 0.05

def chain(maxm):
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    # cutoff=0 keeps every Schmidt state, so maxm=32 really holds the whole
    # Hilbert space of 10 sites
    sc.maxm = maxm ; sc.nsweeps = 20 ; sc.cutoff = 0.
    return sc

fig,axs = plt.subplots(2,1,figsize=(6,6),sharex=True)
for ax,maxm in zip(axs,[32,6]):
    sc = chain(maxm)
    pair = (sc.Sz[n//2],sc.Sz[n//2])
    (x,y,info) = sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,
            delta=delta,nkry=200,return_poles=True)
    (x,yed) = sc.get_dynamical_correlator(mode="ED",submode="ED",name=pair,
            es=es,delta=delta)
    err = np.max(np.abs(y-yed))/np.max(np.abs(yed))
    print("maxm =",maxm,"tangent-space dimension =",info["dimension"],
          "largest deviation from ED (relative to peak) =",err)
    if maxm==32: assert err < 1e-8 # saturated: exact
    ax.plot(x,yed.real,c="k",lw=1,label="ED")
    ax.plot(x,y.real,c="C0",ls="--",label="TaSK, maxm=%d"%maxm)
    w,wt,var = info["omega"],info["weight"].real,info["variance"]
    keep = wt > 1e-3*np.max(wt)
    ax2 = ax.twinx()
    ax2.errorbar(w[keep],wt[keep],xerr=np.sqrt(np.abs(var[keep])),fmt="o",
            c="C3",ms=3,label="poles, error bar = sqrt(variance)")
    ax2.set_ylabel("pole weight")
    ax.set_ylabel(r"$S^{zz}(\omega)$")
    ax.legend(loc="upper left") ; ax2.legend(loc="upper right")
axs[1].set_xlabel(r"$\omega/J$")
plt.tight_layout()
plt.savefig("task_VS_ED.png",dpi=150)
plt.show()
