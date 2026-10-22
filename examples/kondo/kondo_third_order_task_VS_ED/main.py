# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order STM/Kondo spectrum (Ternes, arXiv:1505.04430) at T=0 of an
# S=1 impurity with easy-plane anisotropy in a field, coupled to a chain of
# S=1/2 spins, computed with mode="DMRG", submode="TASK": the exact sums of
# the ED route evaluated on the tangent-space Krylov spectrum of the
# ground-state MPS (arXiv:2510.07279), with no frequency or time grid.
# With cutoff=0 and maxm large enough to hold the whole Hilbert space the
# result is the ED one (asserted below); at maxm=4 the one-site tangent
# space misses part of the intermediate states, and the curve moves.

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

G, MUB = 2.0, 5.7883818066e-5 # g factor, Bohr magneton in eV/T
nb = 6 # spectator S=1/2 spins

def chain(maxm):
    sc = spinchain.Spin_Chain(["S=1"]+["S=1/2"]*nb,itensor_version="python")
    h = 1e-3*sc.Sz[0]*sc.Sz[0] + G*MUB*2.0*sc.Sz[0] # D=1 meV, B=2 T
    for i in range(nb):
        J = 0.6e-3 if i==0 else 1e-3
        h = h + J*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    sc.maxm = maxm ; sc.nsweeps = 20 ; sc.cutoff = 0.
    return sc

eV = np.linspace(-5e-3,5e-3,201)
kw = dict(site=0,Jrho_s=-0.1,U=0.3,T=0.,omega0=20e-3,Gamma0=5e-5,order=3)

sc = chain(64)
t0 = time.time()
_,ed = sc.get_kondo_spectrum(eV,mode="ED",**kw)
t1 = time.time()
_,tk = sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=200,**kw)
t2 = time.time()
err = np.max(np.abs(tk-ed))/np.max(np.abs(ed))
print("saturated: ED %.1f s, TASK %.1f s, largest deviation %.2e"%(t1-t0,t2-t1,err))
assert err < 1e-8

sc4 = chain(4)
_,tk4 = sc4.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=200,**kw)
err4 = np.max(np.abs(tk4-ed))/np.max(np.abs(ed))
print("maxm=4: largest deviation %.2e"%err4)

plt.plot(eV*1e3,ed,c="k",lw=1,label="ED")
plt.plot(eV*1e3,tk,c="C0",ls="--",label="TASK, maxm=64 (saturated)")
plt.plot(eV*1e3,tk4,c="C3",ls=":",label="TASK, maxm=4")
plt.xlabel("eV (meV)") ; plt.ylabel("dI/dV")
plt.legend()
plt.tight_layout()
plt.savefig("kondo_third_order_task_VS_ED.png",dpi=150)
plt.show()
