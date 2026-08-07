# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

def get_traces(n):
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0 # generate a Heisenberg Hamiltonian
    for i in range(n-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]

    h = h*h + 1j

    tr_ed = sc.trace(h,mode="ED")
    tr_mps = sc.trace(h,mode="DMRG")
    itr_ed = sc.inverse_trace(h,mode="ED")
    itr_mps = sc.inverse_trace(h,mode="DMRG")
    return tr_ed,tr_mps,itr_ed,itr_mps

ns = [2,4,6] # sweep the chain length (ED cost grows exponentially, keep small)
tr_eds,tr_mpss,itr_eds,itr_mpss = [],[],[],[]
for n in ns:
    tr_ed,tr_mps,itr_ed,itr_mps = get_traces(n)
    print("n =",n)
    print("Trace ED",tr_ed)
    print("Trace MPS",tr_mps)
    print("Inverse trace ED",itr_ed)
    print("Inverse trace MPS",itr_mps)
    tr_eds.append(tr_ed)
    tr_mpss.append(tr_mps)
    itr_eds.append(itr_ed)
    itr_mpss.append(itr_mps)

import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(10,4))
axes[0].plot(ns,[abs(x) for x in tr_eds],marker="o",label="ED")
axes[0].plot(ns,[abs(x) for x in tr_mpss],marker="s",linestyle="--",label="MPS")
axes[0].set_xlabel("Chain length n")
axes[0].set_ylabel("|Trace(H^2 + i)|")
axes[0].set_yscale("log")
axes[0].legend()

axes[1].plot(ns,[abs(x) for x in itr_eds],marker="o",label="ED")
axes[1].plot(ns,[abs(x) for x in itr_mpss],marker="s",linestyle="--",label="MPS")
axes[1].set_xlabel("Chain length n")
axes[1].set_ylabel("|Inverse trace(H^2 + i)|")
axes[1].set_yscale("log")
axes[1].legend()

plt.tight_layout()
plt.show()
