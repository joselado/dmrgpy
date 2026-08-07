# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain

def gs_energy(mode="ED",gamma=0.3):
    n = 4
    if mode=="DMRGJ":
        # create the fermion chain
        fc = fermionchain.Fermionic_Chain(n,itensor_version="julia_live")
    else:
        fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
    if mode=="ED": mode0="ED"
    else: mode0 = "DMRG"
    mh = np.zeros((n,n),dtype=np.complex128) # TB matrix
    for i in range(n-1):
        mh[i,i+1] = 1.0
        mh[i+1,i] = 1.0
    for i in range(n):
        mh[i,i] = 1j*gamma*np.cos(2*np.pi*np.sqrt(2)*i) # imaginary AAH
    h = 0 # initialize Hamiltonian
    for i in range(n):
        for j in range(n):
            h = h + mh[i,j]*fc.Cdag[i]*fc.C[j]
    for i in range(n-1):
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    e = fc.gs_energy(mode=mode0)
    return e

for mode in ["ED","DMRGJ"]:
    e = gs_energy(mode=mode)
    print(mode,e)

# sweep the imaginary AAH potential strength, comparing ED VS the live
# Julia DMRG backend (itensor_version="julia_live") over the whole range
gammas = np.linspace(0.,0.6,5)
es_ed = [gs_energy(mode="ED",gamma=g) for g in gammas]
es_dmrgj = [gs_energy(mode="DMRGJ",gamma=g) for g in gammas]
for g,eed,ej in zip(gammas,es_ed,es_dmrgj):
    print("gamma =",round(g,3),"ED =",eed,"DMRGJ =",ej)

plt.plot(gammas,[e.real for e in es_ed],marker="o",label="ED")
plt.scatter(gammas,[e.real for e in es_dmrgj],marker="x",c="red",label="DMRG (Julia)")
plt.xlabel("Imaginary AAH potential strength $\\gamma$")
plt.ylabel("Ground state energy (Re)")
plt.legend()
plt.show()








