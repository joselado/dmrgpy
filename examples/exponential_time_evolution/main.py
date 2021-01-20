# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 5 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

Mz = sum(sc.Sz) # total magnetization in Z
Mx = sum(sc.Sx) # total magnetization in X

h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
#sc.itensor_version = "julia"

h = h + h.get_dagger()
sc.set_hamiltonian(h) # set Hamiltonian
h = h/2.
e = sc.gs_energy()
h = h - e

def get(ts,mode):
    sc.set_hamiltonian(h) # set Hamiltonian
    if mode=="DMRG": wf = sc.get_gs(mode="DMRG")
    if mode=="ED": wf = sc.get_gs(mode="ED",array_mode=False)
    wf = sc.Sx[0]*wf # apply the lowering operator
    wf = wf.normalize()
    wf0 = wf.copy()
    # The operator in the exponent MUST be Hermitian
    from dmrgpy.timeevolution import evolve_WF
    wfs = evolve_WF(h,wf,ts=ts)
    return [wf0.dot(wfi).real for wfi in wfs] 



ts = np.linspace(0.0,5.0,50)
sc.set_hamiltonian(h) # set Hamiltonian
# do DMRG
us0 = get(ts,"DMRG")
# do ED
us1 = get(ts,"ED")

import matplotlib.pyplot as plt


plt.scatter(ts,us0,label="DMRG",marker="o",s=200,c="blue")
plt.scatter(ts,us1,label="ED",marker="o",s=60,c="red")
plt.legend()
plt.show()


