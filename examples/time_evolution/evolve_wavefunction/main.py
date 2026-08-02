# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
#------------------------------------------------------------------
import numpy as np
from dmrgpy import spinchain
from dmrgpy.timeevolution import imaginary_exponential # function to perform t-evol

# Heisenberg Hamiltonian S12
n = 4 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
h = h + h.get_dagger()
sc.set_hamiltonian(h) # set Hamiltonian

# perturb the ground state and use it as the initial wavefunction
wf0 = sc.get_gs() # get the ground state
wf0 = sc.Sx[0]*wf0
wf0 = wf0.normalize()

# evolve the perturbed state e^{iHt}|wf0> at a range of times
ts = np.linspace(0.0,3.0,30)
wfs = imaginary_exponential(h,wf0,ts=ts)

overlaps = np.array([abs(wf0.dot(wf)) for wf in wfs]) # |<wf0|wf(t)>|
sz0 = np.array([wf.dot(sc.Sz[0]*wf).real for wf in wfs]) # <Sz_0>(t)

import matplotlib.pyplot as plt

fig,axes = plt.subplots(1,2,figsize=(10,4))
axes[0].plot(ts,overlaps,marker="o")
axes[0].set_xlabel("Time") ; axes[0].set_ylabel(r"$|\langle\Psi_0|\Psi(t)\rangle|$")
axes[1].plot(ts,sz0,marker="o",color="red")
axes[1].set_xlabel("Time") ; axes[1].set_ylabel(r"$\langle S_0^z\rangle(t)$")
plt.tight_layout()
plt.show()
