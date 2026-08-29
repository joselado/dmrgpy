# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import infinitechain
import numpy as np
# a one-site unit cell = a uniform, infinite Heisenberg chain
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg" # growing algorithm ("vumps" is the default)
# C = the central cell, R = the next cell to the right
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h) # set the Hamiltonian of the infinite chain
ic.maxm = 30 # bond dimension
print("Energy per site",ic.gs_energy())
print("Exact (Bethe ansatz)",0.25-np.log(2))
rs = range(1,11) # distances in units of sites
cs = [ic.correlator("Sz",0,"Sz",r).real for r in rs] # Sz-Sz correlator
print("<Sz(0)Sz(r)>",cs)


import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 14})

plt.plot(list(rs),cs,marker="o",c="blue")
plt.axhline(0,c="black",ls="--",lw=0.8)
plt.ylabel("$\\langle S_z(0) S_z(r)\\rangle$")
plt.xlabel("Distance $r$")

plt.tight_layout()
plt.show()
