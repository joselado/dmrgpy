# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h1,h2 = 0,0

h1 = sc.SS(0,1) + sc.SS(1,2) + sc.SS(2,3) + sc.SS(3,0)
h2 = sc.SS(0,2) + sc.SS(1,3) 

j2s = np.linspace(0.,2.,10)
wfs = []
for j2 in j2s:
    h = h1 + j2*h2
    sc.set_hamiltonian(h)
    wf = sc.get_gs(mode="DMRG") # get ground state
    wfs.append(wf) # store

overlaps = []
for i in range(1,len(j2s)):
    ov = np.abs(wfs[i-1].dot(wfs[i]))
    print(ov)
    overlaps.append(ov)

# plot the overlap between consecutive ground states as J2 (the diagonal
# frustrating coupling) is swept, which drops near a level crossing/phase
# transition of the frustrated square plaquette
j2mid = (j2s[1:]+j2s[:-1])/2.
plt.plot(j2mid,overlaps,marker="o")
plt.xlabel("J2")
plt.ylabel("|<GS(J2_i)|GS(J2_i+1)>|")
plt.show()









