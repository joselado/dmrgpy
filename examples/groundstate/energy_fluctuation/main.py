# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np # conventional numpy library
import matplotlib.pyplot as plt # library to plot the results

from dmrgpy import spinchain
spins = ["S=1" for i in range(10)] # 2*S+1=2 for S=1/2
sc = spinchain.Spin_Chain(spins) # create spin chain object


h = 0 # initialize
for i in range(len(spins)-1): 
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]


maxms = [1,2,5,10,20,30,40] # bond dimensions to sweep over
es,des = [],[] # storage for the energies and their fluctuations
for maxm in maxms: # loop over bond dimension
    sc.set_hamiltonian(h)
    sc.maxm = maxm # set the bond dimension
    e = sc.gs_energy() # get the ground state energy
    de = sc.gs_energy_fluctuation() # fluctuation
    print("Energy",e,"fluctuation",de,"for bond dimension",maxm)
    es.append(e) ; des.append(de)


# The point of the sweep is that the energy fluctuation
# sqrt(<H^2>-<H>^2) measures how far the variational state still is from
# an exact eigenstate, and collapses as the bond dimension grows -- so
# plot it on a log scale next to the energy it certifies.
fig,ax = plt.subplots()
ax.plot(maxms,es,marker="o",c="red")
ax.set_xlabel("Bond dimension")
ax.set_ylabel("Energy",color="red")
ax.tick_params(axis="y",labelcolor="red")

ax2 = ax.twinx() # the fluctuation spans decades, the energy does not
ax2.semilogy(maxms,des,marker="s",ls="--",c="blue")
ax2.set_ylabel("Energy fluctuation $\\sqrt{\\langle H^2\\rangle-\\langle H\\rangle^2}$",
               color="blue")
ax2.tick_params(axis="y",labelcolor="blue")

fig.tight_layout()
plt.show()










