# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(2)] # or just a dimer with S=1/2
sc = spinchain.Spin_Chain(spins) # create the spin chain object

H_I = sc.Sz[0]*sc.Sz[1] # this would be an Ising Hamiltonian
H_XY = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] # this would be an XY Hamiltonian
H_Hei = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1] # this would be a Heisenberg Hamiltonian


sc.set_hamiltonian(H_I) # set the Hamiltonian
E0_I = sc.gs_energy(mode="ED") # compute the energy for the Ising model

sc.set_hamiltonian(H_XY) # set the Hamiltonian
E0_XY = sc.gs_energy(mode="ED") # compute the energy for the Ising model

sc.set_hamiltonian(H_Hei) # set the Hamiltonian
E0_Hei = sc.gs_energy(mode="ED") # compute the energy for the Ising model
print("Energy of the Ising model",E0_I)
print("Energy of the XY model",E0_XY)
print("Energy of the Heisenberg model",E0_Hei)

# The three models above are special cases of the one-parameter family
#   H(alpha) = Sz*Sz + alpha*(Sx*Sx + Sy*Sy)
# with alpha=0 -> Ising, alpha=1 -> Heisenberg (XY has no Sz*Sz term, so
# it isn't exactly on this line, but sweeping alpha still shows how the
# ground state energy interpolates between the Ising and Heisenberg
# limits already computed above).
alphas = np.linspace(0.,1.,6)
Es = []
for alpha in alphas:
    H_alpha = sc.Sz[0]*sc.Sz[1] + alpha*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1])
    sc.set_hamiltonian(H_alpha)
    E_alpha = sc.gs_energy(mode="ED")
    print("alpha =",alpha,"  Energy =",E_alpha)
    Es.append(E_alpha)

import matplotlib.pyplot as plt
plt.plot(alphas,Es,marker="o",label="H = Sz.Sz + alpha*(Sx.Sx+Sy.Sy)")
plt.scatter([0.],[E0_I],color="red",zorder=5,label="Ising (alpha=0)")
plt.scatter([1.],[E0_Hei],color="green",zorder=5,label="Heisenberg (alpha=1)")
plt.axhline(E0_XY,color="gray",linestyle="--",label="XY model")
plt.xlabel("alpha")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
