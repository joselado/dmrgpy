# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin structure factor S(q,w) of the Heisenberg chain, from the states of
# the Sz=+1 sector (submode SECTOR). S+ raises Sz by one, so all the weight
# lives in that sector, and the same states give every site. We add the
# sites with a sine transform (open chain). The weight sits above the
# des Cloizeaux-Pearson edge (pi/2)|sin q| and is largest near q=pi.
import numpy as np
from dmrgpy import spinchain, sectordc
from dmrgpy.operatornames import name2MO

L = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(L)])
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 60 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)

# matrix elements <n|S+_i|0> of every excited state n, for every site i
mel = []
for i in range(L):
    e, w, info = sectordc.sector_poles(sc, name=(Sm[0], Sp[i]), nex=30, quiet=True)
    mel.append(info["matrix_elements"])
mel = np.array(mel).T # (state, site)

qs = np.pi*np.arange(1, L+1)/(L+1) # momenta of the open chain
sine = np.sqrt(2/(L+1))*np.sin(np.outer(qs, np.arange(1, L+1)))
mq = mel@sine.T # matrix elements for every q
es = np.linspace(0, 4, 400) # energies
delta = 0.12 # broadening of the peaks
S = np.zeros((len(qs), len(es)))
for iq in range(len(qs)):
    for k in range(len(e)):
        S[iq] += abs(mq[k, iq])**2*delta/np.pi/((es-e[k])**2+delta**2)
assert np.argmax(S.max(axis=1)) >= L-3 # strongest response near q=pi

import matplotlib.pyplot as plt
plt.pcolormesh(qs/np.pi, es, S.T, shading="auto")
plt.plot(qs/np.pi, 0.5*np.pi*np.sin(qs), "w--", label="lower edge")
plt.xlabel("q/pi")
plt.ylabel("energy w")
plt.legend()
plt.show()
