# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
import numpy as np

from dmrgpy import spinchain
n = 100 # number of sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # initialize
for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1] # Ising coupling
for i in range(n): h = h + 0.5*sc.Sx[i] # transverse field
sc.set_hamiltonian(h) # set the Hamiltonian
sc.maxm = 200 # increase bond dimension for a critical system
wf = sc.get_gs() # compute ground state
c = wf.get_CFT_central_charge() # compute central charge
print("Central charge",c)

# get_CFT_central_charge() fits c/6*log(2L/(pi)*sin(pi*l/L)) + const to the
# bond entanglement entropy profile (see entropy.py's central_charge()) --
# reproduce that same profile here and plot it against the fitted curve,
# rather than reporting only the single fitted scalar.
L = len(wf.MBO.sites)
ls = list(range(1,L))
sr = [wf.get_bond_entropy(i-1,i) for i in ls]

a = 1.0
fit = [c/6.*np.log(2*L/np.pi*np.sin(np.pi*l/L)) for l in ls]
offset = np.mean(np.array(sr)-np.array(fit)) # const term from the fit

import matplotlib.pyplot as plt

fig = plt.figure()
plt.plot(ls,sr,marker="o",ls="none",label="DMRG bond entropy")
plt.plot(ls,np.array(fit)+offset,c="red",
         label="CFT fit, c=%.4f"%c)
plt.xlabel("bond index")
plt.ylabel("Entanglement entropy")
plt.legend()
plt.title("Central charge from CFT scaling of the entropy profile (n=%d)"%n)
plt.savefig("central_charge_critical.png",dpi=150)
print("Plot saved to central_charge_critical.png")
plt.show()
