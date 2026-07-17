# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time quench evolution of an MPS, contrasting the two available
# itensor_version=3 time-evolution engines: TDVP (mpscpp3/TDVP/, the
# default, sc.tevol_method=="TDVP") against the legacy MPO-Taylor method
# kept as a backup (sc.tevol_method="MPO", also the only option for
# itensor_version=2). See timedependent.py's evolve_and_measure_dmrg().
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 6 # number of sites in the chain
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins) # create the chain

print("Default time-evolution method:",sc.tevol_method) # TDVP by default

# create two Hamiltonians
h0 = 0
h1 = 0
for i in range(n):
    h0 = h0+(-1)**i*sc.Sz[i]
for i in range(n-1):
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h0) # set AF hamiltonian
wf = sc.get_gs() # compute ground state
sc.set_hamiltonian(h1) # set the (new) Heisenberg Hamiltonian

op = sc.Sz[0] # operator to compute

nt = 300 # number of time steps
dt = 1e-2 # time step

# TDVP (default)
sc.tevol_method = "TDVP"
(ts_tdvp,sz_tdvp) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)

# legacy MPO-Taylor method, kept as a backup
sc.tevol_method = "MPO"
(ts_mpo,sz_mpo) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)

print("<Sz_0>(t) sample (TDVP):",sz_tdvp[:5].real)
print("<Sz_0>(t) sample (MPO Taylor backup):",sz_mpo[:5].real)
print("Max abs difference between the two methods =",np.max(np.abs(sz_tdvp-sz_mpo)))

# now plot the result
import matplotlib.pyplot as plt
plt.plot(ts_tdvp,sz_tdvp.real,label="TDVP (default)",c="blue")
plt.plot(ts_mpo,sz_mpo.real,label="MPO Taylor (backup)",c="red",linestyle="--")
plt.legend()
plt.xlabel("time")
plt.ylabel(r"$\langle S_0^z\rangle(t)$")
plt.show()
