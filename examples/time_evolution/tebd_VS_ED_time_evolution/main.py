# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test: real-time quench evolution of an MPS via TEBD
# (pyitensor/tebd.py's TEBDEvolver, self.tevol_method="TEBD",
# itensor_version="python" only -- see CLAUDE.md and
# timedependent.py's evolve_and_measure_dmrg()) must agree with exact
# diagonalization on a small nearest-neighbor system. Modeled directly on
# examples/time_evolution/tdvp_VS_ED_time_evolution, with the same
# structure and tolerance, but exercising TEBD instead of TDVP.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 4 # small enough for exact diagonalization
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins)
sc.setup_python() # TEBD only exists on the pure-Python backend
sc.tevol_method = "TEBD"
assert sc.itensor_version=="python" and sc.tevol_method=="TEBD"

# create two Hamiltonians: prepare the GS of h0, then quench to h1
h0 = 0
h1 = 0
for i in range(n):
    h0 = h0+(-1)**i*sc.Sz[i] # Neel-favoring field (onsite -> split across bonds)
for i in range(n-1):
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1] # strictly NN

sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG (TEBD) ground state
wfED = sc.get_gs(mode="ED") # ED ground state, same Hamiltonian
sc.set_hamiltonian(h1) # quench to the Heisenberg Hamiltonian

op = sc.Sz[0] # operator to measure

nt = 50 # number of time steps
dt = 0.05 # time step

(ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
(tsED,szED) = timedependent.evolve_and_measure(sc,
        operator=op,nt=nt,dt=dt,wf=wfED,mode="ED")

diff = np.max(np.abs(sz-szED))
print("<Sz_0>(t) sample (DMRG, TEBD):",sz[:5].real)
print("<Sz_0>(t) sample (ED):",szED[:5].real)
print("Max abs difference between TEBD and ED =",diff)

tol = 1e-4
assert diff<tol, "TEBD time evolution disagrees with ED by %g (tol=%g)"%(diff,tol)

print("TEST PASSED")

# visualize the comparison behind the assert above
import matplotlib.pyplot as plt
plt.plot(ts,sz.real,label="DMRG (TEBD)",c="blue")
plt.scatter(tsED,szED.real,label="ED",c="red")
plt.xlabel("time") ; plt.ylabel(r"$\langle S_0^z\rangle(t)$")
plt.legend()
plt.show()
