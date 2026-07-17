# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: real-time quench evolution of an MPS via the new TDVP
# engine (mpscpp3/TDVP/, the itensor_version=3 default -- see
# CLAUDE.md's "real-time MPS evolution defaults to TDVP" and
# timedependent.py's evolve_and_measure_dmrg()) must agree with exact
# diagonalization on a small system. Modeled on
# examples/time_evolution/time_evolution_quench and
# examples/v2_VS_v3_time_evolution_quench, but with actual asserts instead
# of just printing/plotting, since this is specifically meant as a test.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 4 # small enough for exact diagonalization
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins)
assert sc.itensor_version==3 and sc.tevol_method=="TDVP" # confirm the default

# create two Hamiltonians: prepare the GS of h0, then quench to h1
h0 = 0
h1 = 0
for i in range(n):
    h0 = h0+(-1)**i*sc.Sz[i] # Neel-favoring field
for i in range(n-1):
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG (TDVP) ground state
wfED = sc.get_gs(mode="ED") # ED ground state, same Hamiltonian
sc.set_hamiltonian(h1) # quench to the Heisenberg Hamiltonian

op = sc.Sz[0] # operator to measure

nt = 50 # number of time steps
dt = 0.05 # time step

(ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
(tsED,szED) = timedependent.evolve_and_measure(sc,
        operator=op,nt=nt,dt=dt,wf=wfED,mode="ED")

diff = np.max(np.abs(sz-szED))
print("<Sz_0>(t) sample (DMRG, TDVP):",sz[:5].real)
print("<Sz_0>(t) sample (ED):",szED[:5].real)
print("Max abs difference between TDVP and ED =",diff)

tol = 1e-4
assert diff<tol, "TDVP time evolution disagrees with ED by %g (tol=%g)"%(diff,tol)

print("TEST PASSED")
