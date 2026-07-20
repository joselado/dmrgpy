# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Same regression test as main.py, but for the pure-Python backend
# (itensor_version="python", pyitensor/gse.py's global_subspace_expand() +
# pyitensor/tdvp.py's one-site TDVP path) instead of mpscpp3.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 4 # small enough for exact diagonalization
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins)
sc.setup_python()
sc.tevol_method = "TDVP_GSE"
sc.tdvp_gse_sweeps = 3
sc.tdvp_gse_krylov_order = 3
sc.tdvp_gse_cutoff = 1e-8

# create two Hamiltonians: prepare the GS of h0, then quench to h1
h0 = 0
h1 = 0
for i in range(n):
    h0 = h0+(-1)**i*sc.Sz[i] # Neel-favoring field
# Mixed in a small XX+YY coupling so h0's ground state isn't an exact
# product state (bond dimension 1) -- see main.py's matching comment: a
# pure-field h0 triggers a real bond-dim-1 edge case in mpscpp3's C++
# global_subspace_expand()/one-site TDVP (not reproduced on this
# itensor_version="python" backend, but kept identical to main.py so
# both scripts test the exact same setup).
for i in range(n-1):
    h0 = h0+0.3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1])
for i in range(n-1):
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG ground state
wfED = sc.get_gs(mode="ED") # ED ground state, same Hamiltonian
sc.set_hamiltonian(h1) # quench to the Heisenberg Hamiltonian

op = sc.Sz[0] # operator to measure

nt = 50 # number of time steps
dt = 0.05 # time step

(ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
(tsED,szED) = timedependent.evolve_and_measure(sc,
        operator=op,nt=nt,dt=dt,wf=wfED,mode="ED")

diff = np.max(np.abs(sz-szED))
print("<Sz_0>(t) sample (python, TDVP_GSE):",sz[:5].real)
print("<Sz_0>(t) sample (ED):",szED[:5].real)
print("Max abs difference between TDVP_GSE (python) and ED =",diff)

tol = 1e-4
assert diff<tol, "TDVP_GSE (python) time evolution disagrees with ED by %g (tol=%g)"%(diff,tol)

print("TEST PASSED")
