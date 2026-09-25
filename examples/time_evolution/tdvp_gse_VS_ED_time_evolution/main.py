# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test: real-time quench evolution of an MPS via one-site TDVP
# with Krylov global subspace expansion (tevol_method="TDVP_GSE",
# mpscpp3/TDVP/basisextension.h's addBasis(), arXiv:2005.06104) must agree
# with exact diagonalization on a small system, same as the plain two-site
# tevol_method="TDVP" already checked by
# examples/tdvp_VS_ED_time_evolution/main.py (this script mirrors it).
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 4 # small enough for exact diagonalization
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins)
sc.tevol_method = "TDVP_GSE"
sc.tdvp_gse_sweeps = 3
sc.tdvp_gse_krylov_order = 3
sc.tdvp_gse_cutoff = 1e-8
assert sc.itensor_version==3

# create two Hamiltonians: prepare the GS of h0, then quench to h1
h0 = 0
h1 = 0
for i in range(n):
    h0 = h0+(-1)**i*sc.Sz[i] # Neel-favoring field
# A field-only h0 has an exact *product-state* ground state (bond
# dimension 1). That used to trip itensor_version=3's TDVP_GSE (up to 5 of
# 8 trials off by 0.02-0.04), and the small XX+YY coupling below was mixed
# in to avoid it. The cause was general, not specific to bond dimension 1:
# any start whose site 0 is one local basis vector lost the expansion at
# the left edge, because Chain::global_subspace_expand truncated with
# Cutoff=0, which in ITensor v3 discards exactly-zero weights. It is fixed
# (docs/audit_2026_09_24c_hole_hunt.md, finding 14; a pure Neel start is
# now exact in 8 of 8 trials), and the coupling is kept only so this
# script's numbers stay what they were.
for i in range(n-1):
    h0 = h0+0.3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1])
for i in range(n-1):
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h0)
wf = sc.get_gs() # DMRG ground state (still plain DMRG, not TDVP)
wfED = sc.get_gs(mode="ED") # ED ground state, same Hamiltonian
sc.set_hamiltonian(h1) # quench to the Heisenberg Hamiltonian

op = sc.Sz[0] # operator to measure

nt = 50 # number of time steps
dt = 0.05 # time step

(ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
(tsED,szED) = timedependent.evolve_and_measure(sc,
        operator=op,nt=nt,dt=dt,wf=wfED,mode="ED")

diff = np.max(np.abs(sz-szED))
print("<Sz_0>(t) sample (DMRG, TDVP_GSE):",sz[:5].real)
print("<Sz_0>(t) sample (ED):",szED[:5].real)
print("Max abs difference between TDVP_GSE and ED =",diff)

tol = 1e-4
assert diff<tol, "TDVP_GSE time evolution disagrees with ED by %g (tol=%g)"%(diff,tol)

print("TEST PASSED")

# visualize the comparison behind the assert above
import matplotlib.pyplot as plt
plt.plot(ts,sz.real,label="DMRG (TDVP_GSE)",c="blue")
plt.scatter(tsED,szED.real,label="ED",c="red")
plt.xlabel("time") ; plt.ylabel(r"$\langle S_0^z\rangle(t)$")
plt.legend()
plt.show()
