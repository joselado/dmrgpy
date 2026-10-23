# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Overlap between two infinite matrix product states (iMPS). An overlap
# over an infinite chain is 0 or infinite, so we compute the fidelity per
# site, which has modulus 1 only if the two states are the same state.
# Here we compare the ground state of the Heisenberg chain with the ground
# state of the XXZ chain, H = sum_i Sx Sx + Sy Sy + Delta Sz Sz.
import numpy as np
from dmrgpy import infinitechain

def get_state(delta):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]
                       + delta*ic.SzC[0]*ic.SzR[0])
    ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9 # bond dimension, iterations
    ic.gs_energy()
    return ic

ref = get_state(1.0) # Heisenberg chain
print("Self overlap =", ref.overlap(ref))
deltas = np.linspace(0.6, 1.4, 9)
fidelity = [abs(ref.overlap(get_state(d))) for d in deltas]

import matplotlib.pyplot as plt
plt.plot(deltas, fidelity, marker="o")
plt.xlabel("anisotropy Delta of the XXZ chain")
plt.ylabel("|overlap| per site with Delta=1")
plt.show()
