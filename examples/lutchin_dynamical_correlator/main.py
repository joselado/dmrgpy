import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
sys.path.append(os.environ["PYGRAROOT"]) # root for pygra
import matplotlib.pyplot as plt
import fermionchain
import geometry

n = 4
g = geometry.chain()
g = g.supercell(n)
h = g.get_hamiltonian()
h.add_rashba(0.5)
h.shift_fermi(1.7)
h.add_zeeman(0.6)

def fd(i,j):
    if i==j: return 0.1
    else: return 0.0

fc = fermionchain.Fermionic_Hamiltonian(len(g.r))
fc.set_spinful_hoppings(h.intra)
fc.set_pairing(fd)
fc.get_dynamical_correlator(mode="DMRG",i=0,j=0,name="cdc",delta=0.03)

