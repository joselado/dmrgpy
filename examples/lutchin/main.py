import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
sys.path.append(os.environ["PYGRAROOT"]) # root for pygra
import matplotlib.pyplot as plt
import fermionchain
import geometry

g = geometry.chain()
g = g.supercell(10)
h = g.get_hamiltonian()
h.add_rashba(0.5)
h.shift_fermi(-1.5)
h.add_zeeman(0.3)

fc = fermionchain.Fermionic_Hamiltonian(len(g.r))
fc.set_spinful_hoppings(h.intra)

print("DMRG",fc.gs_energy())
print("ED",fc.gs_energy_free())


