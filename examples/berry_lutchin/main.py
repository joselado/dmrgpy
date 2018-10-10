import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
sys.path.append(os.environ["PYGRAROOT"]) # root for pygra
import matplotlib.pyplot as plt
import fermionchain
import geometry


# This is not working.....
raise

g = geometry.chain()
g = g.supercell(2)
h = g.get_hamiltonian()
h.add_rashba(0.5)
h.shift_fermi(1.8)
h.add_zeeman(0.6)
#h.add_swave(0.2)
#h.get_bands()
#exit()
hk = h.get_hk_gen()

fc = fermionchain.Fermionic_Hamiltonian(len(g.r))

def getfc(k):
  fc.set_spinful_hoppings(hk(k)) # add hoppings
#  print(k,fc.gs_energy(),fc.gs_energy_free())
  def fd(i,j):
      if i==j: return 0.1
      return 0.0
  fc.set_pairing(fd)
#  print(fc.get_excited(3))
  return fc

import topology
print(topology.berry_phase(getfc))
#h.get_bands()
exit()

fc = fermionchain.Fermionic_Hamiltonian(len(g.r))
fc.set_spinful_hoppings(h.intra)

print("DMRG",fc.gs_energy())
print("ED",fc.gs_energy_free())


