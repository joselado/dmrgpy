import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
sys.path.append(os.environ["PYGRAROOT"]) # root for pygra
import matplotlib.pyplot as plt
import fermionchain
import geometry

def getfc(n=10):
  g = geometry.chain()
  g = g.supercell(n)
  h = g.get_hamiltonian()
  h.add_rashba(0.5)
  h.shift_fermi(1.8)
  h.add_zeeman(0.3)
  hk = h.get_hk_gen()
  
  def fd(i,j):
      if i==j: return 0.2
      else: return 0.0
  
  fc = fermionchain.Fermionic_Hamiltonian(len(g.r))
  fc.set_spinful_hoppings(h.intra)
  fc.set_pairing(fd)
  return fc

import matplotlib.pyplot as plt
ns = range(4,44,4)
f = open("GAP.OUT","w")
for n in ns:
    print(n)
    g = getfc(n).get_gap()
    f.write(str(n)+"  ")
    f.write(str(g)+"\n")
    f.flush()
f.close()

