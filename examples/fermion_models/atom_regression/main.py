# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of an atom with a 5-orbital d shell and 5 electrons.
# The model has hopping between orbitals, Hubbard U, Hund coupling J,
# spin-orbit coupling and a magnetic field. We compute the energy and
# the total spin <S^2> with ED, as a function of U.
import numpy as np
from dmrgpy.atom import generate_atom

orbs = ["dz2", "dxz", "dyz", "dx2y2", "dxy"]
tij = np.array([[ 0.20,  0.10,  0.00,  0.00,  0.00],
                [ 0.10, -0.10,  0.05,  0.00,  0.00],
                [ 0.00,  0.05,  0.05,  0.00,  0.00],
                [ 0.00,  0.00,  0.00,  0.15,  0.02],
                [ 0.00,  0.00,  0.00,  0.02, -0.05]]) # hopping matrix

Us = [2., 3., 4., 5., 6.] # Hubbard interactions
es, s2s = [], []
for U in Us:
    fc = generate_atom(orbs=orbs, tij=tij, U=U, B=[0.05, 0., 0.1],
                       Js=[0., 0.1, 0.], soc=0.2, J=0.5, Ne=5, lamb_Ne=200)
    es.append(fc.gs_energy(mode="ED").real)
    s2s.append(fc.vev(fc.ST2, mode="ED").real)

# reference value of the energy at U=4
print("E0(U=4) =", es[2], " reference = -35.37282567472432")
assert abs(es[2] + 35.37282567472432) < 1e-6

import matplotlib.pyplot as plt
plt.plot(Us, s2s, "o-")
plt.xlabel("Hubbard U")
plt.ylabel("total spin <S^2>")
plt.show()
