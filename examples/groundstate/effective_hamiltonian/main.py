# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Effective spin Hamiltonian of a Hubbard dimer with a Peierls phase phi
# on the hopping. At strong coupling H = J[cos(2phi)(SxSx+SySy)
# + sin(2phi)(SxSy-SySx) + SzSz], with J = 4t^2/U. We fit the couplings
# from the low energy states and compare with this formula.
import numpy as np
from dmrgpy import fermionchain, effectivehamiltonian

def dimer(phi): # Hubbard dimer, hopping phase phi in units of pi
    fc = fermionchain.Spinful_Fermionic_Chain(2)
    t = np.exp(1j*phi*np.pi)
    h = t*fc.Cdagup[0]*fc.Cup[1] + np.conjugate(t)*fc.Cdagdn[0]*fc.Cdn[1]
    h = h + 20*(fc.Nup[0]-.5)*(fc.Ndn[0]-.5) + 20*(fc.Nup[1]-.5)*(fc.Ndn[1]-.5)
    fc.set_hamiltonian(h + h.get_dagger()) # U=40 after adding the dagger
    return fc

J = 4./40. # J = 4 t^2/U
phis = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
xx, xy, zz = [], [], []
for phi in phis:
    fc = dimer(phi)
    Q = effectivehamiltonian.single_occupancy_projector(fc)
    c = effectivehamiltonian.get_effective_hamiltonian_couplings(
            fc, method="single", mode="ED", projector=Q)
    xx.append(np.real(c.get(("S^x_1", "S^x_2"), 0.0))) # absent if zero
    xy.append(np.real(c.get(("S^x_1", "S^y_2"), 0.0)))
    zz.append(np.real(c.get(("S^z_1", "S^z_2"), 0.0)))
assert np.allclose(np.sqrt(np.array(xx)**2+np.array(xy)**2), J, atol=1e-3)

import matplotlib.pyplot as plt
g = np.linspace(0, 0.5, 100)
plt.plot(g, J*np.cos(2*np.pi*g), "C0", label="J cos(2phi)")
plt.plot(g, J*np.sin(2*np.pi*g), "C1", label="J sin(2phi)")
plt.plot(phis, xx, "C0o", label="fitted SxSx")
plt.plot(phis, xy, "C1s", label="fitted SxSy")
plt.plot(phis, zz, "C2^", label="fitted SzSz")
plt.xlabel("phi/pi")
plt.ylabel("effective coupling")
plt.legend()
plt.show()
