# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin excitations of a S=5/2 molecule (Fe-OEP-Cl) on superconducting Pb,
# measured with a Pb tip. H = D Sz^2, so the first excitation is 2D.
# On a normal metal each excitation gives a step at eV = 2D. With a gap in
# tip and sample the step becomes a peak at eV = Delta_t + Delta_s + 2D,
# next to the elastic BCS peak at Delta_t + Delta_s.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.ysrtk import bcsperturbative

D, Delta = 0.7e-3, 1.35e-3 # anisotropy and gap of Pb, in eV
sc = spinchain.Spin_Chain(["S=5/2"]) # one spin 5/2
sc.set_hamiltonian(D*sc.Sz[0]*sc.Sz[0])

eV = np.linspace(-6e-3, 6e-3, 2401) # bias
U = 1.15 # direct tunneling, sets the elastic weight
# normal tip and normal substrate
(e, G0) = bcsperturbative.get_spin_excitation_spectrum(sc, eV, T=1.2, U=U)
# Pb tip on Pb, Dynes broadening gamma = Delta/100
(e, G1) = bcsperturbative.get_spin_excitation_spectrum(sc, eV, T=1.2, U=U,
            Delta_s=Delta, Delta_t=Delta, gamma=1e-2*Delta)
# the same at 4.8 K, where the state Sz=3/2 is thermally occupied
(e, G2) = bcsperturbative.get_spin_excitation_spectrum(sc, eV, T=4.8, U=U,
            Delta_s=Delta, Delta_t=Delta, gamma=1e-2*Delta)

big = (e > 0) & (G1 > 0.01*np.max(G1)) # positive bias, the main peaks
peaks = [e[i] for i in range(1, len(e)-1)
         if G1[i] > G1[i-1] and G1[i] > G1[i+1] and big[i]]
print("Peaks at 1.2 K (meV):", np.round(np.array(peaks)*1e3, 3))
print("Expected 2 Delta =", round(2*Delta*1e3, 3),
      " 2 Delta + 2D =", round((2*Delta + 2*D)*1e3, 3))
k1, k2 = np.argmin(np.abs(e - 2*Delta)), np.argmin(np.abs(e - 2*Delta - 2*D))
print("Height ratio inelastic/elastic:", round(G1[k2]/G1[k1], 3))

import matplotlib.pyplot as plt
plt.plot(e*1e3, G0/np.max(G0), linewidth=3, label="normal metal")
plt.plot(e*1e3, G1/np.max(G1), label="superconductor, 1.2 K")
plt.plot(e*1e3, G2/np.max(G2), "--", label="superconductor, 4.8 K")
plt.xlabel("bias (mV)")
plt.ylabel("dI/dV (normalized)")
plt.legend()
plt.show()
