# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order STM/Kondo spectrum dI/dV at T>0 of a Heisenberg chain in a
# field that turns from site to site in the xz plane (no Sz conservation).
# submode="KRYLOV" takes the thermal states from a sparse solver and the
# other states from one Krylov space per thermal state, without diagonalizing
# the whole Hamiltonian. We compare it with full exact diagonalization.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

n = 8 # number of S=1/2 spins
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
for i in range(n): # a 120 degree field pattern, 0.2 meV
    h = h + 2e-4*(np.cos(2*np.pi*i/3)*sc.Sx[i] + np.sin(2*np.pi*i/3)*sc.Sz[i])
sc.set_hamiltonian(h)

eV = np.linspace(-4e-3, 4e-3, 101) + 1e-5 # bias voltage
kw = dict(site=1, Jrho_s=-0.05, U=0.2, omega0=20e-3, Gamma0=5e-5, order=3)
for T in [0.3, 1.0]: # temperature in K
    _, ed = sc.get_kondo_spectrum(eV, T=T, **kw) # full diagonalization
    _, kr = sc.get_kondo_spectrum(eV, T=T, submode="KRYLOV", nkry=100, **kw)
    print("T =", T, "K, max difference =", np.max(np.abs(kr-ed))/np.max(ed))
    assert np.max(np.abs(kr-ed)) < 1e-6*np.max(ed)
    plt.plot(eV*1e3, ed, linewidth=4, label="ED, T="+str(T)+" K")
    plt.plot(eV*1e3, kr, "k--", linewidth=2, label="Krylov, T="+str(T)+" K")
plt.xlabel("eV (meV)")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
