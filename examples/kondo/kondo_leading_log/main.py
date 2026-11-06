# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# STM spectrum dI/dV of a magnetic atom beyond third order.
# Each order adds one more power of the log Jrho_s*ln(omega0/|eV|),
# which is large near zero bias and near the spin-excitation steps.
# order="LL" sums these leading logs to all orders (renormalization group),
# and its first order in Jrho_s is exactly order=3.
import numpy as np
from dmrgpy import spinchain

muB = 5.7883818066e-5 # Bohr magneton in eV/T
eV = np.linspace(-3e-3, 3e-3, 241) + 1e-5 # bias in eV

# a spin-1/2 at zero field, and a spin-1 with anisotropy in a field
s1 = spinchain.Spin_Chain(["S=1/2"])
s1.set_hamiltonian(0*s1.Sz[0])
s2 = spinchain.Spin_Chain(["S=1"])
s2.set_hamiltonian(-1e-3*s2.Sz[0]*s2.Sz[0] + 2*muB*2.0*s2.Sz[0]) # D=-1 meV, B=2 T

# at weak coupling the leading logs and the third order must agree
k = dict(Jrho_s=-0.002, U=0.25, T=0.0)
d2 = s2.get_kondo_spectrum(eV, order=2, **k)[1]
d3 = s2.get_kondo_spectrum(eV, order=3, **k)[1]
dl = s2.get_kondo_spectrum(eV, order="LL", **k)[1]
diff = np.max(np.abs(dl - d3))/np.max(np.abs(d3 - d2))
print("Max difference at weak coupling, relative to the 3rd order term", diff)
assert diff < 0.05

import matplotlib.pyplot as plt
for i, sc in enumerate([s1, s2]):
    plt.subplot(1, 2, i+1)
    for order in [2, 3, "LL"]:
        e, d = sc.get_kondo_spectrum(eV, order=order, Jrho_s=-0.05, U=0.25, T=0.0)
        plt.plot(e*1e3, d/(2*np.pi), label="order="+str(order))
    plt.xlabel("bias eV (meV)")
    plt.ylabel("dI/dV (e^2 T0^2/h)")
    plt.legend()
plt.show()
