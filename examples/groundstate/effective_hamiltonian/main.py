# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import effectivehamiltonian
n = 2 # number of spinful fermionic sites

def build_chain(phi):
    """Spinful fermionic dimer with a Peierls hopping phase phi (in units
    of pi) and a strong Hubbard repulsion, same construction used below
    for the single latex printout."""
    fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
    h = 0
    t = 1.0*np.exp(1j*phi*np.pi)
    for i in range(n-1): # hopping
        h = h + t*fc.Cdagup[i]*fc.Cup[i+1]
        h = h + np.conjugate(t)*fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n): # Hubbard
        h = h + 20*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    return fc

phi0 = 0.1
fc = build_chain(phi0)

# print the effective Hamiltonian in latex form
l = effectivehamiltonian.get_effective_hamiltonian(fc,method="single",
        mode="ED")
print("Effective Hamiltonian in latex form")
print(l) # write the Hamiltonian in latex

# Sweep the hopping phase and check the fitted couplings against the
# analytic answer. At strong coupling the Hubbard dimer maps onto a
# two-spin exchange model, and a Peierls phase phi twists the two spin
# species oppositely, so the exchange stays isotropic in magnitude while
# its XY part is rotated by 2*phi:
#
#   H_eff = J [ cos(2 phi)(SxSx + SySy) + sin(2 phi)(SxSy - SySx) + SzSz ]
#
# with J = 4 t^2 / U. Note h = h + h.get_dagger() below doubles the
# Hubbard term (it is already Hermitian) while leaving |t| = 1, so the
# effective U is 40 and J = 4/40 = 0.1.
J = 4.0 * 1.0**2 / 40.0

phis = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
xx, xy, zz = [], [], []
for phi in phis:
    coef = effectivehamiltonian.get_effective_hamiltonian_couplings(
            build_chain(phi), method="single", mode="ED")
    # get() rather than [] : a coupling that is exactly zero at some phi
    # (e.g. the XY rotation at phi=0) falls below the fit's cutoff and is
    # simply absent from the returned dict
    xx.append(np.real(coef.get(("S^x_1", "S^x_2"), 0.0)))
    xy.append(np.real(coef.get(("S^x_1", "S^y_2"), 0.0)))
    zz.append(np.real(coef.get(("S^z_1", "S^z_2"), 0.0)))
xx, xy, zz = np.array(xx), np.array(xy), np.array(zz)

print("phi/pi   J_xx      J_xy      J_zz      |J|")
for k, phi in enumerate(phis):
    print("%4.1f   %8.5f  %8.5f  %8.5f  %8.5f"
          % (phi, xx[k], xy[k], zz[k], np.sqrt(xx[k]**2 + xy[k]**2)))
# the magnitude must stay at J for every phi -- the phase only rotates it
assert np.allclose(np.sqrt(xx**2 + xy**2), J, atol=1e-3)
assert np.allclose(zz, J, atol=1e-3)

grid = np.linspace(0, 0.5, 200)
plt.figure(figsize=(7, 4.5))
plt.plot(grid, J*np.cos(2*np.pi*grid), "C0-", lw=1, label=r"$J\cos 2\phi$ (exact)")
plt.plot(grid, J*np.sin(2*np.pi*grid), "C1-", lw=1, label=r"$J\sin 2\phi$ (exact)")
plt.axhline(J, color="C2", lw=1, label=r"$J$ (exact)")
plt.plot(phis, xx, "C0o", label=r"fitted $S^x_1S^x_2$")
plt.plot(phis, xy, "C1s", label=r"fitted $S^x_1S^y_2$")
plt.plot(phis, zz, "C2^", label=r"fitted $S^z_1S^z_2$")
plt.xlabel(r"hopping phase $\phi/\pi$")
plt.ylabel("effective coupling")
plt.title("Effective spin couplings of a Hubbard dimer vs Peierls phase")
plt.legend(fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
