# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test for dmrgpy.atom.generate_atom() (a 5-orbital d-shell
# spinful-fermion atom: hopping + Hund's coupling via S^2/L^2 + SOC +
# external/exchange fields + a Lagrange-multiplier filling constraint).
# Until now atomtk/ had no test or example coverage at all, so a change
# there (e.g. reordering how the internal one2many() helper accumulates
# MultiOperator terms) had nothing locking in correctness. Golden values
# below were computed with mode="ED" on a small, hand-written (not
# random, for reproducibility) hopping matrix.
import numpy as np
from dmrgpy.atom import generate_atom

tol = 1e-6

orbs = ["dz2", "dxz", "dyz", "dx2y2", "dxy"]
tij = np.array([
    [ 0.20,  0.10,  0.00,  0.00,  0.00],
    [ 0.10, -0.10,  0.05,  0.00,  0.00],
    [ 0.00,  0.05,  0.05,  0.00,  0.00],
    [ 0.00,  0.00,  0.00,  0.15,  0.02],
    [ 0.00,  0.00,  0.00,  0.02, -0.05],
])

fc = generate_atom(orbs=orbs, tij=tij, U=4., B=[0.05, 0., 0.1],
                    Js=[0., 0.1, 0.], soc=0.2, J=0.5, Ne=5, lamb_Ne=200)

e0 = fc.gs_energy(mode="ED")
s2 = fc.vev(fc.ST2, mode="ED")
l2 = fc.vev(fc.LT2, mode="ED")
ntot = sum(fc.vev(Ni, mode="ED") for Ni in fc.N)

e0_ref = -35.37282567472432
s2_ref = 8.747265937476783
l2_ref = 0.001100335706273672
ntot_ref = 5.0 # enforced by the lamb_Ne filling constraint

print("E0   = %.12f (ref %.12f)"%(e0.real, e0_ref))
print("<S2> = %.12f (ref %.12f)"%(s2.real, s2_ref))
print("<L2> = %.12f (ref %.12f)"%(l2.real, l2_ref))
print("<N>  = %.12f (ref %.12f)"%(ntot.real, ntot_ref))

assert abs(e0.real-e0_ref)<tol, "ground state energy drifted from golden value"
assert abs(s2.real-s2_ref)<tol, "<S^2> drifted from golden value"
assert abs(l2.real-l2_ref)<tol, "<L^2> drifted from golden value"
assert abs(ntot.real-ntot_ref)<1e-3, "filling constraint not enforced"

print("TEST PASSED")

# Sweep the Hubbard U to see how the ground-state energy and total spin
# respond, reusing the same hopping matrix/orbitals/reference physics
# above (cheap: same small 5-orbital atom, ED only).
Us = [2.,3.,4.,5.,6.]
e0s = []
s2s = []
for U in Us:
    fcU = generate_atom(orbs=orbs, tij=tij, U=U, B=[0.05, 0., 0.1],
                         Js=[0., 0.1, 0.], soc=0.2, J=0.5, Ne=5, lamb_Ne=200)
    e0s.append(fcU.gs_energy(mode="ED").real)
    s2s.append(fcU.vev(fcU.ST2, mode="ED").real)
    print("U =",U,"E0 =",e0s[-1],"<S2> =",s2s[-1])

import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(Us, e0s, marker="o")
ax1.set_xlabel("U")
ax1.set_ylabel("Ground-state energy")
ax2.plot(Us, s2s, marker="o", color="tab:orange")
ax2.set_xlabel("U")
ax2.set_ylabel(r"$\langle S^2\rangle$")
plt.tight_layout()
plt.show()
