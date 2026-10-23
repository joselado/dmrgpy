# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Rotation of a fully polarized spin chain around the x axis,
# |psi(theta)> = exp(i theta Sx_total)|up...up>. The generators Sx_i of
# different sites commute, so applying the n single-site exponentials one
# after the other must give the same state as one exponential of the sum.
# We check it with the overlap <psi(0)|psi(theta)>.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.set_hamiltonian(-sum(sc.Sz))
wf0 = sc.get_gs() # all spins up
Mx = sum(sc.Sx) # total Sx

thetas = np.linspace(0., np.pi, 10)
seq = [] ; one = []
for theta in thetas:
    wf = wf0.copy()
    for i in range(n): # one site at a time
        wf = sc.exponential(1j*theta*sc.Sx[i], wf)
    seq.append(abs(wf0.overlap(wf)))
    one.append(abs(wf0.overlap(sc.exponential(1j*theta*Mx, wf0)))) # all at once
print("Max difference =", np.max(np.abs(np.array(seq)-one)))

import matplotlib.pyplot as plt
plt.plot(thetas, one, linewidth=4, label="one exponential")
plt.plot(thetas, seq, "--", linewidth=2, label="product of exponentials")
plt.xlabel("angle theta")
plt.ylabel("|<psi(0)|psi(theta)>|")
plt.legend()
plt.show()
