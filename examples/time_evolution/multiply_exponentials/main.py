# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

hfe = sum(sc.Sz) # Fully ferromagnetic
sc.set_hamiltonian(-hfe) ; wf0 = sc.get_gs() # Fully up wavefunction

Mx = sum(sc.Sx) # total Sx, the sum of the same on-site generators

# Consistency check: rotating every site one at a time by the same angle
# (multiplying n single-site exponentials) must give the same state as
# applying a single exponential of the summed generator, since the local
# Sx[i] all commute with each other. Sweep theta up to pi, the angle at
# which every S=1/2 site is fully flipped and the state becomes orthogonal
# to the starting one.
thetas = np.linspace(0.,np.pi,10)
overlap_seq = [] # apply the n single-site rotations sequentially
overlap_sum = [] # apply exp(i*theta*Mx) directly, in one shot
for theta in thetas:
    wf_seq = wf0.copy()
    for i in range(n):
        wf_seq = sc.exponential(1j*theta*sc.Sx[i],wf_seq)
    wf_sum = sc.exponential(1j*theta*Mx,wf0)
    overlap_seq.append(wf0.overlap(wf_seq))
    overlap_sum.append(wf0.overlap(wf_sum))

overlap_seq = np.abs(overlap_seq)
overlap_sum = np.abs(overlap_sum)
print("Max |sequential - single-shot| overlap difference:",
        np.max(np.abs(overlap_seq-overlap_sum)))

import matplotlib.pyplot as plt

plt.plot(thetas,overlap_seq,marker="o",label="sequential single-site exponentials")
plt.plot(thetas,overlap_sum,marker="x",linestyle="--",
        label="exponential of the summed generator")
plt.xlabel(r"Rotation angle $\theta$")
plt.ylabel(r"$|\langle\Psi_0|\Psi(\theta)\rangle|$")
plt.legend()
plt.show()
