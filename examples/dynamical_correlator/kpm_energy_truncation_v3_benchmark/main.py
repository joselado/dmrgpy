# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0) of
# the S=1/2 Heisenberg chain with the kernel polynomial method (KPM) of
# ITensor v3. Energy truncation lets the KPM window be narrower than the
# full bandwidth, which gives more resolution per Chebyshev moment. The
# resonance must be at the same energy with and without truncation.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
es = np.linspace(0.3, 1.2, 41) # energies
spectra = []
for truncate in [False, True]:
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    sc.setup_cpp(3) # ITensor v3
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    if truncate:
        sc.kpm_scale = 0.65 # narrower energy window
        sc.kpm_energy_truncate = True
    name = (sc.Sz[0], sc.Sz[0])
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                         es=es, delta=0.08)
    spectra.append(np.array(s).real)
print("Max difference =", np.max(np.abs(spectra[0]-spectra[1])))

import matplotlib.pyplot as plt
plt.plot(es, spectra[0], linewidth=4, label="full window")
plt.plot(es, spectra[1], "--", linewidth=2, label="energy truncation")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
