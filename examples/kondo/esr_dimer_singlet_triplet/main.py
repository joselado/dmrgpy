# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# ESR-STM of the singlet to T0 transition of two coupled spins 1/2, as seen
# by Bae et al. (arXiv:1807.04322). Without a tip field, S and T0 have the
# same spin on the atom under the tip, so the DC readout does not see the
# line and the homodyne part dominates: the line is antisymmetric (a small
# symmetric part comes from the relaxation of S and T0 to T- and T+). A tip
# field mixes S and T0 into |01> and |10>, and the line becomes symmetric.
import numpy as np
from dmrgpy import spinchain

h = 4.135667696e-6 # Planck constant in eV/GHz
J, B = 15.0, 25.0 # exchange and Zeeman energy in GHz

import matplotlib.pyplot as plt
for btip in [0.0, 4.0, 16.0]: # tip field on the first spin, in GHz
    sc = spinchain.Spin_Chain(["S=1/2", "S=1/2"])
    hh = J*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    hh = hh + B*(sc.Sz[0] + sc.Sz[1]) + btip*sc.Sz[0]
    sc.set_hamiltonian(h*hh) # energies in eV
    f0 = np.sqrt(J**2 + btip**2) # frequency of the S-T0 line
    f = f0 + np.linspace(-0.1, 0.1, 801) # RF frequency in GHz
    f, dI, poles = sc.get_esr_spectrum(f, tip=[0.5, 0., 0.87], T=1.2,
                                       T1=100., T2=20., f_ref=B, rabi=5e-4,
                                       vrf_over_vdc=0.2, return_poles=True)
    k = np.argmin(np.abs(poles["frequency"] - f0)) # the S-T0 line
    print("tip field", btip, "GHz: form factor of the S-T0 line =",
          poles["form_factor"][k])
    plt.plot(f - f0, dI/np.max(np.abs(dI)), linewidth=3,
             label="tip field = "+str(btip)+" GHz")
plt.xlabel("f - f_ST0 (GHz)")
plt.ylabel("dI (normalized)")
plt.legend()
plt.show()
