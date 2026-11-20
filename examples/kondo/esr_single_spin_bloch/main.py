# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# ESR-STM line of one spin 1/2 in a field. The RF drive moves population
# from the lower to the upper Zeeman state, and the tip reads it through the
# tunnel magnetoresistance. At weak drive the line is a Lorentzian of full
# width 1/(pi T2); a stronger drive saturates it and makes it broader, as in
# the Bloch formula Omega^2 T1 T2/(1 + (2 pi df T2)^2 + Omega^2 T1 T2).
import numpy as np
from dmrgpy import spinchain

h = 4.135667696e-6 # Planck constant in eV/GHz
kB = 8.617333262e-5 # Boltzmann constant in eV/K
sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(h*20.0*sc.Sz[0]) # Zeeman splitting of 20 GHz

th = 0.7 # angle between the field and the tip moment
tip = [np.sin(th), 0., np.cos(th)]
f = np.linspace(19.9, 20.1, 801) # RF frequency in GHz
T, T1, T2 = 1.2, 100., 50. # temperature in K, times in ns at 20 GHz

import matplotlib.pyplot as plt
for rabi in [1e-3, 3e-3, 1e-2]: # Rabi frequency in GHz
    f, dI = sc.get_esr_spectrum(f, tip=tip, T=T, T1=T1, T2=T2, f_ref=20.0,
                                rabi=rabi)
    # Bloch formula, with the Rabi frequency of this line, rabi*sin(th)
    xi2 = (2*np.pi*rabi*np.sin(th))**2*T1*T2
    x2 = (2*np.pi*(f - 20.0)*T2)**2
    bloch = 0.5*np.cos(th)*np.tanh(h*20.0/(2*kB*T))*xi2/(1 + x2 + xi2)
    print("rabi =", rabi, "GHz, max difference with Bloch =",
          np.max(np.abs(dI - bloch)))
    assert np.max(np.abs(dI - bloch)) < 1e-10
    plt.plot(f, bloch, linewidth=4, label="Bloch, rabi = "+str(rabi)+" GHz")
    plt.plot(f, dI, "--", linewidth=2, color="black")
plt.xlabel("RF frequency (GHz)")
plt.ylabel("dI (units of G a V_DC)")
plt.legend()
plt.show()
