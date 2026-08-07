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

# Sweep the hopping phase and track how the magnitude of each fitted
# effective coupling (the same coefficients that get formatted into the
# latex string above) evolves with it -- a cheap way to visualize the
# effective-Hamiltonian fit instead of only printing one phi's latex form.
phis = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
coupling_series = {}
for idx, phi in enumerate(phis):
    fci = build_chain(phi)
    coef = effectivehamiltonian.get_effective_hamiltonian_couplings(fci,
            method="single", mode="ED")
    # acceptable_matrix() (effectivehamiltonian.py) can drop a coupling at
    # a given phi if its matrix representation vanishes there (e.g. some
    # terms are exactly zero at phi=0) or is too collinear with another
    # already-kept one -- so the set of returned keys isn't the same at
    # every phi; pad with 0 for phis where a given key wasn't fitted.
    for key, val in coef.items():
        coupling_series.setdefault(key, [0.0]*len(phis))
        coupling_series[key][idx] = abs(val)

plt.figure(figsize=(7, 4.5))
for key, vals in coupling_series.items():
    if max(vals) < 1e-3: continue # skip negligible couplings from the plot
    plt.plot(phis, vals, "o-", label=" ".join(key) if isinstance(key, tuple) else str(key))
plt.xlabel(r"hopping phase $\phi/\pi$")
plt.ylabel("|fitted effective coupling|")
plt.title("Effective Hamiltonian couplings vs hopping phase (n=%d dimer)" % n)
plt.legend(fontsize=7)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()












