# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import time
import numpy as np
from dmrgpy import fermionchain

# Demonstrates fermionchain.Spinful_Fermionic_Chain_Native: a Hubbard-model
# chain built directly on native spinful (Electron/Hubbard) sites -- one
# dimension-4 tensor-network site per orbital (Empty/Up/Down/UpDn) -- instead
# of Spinful_Fermionic_Chain.py's two interleaved dimension-2 sites per
# orbital. Only itensor_version=3 wires up this site type on the DMRG side.
#
# Part 1 checks the two classes agree exactly (same physics, same
# Jordan-Wigner sign convention, see multioperatortk/jordanwigner_spinful.py)
# for both the ground-state energy and a KPM dynamical correlator.
#
# Part 2 is the actual performance comparison this class was built to
# answer: does packing both spin flavors into one site (fewer, "fatter"
# tensor-network sites) pay off for two-site DMRG? Measured answer: no --
# see fermionchain.Spinful_Fermionic_Chain_Native's docstring for why.

n = 4
U = 2.0


def build(chain_cls, maxm=60, nsweeps=14):
    fc = chain_cls(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    fc.set_hamiltonian(h)
    fc.setup_cpp(3)
    return fc


print("### Part 1: correctness cross-check ###")
fc_native = build(fermionchain.Spinful_Fermionic_Chain_Native)
fc_double = build(fermionchain.Spinful_Fermionic_Chain)

e_native = fc_native.gs_energy(mode="DMRG")
e_double = fc_double.gs_energy(mode="DMRG")
print("Native  E0 =", e_native)
print("Doubled E0 =", e_double)
assert abs(e_native - e_double) < 1e-5

es = np.linspace(-4, 4, 60)
(_, y_native) = fc_native.get_dynamical_correlator(
        name=(fc_native.Cup[0], fc_native.Cdagup[0]), es=es, delta=0.3)
(_, y_double) = fc_double.get_dynamical_correlator(
        name=(fc_double.Cup[0], fc_double.Cdagup[0]), es=es, delta=0.3)
assert np.max(np.abs(y_native - y_double)) < 5e-3
print("KPM dynamical correlator matches between the two site conventions")

print()
print("### Part 2: performance comparison (fixed maxm, growing n) ###")
for n_size in [6, 10, 14]:
    n = n_size
    for chain_cls, label in [(fermionchain.Spinful_Fermionic_Chain_Native, "native (n sites, dim 4)"),
                              (fermionchain.Spinful_Fermionic_Chain, "doubled (2n sites, dim 2)")]:
        fc = build(chain_cls, maxm=80, nsweeps=16)
        t0 = time.time()
        e0 = fc.gs_energy(mode="DMRG")
        dt = time.time() - t0
        print(f"  n={n_size:2d}  {label:26s}  n_sites={fc.ns:2d}  E0={e0:.6f}  wall={dt:.2f}s")
