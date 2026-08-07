# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") for a (spinful) Hubbard chain, built the same way as
# examples/hubbard_gap/main.py: two interleaved spinless-fermion sites
# (up/down) per physical site, plus an on-site U(n_up-1/2)(n_dn-1/2)
# interaction. Exercises Fermionic_Chain (not just Spin_Chain) and
# pyitensor's ElectronSite/Jordan-Wigner threading on the python backend.
import numpy as np
from dmrgpy import fermionchain

ns = 3 # physical sites (2*ns fermionic sites total)
U = 4.0

def get_energy(itensor_version,U):
    n = ns*2
    fc = fermionchain.Fermionic_Chain(n,itensor_version=itensor_version)
    C = fc.C ; Cdag = fc.Cdag ; N = fc.N
    h = 0
    for i in range(ns-1):
        for j in range(2):
            h = h + Cdag[2*i+j]*C[2*(i+1)+j]
    h = h + h.get_dagger()
    for i in range(ns):
        h = h + U*(N[2*i]-0.5)*(N[2*i+1]-0.5)
    fc.set_hamiltonian(h)
    return fc.gs_energy()

e2 = get_energy(2,U)
e3 = get_energy(3,U)
epy = get_energy("python",U)

print("Ground state energy (ITensor v2)     =",e2)
print("Ground state energy (ITensor v3)     =",e3)
print("Ground state energy (pure Python)    =",epy)
print("Difference v2 vs v3                  =",abs(e2-e3))
print("Difference v3 vs pure Python         =",abs(e3-epy))

# sweep the Hubbard U (cheap for this small ns=3 chain), comparing all
# three backends across the interaction strength
Us = [0.0,1.0,2.0,4.0,8.0]
e2s = [get_energy(2,Ui) for Ui in Us]
e3s = [get_energy(3,Ui) for Ui in Us]
epys = [get_energy("python",Ui) for Ui in Us]
for Ui,e2i,e3i,epyi in zip(Us,e2s,e3s,epys):
    print("U =",Ui,"v2 =",e2i,"v3 =",e3i,"python =",epyi)

import matplotlib.pyplot as plt
plt.plot(Us,e2s,marker="o",label="ITensor v2")
plt.plot(Us,e3s,marker="x",linestyle="--",label="ITensor v3")
plt.plot(Us,epys,marker="s",linestyle=":",label="pure Python")
plt.xlabel("U")
plt.ylabel("Ground-state energy")
plt.legend()
plt.show()
