# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") for a spin-1 Heisenberg chain (a different local Hilbert space
# than spin-1/2, exercising the stock SpinOne site type on all three
# backends).
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

n = 6

def get_energy(itensor_version,n):
    spins = ["S=1" for i in range(n)]
    sc = spinchain.Spin_Chain(spins,itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc.gs_energy()

e2 = get_energy(2,n)
e3 = get_energy(3,n)
epy = get_energy("python",n)

print("Ground state energy (ITensor v2)     =",e2)
print("Ground state energy (ITensor v3)     =",e3)
print("Ground state energy (pure Python)    =",epy)
print("Difference v2 vs v3                  =",abs(e2-e3))
print("Difference v3 vs pure Python         =",abs(e3-epy))

# sweep the chain length and overlay all three backends
ns = [3,4,5,6]
e2s = [get_energy(2,ni) for ni in ns]
e3s = [get_energy(3,ni) for ni in ns]
epys = [get_energy("python",ni) for ni in ns]

plt.plot(ns,e2s,marker="o",label="ITensor v2")
plt.plot(ns,e3s,marker="s",label="ITensor v3")
plt.plot(ns,epys,marker="x",label="pure Python")
plt.xlabel("Chain length")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
