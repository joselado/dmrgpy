# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Time needed for the ground state energy of the S=1/2 Heisenberg chain
# as a function of the chain length, for the three DMRG backends: ITensor
# v2, ITensor v3 and the pure Python one. The result shows how much
# slower the Python backend is, which does not need any compiler.
import time
from dmrgpy import spinchain

def run(version, n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    t0 = time.time()
    sc.gs_energy()
    return time.time() - t0

sizes = [8, 12, 16, 20]
import matplotlib.pyplot as plt
for version in [2, 3, "python"]:
    ts = [run(version, n) for n in sizes]
    print(version, ts)
    plt.plot(sizes, ts, "o-", label=str(version))
plt.xlabel("chain length n")
plt.ylabel("time (s)")
plt.yscale("log")
plt.legend()
plt.show()
