# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Time needed for the finite-temperature correlator <Sz_0(t) Sz_0(0)> of
# the Heisenberg chain with METTS (minimally entangled typical thermal
# states), as a function of the number of samples. It is computed
# with ITensor v3 and with the pure Python backend. The cost grows
# linearly with the number of samples, since each one is a time evolution.
import time
from dmrgpy import spinchain

def run(version, nsamples):
    n = 8 # number of sites
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    t0 = time.time()
    sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), T=1.0, nt=10, dt=0.1,
            nsamples=nsamples, nwarmup=10, dbeta_half_step=0.08, seed=42)
    return time.time() - t0

nsamples = [5, 10, 20]
import matplotlib.pyplot as plt
for version in [3, "python"]:
    ts = [run(version, ns) for ns in nsamples]
    print(version, ts)
    plt.plot(nsamples, ts, "o-", label=str(version))
plt.xlabel("number of samples")
plt.ylabel("time (s)")
plt.legend()
plt.show()
