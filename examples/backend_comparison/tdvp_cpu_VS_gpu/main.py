# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time TDVP evolution of <Sz_0>(t) after a quench of the S=1/2 chain.
# The state starts as the ground state of a staggered field and evolves
# with the Heisenberg Hamiltonian. We run it with NumPy (CPU) and with
# JAX (GPU if there is one), and check that both give the same curve.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent
from dmrgpy.pyitensor import backend

n = 6 # number of sites
def evolve(name):
    backend.set_backend(name) # "numpy" is the CPU, "jax" the GPU
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)],itensor_version="python")
    sc.maxm = 16 # bond dimension
    sc.nsweeps = 4 # number of DMRG sweeps
    sc.set_hamiltonian(sum([(-1)**i*sc.Sz[i] for i in range(n)])) # staggered field
    wf = sc.get_gs()
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h) # sudden quench to the Heisenberg chain
    (ts,sz) = timedependent.evolve_and_measure(sc,operator=sc.Sz[0],nt=12,dt=0.05,wf=wf)
    return np.array(ts).real, np.array(sz).real

(t1,cpu) = evolve("numpy")
(t2,gpu) = evolve("jax")
print("Max difference CPU-GPU =", np.max(np.abs(cpu-gpu)))

import matplotlib.pyplot as plt
plt.plot(t1, cpu, linewidth=4, label="CPU (NumPy)")
plt.plot(t2, gpu, "--", linewidth=2, label="GPU (JAX)")
plt.xlabel("time")
plt.ylabel("<Sz_0>(t)")
plt.legend()
plt.show()
