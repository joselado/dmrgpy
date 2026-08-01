# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator computed via real-time evolution of an MPS
# (submode="TD", dynamics.py -> timedependent.dynamical_correlator() ->
# evolution_dmrg_DC()), here backed by TEBD (self.tevol_method="TEBD",
# pyitensor/tebd.py's TEBDEvolver -- itensor_version="python" only, and
# only valid for a strictly nearest-neighbor Hamiltonian, see
# tebd.py/CLAUDE.md) instead of the default TDVP engine. Modeled directly
# on dynamical_correlator_time_evolution_tdvp, with TDVP added alongside
# TEBD as a third curve so all three independent methods (KPM, TD/TDVP,
# TD/TEBD) can be compared on the same plot -- since TDVP's own dynamical
# correlator is already validated elsewhere, close TEBD-vs-TDVP agreement
# here is itself evidence TEBD's quench_tebd() is computing the same
# physical quantity, on top of the ED cross-check done directly in
# tests/test_time_evolution.py (which compares raw C(t), not the
# Fourier-transformed spectrum plotted here).
#
# maxm=30/delta=0.15 (rather than a smaller maxm/delta) are chosen
# deliberately: at a small bond-dimension cap, two independently-
# truncated long real-time evolutions (TDVP's per-step Krylov SVD vs
# TEBD's per-gate SVD) accumulate different truncation error and their
# trajectories diverge over long times -- confirmed directly, maxm=10/
# delta=0.05 (the default used elsewhere in this directory) pushes the
# required simulated time (damping_periods/delta) long enough that the
# two methods' peaks visibly shift apart, even though their raw C(t)
# time series still agree to ~1e-3 over the first few hundred steps.
# This is truncation sensitivity, not a TEBD correctness bug (the ED
# cross-check in tests/test_time_evolution.py already rules that out) --
# raising maxm and delta here keeps this demo in the well-converged
# regime where TDVP and TEBD agree closely.
import time
import numpy as np
from dmrgpy import spinchain

n = 5
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
sc.setup_python() # TEBD only exists on the pure-Python backend
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)

sc.maxm = 30
es = np.linspace(-1.0,10.0,2000)
name = (sc.Sx[0],sc.Sx[0]) # correlator to compute
delta = 0.15

print("Starting")
t0 = time.time()
(x1,y1) = sc.get_dynamical_correlator(es=es,name=name,submode="KPM",delta=delta)
t1 = time.time()

sc.tevol_method = "TDVP"
(x2,y2) = sc.get_dynamical_correlator(es=es,name=name,submode="TD",delta=delta)
t2 = time.time()

sc.tevol_method = "TEBD"
(x3,y3) = sc.get_dynamical_correlator(es=es,name=name,submode="TD",delta=delta)
t3 = time.time()

print("Time in KPM",t1-t0)
print("Time in TD (TDVP)",t2-t1)
print("Time in TD (TEBD)",t3-t2)
print("Max |TD(TDVP) - TD(TEBD)| =",np.max(np.abs(y2-y3)))


### Plot the result ###

import matplotlib.pyplot as plt
plt.plot(x1,y1.real,label="KPM")
plt.plot(x2,np.abs(y2),label="TD (TDVP)")
plt.plot(x3,np.abs(y3),label="TD (TEBD)",linestyle="--")
plt.legend()
plt.xlabel("frequency")
plt.ylabel("Dynamical correlator")
plt.show()
