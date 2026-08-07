# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

n = 60 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

wf0 = sc.random_mps()

h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

Op = h # operator
from dmrgpy.multioperatortk.staticoperator import StaticOperator
SOp = StaticOperator(Op,sc) # create a static operator

import time
dots = {} # store the trajectory of <wf|wf> for each mode
times = {} # store the total time spent for each mode
for (O,mode) in [(SOp,"static"),(Op,"default")]:
    wf = wf0.copy()
    ds = []
    t0 = time.time()
    for i in range(50):
        d = wf.dot(wf)
        print(d)
        ds.append(d.real)
        wf = O*wf
    t1 = time.time()
    print("Time in",mode,t1-t0)
    dots[mode] = ds
    times[mode] = t1-t0

import matplotlib.pyplot as plt

plt.subplot(1,2,1)
for mode in dots:
    plt.plot(dots[mode],marker="o",label=mode)
plt.xlabel("Iteration")
plt.ylabel(r"$\langle \psi | \psi \rangle$")
plt.legend()

plt.subplot(1,2,2)
plt.bar(list(times.keys()),list(times.values()))
plt.ylabel("Time (s)")

plt.tight_layout()
plt.show()
