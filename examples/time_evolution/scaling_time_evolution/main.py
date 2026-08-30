# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain  # was a bare "import spinchain", which has
# not resolved since the library became a package -- this script could not
# run at all
import time

# How the cost of a real-time-evolution ("TD") dynamical correlator grows
# with chain length. Everything below the Hamiltonian was dead: the script
# imported `spinchain` bare (unresolvable since the library became a
# package), passed an empty energy grid, and asked for the correlator with
# `use_kpm=False, dt=0.1` -- neither of which the KPM path has ever read,
# so the "time evolution" it claims to time was in fact plain KPM. It also
# overwrote its own size sweep with `ns=[6]` and then fitted a power law
# through that single point.
es = np.linspace(-1.0, 6.0, 400)


def compute(n):
  spins = [2 for i in range(n)]
  sc = spinchain.Spin_Chain(spins) # create the chain
  def fj(i,j):
      if 0.9<abs(i-j)<1.1: return 1.0
      return 0.0
  # set_exchange() was removed; this is the Hamiltonian it built --
  # the coupling summed over BOTH orderings of every pair, so a
  # nearest-neighbour fj=1 gives 2*sum_<ij> S_i.S_j, not sum_<ij>.
  h = 0
  for i in range(n):
    for j in range(n):
      if fj(i,j)!=0.0: h = h + fj(i,j)*sc.SS(i,j)
  sc.set_hamiltonian(h) # set exchange couplings
  #sc.set_fields(lambda x: [0.2,0.2,0.2]) # set exchange couplings
  sc.maxm = 10
  sc.nsweeps = 2
  sc.get_gs()
  t0 = time.time()
  (x2,y2) = sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[1]),
          submode="TD",es=es,dt=0.1,nt=200)
  t1 = time.time()
  print("n =",n," time =",round(t1-t0,3))
  return t1-t0,np.array(x2).real,np.array(y2).real


ns = [4,6,8,10]
results = [compute(n) for n in ns]
ts = np.array([r[0] for r in results])
out = np.polyfit(np.log(ns),np.log(ts),deg=1)
print("Power",out[0])

import matplotlib.pyplot as plt
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4))
ax1.plot(np.log(ns),np.log(ts),marker="o")
ax1.plot(np.log(ns),np.polyval(out,np.log(ns)),ls="--",
        label="slope %.2f"%out[0])
ax1.set_xlabel("log size") ; ax1.set_ylabel("log time") ; ax1.legend()
ax1.set_title("cost of the TD dynamical correlator")
for n,(_t,x,y) in zip(ns,results):
    ax2.plot(x,y,label="n = %d"%n)
ax2.set_xlabel("energy") ; ax2.set_ylabel(r"$S_{zz}(\omega)$") ; ax2.legend()
ax2.set_title("the spectra those timings came from")
plt.tight_layout()
plt.show()











