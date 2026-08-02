# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
def get(n):
  spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
  sc = spinchain.Spin_Chain(spins) # create the spin chain
  h = 0
  for i in range(n-1):
      h = h +sc.Sx[i]*sc.Sx[i+1]
      h = h +sc.Sy[i]*sc.Sy[i+1]
      h = h +sc.Sz[i]*sc.Sz[i+1]
  sc.set_hamiltonian(h)
  return sc


import time

def compare(n):
  t0 = time.time()
  sc = get(n)
  e0 = sc.gs_energy() # compute the ground state energy
  t1 = time.time()
  dt_julia = None
  try:
      sc = get(n)
      sc.setup_julia() # switch to the live Julia/ITensors.jl backend (mpsjulialive/)
      e1 = sc.gs_energy() # compute the ground state energy
      t2 = time.time()
      dt_julia = t2-t1
      print("Energy with C++",e0)
      print("Energy with Julia",e1)
  except Exception as e:
      # pyjulia/ITensors.jl may simply not be installed in this environment
      # (see install_julia.py); still report the C++ timing in that case.
      print("Julia backend unavailable (%s), skipping Julia timing for n=%d"%(e,n))
  return t1-t0,dt_julia

ns = [10,20,40,80,160,320,640]
ts_cpp,ts_julia = [],[]
f = open("TIMES.OUT","w")
for n in ns:
    dt0,dt1 = compare(n)
    ts_cpp.append(dt0)
    ts_julia.append(dt1)
    f.write(str(n)+"  "+str(dt0)+"  "+str(dt1)+"\n")
    f.flush()
    print("n =",n," Time with C++",dt0," Time with Julia",dt1)
    print()
f.close()

import matplotlib.pyplot as plt

plt.plot(ns,ts_cpp,marker="o",label="C++ (ITensor)")
ns_julia = [n for n,dt in zip(ns,ts_julia) if dt is not None]
dt_julia = [dt for dt in ts_julia if dt is not None]
if len(dt_julia)>0:
    plt.plot(ns_julia,dt_julia,marker="o",label="Julia (ITensors.jl)")
plt.xlabel("Number of sites")
plt.ylabel("Time (s)")
plt.xscale("log") ; plt.yscale("log")
plt.legend()
plt.show()
