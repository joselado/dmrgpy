from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 4
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

ps = np.linspace(0.,2.0,40) # array
wfs = [] # wavefunctions

eout = []

fo = open("SWEEP.OUT","w")

for p in ps:
  # dimerized coupling
  def fj(i,j):
    out = 0.0
    ij = (i+j)%4
    ddj = -0.5
    if ij==1: dj = -ddj
    else: dj = ddj
    if i==j+1: out = 1.0 + dj
    # rotate in the case of the last link
    if i==0 and j==(n-1):
      cp = np.cos(p*np.pi)
      sp = np.sin(p*np.pi)
      out = np.matrix([[cp,sp,0.],[-sp,cp,0.],[0.,0.,1.]])*(1.0+dj)
#      print(out,"Here")
    return out
  
  sc.set_exchange(fj) # set those exchange couplings
#  es = sc.get_excited(n=10) ; eout.append(es)
  (x,y) = sc.get_dynamical_correlator(n=300,mode="DMRG",i=0,j=0)
  for (ix,iy) in zip(x,y):
    fo.write(str(p)+"  ")
    fo.write(str(ix)+"  ")
    fo.write(str(iy)+"\n")
  fo.flush()
  print(p)

fo.close()
