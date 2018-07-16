from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 10
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

ps = np.linspace(0.,1.0,20) # array
wfs = [] # wavefunctions

eout = []


for p in ps:
  # dimerized coupling
  def fj(i,j):
    out = 0.0
    ij = (i+j)%4
    ddj = -0.2
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
  e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
  print(e0,p)
  es = sc.get_excited(n=2) ; eout.append(es)
  wf = sc.get_gs() # get the ground state as an MPS object
  wfs.append(wf.copy()) # store wavefunction

import matplotlib.pyplot as plt

#exit()


# now perform the product
berry = 1.0+0.0j

for i in range(len(wfs)-1): # loop
  fac = wfs[i].dot(wfs[i+1])
  berry *= fac
#  print(fac)

berry *= wfs[-1].dot(wfs[0])

print("Berry phase",berry)



for (p,e) in zip(ps,eout):
  plt.scatter(e*0.+p,e)

plt.show()


