from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain

n = 20
sc = fermionchain.Fermionic_Hamiltonian(n) # create the fermionic chain

ps = np.linspace(0.,1.0,10) # array
wfs = [] # wavefunctions

eout = []





def get_sc(p):
  def ft(i,j):
      dt = -0.3
      if j==(i+1): return 1.0 + dt*(-1)**i # NN
      if i==(j+1): return 1.0 +dt*(-1)**j # NN
      if j==0 and i==(n-1): return np.exp(1j*p*np.pi*2)*(1.+dt*(-1)**i) # TBC
      if i==0 and j==(n-1): return np.exp(-1j*p*np.pi*2)*(1.+dt*(-1)**j) # TBC
      if i==j: return 5.0
      return 0.0
  
  sc.set_hoppings(ft) # set hopping
#  sc.nsweeps = 30
  sc.set_fields(lambda i: [0.,0.,10.0]) # set zeeman
  return sc


import topology

print("phi",topology.berry_phase(get_sc,nk=10))
exit()






for p in ps:
  # dimerized coupling
#  exit()
#  print(e0)
  sc = get_sc(p)
  es = sc.get_excited(n=2) ; eout.append(es)
  print(es)
  wf = sc.get_gs() # get the ground state as an MPS object
  wfs.append(wf.copy()) # store wavefunction

import matplotlib.pyplot as plt

#exit()


# now perform the product
berry = 1.0+0.0j

#for i in range(1,len(wfs)): # gauge to the first one
#  wfs[i] = wfs[i]*np.exp(-1j*np.angle(wfs[0].dot(wfs[i])))

facs = []

for i in range(len(wfs)-1): # loop
  fac = wfs[i].dot(wfs[i+1])
  berry *= fac
  print(np.abs(fac),np.angle(fac))
  facs.append(fac)

facs = np.array(facs)

#plt.scatter(range(len(facs)),facs.imag)
#plt.show()
#exit()
#berry *= wfs[-1].dot(wfs[0])

print("Berry phase",np.arctan2(berry.imag,berry.real)/np.pi)
print(berry)



for (p,e) in zip(ps,eout):
  plt.scatter(e*0.+p,e)

plt.show()

sc.clean()
