# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)
wf = sc.get_gs()
def f(e):
  wf1 = sc.applyinverse(h -e +1j*1e-1,wf)
  return wf.dot(wf1)
es = np.linspace(-4,4,100)
ds = [-f(e).imag for e in es]
dsz = [f(1j*e) for e in es]
from dmrgpy.analyticcontinuation import imag2real
es1,ds1 = imag2real(1j*es+1j*1e-1,dsz,x=es+1j*1e-1)

e = sc.gs_energy() # compute the ground state energy
print("Energy",e)
import matplotlib.pyplot as plt
plt.plot(es,ds,marker="x")
plt.plot(es1.real,-ds1.imag,marker="o")
plt.show()
