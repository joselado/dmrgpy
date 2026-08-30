# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 12
spins = [3 for i in range(n)] # spin 1 Haldane chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1): # nearest-neighbor Heisenberg exchange
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 30
sc.nsweeps = 10
sc.kpmmaxm = 10 # KPM max m
fo = open("DCF.OUT","w") # dynamical correlation function
data = [] # (site, xs, ys) for plotting
for i in range(n): # loop over sites
  name = (sc.Sz[i], sc.Sz[i])
  # n= removed: it never reached the KPM recursion (see the audit doc)
  (xs,ys) = sc.get_dynamical_correlator(mode="DMRG",name=name)
  print("Doing",i)
  data.append((xs,ys))
  for (x,y) in zip(xs,ys):
    fo.write(str(i)+"  ")
    fo.write(str(x)+"  ")
    fo.write(str(y)+"\n")
  fo.flush()
fo.close()

# dynamical structure factor <Sz_i(t) Sz_i(0)> as a function of site and
# frequency: a natural site-vs-frequency heatmap
import matplotlib.pyplot as plt
xs0 = data[0][0]
cmap = np.array([np.abs(ys) for (xs,ys) in data])
plt.imshow(cmap,aspect="auto",origin="lower",
        extent=[xs0[0],xs0[-1],0,n-1])
plt.colorbar(label="|DCF(i,w)|")
plt.xlabel("frequency [J]")
plt.ylabel("site i")
plt.show()
