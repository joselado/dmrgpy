# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
import numpy as np

### This script will compute the dynamical correlator for a
# S=1/2 chain with an S=1 impurity in the middle
# and save it to the file CORRELATOR_VS_SITE.OUT


from dmrgpy import spinchain
n = 9 # total number of sites (select an odd number)
# create an S=1/2 spin chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain (2 means S=1/2)

# Comment the following line to have a pristine system
spins[n//2] = 3 # S=1 in the middle, (3 means S=1)


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
def fj(i,j): # first neighbor coupling
  if abs(i-j)==1: # first neighbors
      return 1.0
  else: return 0.0


sc.set_exchange(fj) # add the exchange couplings to the Hamiltonian

sc.kpmmaxm = 20 # KPM maxm
xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
for i in range(n): # loop over sites
  print("Doing ",i)
  sc.kpmmaxm = 10 # bond dimension of the DMRG-KPM
  # compute the dynamical correlator
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",i=i,j=i,name="ZZ",
          es=np.linspace(-0.5,4.0,200),delta=0.1)
  zs.append(s.real) # store data


# prepare the data for saving and plotting
xs = range(n)
ys = e
zs = np.array(zs).transpose()

# write result in a file
f = open("CORRELATOR_VS_SITE.OUT","w")
f.write("# site, energy, correlator\n")
for i in range(len(xs)):
  for j in range(len(ys)):
      f.write(str(xs[i])+"  ")
      f.write(str(ys[j])+"  ")
      f.write(str(zs[j,i])+"\n")
f.close()

# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
matplotlib.rcParams.update({'font.size': 18})
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.contourf(xs,ys,zs,500,cmap = plt.get_cmap("hot"),vmax=np.percentile(zs,95))
plt.ylabel("frequency [J]")
plt.xlabel("Site")
plt.show()












