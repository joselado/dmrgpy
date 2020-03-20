# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results

n = 8 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Chain(spins) # create the spin chain object
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0
    else: return 0.0 # otherwise
sc.set_exchange(fj) # add the exchange couplings
pairs = [(0,i) for i in range(n)] # between the edge and the rest
# compute with the method
cs = sc.get_correlator(pairs=pairs,mode="DMRG",name="XX").real # get the static correlators

# compute by hand
sxs = [sc.get_operator("Sx",i) for i in range(n)]
cs1 = [sc.vev(sxs[p[0]]*sxs[p[1]],mode="ED") for p in pairs]



########################
# Now plot the results #
########################
# now plot the result
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(range(n),cs,marker="o",c="blue",label="DMRG")
plt.scatter(range(n),cs1,marker="o",c="red",s=100,label="ED")
plt.legend()
plt.xlabel("N")
plt.ylabel("<GS|S_0 S_N |GS>")
plt.show()


