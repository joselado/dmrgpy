# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')
#------------------------------------------------------------------
import numpy as np
from dmrgpy import spinchain

# Heisenberg Hamiltonian S12
n = 4 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
h = h + h.get_dagger()
sc.set_hamiltonian(h) # set Hamiltonian

# set default initial state as ground state
wf = sc.get_gs() # get the ground state
wf = sc.Sx[0]*wf 
wf = wf.normalize() 
wf0 = wf.copy()

# change second value for different times

def get_wf(h, ts, wf_initial=wf0):
    ts = np.linspace(0.0, 1, 1) 
    from dmrgpy.timeevolution import imaginary_exponential # function to perform t-evol
    wfs = imaginary_exponential(h,wf_initial,ts=ts) # this computes e^{iht}*wf for each t 
    wf_final = wfs[1]
    return wf_final 

wf_final = get_wf(h, ts, wf_initial=wf0)










