# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
#sc.itensor_version = "julia"

h = h + h.get_dagger()
sc.set_hamiltonian(h) # set Hamiltonian

def get(ts,mode):
    """Function to evolve and perform overlap
    this will compute
    <WF|e^{iHt}|WF>"""
    sc.set_hamiltonian(h) # set Hamiltonian
    wf = sc.get_gs(mode=mode) # get the ground state
    wf = sc.Sx[0]*wf # apply a flip operator to not have an eigenstate
    wf = wf.normalize() # normalize the initial wavefunction
    wf0 = wf.copy() # copy the initial state
    from dmrgpy.timeevolution import imaginary_exponential # function to perform t-evol
    wfs = imaginary_exponential(h,wf,ts=ts) # this computes e^{iht}*wf for each t 
    # if wf was a conventional array and h a matrix, the previous line would be like
    # wfs = [scipy.linalg.expm(1j*t*h)@wf for t in ts] 
    out = [wf0.dot(wfi) for wfi in wfs] # compute overlap
    out = np.array(out).real # take the real part of the time evolution
    return out # return result



ts = np.linspace(0.0,10.0,100) # times to compute
# do DMRG
us0 = get(ts,"DMRG")
# do ED
us1 = get(ts,"ED")

import matplotlib.pyplot as plt # and plot a comparison
plt.scatter(ts,us0,label="DMRG",marker="o",s=200,c="blue")
plt.scatter(ts,us1,label="ED",marker="o",s=60,c="red")
plt.legend()
plt.show()











