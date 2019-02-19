from .. import pychainwrapper
import numpy as np
from ..algebra import algebra

def get_magnetization(sc):
    """Compute a static correlator"""
    scp = sc.get_pychain() # get pychain spinchain object
    h = pychainwrapper.get_full_hamiltonian(sc) # get Hamiltonian
    wf = algebra.ground_state(h)[1] # get GS wavefunction 
    mx,my,mz = [],[],[]
    for i in range(sc.ns): # loop over sites
        mx.append(algebra.braket_wAw(wf,scp.sxi[i]).real)
        my.append(algebra.braket_wAw(wf,scp.syi[i]).real)
        mz.append(algebra.braket_wAw(wf,scp.szi[i]).real)
    mx = np.array(mx)
    my = np.array(my)
    mz = np.array(mz)
    return (mx,my,mz)
