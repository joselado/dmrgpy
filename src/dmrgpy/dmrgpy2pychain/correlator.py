from .. import pychainwrapper
from ..pychain import correlator as pychaincorrelator
import numpy as np

def correlator(sc,pairs=[[]]):
    """Compute a static correlator"""
    scp = pychainwrapper.get_pychain(sc) # get pychain spinchain object
    h = pychainwrapper.get_full_hamiltonian(sc) # get Hamiltonian
    cs = [] # initialize
    for (i,j) in pairs:
      c = 0.0
      c += pychaincorrelator.static(scp,h,namei="X",namej="X",i=i,j=j)
      c += pychaincorrelator.static(scp,h,namei="Y",namej="Y",i=i,j=j)
      c += pychaincorrelator.static(scp,h,namei="Z",namej="Z",i=i,j=j)
      cs.append(c)
    return np.array(cs).real # return
