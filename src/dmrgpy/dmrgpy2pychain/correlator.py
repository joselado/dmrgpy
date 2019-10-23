from .. import pychainwrapper
from .. import operatornames
from ..pychain import correlator as pychaincorrelator
import numpy as np

def correlator(sc,pairs=[[]],name="SS",**kwargs):
    """Compute a static correlator"""
    if name=="SS": # total correlator
        f = lambda n: correlator(sc,pairs=pairs,name=n,**kwargs)
        return f("XX") + f("YY") + f("ZZ")
    scp = pychainwrapper.get_pychain(sc) # get pychain spinchain object
    h = pychainwrapper.get_full_hamiltonian(sc) # get Hamiltonian
    cs = [] # initialize
    namei,namej = operatornames.recognize(name) # get the two operators
    print(namei,namej)
    for (i,j) in pairs:
      c = pychaincorrelator.static(scp,h,namei=namei,namej=namej,
              i=i,j=j,**kwargs)
      cs.append(c)
    cs = np.array(cs) # convert to array
    if np.max(np.abs(cs.imag))<1e-5: cs = cs.real # real part
    return cs # return
