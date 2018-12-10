from __future__ import print_function
import operatornames
import numpy as np
from scipy.sparse import csc_matrix


def evolution(self,name="XX",i=0,j=0,nt=100,dt=0.01):
    """Perform time evolution exactly"""
    h = self.get_full_hamiltonian() # get full Hamiltonian
    from ..pychainwrapper import get_pychain
    sc = get_pychain(self) # the pychain spin object
    namei,namej = operatornames.recognize(self,name) # name of the operator
    nd = {"Sx":0,"Sy":1,"Sz":2}
    ni = nd[namei] # index
    nj = nd[namej]
    Ai = sc.ski[ni][i] # operator
    Aj = sc.ski[nj][j] # operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    from ..pychain.evolution import evolve # evolve the wavefunction
    wf0 = self.get_gs(mode="ED") # get wavefunction
    wf0 = csc_matrix(wf0).T # turn into a matrix
    wf = wf0.copy() # copy wavefunction
    from ..pychain.algebra import Av,braket
    wf = Av(Ai,wf) # apply operator
    cs = [] # empty list
    for it in range(nt): # loop
        wf = evolve([wf],h,t=dt,dt=dt)[0] # evolve wavefunction
        wf = csc_matrix(wf).T # turn into a matrix
        c = braket(wf0,Av(Aj,wf)) # store
        cs.append(c) # store value
    cs = np.array(cs) # to array
    return ts,cs # return


