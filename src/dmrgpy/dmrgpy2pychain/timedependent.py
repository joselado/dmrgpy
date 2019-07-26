from __future__ import print_function
from .. import operatornames
import numpy as np
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as slg
from scipy.sparse import identity

def evolution(self,name="XX",i=0,j=0,nt=100,dt=0.01):
    """Perform time evolution exactly"""
    h = self.get_full_hamiltonian() # get full Hamiltonian
    from ..pychainwrapper import get_pychain
    sc = get_pychain(self) # the pychain spin object
    namei,namej = operatornames.recognize(name) # name of the operator
    nd = {"Sx":0,"Sy":1,"Sz":2}
    ni = nd[namei] # index
    nj = nd[namej]
    Ai = sc.ski[ni][i] # operator
    Aj = sc.ski[nj][j] # operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    from ..pychain.evolution import evolve # evolve the wavefunction
    h = self.get_full_hamiltonian() # get Hamiltonian
    e0,wf0 = slg.eigsh(-h,k=1,ncv=20,which="LA")
    wf0 = wf0.reshape(wf0.shape[0])
    wf = wf0.copy() # copy wavefunction
    wf = Ai@wf # apply operator
    cs = [] # empty list
    ht = h + e0[0]*identity(h.shape[0],dtype=np.complex)
    from ..algebra.algebra import braket_wAw
    for it in range(nt): # loop
        wf = evolve([wf],-ht,t=dt,dt=dt)[0] # evolve wavefunction
        wf = wf.reshape((wf.shape[0]))
        c = braket_wAw(wf0,Aj,wi=wf) # store
        cs.append(c) # store value
    cs = np.array(cs) # to array
    return ts,cs # return


