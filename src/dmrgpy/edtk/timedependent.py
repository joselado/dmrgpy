import numpy as np
from scipy.sparse import linalg as slg
from scipy.sparse import identity
from .tdtk import evolve # evolve the wavefunction
from .. import multioperator

def evolution_ABC(self,h,A=None,B=None,C=None,wf=None,nt=100,dt=0.01):
    """Aply operator C, evolve, apply operator B, evolve back,
    apply operator A <AU-1BUC>"""
    nt = int(nt)
    Aop = self.get_operator(A) # get operator
    Bop = self.get_operator(B) # get operator
    Cop = self.get_operator(C) # get operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    Hop = self.get_operator(h) # get Hamiltonian
    if wf is None: 
        e0,wf = slg.eigsh(-Hop,k=1,ncv=20,which="LA")
        wf = wf.reshape(wf.shape[0])
    wfA = Aop@wf # apply operator
    wfC = Cop@wf # apply operator
    cs = [] # empty list
    for it in range(nt): # loop
        wfA = evolve(wfA,Hop,t=dt,dt=dt) # evolve wavefunction
        wfC = evolve(wfC,Hop,t=dt,dt=dt) # evolve wavefunction
        c = np.conjugate(wfC)@Bop@wfA # compute braket
        cs.append(c) # store value
    cs = np.array(cs) # to array
    return ts,cs # return

def evolve_and_measure(self,h,operator=None,**kwargs):
    """Evolve and measure"""
    one = multioperator.identity()
    return evolution_ABC(self,h,A=one,B=operator,C=one,**kwargs)

def evolution_ABA(self,h=None,A=None,B=None,**kwargs):
    """Evolve and measure"""
    return evolution_ABC(self,h,A=A,B=B,C=A,**kwargs)




def evolution_DC(self,h=None,name=None,nt=100,dt=0.01,**kwargs):
    """Special time evolution for the dynamical correlator"""
    (A,B) = name[0],name[1] # get the operators
    Hop = self.get_operator(h) # return Hamiltonian
    Aop = self.get_operator(A) # return operator
    Bop = self.get_operator(B) # return operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    e0,wf0 = slg.eigsh(-Hop,k=1,ncv=20,which="LA")
    wf0 = wf0.reshape(wf0.shape[0])
    wf = wf0.copy() # copy wavefunction
    wf = Aop@wf # apply operator
    wfc = np.conjugate(wf0) # conjugate wavefunction
    cs = [] # empty list
    ht = Hop + e0[0]*identity(Hop.shape[0],dtype=np.complex)
    for it in range(nt): # loop
        wf = evolve(wf,ht,t=dt,dt=dt) # evolve wavefunction
        wf = wf.reshape((wf.shape[0]))
        c = wfc@Bop@wf # store
        cs.append(c) # store value
    cs = np.array(cs) # to array
    return ts,cs # return





