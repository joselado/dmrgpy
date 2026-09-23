import numpy as np
from scipy.sparse import linalg as slg
from scipy.sparse import identity
from .tdtk import evolve # evolve the wavefunction
from .. import multioperator
from .edchain import State

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
    else: # use the wavefunction provided
        if type(wf)==State: wf = wf.v # get the vector
    wfA = Aop@wf # apply operator
    wfC = Cop@wf # apply operator
    cs = [] # empty list
    for it in range(nt): # loop
        # Measure *before* evolving, so cs[k] is the t=k*dt value and lines
        # up with the ts grid built above -- same fix, and same reason, as
        # evolution_DC() below and the DMRG backends' own copies of this
        # loop (see timedependent.py's evolution_dmrg_DC()). This function
        # backs evolve_and_measure()/evolution_ABA() on the ED side, which
        # tests/test_time_evolution.py compares directly against the DMRG
        # backends step by step -- so it has to use the identical
        # convention or every such comparison is off by one step.
        c = np.conjugate(wfC)@Bop@wfA # compute braket
        cs.append(c) # store value
        wfA = evolve(wfA,Hop,t=dt,dt=dt) # evolve wavefunction
        wfC = evolve(wfC,Hop,t=dt,dt=dt) # evolve wavefunction
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
    """Special time evolution for the dynamical correlator.

    The operator convention is the DMRG backends' one (see
    timedependent.evolution_dmrg_DC): the pair (A,B) the caller wrote
    means the density sum_n <GS|A|n><n|B|GS>, so B is the operator that
    acts on the ket and A the one that acts on the bra. This route used
    to read the two the other way round and so returned the correlator
    of the *swapped* pair, C[B,A) -- invisible whenever the two operators
    are the same one, which is every example in the documentation, and
    measured on the 2026-09 audit's own seeded 4-site complex-hopping
    chain (A=Cdag_0, B=C_2) as agreeing with C[B,A] to 3.3e-04 while
    mode="DMRG" agrees with C[A,B] to 2.6e-04, i.e. the two solvers
    computed different quantities under one submode name."""
    (A,B) = name[1],name[0] # get the operators, bra first
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
    ht = Hop + e0[0]*identity(Hop.shape[0],dtype=np.complex128)
    for it in range(nt): # loop
        # Measure *before* evolving, so cs[k] is C(k*dt) and lines up with
        # the ts=[0,dt,...,(nt-1)*dt] grid built above -- see the
        # "measure before evolving" comment on timedependent.py's
        # evolution_dmrg_DC() for why the old evolve-then-measure order was
        # a real bug (it dropped C(0), the largest sample, from the Fourier
        # sum, leaving an O(dt) error in every submode="TD" spectrum). Kept
        # in lockstep with the DMRG backends' own copies of this loop so
        # ED stays a valid cross-check of them.
        c = wfc@Bop@wf # store
        cs.append(c) # store value
        wf = evolve(wf,ht,t=dt,dt=dt) # evolve wavefunction
        wf = wf.reshape((wf.shape[0]))
    cs = np.array(cs) # to array
    return ts,cs # return





