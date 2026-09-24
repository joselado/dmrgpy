import numpy as np
from scipy.sparse import linalg as slg
from scipy.sparse import identity
from .tdtk import evolve # evolve the wavefunction
from .. import multioperator
from .edchain import State

def evolution_ABC(self,h,A=None,B=None,C=None,wf=None,nt=100,dt=0.01):
    """<wf|C^dagger U(t)^dagger B U(t) A|wf> on the ts grid, with the
    Schrodinger propagator U(t) = e^{-iHt}: A|wf> and C|wf> are both
    evolved forward in time and B is measured between them, so with
    A = C = identity this is <psi(t)|B|psi(t)>, psi(t) = e^{-iHt}|wf>,
    what every DMRG backend's evolve_and_measure returns.

    Both states are advanced with evolve(..., -Hop, ...), because
    tdtk.evolve integrates dpsi/dt = +i*h@psi, i.e. it applies e^{+iht}.
    That sign is what evolution_DC below is built on (its e^{-i w t}
    Fourier kernel needs the e^{+iHt} series), so it stays as it is and
    the minus goes here. This function used to pass Hop itself, so it
    returned the time-reversed trajectory <psi(-t)|B|psi(-t)>, which a
    real Hamiltonian, a real start and a real observable cannot tell
    apart from the forward one, and which is why no ED-versus-DMRG test
    saw it: Larmor precession under B*sum Sz from a +x start came out as
    <Sy_0>(t=1) = -0.4207 against the closed form +sin(1)/2 = +0.4207,
    and one fermion on a 3-site ring with flux pi/6 went round the wrong
    way, <N_1>(t=1.5) = 0.1024 against 0.8413 (2026-09-24b audit,
    finding 9)."""
    nt = int(nt)
    Aop = self.get_operator(A) # get operator
    Bop = self.get_operator(B) # get operator
    Cop = self.get_operator(C) # get operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    Hop = self.get_operator(h) # get Hamiltonian
    if wf is None: # the cached ground state, see evolution_DC below
        wf = self.get_gs_array()
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
        # -Hop: evolve() applies e^{+iht}, so this is e^{-iHt}, forward in
        # time (see the docstring)
        wfA = evolve(wfA,-Hop,t=dt,dt=dt) # evolve wavefunction
        wfC = evolve(wfC,-Hop,t=dt,dt=dt) # evolve wavefunction
    cs = np.array(cs) # to array
    return ts,cs # return

def evolve_and_measure(self,h,operator=None,**kwargs):
    """Evolve and measure"""
    one = multioperator.identity()
    return evolution_ABC(self,h,A=one,B=operator,C=one,**kwargs)

def evolution_ABA(self,h=None,A=None,B=None,**kwargs):
    """Evolve and measure"""
    return evolution_ABC(self,h,A=A,B=B,C=A,**kwargs)




def evolution_DC(self,h=None,name=None,nt=100,dt=0.01,wf0=None,**kwargs):
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
    computed different quantities under one submode name.

    wf0 is the state the correlator is measured in, an array or a State,
    and defaults to the EDchain's cached ground state, the one the KPM,
    INV, CVM and ROOTN submodes use. The energy origin is the cached
    ground-state energy self.e0 either way, as in those submodes. This
    used to re-solve the ground state with a randomly started eigsh on
    every call and ignore both states it had access to, so on a
    degenerate ground state each call measured a different member of the
    manifold, and a pair that is not provably self-adjoint, whose density
    takes two calls (timedependent.lehmann_density_from_one_sided), got
    its two halves from two different members: the density of no state
    and of no mixture, off the cached-state density by up to 3.4e-01 on a
    peak of 1.8e-01 for (Sp_0,Sm_2) on a 3-site Heisenberg chain
    (2026-09-24 audit, finding 8). On a non-degenerate ground state the
    two constructions agree to ~1e-11. Note the sign of the shift below:
    that eigsh ran on -H, so its eigenvalue was -E_0 and the line read
    H + e0*I; with the cached +E_0 it is H - e0*I."""
    (A,B) = name[1],name[0] # get the operators, bra first
    Hop = self.get_operator(h) # return Hamiltonian
    Aop = self.get_operator(A) # return operator
    Bop = self.get_operator(B) # return operator
    ts = np.array([dt*ii for ii in range(nt)]) # times
    gs = self.get_gs_array() # cached; this is also what sets self.e0
    if wf0 is None: wf0 = gs
    elif type(wf0)==State: wf0 = wf0.v # get the vector
    wf0 = np.asarray(wf0).reshape(-1)
    e0 = self.e0 # +E_0, the ground-state energy
    wf = wf0.copy() # copy wavefunction
    wf = Aop@wf # apply operator
    wfc = np.conjugate(wf0) # conjugate wavefunction
    cs = [] # empty list
    ht = Hop - e0*identity(Hop.shape[0],dtype=np.complex128)
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





