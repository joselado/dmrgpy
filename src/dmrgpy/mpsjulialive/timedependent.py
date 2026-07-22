import numpy as np

from .. import operatornames
from .. import multioperator
from .mpo import MPO


def evolution_dmrg_DC(self,name="XX",nt=10000,dt=0.1,**kwargs):
    """Real-time quench dynamical correlator on the Julia backend, via
    native TDVP (mpsjulialive/tdvp.jl's quench_tdvp -- the whole nt-step
    trajectory runs in one Julia call, same design as kpm.jl). Mirrors
    timedependent.py::evolution_dmrg_DC's TDVP branch (itensor_version=3/
    "python"); only TDVP is implemented here, there is no MPO-Taylor
    fallback for julia_live."""
    name = operatornames.str2MO(self,name,**kwargs)
    name[0] = name[0].get_dagger()
    A,B = name[0],name[1]
    wf0 = self.get_gs() # also sets self.e0
    H = self.hamiltonian
    # shift H by the ground-state energy so psi1's trivial overall phase
    # exp(-i*e0*t) cancels in the <psi2|psi1(t)> correlator, mirroring
    # pyitensor.chain.Chain.quench_tdvp()/quench()
    Hshift_MO = H-self.e0*multioperator.identity()
    Hshift = MPO(Hshift_MO,MBO=self)
    A1 = MPO(A,MBO=self)
    A2 = MPO(B,MBO=self)
    from .juliasession import Main as Mainjl
    correlator,_wf = Mainjl.quench_tdvp(Hshift.jlmpo,
            A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    return ts,cs.real-1j*cs.imag


def evolve_and_measure_dmrg(self,operator=None,nt=1000,h=None,
        dt=1e-2,wf=None,return_wf=False,**kwargs):
    """Real-time evolution + measurement on the Julia backend, via native
    TDVP (mpsjulialive/tdvp.jl's evolve_and_measure_tdvp). Mirrors
    timedependent.py::evolve_and_measure_dmrg's TDVP branch."""
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    Hop = MPO(h,MBO=self)
    Aop = MPO(operator,MBO=self)
    from .juliasession import Main as Mainjl
    correlator,wf_final_jl = Mainjl.evolve_and_measure_tdvp(Hop.jlmpo,
            Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,self.maxm)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    if return_wf:
        from .mps import MPS
        wf_final = MPS(wf_final_jl,MBO=self)
        return ts,cs.real-1j*cs.imag,wf_final
    return ts,cs.real-1j*cs.imag
