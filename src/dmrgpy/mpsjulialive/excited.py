import numpy as np

from .mpo import MPO


def get_excited_states_dmrg(self,n=2,noise=0.0,scale=10.0):
    """Excited states on the Julia backend, via ITensorMPS.jl's own
    orthogonality-penalty dmrg() (mpsjulialive/excited.jl's
    excited_states_dmrg -- the whole n-state loop runs in one Julia
    call). Mirrors excited.py::get_excited_states_dmrg's C++ path
    (mpscpp2/mpscpp3's Chain::excited_states) and
    pyitensor.chain.Chain.excited_states()."""
    wf0 = self.get_gs()
    H = self.hamiltonian
    e0 = self.e0
    from .dynamics import _max_energy_bound
    emax = _max_energy_bound(self,H) # mutates self.wf0/self.e0, see below
    self.wf0,self.e0,self.computed_gs = wf0,e0,True # restore GS cache
    weight = (emax-e0)*scale
    Hop = MPO(H,MBO=self)
    from .juliasession import Main as Mainjl
    energies_jl,wfs_jl = Mainjl.excited_states_dmrg(Hop.jlmpo,wf0.jlmps,
            int(n),weight,self.jlsites,self.nsweeps,self.cutoff,self.maxm)
    from .mps import MPS
    wfs = [MPS(w,MBO=self) for w in wfs_jl]
    energies = np.array([complex(e).real for e in energies_jl])
    return energies,wfs
