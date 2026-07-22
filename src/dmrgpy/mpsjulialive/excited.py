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
    if self.excited_gram_schmidt:
        # mirrors mpscpp3/chain_session.h's Chain::excited_states and
        # pyitensor.chain.Chain.excited_states, both of which apply
        # gram-schmidt to wfs *before* computing the returned energies
        # (previously silently skipped here -- excited.jl's own
        # excited_states_dmrg has no gram-schmidt step at all, unlike the
        # C++/pyitensor backends' self._session.excited_states(n,scale,
        # self.excited_gram_schmidt)). gram_smith() only uses generic
        # MPS algebra (.copy()/.normalize()/.dot()), already julia_live-safe.
        from ..algebra.arnolditk import gram_smith
        wfs = gram_smith(wfs)
        energies = np.array([wfi.dot(Hop*wfi).real for wfi in wfs])
    else:
        energies = np.array([complex(e).real for e in energies_jl])
    return energies,wfs
