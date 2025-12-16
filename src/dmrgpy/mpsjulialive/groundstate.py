from .juliasession import Main as Mainjl

def get_gs_dmrg(self,ishermitian=True,wf0=None):
    """Compute the ground state using DMRG"""
    if self.computed_gs: # if stored, just return
        return self.wf0 # return wavefunction
    else: # if not computed
        H = self.toMPO(self.hamiltonian) # get the Hamiltonian
        if wf0 is None:
            psi0 = self.random_state() # random state
        else: psi0 = wf0 # from input
        e0,wf0 = Mainjl.get_gs_dmrg(H.jlmpo,psi0.jlmps,nsweeps=self.nsweeps,
                cutoff=self.cutoff,maxm=self.maxm,ishermitian=ishermitian)
        from .mps import MPS
        WF = MPS(wf0,MBO=self) 
        self.wf0 = WF # store wavefunction
        self.e0 = e0 # store energy
        self.computed_gs = True # assume that the ground state is computed
        return e0,WF # return energy and wavefunction




