from .juliasession import Main as Mainjl

#NH_biorthoalg = "biorthoblock"
NH_biorthoalg = "fidelity"
NH_alg="onesided"
NH_dmrg = True

def get_gs_dmrg(self,ishermitian=True,wf0=None):
    """Compute the ground state using DMRG"""
    if self.computed_gs: # if stored, just return
        return self.wf0 # return wavefunction
    else: # if not computed
        H = self.toMPO(self.hamiltonian) # get the Hamiltonian
        if wf0 is None:
            psi0 = self.random_state() # random state
        else: psi0 = wf0 # from input
        # technically itensor can also run with nonhermitian=False,
        # but by default we will use the NH routine
        if ishermitian: use_dmrg = True # use conventional DMRG
        else: # non-Hermitian Hamiltonian
            if NH_dmrg: use_dmrg = False # use the NH version
            else: use_dmrg = True # use conventional DMRG
        if use_dmrg: # for Hermitian Hamiltonians, usual DMRG
            e0,wf0 = Mainjl.get_gs_dmrg(H.jlmpo,psi0.jlmps,nsweeps=self.nsweeps,
                cutoff=self.cutoff,maxm=self.maxm,ishermitian=ishermitian)
        else: # for non-Hermitian, the specialized routine
            e0,wfl0,wfr0 = Mainjl.get_gs_nhdmrg(H.jlmpo,psi0.jlmps,
                nsweeps=self.nsweeps,
                biorthoalg = NH_biorthoalg,
                alg = NH_alg,
                cutoff=self.cutoff,maxm=self.maxm)
            wf0 = wfr0 # just the right one
        from .mps import MPS
        WF = MPS(wf0,MBO=self) 
        self.wf0 = WF # store wavefunction
        self.e0 = e0 # store energy
        self.computed_gs = True # assume that the ground state is computed
        return e0,WF # return energy and wavefunction




