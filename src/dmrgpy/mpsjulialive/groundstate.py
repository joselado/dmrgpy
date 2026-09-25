from .juliasession import Main as Mainjl

#NH_biorthoalg = "biorthoblock"
NH_biorthoalg = "fidelity"
NH_alg="onesided"
NH_dmrg = True

def get_gs_dmrg(self,ishermitian=True,wf0=None):
    """Compute the ground state using DMRG, from wf0 when given and from a
    random state otherwise. Always solves: whether a stored state answers
    the call is decided by the caller (Many_Body_Chain.gs_energy() through
    groundstate.gs_is_current(), then groundstate._gs_energy_julia()).
    It used to return the stored wf0 alone when computed_gs was set, which
    its only caller unpacked as (e0,wf0), so gs_energy(wf0=x) on a solved
    chain raised instead of sweeping from x."""
    H = self.toMPO(self.hamiltonian) # get the Hamiltonian
    if wf0 is None:
        psi0 = self.random_state().jlmps # random state
    else: # from input
        # An MPS made by operator algebra (MPO.__mul__, i.e. any state a
        # caller builds as A*psi) carries a stale prime level on its Link
        # indices, which dmrg()'s environments cannot take (a
        # DimensionMismatch inside eigsolve, for a doublet member built as
        # (Sz_tot+1/2)*psi); tdvp.jl's tdvp_step clears it the same way.
        # On a copy, so the caller's MPS is left as it is.
        psi0 = Mainjl.noprime(Mainjl.copy(wf0.jlmps),"Link")
    # technically itensor can also run with nonhermitian=False,
    # but by default we will use the NH routine
    if ishermitian: use_dmrg = True # use conventional DMRG
    else: # non-Hermitian Hamiltonian
        if NH_dmrg: use_dmrg = False # use the NH version
        else: use_dmrg = True # use conventional DMRG
    if use_dmrg: # for Hermitian Hamiltonians, usual DMRG
        e0,wf0 = Mainjl.get_gs_dmrg(H.jlmpo,psi0,nsweeps=self.nsweeps,
            cutoff=self.cutoff,maxm=self.maxm,ishermitian=ishermitian)
    else: # for non-Hermitian, the specialized routine
        e0,wfl0,wfr0 = Mainjl.get_gs_nhdmrg(H.jlmpo,psi0,
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




