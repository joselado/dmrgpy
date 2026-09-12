from .manybodychain import Many_Body_Chain
import numpy as np


class Parafermionic_Chain(Many_Body_Chain):
    """Class for a parafermionic chain"""
    def __init__(self,n,Z=3,**kwargs):
        self.Z = Z # type of parafermion
        self.N = [self.get_operator("N",i) for i in range(n)]
        self.Sig = [self.get_operator("Sig",i) for i in range(n)]
        self.Sigd = [self.get_operator("SigDag",i) for i in range(n)]
        self.Tau = [self.get_operator("Tau",i) for i in range(n)]
        self.Taud = [self.get_operator("TauDag",i) for i in range(n)]
        self.Id = self.get_operator("Id",1)
        # **kwargs (itensor_version=, ...) is forwarded at the
        # Many_Body_Chain.__init__ call site rather than at the top of
        # this method, because the get_operator() calls above have to run
        # before it
        if Z==3: Many_Body_Chain.__init__(self,[-2 for i in range(n)],**kwargs)
        elif Z==4: Many_Body_Chain.__init__(self,[-3 for i in range(n)],**kwargs)
        elif Z==2: 
             print("Not with DMRG!")
             Many_Body_Chain.__init__(self,[-1 for i in range(n)],**kwargs)
        else: raise ValueError("only Z=2,3,4 parafermions are implemented, got "+repr(Z))
        self.use_ampo_hamiltonian = True # use ampo
        self.Chi = []
        self.Dis = []
        self.Psi = []
        for i in range(n): 
            t = 1
            for j in range(i): t = t*self.Tau[j]
            self.Chi.append(t*self.Sig[i])
            self.Psi.append(t*self.Sig[i]*self.Tau[i])
        for i in range(n): 
            t = 1
            for j in range(i+1): t = t*self.Tau[j]
            self.Dis.append(t)
        self.Psid = [o.get_dagger() for o in self.Psi]
        self.Chid = [o.get_dagger() for o in self.Chi]
    def get_ED_obj(self):
        """Return the associated ED object"""
        # Same has_ED_obj/ED_obj caching protocol every other chain class
        # uses (see Fermionic_Chain.get_ED_obj): restart() clears the
        # flag, and set_hamiltonian() calls restart(), so the cached
        # object can never outlive the Hamiltonian it was built for.
        # Without it every single mode="ED" call rebuilt the full sparse
        # operator dictionary and redid the ground-state solve.
        if self.has_ED_obj: # if the ED object has been computed
            return self.ED_obj # return the stored object
        from .pyparafermion import parafermion
        obj = parafermion.Parafermion_Chain(self)
        self._apply_sector_to_ed(obj) # conserved sector, if any (there is
            # none for parafermions -- get_sector_charge_operators raises
            # -- so this is a no-op kept for consistency with the others)
        self.ED_obj = obj # store the object
        self.has_ED_obj = True # set to True
        return self.ED_obj
    # get_dynamical_correlator used to be overridden here, as
    #   def get_dynamical_correlator(self,mode="DMRG",**kwargs):
    #       if mode=="DMRG": return super().get_dynamical_correlator_MB(**kwargs)
    #       elif mode=="ED": return self.get_ED_obj().get_dynamical_correlator(**kwargs)
    # i.e. the two branches of Many_Body_Chain.get_dynamical_correlator
    # verbatim (get_dynamical_correlator_MB *is* dynamics.
    # get_dynamical_correlator, which is what the base class's DMRG branch
    # calls) minus the two things that method does before them, and it was
    # the only model-class override in the tree skipping either:
    #   - `mode = self.get_mode(mode=mode)`, so an enforced self.mode="ED"
    #     -- and every DMRG->ED fallback mode.py makes on its own (no
    #     compiled extension, itensor_version=3 with ns<3) -- was ignored
    #     and the call went to DMRG anyway. On itensor_version=3 that
    #     reached Chain::kpm_dynamical_correlator with have_H_ false,
    #     where ITensor's Error() calls abort(): the user's whole process
    #     died with SIGABRT, uncatchable from Python (2026-09 audit #12).
    #   - resolution of the documented *string* form of name= (via
    #     operatornames.str2MO), so name="ZZ" died several frames deep
    #     inside EDOperator with "takes a MultiOperator or another
    #     EDOperator, got str" instead of being resolved, or -- for an
    #     operator family a parafermion chain genuinely does not have --
    #     refused by name with a ValueError.
    # Deleting it fixes both at once; there was nothing else in it.
    def test(self,**kwargs):
        return test_commutation(self,**kwargs)








def test_commutation(self,ntries=3):
    """Perform a test of the commutation relations.

    ntries random site pairs are drawn (pairs with i>=j are skipped, so
    the effective number of checks is smaller); it used to be hardcoded to
    3, while Parafermionic_Chain.test() forwarded **kwargs here and so
    died with a TypeError on the ntries= the sibling Spin_Chain.test()
    accepts."""
    Chi = self.Chi
    Chid = self.Chid
    Tau = self.Tau
    Taud = self.Taud
    Psi = self.Psi
    Psid = self.Psid
    n = len(Chi) # number of sites
    omega = np.exp(1j*2*np.pi/self.Z)
    for ii in range(ntries):
        i = np.random.randint(n)
        j = np.random.randint(n)
        if i>=j: continue
        # each bare `raise` below used to report "RuntimeError: No active
        # exception to reraise", so the diagnostic existed only in the
        # print above it and never in the exception itself
        def fail(what):
            raise AssertionError("parafermion commutation test failed: "
                    +what+" do not satisfy the Z%d relation at sites "
                    "(%d,%d)"%(self.Z,i,j))
        d = Chi[i]*Chi[j] - omega*Chi[j]*Chi[i]
        if not self.is_zero_operator(d): fail("Chi,Chi")
        d = Psi[i]*Psi[j] - omega*Psi[j]*Psi[i]
        if not self.is_zero_operator(d): fail("Psi,Psi")
        d = Chi[i]*Psi[j] - omega*Psi[j]*Chi[i]
        if not self.is_zero_operator(d): fail("Psi,Chi")
    print("Commutation test passed")

