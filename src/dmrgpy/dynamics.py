from . import kpmdmrg
from . import timedependent
from . import cvm
from . import dcex
from . import tdz

def get_dynamical_correlator(self,submode="KPM",**kwargs):
    if self.itensor_version in (2,3,"python"): # C++ or pure-Python
        self.set_initial_wf(self.wf0) # set the initial wavefunction
        if not self.is_hermitian(self.hamiltonian): # non Hermitian Hamiltonian
            from .nonhermitian.dynamics import dynamical_correlator_non_hermitian
            return dynamical_correlator_non_hermitian(self,**kwargs)
    #        submode = "CVM_explicit" # only mode that works with non-Hemrmitian
        if submode=="KPM": # KPM method
            return kpmdmrg.get_dynamical_correlator(self,**kwargs)
        elif submode=="TD": # time dependent
            return timedependent.dynamical_correlator(self,**kwargs)
        elif submode=="TDZ": # complex-time evolution (arXiv:2311.10909)
            return tdz.dynamical_correlator_tdz(self,**kwargs)
        elif submode=="CVM_explicit": # CVM mode
            return cvm.dynamical_correlator_cvm_explicit(self,**kwargs)
        elif submode=="CVM": # CVM mode
            return cvm.dynamical_correlator(self,**kwargs)
        elif submode=="CVMimag": # CVM mode
            return cvm.dynamical_correlator_analytic_continuation(self,**kwargs)
        elif submode=="EX": # EX mode
            return dcex.dynamical_correlator(self,**kwargs)
        elif submode=="maxent": # Max ent mode
            from .distribution import dynamical_correlator_positive_defined
            return dynamical_correlator_positive_defined(self,**kwargs)
        else: raise
    elif self.itensor_version=="julia_live": # Julia version
        # Only KPM/CVM/TDZ assume a Hermitian Hamiltonian (Chebyshev
        # spectrum-in-[-1,1], resolvent CG solve, and the TDZ damping
        # mechanism respectively) -- EX and maxent are already
        # backend-agnostic MultiOperator/MPS algebra with their own
        # working non-Hermitian path (dcex.py -> excited.py's
        # excited_states_non_hermitian, not itensor_version-gated) and
        # must not be blocked here. This check used to run before submode
        # dispatch entirely, which also rejected EX/maxent even though
        # they work fine -- confirmed directly, this was blocking a
        # working code path.
        if submode in ("KPM","CVM","TDZ") and not self.is_hermitian(self.hamiltonian):
            # unlike the (2,3,"python") branch above, there is no
            # non-Hermitian route to fall back to here for these
            # submodes: dynamical_correlator_non_hermitian ultimately
            # needs applyinverse_dmrg(), which is self._session-only
            # (mpsalgebra.py) and also dispatches on type(wf)==mps.MPS --
            # the *top-level* MPS class, not mpsjulialive.mps.MPS -- so it
            # would fail regardless. Silently running the Hermitian-only
            # KPM/CVM/TDZ math on a non-Hermitian Hamiltonian produces
            # numerically wrong output with no error; raise instead.
            raise NotImplementedError(
                "get_dynamical_correlator: itensor_version='julia_live' "
                "does not implement non-Hermitian Hamiltonians for "
                "submode=%r (KPM/CVM/TDZ all assume a Hermitian one); use "
                "submode='EX'/'maxent', or itensor_version in "
                "(2,3,'python') instead"%submode)
        from .mpsjulialive import dynamics as dynamicsjl
        return dynamicsjl.get_dynamical_correlator(self,submode=submode,**kwargs)
    else: raise



