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
        if not self.is_hermitian(self.hamiltonian): # non-Hermitian Hamiltonian
            # unlike the (2,3,"python") branch above, there is no
            # non-Hermitian route to fall back to here: dynamical_correlator_
            # non_hermitian ultimately needs applyinverse_dmrg(), which is
            # self._session-only (mpsalgebra.py) and also dispatches on
            # type(wf)==mps.MPS -- the *top-level* MPS class, not
            # mpsjulialive.mps.MPS -- so it would fail regardless. Silently
            # running the Hermitian-only KPM/CVM/TDZ math on a non-Hermitian
            # Hamiltonian (the previous behavior here) produces numerically
            # wrong output with no error; raise instead.
            raise NotImplementedError(
                "get_dynamical_correlator: itensor_version='julia_live' "
                "does not implement non-Hermitian Hamiltonians (the KPM/"
                "CVM/TDZ submodes all assume a Hermitian one); use "
                "itensor_version in (2,3,'python') instead")
        from .mpsjulialive import dynamics as dynamicsjl
        return dynamicsjl.get_dynamical_correlator(self,submode=submode,**kwargs)
    else: raise



