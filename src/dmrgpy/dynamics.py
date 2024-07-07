from . import kpmdmrg
from . import timedependent
from . import cvm
from . import dcex

def get_dynamical_correlator(self,submode="KPM",**kwargs):
    if self.itensor_version==2: # C++
        self.set_initial_wf(self.wf0) # set the initial wavefunction
        if not self.is_hermitian(self.hamiltonian): # non Hermitian Hamiltonian
            from .nonhermitian.dynamics import dynamical_correlator_non_hermitian
            return dynamical_correlator_non_hermitian(self,**kwargs)
    #        submode = "CVM_explicit" # only mode that works with non-Hemrmitian
        if submode=="KPM": # KPM method
            return kpmdmrg.get_dynamical_correlator(self,**kwargs)
        elif submode=="TD": # time dependent 
            return timedependent.dynamical_correlator(self,**kwargs)
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
        from .mpsjulialive import dynamics as dynamicsjl
        return dynamicsjl.get_dynamical_correlator(self,submode=submode,**kwargs)
    else: raise



