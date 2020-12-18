from . import kpmdmrg
from . import timedependent
from . import cvm
from . import dcex

def get_dynamical_correlator(self,submode="KPM",**kwargs):
    self.set_initial_wf(self.wf0) # set the initial wavefunction
    if submode=="KPM": # KPM method
        return kpmdmrg.get_dynamical_correlator(self,**kwargs)
    elif submode=="TD": # time dependent 
        return timedependent.dynamical_correlator(self,**kwargs)
    elif submode=="CVM": # CVM mode
        return cvm.dynamical_correlator(self,**kwargs)
    elif submode=="CVMimag": # CVM mode
        return cvm.dynamical_correlator_analytic_continuation(self,**kwargs)
    elif submode=="EX": # EX mode
        return dcex.dynamical_correlator(self,**kwargs)
    else: raise



