import numpy as np
from . import operatornames
from . import taskdmrg
from . import multioperator

def dynamical_correlator(self,es=np.linspace(0.,10.0,100),
        delta=1e-1,name="XX",i=0,j=0):
    """
    Compute the dynamical correlator using CVM method in DMRG
    """
    if not self.computed_gs: self.get_gs() # compute ground state
    out = [] # empty list
    for e in es: # loop over energies
        print("CVM in E = ",e)
        o = cvm_dmrg(self,name=name,i=i,j=j,delta=delta,e=e)
        out.append(o) # store
    return (es,np.array(out)) # return result







def cvm_dmrg(self,name="XX",delta=1e-1,e=0.0,**kwargs):
    """
    Return the dynamical correlator for a single energy
    """
    name = operatornames.str2MO(self,name,**kwargs)
    if self.fit_td: fittd = "true"
    else: fittd = "false"
    fittd = "true"
    task = {"cvm":"true",
            "cvm_delta":str(delta),
            "cvm_energy":str(e),
            "cvm_e0":str(self.e0),
            "cvm_nit":str(int(self.cvm_nit)),
            "cvm_tol":str(self.cvm_tol),
            }
    self.task = task # override tasks
#    name[0] = name[0].get_dagger()
    self.execute(lambda: name[0].write(name="dc_multioperator_i.in"))
    self.execute(lambda: name[1].write(name="dc_multioperator_j.in"))
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("CVM.OUT") # read the correlator
    return cs[0] + 1j*cs[1] # return the correlator

