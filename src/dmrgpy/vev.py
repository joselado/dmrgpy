# compute vacuum expectation values using multioperators

from . import multioperator

import numpy as np

def multi_vev(self,MO,excited=False,n=4,scale=10.0):
    """
    Compute a VEV using multioperators
    """
    MO = multioperator.obj2MO(MO,name="vev_multioperator")
    if MO.name!="vev_multioperator": raise
    self.get_gs()
    taskd = MO.get_dict() # get the dictionary
    if excited: 
        self.task["excited_vev"] = "true" # do a VEV
        self.task["nexcited"] = n # do a VEV
        self.task["scale_lagrange_excited"] = scale
    else: self.task["vev"] = "true" # do a VEV
    self.task["GS"] = "false" # do a VEV
    self.write_task() # write the tasks in a file
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.execute(lambda: MO.write()) # write multioperator
    self.run() # perform the calculation
    m = self.execute(lambda: np.genfromtxt("VEV.OUT"))
    if excited: m = m.T
    return m[0]+1j*m[1] # return result


def vev(*args):
    return multi_vev(*args,excited=False)


def excited_vev(*args,**kwargs):
    return multi_vev(*args,excited=True,**kwargs)












