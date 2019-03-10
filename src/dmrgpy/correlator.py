from . import operatornames
from . import writemps
import numpy as np


def get_correlator(self,pairs=[[]],name="SS"):
    """Compute a certain static correlator in the spin chain"""
    ########################################
    # workaround for total spin correlator #
    ########################################
    if name=="SS":
        m0 = get_correlator(self,pairs=pairs,name="XX")
        m1 = get_correlator(self,pairs=pairs,name="YY")
        m2 = get_correlator(self,pairs=pairs,name="ZZ")
        return m0+m1+m2
    ###################
    # normal workflow #
    ###################
    self.to_folder() # go to temporal folder
    namei,namej = operatornames.recognize(name) # return that one
    task = {"correlator_operator_i":namei,"correlator_operator_j":namej}
    self.setup_task("correlator",task=task)
    self.write_hamiltonian() # write the Hamiltonian to a file
    writemps.write_correlators(pairs) # write the input file
    self.run() # perform the calculation
    # return the correlators
    m = self.execute(lambda: np.genfromtxt("CORRELATORS.OUT").transpose())
    self.to_origin() # go to main folder
    return m[1]+1j*m[2] # return correlator

