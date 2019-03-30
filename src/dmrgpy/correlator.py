from . import operatornames
from . import writemps
import numpy as np
from . import taskdmrg


def get_correlator(self,pairs=[[]],name="SS"):
    """Compute a certain static correlator in the spin chain"""
    self.get_gs() # compute ground state
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
    namei,namej = operatornames.recognize(name) # return that one
    task = {"correlator":"true",
            "correlator_operator_i":namei,
            "correlator_operator_j":namej
            }
    self.task = task # override tasks
    self.execute( lambda : taskdmrg.write_tasks(self))
    # write the Hamiltonian to a file
    self.execute( lambda: self.write_hamiltonian())
    self.execute( lambda: write_correlators(pairs)) # write the input file
    self.run() # perform the calculation
    # return the correlators
    m = self.execute(lambda: np.genfromtxt("CORRELATORS.OUT").transpose())
    return m[1]+1j*m[2] # return correlator


def write_correlators(pairs):
  """Write the pairs of correlators in a file"""
  fo = open("correlators.in","w") # open the file
  fo.write(str(len(pairs))+"\n") # write in a file
  for p in pairs:
    # the first has to be smaller than the second
    fo.write(str(p[0])+"  "+str(p[1])+"\n")
  fo.close()

