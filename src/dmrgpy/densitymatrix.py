# routines to compute density matrices
from . import taskdmrg
import numpy as np

def reduced_dm(self,i=0):
    """
    Compute the reduced density matrix
    """
    self.get_gs() # compute ground state
    task = {"density_matrix":"true",
            "index_i_DM":str(i+1)
#            "index_j_DM":str(j+1)
            }
    self.task = task # add the tasks
    self.write_task() # write the tasks
    self.run() # run the calculation
    out = self.execute(lambda: np.genfromtxt("DM.OUT").T) # read density matrix
    m = out[0] + out[1]*1j # get matrix
    n = int(np.sqrt(m.shape[0]))
    m = m.reshape((n,n)) # transform to matrix
    return m # return the density matrix



