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



def reduced_dm_projective(self,wf,i=0,j=None):
    """Compute the reduced density matrix using a brute force approach"""
    def projectors(k):
        """Build the projectors onto the different single-site components"""
        site = self.sites[k] # take this site
        if site==2: # S=1/2 site
            Szk = self.Sz[k] # Sz operator
            P01 = self.Sx[k]+1j*self.Sy[k] # project and rotate
            # we need to project on up/dn and rotate to up!
            return [Szk+0.5,P01] 
        elif site==3: # S=1 site
            raise # not finished
            Szk = self.Sz[k]
            return [Szk*(Szk+1)/2.,-(Szk-1)*(Szk+1),Szk*(1-Szk)/2.] 
        else: raise # not implemented
    Pi = projectors(i) # projectors for site i
    if j is not None: Pj = projectors(j) # projectors for site j
    else: Pj = [1] # workaorund for a single site
    Pk = [] # projectors in the ij subspace
    for pi in Pi:
        for pj in Pj: Pk.append(pi*pj) # store the projector in this subspace
    # now compute the density matrix in this subspace
    n = len(Pk) # number of components
    dm = np.zeros((n,n),dtype=np.complex) # initialize
    for i in range(n):
        wfi = Pk[i]*wf # projector
        for j in range(n):
            dm[i,j] = wfi.aMb(Pk[j],wf) # project
    return dm # return density matrix




