import numpy as np



def dynamical_correlator_kpm(self):
    """Return the entropies of the dynamical correlator"""
    s = self.execute(lambda: np.genfromtxt("KPM_ENTROPY.OUT"))
    return s

from statistics import geometric_mean

def gmean(x):
    if np.max(x)<1e-10: return 0.0
    else: return geometric_mean(x)


def compute_entropy(self,psi,b=1):
    if b is None:
        out = np.array([compute_entropy_single(self,psi,b=i) 
            for i in range(1,self.ns)])
        return gmean(out)
    else: return compute_entropy_single(self,psi,b=b)


def compute_entropy_single(self,psi,b=1):
    """Compute entanglement entropy in a bond"""
    if b<1 or b>self.ns: raise
    self.execute(lambda: psi.write(name="wavefunction.mps"))
    # write the task
    task = {    "entropy": "true",
                "bond_entropy":str(b),
                }
    self.task = task # assign tasks
    self.run() # perform the calculation
    entr = self.execute(lambda: np.genfromtxt("ENTROPY.OUT"))
    return np.abs(entr)



def bond_entropy(self,wf,i,j):
    """Compute the entropy of a state in bond i"""
    if abs(i-j)==1: # use the DMRG approach
        return compute_entropy_single(self,wf,b=max([i,j]))
    else: raise
#    from .densitymatrix import reduced_dm_projective
#    dm = reduced_dm_projective(self,wf,i=i,j=j) # compute density matrix
#    return entropy_dm(dm,normalize=True) # return the entropy


def pair_entropy(self,wf,i,j):
    """Compute the entropy of a state in bond i"""
    from .densitymatrix import reduced_dm_projective
    dm = reduced_dm_projective(self,wf,i=i,j=j) # compute density matrix
    return entropy_dm(dm,normalize=False) # return the entropy


def site_entropy(self,wf,i):
    """Compute the entropy of a state in bond i"""
    from .densitymatrix import reduced_dm_projective
    dm = reduced_dm_projective(self,wf,i=i,j=None) # compute density matrix
    return entropy_dm(dm,normalize=True) # return the entropy


def entropy_dm(dm,normalize=False):
    from scipy.linalg import eigvalsh
    if np.abs(1.-np.trace(dm))>1e-3: raise
    ds = eigvalsh(dm) # compute eigenvalues
    ds = ds[ds>1e-6]
    if normalize: 
        n = dm.shape[0]
        norm = np.log(n) # normalize
    else: norm = 1.
    return -np.sum(ds*np.log(ds))/norm # return the entropy


    



