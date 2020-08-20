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



