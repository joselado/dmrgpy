import numpy as np



def dynamical_correlator_kpm(self):
    """Return the entropies of the dynamical correlator"""
    s = self.execute(lambda: np.genfromtxt("KPM_ENTROPY.OUT"))
    return s



def compute_entropy(self,psi,b=1):
    """Compute entanglement entropy in a bond"""
    if i<1 or i>self.ns: raise
    self.execute(lambda: psi.write(name="wavefunction.mps"))
    # write the task
    task = {    "entropy": "true",
                "bond_entropy":str(b),
                }
    self.task = task # assign tasks
    self.run() # perform the calculation
    entr = self.execute(lambda: np.genfromtxt("ENTROPY.OUT"))
    return entr



