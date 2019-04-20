import numpy as np



def dynamical_correlator_kpm(self):
    """Return the entropies of the dynamical correlator"""
    s = self.execute(lambda: np.genfromtxt("KPM_ENTROPY.OUT"))
    return s
