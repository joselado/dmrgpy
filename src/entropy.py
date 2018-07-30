import numpy as np



def dynamical_correlator_kpm(self):
    """Return the entropies of the dynamical correlator"""
    self.to_folder()
    s = np.genfromtxt("KPM_ENTROPY.OUT")
    self.to_origin()
    return s
