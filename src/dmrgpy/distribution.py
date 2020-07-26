from .algebra.kpm import generate_profile
import numpy as np



def get_moments_distribution(self,X=None,num_p=None):
  """Get the mments with DMRG"""
  # define the dictionary
  self.gs_energy()
  task = {      "distribution": "true",
                "kpmmaxm":str(self.kpmmaxm),
                "kpm_scale":str(self.kpm_scale),
                "kpm_accelerate":self.kpm_accelerate,
                "kpm_n_scale":str(self.kpm_n_scale),
                "kpm_num_polynomials":str(num_p),
                "kpm_cutoff":str(self.kpmcutoff),
                }
  self.execute(lambda: X.write(name="kpm_distribution_multioperator.in")) 
  self.task = task # assign tasks
  self.write_task()
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  m = self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose())
  return m[0]+1j*m[1]





def get_distribution(self,**kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
    from .kpmdmrg import general_kpm
    return general_kpm(self,**kwargs)

