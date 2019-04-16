import numpy as np
from .algebra.kpm import generate_profile

# library to compute the DOS

def get_moments_dos_dmrg(self,delta=1e-1):

  """Get the moments with DMRG"""
  task= {       "kpmmaxm":str(self.kpmmaxm),
                "kpm_scale":str(self.kpm_scale),
                "kpm_n_scale":str(self.kpm_n_scale),
                "kpm_delta":str(delta),
                "kpm_cutoff":str(self.kpmcutoff)}
  self.setup_task("dos",task=task)
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  m = self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose())
  return m[0]+1j*m[1] # return the moments




def get_dos(self,delta=2e-2,**kwargs):
  """
  Compute the many body DOS using the KPM-DMRG method
  """
  self.to_folder() # go to temporal folder
  mus = get_moments_dos_dmrg(self,delta=delta,**kwargs)
  # scale of the dos
  kpmscales = self.execute(lambda: np.genfromtxt("KPM_SCALE.OUT"))
  emin = kpmscales[0] # minimum energy
  emax = kpmscales[1] # maximum energy
  scale = kpmscales[2] # scaling of the energies
  n = self.execute(lambda: np.genfromtxt("KPM_NUM_POLYNOMIALS.OUT"))
  xs = 0.99*np.linspace(-1.0,1.0,n*10,endpoint=True) # energies
  # generate the DOS
  ys = generate_profile(mus,xs,use_fortran=False,kernel="jackson").real
  xs /= scale # scale back the energies
  xs += emin # shift energy
  xs += (emin+emax)/2.-emin # shift the energies
  ys *= scale # renormalize the y positions
  self.to_origin() # go to origin folder
  return (xs,ys) # return the DOS











