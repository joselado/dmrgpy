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





def get_distribution(self,X=None,scale=None,delta=1e-1,xs=None,**kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
    if scale is None: raise
    X = X/scale # renormalize the operator for KPM
    num_p = int(3*scale/delta) # number of polynomial
    mus = get_moments_distribution(self,X=X,num_p=num_p)
    # scale of the dos
    kpmscales = scale
    n = self.execute(lambda: np.genfromtxt("KPM_NUM_POLYNOMIALS.OUT"))
    xs2 = 0.99*np.linspace(-1.0,1.0,int(n*10),endpoint=False) # energies
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel="lorentz") # generate the DOS
    xs2 *= scale
    ys2 /= scale
    return xs2,ys2
#    xs /= scale # scale back the energies
#    xs += (emin+emax)/2. -emin # shift the energies
#    ys *= scale # renormalize the y values
#    from scipy.interpolate import interp1d
#    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
#    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
#    return (es,fr(es)+1j*fi(es))

