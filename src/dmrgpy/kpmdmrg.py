from __future__ import print_function
import numpy as np
from . import multioperator
from . import operatornames
from .algebra import kpm

def get_moments_dmrg(self,n=1000):
  """Get the moments with DMRG"""
  self.setup_task("dos",task={"nkpm":str(n)})
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  return self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose()[0])









def get_moments_dynamical_correlator_dmrg(self,name=None,delta=1e-1):
  """Get the moments with DMRG"""
  # do some sanity checks
  if delta<0.0: raise 
  if self.kpm_extrapolate: delta = delta*self.kpm_extrapolate_factor
  if self.itensor_version!="julia": self.get_gs() # compute ground state
  # define the dictionary
  task = {      "dynamical_correlator": "true",
                "kpmmaxm":str(self.kpmmaxm),
                "kpm_scale":str(self.kpm_scale),
                "kpm_accelerate":self.kpm_accelerate,
                "kpm_n_scale":str(self.kpm_n_scale),
                "kpm_delta":str(delta),
                "kpm_cutoff":str(self.kpmcutoff),
                }
  # go on, check the kind of input used to define the correlator
  if type(name[0])==multioperator.MultiOperator: 
      task["kpm_multioperator_i"] = "true"
      task["kpm_multioperator_j"] = "true"
      mi = name[1] # first operator
      mj = name[0] # second operator
      mj = mj.get_dagger()
      self.execute(lambda: mi.write(name="kpm_multioperator_i.in")) # write
      self.execute(lambda: mj.write(name="kpm_multioperator_j.in")) # write
  else: raise
  self.task = task # assign tasks
  self.write_task() 
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  m = self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose())
#  return m[1]
  mus = m[0]+1j*m[1]
  from .algebra import kpm
  if self.kpm_extrapolate: 
      return kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
              extrapolation_mode=self.kpm_extrapolate_mode)
  else: return mus






from . import pychain
from .algebra.kpm import generate_profile


def restrict_interval(x,y,window):
  """Restrict the result to a certain energy window"""
  if window is None: return (x,y)
  i = np.argwhere(x<window[0]) # last one
  j = np.argwhere(x>window[1]) # last one
  if len(i)==0: i = 0
  else: i = i[0][-1]
  if len(j)==0: j = len(x)
  else: j = j[0][0]
  return x[i:j].real,y[i:j]







def get_dynamical_correlator(self,n=1000,
             name=None,delta=1e-1,kernel="jackson",
             es=np.linspace(-1.,10,500),deconvolve=None,
             **kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
# get the moments
    mus = get_moments_dynamical_correlator_dmrg(self,delta=delta,
            name=name,**kwargs) 
    # scale of the dos
    kpmscales = self.execute(lambda: np.genfromtxt("KPM_SCALE.OUT"))
    emin = kpmscales[0] # minimum energy
    emax = kpmscales[1] # maximum energy
    scale = kpmscales[2] # scaling of the energies
    # ground state energy
    e0 = self.execute(lambda: np.genfromtxt("GS_ENERGY.OUT"))
    self.e0 = e0 # add this quantity
    n = self.execute(lambda: np.genfromtxt("KPM_NUM_POLYNOMIALS.OUT"))
    xs = 0.99*np.linspace(-1.0,1.0,int(n*10),endpoint=False) # energies
#    if self.kpm_extrapolate: kernel = None # no kernel
    ys = generate_profile(mus,xs,use_fortran=False,kernel=kernel) # generate the DOS
    xs /= scale # scale back the energies
    xs += (emin+emax)/2. -emin # shift the energies
    ys *= scale # renormalize the y values
    from scipy.interpolate import interp1d
    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
    (es,z) = (es,fr(es)+1j*fi(es)) # interpolate
    return (es,z)
#    from .algebra import kpm
#    (es,z) = kpm.deconvolution(es,z,mode=deconvolve,delta=delta)
#    return es,z





def general_kpm_moments(self,X=None,A=None,B=None,
        scale=None,wf=None,a=-0.8,b=0.8,
        delta=1e-1,kernel="jackson",xs=None,**kwargs):
    """
    Compute a dynamical correlator of Bdelta(X)A using the KPM-DMRG method
    """
    if X is None: raise
    # extrapolate
    if self.kpm_extrapolate: delta = delta*self.kpm_extrapolate_factor
    if scale is not None: 
        X = X/scale # renormalize the operator for KPM
        shift = 0.0 # no additional shift
    else: # no scale provided
        X,scale,shift = scale_operator(self,X,a=a,b=b)
    num_p = int(3*scale/delta) # number of polynomials
    if wf is None: wf = self.get_gs() # no wavefunction provided
    # compute the wavefunctions
    if A is not None: wfa = self.applyoperator(A,wf)
    else: wfa = wf
    if B is not None: wfb = self.applyoperator(B,wf)
    else: wfb = wf
    # write the wavefunctions
    self.execute(lambda: wfa.write(name="wfa.mps"))
    self.execute(lambda: wfb.write(name="wfb.mps"))
    # write the task
    task = {    "general_kpm": "true",
                "kpmmaxm":str(self.maxm),
                "kpm_accelerate":self.kpm_accelerate,
                "kpm_num_polynomials":str(num_p),
                "kpm_cutoff":str(self.cutoff),
                }
    self.execute(lambda: X.write(name="kpm_operator.in"))
    self.task = task # assign tasks
    self.run() # perform the calculation
    m = self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose())
    mus = m[0]+1j*m[1]
    # perform extrapolation if 
    if self.kpm_extrapolate: 
      mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
              extrapolation_mode=self.kpm_extrapolate_mode)
    return mus,shift,scale

def general_kpm(self,kernel="jackson",xs=None,**kwargs):
    """
    Compute a dynamical correlator of Bdelta(X)A using the KPM-DMRG method
    """
    mus,shift,scale = general_kpm_moments(self,**kwargs)
    # scale of the distribution
    kpmscales = scale
    num_p = len(mus)
    xs2 = 0.99*np.linspace(-1.0,1.0,int(num_p*10),endpoint=False) # energies
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel=kernel) # generate the DOS
    xs2 += shift # add the shift
    xs2 *= scale # scale
    ys2 /= scale # scale
    if xs is None: return xs2,ys2
    else:
      from scipy.interpolate import interp1d
      fr = interp1d(xs2, ys2.real,fill_value=0.0,bounds_error=False)
      fi = interp1d(xs2, ys2.imag,fill_value=0.0,bounds_error=False)
      return xs,fr(xs)+1j*fi(xs)


def scale_operator(self,X,a=-0.9,b=0.9):
    """Scale an operator so its spectra falls in the interval [a,b]"""
    A = self.lowest_eigenvalue(X) # compute the lowest eigenvalue
    B = -self.lowest_eigenvalue(-X) # compute the highest eigenvalue
    ab = b - a  # width of the output
    AB = B - A # width of the input
    X = X - A # shift to zero
    X = X/AB*ab # redefine width
    X = X + a # shift to the right point
    X = X.simplify() # simplify
    scale = AB/ab # scaling of the operator
    shift = (a*B - b*A)/(A-B) # required shift
    return X,scale,shift # return the new operator



