from __future__ import print_function
import numpy as np
from . import multioperator
from . import operatornames

def get_moments_dmrg(self,n=1000):
  """Get the moments with DMRG"""
  self.setup_task("dos",task={"nkpm":str(n)})
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  return self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose()[0])









def get_moments_dynamical_correlator_dmrg(self,name=None,delta=1e-1):
  """Get the moments with DMRG"""
#  if type(name)==str: # string input
#        namei,namej = operatornames.recognize(name)
#        namei = self.get_operator(namei,i)
#        namej = self.get_operator(namej,j)
#        return get_moments_dynamical_correlator_dmrg(self,
#                name=(namei,namej),delta=delta)
  # do some sanity checks
  if delta<0.0: raise 
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
  return m[0]+1j*m[1]





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
             name=None,
             es=np.linspace(-1.,10,500),
             **kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
# get the moments
    mus = get_moments_dynamical_correlator_dmrg(self,
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
    ys = generate_profile(mus,xs,use_fortran=False,kernel="lorentz") # generate the DOS
    xs /= scale # scale back the energies
    xs += (emin+emax)/2. -emin # shift the energies
    ys *= scale # renormalize the y values
    from scipy.interpolate import interp1d
    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
    return (es,fr(es)+1j*fi(es))





def general_kpm(self,X=None,A=None,B=None,scale=None,wf=None,
        delta=1e-1,xs=None,**kwargs):
    """
    Compute a dynamical correlator of Bdelta(X)A using the KPM-DMRG method
    """
    if scale is None: scale = self.bandwidth(X)*1.1
    if X is None: raise
    X = X/scale # renormalize the operator for KPM
    num_p = int(3*scale/delta) # number of polynomial
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
    # scale of the dos
    kpmscales = scale
    xs2 = 0.99*np.linspace(-1.0,1.0,int(num_p*10),endpoint=False) # energies
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel="lorentz") # generate the DOS
    xs2 *= scale
    ys2 /= scale
    return xs2,ys2


