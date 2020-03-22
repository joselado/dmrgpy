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
  self.get_gs() # compute ground state
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
      mi = name[0] # first operator
      mj = name[1] # second operator
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
    xs = 0.99*np.linspace(-1.0,1.0,n*10,endpoint=False) # energies
    ys = generate_profile(mus,xs,use_fortran=False,kernel="lorentz") # generate the DOS
    xs /= scale # scale back the energies
    xs += (emin+emax)/2. -emin # shift the energies
    ys *= scale # renormalize the y values
    from scipy.interpolate import interp1d
    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
    return (es,fr(es)+1j*fi(es))
#    e0 = self.gs_energy() # ground state energy
    # now retain only an energy window
#  else: 
#    h = self.get_full_hamiltonian()
#    sc = self.get_pychain()
#    from .pychain import correlator as pychain_correlator
#    if delta is None: delta = float(self.ns)/n*1.5
#    if mode=="fullKPM":
#      (xs,ys) = pychain_correlator.dynamical_correlator_kpm(sc,h,n=n,i=i,j=j,
#                         namei=name[0],namej=name[1])
#    elif mode=="ED":
#      (xs,ys) = pychain_correlator.dynamical_correlator(sc,h,delta=delta,i=i,
#                        j=j,namei=name[0],namej=name[1])
#    else: raise
#  if es is None:
#    (xs,ys) = restrict_interval(xs,ys,window) # restrict the interval
#  else:
#    (xs,ys) = restrict_interval(xs,ys,[min(es),max(es)]) # restrict the interval
#  from scipy.interpolate import interp1d
#  fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
#  fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
#  if es is None: 
#      ne = int(100*(window[1] - window[0])/delta) # number of energies
#      xs = np.linspace(window[0],window[1],ne)
#  else: xs = np.array(es).copy() # copy input array
#  ys = fr(xs) + 1j*fi(xs) # evaluate the interpolator
##  np.savetxt("DYNAMICAL_CORRELATOR.OUT",np.matrix([xs.real,ys.real,ys.imag]).T)
#  return (xs,ys)



