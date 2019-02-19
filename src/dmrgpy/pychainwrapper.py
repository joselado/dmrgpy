import numpy as np
from .kpmdmrg import restrict_interval


# wrapper function for pychain

def get_full_hamiltonian(self):
  sc = get_pychain(self) # get pychain object
  def get_coupling(i,j):
    """Return the coupling between two sites"""
    for c in self.exchange:
      if i==c.i and j==c.j: return c.g
    return np.zeros((3,3))
  h = sc.add_tensor_interaction(get_coupling) # add interaction
  h = h + sc.add_exchange(self.fields) # add magnetic fields
  return h


def get_pychain(self):
  from .pychain import build
  sc = build.Spin_chain()
  # the pychain library assumes that s=1/2 is for spin one-half
  # whereas in DMRG s = 2 is for S=1/2
  sc.build((np.array(self.sites)-1.)/2.) 
  return sc





def get_dynamical_correlator(self,n=1000,mode="DMRG",i=0,j=0,
             window=[-1,10],name="XX",delta=2e-2,es=None):
  """
  Compute a dynamical correlator using the KPM-DMRG method
  """
  if delta is not None: # estimate the number of polynomials
    scale = 0.
    for s in self.sites:
        if s>1: scale += s**2 # spins
        else: scale += 4. # anything else
    n = int(scale/(4.*delta))
  self.to_folder() # go to temporal folder
  h = self.get_full_hamiltonian()
  sc = self.get_pychain()
  from .pychain import correlator as pychain_correlator
  if delta is None: delta = float(self.ns)/n*1.5
  if mode=="fullKPM":
    (xs,ys) = pychain_correlator.dynamical_correlator_kpm(sc,h,n=n,i=i,j=j,
                       namei=name[0],namej=name[1])
  elif mode=="ED":
    (xs,ys) = pychain_correlator.dynamical_correlator(sc,h,delta=delta,i=i,
                      j=j,namei=name[0],namej=name[1])
  else: raise
  self.to_origin() # go to origin folder
  if es is None:
    (xs,ys) = restrict_interval(xs,ys,window) # restrict the interval
  else:
    (xs,ys) = restrict_interval(xs,ys,[min(es),max(es)]) 
  from scipy.interpolate import interp1d
  fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
  fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
  if es is None:
      ne = int(100*(window[1] - window[0])/delta) # number of energies
      xs = np.linspace(window[0],window[1],ne)
  else: xs = np.array(es).copy() # copy input array
  ys = fr(xs) + 1j*fi(xs) # evaluate the interpolator
  np.savetxt("DYNAMICAL_CORRELATOR.OUT",np.matrix([xs.real,ys.real,ys.imag]).T)
  return (xs,ys)


