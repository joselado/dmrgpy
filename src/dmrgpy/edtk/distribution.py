import scipy.sparse.linalg as slg
from ..algebra import kpm
from ..algebra import algebra
import numpy as np
from scipy.interpolate import interp1d
from .edchain import State

def get_distribution(self,X=None,wf=None,method="KPM",**kwargs):
    """Get a certain distribution"""
    if wf is None: wf = self.get_gs_array() # ground state
    elif type(wf)==State: wf = wf.v
    X = self.get_operator(X)
    if method=="KPM":
      return distribution_kpm(wf,X=X,**kwargs)
    if method=="INV":
      return distribution_inversion(wf,X=X,**kwargs)


def distribution_kpm(wf0,X=None,scale=10.0,
        delta=1e-1,xs=None):
    """Compute <0| \delta (m-M) |0> using the KPM"""
    if X is None: raise
    M = X # assign
    vi = wf0 # first wavefunction
    vj = wf0 # second wavefunction
    if scale is None:
        emax = slg.eigsh(M,k=1,ncv=20,which="LA")[0] # upper energy
        emin = slg.eigsh(-M,k=1,ncv=20,which="LA")[0] # upper energy
        scale = np.max(np.abs([emin,emax]))*2. # compute the scale
    n = int(scale/delta) # number of polynomials
    (xs2,ys2) = kpm.dm_vivj_energy(M,vi,vj,scale=scale,
                                npol=n*4,ne=n*10)
    ys2 /= np.pi # normalization from the KPM
    if xs is None: return xs2,ys2
    ys = interp1d(xs2,ys2.real,fill_value=0.,bounds_error=False)(xs) 
    ys = ys+ 1j*interp1d(xs2,ys2.real,fill_value=0.,bounds_error=False)(xs) 
    return xs,ys # return correlator


def distribution_inversion(wf0,X=None,scale=10.0,
        delta=2e-1,mode="full"):
  """Calculate a correlation function SiSj in a frequency window"""
  es=np.linspace(-scale,scale,5*int(scale/delta))
  ## default method
  iden = np.identity(X.shape[0],dtype=np.complex_) # identity
  out = []
  for e in es: # loop over energies
      if mode=="full": # using exact inversion
        g1 = algebra.inv(iden*(e+1j*delta)-X)
        g2 = algebra.inv(iden*(e-1j*delta)-X)
        g = 1j*(g1-g2)/2./np.pi
        op = g # operator
        o = algebra.braket_wAw(wf0,op) # correlator
#      elif mode=="cv": # correction vector algorithm
#          o1 = solve_cv(h0,wf0,A,B,e+e0,delta=delta) # conjugate gradient
#          o2 = solve_cv(h0,wf0,A,B,e+e0,delta=-delta) # conjugate gradient
#          o = 1j*(o1 - o2)/2. # substract
      else: raise # not recognised
      out.append(o)
  return es,np.array(out) # return result

