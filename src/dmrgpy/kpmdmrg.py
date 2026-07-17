from __future__ import print_function
import numpy as np
from . import multioperator
from . import operatornames
from .algebra import kpm





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
    Compute a dynamical correlator using the KPM-DMRG method, via the
    in-process pybind11 extension (mpscpp2/chain_session.h's
    Chain::kpm_dynamical_correlator).
    """
    if delta<0.0: raise
    if self.kpm_extrapolate: delta = delta*self.kpm_extrapolate_factor
    self.get_gs() # compute ground state (also sets self.e0)
    if type(name[0])!=multioperator.MultiOperator: raise
    mi = name[1] # first operator
    mj = name[0].get_dagger() # second operator
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    moments,emin,emax,scale,n = self._session.kpm_dynamical_correlator(
            mi.to_terms(),mj.to_terms(),
            self.kpmmaxm,self.kpm_scale,self.kpm_accelerate,
            self.kpm_n_scale,delta,self.kpmcutoff)
    mus = np.array(moments)
    if self.kpm_extrapolate:
        mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
                extrapolation_mode=self.kpm_extrapolate_mode)
    xs = 0.99*np.linspace(-1.0,1.0,int(n*10),endpoint=False) # energies
    ys = generate_profile(mus,xs,use_fortran=False,kernel=kernel) # generate the DOS
    xs /= scale # scale back the energies
    xs += (emin+emax)/2. -emin # shift the energies
    ys *= scale # renormalize the y values
    from scipy.interpolate import interp1d
    fr = interp1d(xs, ys.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(xs, ys.imag,fill_value=0.0,bounds_error=False)
    return (es,fr(es)+1j*fi(es)) # interpolate





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
    mus = general_kpm_moments_cpp_ext(self,X,wfa,wfb,num_p,self.kpm_accelerate)
    if self.kpm_extrapolate:
        mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
                extrapolation_mode=self.kpm_extrapolate_mode)
    return mus,shift,scale


def general_kpm_moments_cpp_ext(self,X,wfa,wfb,num_p,accelerate):
    """
    KPM moments of an arbitrary operator between two wavefunctions via the
    in-process pybind11 extension (mpscpp2/chain_session.h's
    Chain::general_kpm), mirroring general_kpm_moments()/
    kpm_moments_wfa_wfb()'s DMRG path exactly but with no file I/O. Shared
    by both callers -- the only difference between them is whether
    kpm_accelerate is self.kpm_accelerate or hardcoded false, which the
    caller now passes in directly.
    """
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    mus = self._session.general_kpm(X.to_terms(),wfa.cpp_handle,wfb.cpp_handle,
            self.maxm,accelerate,int(num_p),self.cutoff)
    return np.array(mus)


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



# this function is a bit redundant with general_kpm
def kpm_wfa_wfb(self,kernel="jackson",xs=None,**kwargs):
    """
    Compute a dynamical correlator of wfa and wfb
    """
    mus,shift,scale = kpm_moments_wfa_wfb(self,**kwargs)
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



def kpm_moments_wfa_wfb(self,X=None,wfa=None,wfb=None,
        scale=None,a=-0.8,b=0.8,
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
    if wfa is None or wfb is None:
        print("Wavefunctions must be provided")
        raise
    mus = general_kpm_moments_cpp_ext(self,X,wfa,wfb,num_p,False)
    if self.kpm_extrapolate:
        mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
                extrapolation_mode=self.kpm_extrapolate_mode)
    return mus,shift,scale
















