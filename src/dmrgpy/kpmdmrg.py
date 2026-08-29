from __future__ import print_function
import numpy as np
from . import multioperator
from . import operatornames
from .algebra import kpm





from .algebra.kpm import generate_profile


def _sync_kpm_energy_truncation(self, enabled=None):
    """Push kpm_energy_truncate*'s current values onto self._session, if
    it understands them. `enabled` overrides self.kpm_energy_truncate for
    this push (used to force the feature *off* on code paths where it is
    not supported -- see general_kpm_moments_cpp_ext -- since the session
    otherwise keeps whatever a previous get_dynamical_correlator() call
    left enabled on it). Only pyitensor.chain.Chain does (see its
    set_kpm_energy_truncation()): the pyitensor backend wires energy
    truncation (Holzner et al., PRB 83, 195115 (2011), Sec. III-B) into
    its *existing* kpm_dynamical_correlator() via this setter + an
    internal branch. itensor_version=3 has its own, independent native
    port instead (mpscpp3/chain_session.h's kpm_dynamical_correlator_
    truncated(), a wholly separate method from kpm_dynamical_correlator
    -- see get_dynamical_correlator()'s own dispatch below, which calls
    it directly rather than through this setter). itensor_version=2 has
    neither."""
    if self.itensor_version != "python":
        return
    if enabled is None:
        enabled = self.kpm_energy_truncate
    self._session.set_kpm_energy_truncation(
        enabled, self.kpm_truncate_dK,
        self.kpm_truncate_nsweeps, self.kpm_truncate_threshold)


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
             hodc_order=6,hodc_eta=None,
             **kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method, via the
    in-process pybind11 extension (mpscpp2/chain_session.h's
    Chain::kpm_dynamical_correlator).

    kernel selects the reconstruction applied to the Chebyshev moments:
    "jackson" (default), "lorentz", "plain", or "hodc" for the high-order
    delta-Chebyshev kernel of arXiv:2512.03149 (see algebra/kpm.py).
    hodc_order/hodc_eta are the order m and regularization width of the
    latter and are ignored by every other kernel; hodc_eta is given in
    the same energy units as delta/es and defaults to delta itself.

    The moment computation is identical for every kernel -- see
    dynamical_correlator_moments()/dynamical_correlator_from_moments(),
    which this function just composes, and which let two kernels be
    compared from a single (expensive) moment run.
    """
    mus,emin,emax,scale,n,delta = dynamical_correlator_moments(self,
            name=name,delta=delta,**kwargs)
    return dynamical_correlator_from_moments(mus,emin,emax,scale,n,es,
            kernel=kernel,delta=delta,
            hodc_order=hodc_order,hodc_eta=hodc_eta)


def dynamical_correlator_moments(self,name=None,delta=1e-1,**kwargs):
    """Compute the raw (undamped) KPM-DMRG Chebyshev moments
    <name[0]|T_k(H_scaled)|name[1]> of a dynamical correlator, plus
    everything needed to turn them into a spectrum.

    Returns (mus,emin,emax,scale,n,delta), where scale/emin/emax define
    the rescaling of the spectrum onto [-1,1], n is the number of
    polynomials the backend chose (before any kpm_extrapolate resampling
    of mus) and delta is the effective resolution actually requested.
    Split out of get_dynamical_correlator() so a caller can reconstruct
    the same moments with several kernels without repeating the DMRG
    work, which is all of the cost."""
    if delta<0.0: raise
    if self.kpm_extrapolate: delta = delta*self.kpm_extrapolate_factor
    self.get_gs() # compute ground state (also sets self.e0)
    if type(name[0])!=multioperator.MultiOperator: raise
    mi = name[1] # first operator
    mj = name[0].get_dagger() # second operator
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    if self.kpm_energy_truncate and self.itensor_version==3:
        # v3's own independent method (see _sync_kpm_energy_truncation's
        # docstring) -- not the setter+branch pattern used for "python".
        moments,emin,emax,scale,n = self._session.kpm_dynamical_correlator_truncated(
                mi.to_terms(),mj.to_terms(),
                self.kpmmaxm,self.kpm_scale,self.kpm_accelerate,
                self.kpm_n_scale,delta,self.kpmcutoff,
                self.kpm_truncate_dK,self.kpm_truncate_nsweeps,self.kpm_truncate_threshold)
    else:
        _sync_kpm_energy_truncation(self)
        moments,emin,emax,scale,n = self._session.kpm_dynamical_correlator(
                mi.to_terms(),mj.to_terms(),
                self.kpmmaxm,self.kpm_scale,self.kpm_accelerate,
                self.kpm_n_scale,delta,self.kpmcutoff)
    mus = np.array(moments)
    if self.kpm_extrapolate:
        mus = kpm.extrapolate_moments(mus,fac=self.kpm_extrapolate_factor,
                extrapolation_mode=self.kpm_extrapolate_mode)
    return mus,emin,emax,scale,n,delta


def dynamical_correlator_from_moments(mus,emin,emax,scale,n,es,
        kernel="jackson",delta=None,hodc_order=6,hodc_eta=None):
    """Reconstruct a dynamical correlator on the energy grid es from the
    output of dynamical_correlator_moments(). See get_dynamical_
    correlator() for the meaning of kernel/hodc_order/hodc_eta."""
    xs = 0.99*np.linspace(-1.0,1.0,int(n*10),endpoint=False) # energies
    if kernel=="hodc":
        # hodc_eta is quoted in physical energy units, like delta and es;
        # generate_profile wants it in the rescaled units of xs, i.e.
        # before the "xs/scale" below. Left unset it defaults to the
        # requested resolution delta, which is also what keeps n*eta in
        # the flat minimum of HODC's error-vs-eta curve -- see
        # algebra/kpm.py's HODC_DEFAULT_P_ETA.
        if hodc_eta is not None: hodc_eta = hodc_eta*scale
        elif delta is not None: hodc_eta = delta*scale
    ys = generate_profile(mus,xs,use_fortran=False,kernel=kernel,
            hodc_order=hodc_order,hodc_eta=hodc_eta) # generate the DOS
    xs = xs/scale # scale back the energies
    xs = xs + (emin+emax)/2. -emin # shift the energies
    ys = ys*scale # renormalize the y values
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
    # Energy truncation is deliberately forced off here: it is defined
    # only for the (A,B) dynamical correlator, whose Chebyshev recursion
    # is driven by the ground-state-anchored rescaled *Hamiltonian*. This
    # path expands an arbitrary operator X instead (already rescaled into
    # [a,b] by scale_operator(), so it needs no such safeguard), and
    # filtering against X's own local spectrum is a different operation
    # that no other backend performs -- itensor_version=3's general_kpm
    # has no truncation at all, so leaving it on here would silently make
    # the two backends disagree. See the user guide's "Available only for
    # the (A,B)-operator dynamical correlator" note.
    _sync_kpm_energy_truncation(self, enabled=False)
    mus = self._session.general_kpm(X.to_terms(),wfa.cpp_handle,wfb.cpp_handle,
            self.maxm,accelerate,int(num_p),self.cutoff)
    return np.array(mus)


def general_kpm(self,kernel="jackson",xs=None,hodc_order=6,hodc_eta=None,**kwargs):
    """
    Compute a dynamical correlator of Bdelta(X)A using the KPM-DMRG method
    """
    mus,shift,scale = general_kpm_moments(self,**kwargs)
    # scale of the distribution
    kpmscales = scale
    num_p = len(mus)
    xs2 = 0.99*np.linspace(-1.0,1.0,int(num_p*10),endpoint=False) # energies
    # hodc_eta is quoted in the units of the returned (physical) x axis,
    # as it is in get_dynamical_correlator; generate_profile wants it in
    # the rescaled units of xs2, before the "xs2 *= scale" below.
    if kernel=="hodc" and hodc_eta is not None: hodc_eta = hodc_eta/scale
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel=kernel,
            hodc_order=hodc_order,hodc_eta=hodc_eta) # generate the DOS
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
def kpm_wfa_wfb(self,kernel="jackson",xs=None,hodc_order=6,hodc_eta=None,**kwargs):
    """
    Compute a dynamical correlator of wfa and wfb
    """
    mus,shift,scale = kpm_moments_wfa_wfb(self,**kwargs)
    # scale of the distribution
    kpmscales = scale
    num_p = len(mus)
    xs2 = 0.99*np.linspace(-1.0,1.0,int(num_p*10),endpoint=False) # energies
    # hodc_eta is quoted in the units of the returned (physical) x axis,
    # as it is in get_dynamical_correlator; generate_profile wants it in
    # the rescaled units of xs2, before the "xs2 *= scale" below.
    if kernel=="hodc" and hodc_eta is not None: hodc_eta = hodc_eta/scale
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel=kernel,
            hodc_order=hodc_order,hodc_eta=hodc_eta) # generate the DOS
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
















