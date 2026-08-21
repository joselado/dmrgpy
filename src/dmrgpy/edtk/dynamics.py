from ..algebra import algebra
from .. import multioperator
import scipy.sparse.linalg as slg
from ..algebra import kpm
import numpy as np
from .edchain import EDOperator


#from numba import jit

is_hermitian = algebra.is_hermitian


def get_dynamical_correlator(self,name=None,submode="KPM",
        wf0=None,T=0.,**kwargs):
    """
    Compute the dynamical correlator
    """
    if name is None:
        raise ValueError("get_dynamical_correlator: name (the two-operator "
                          "tuple) must be given")

    A = EDOperator(name[0],self).SO # create first operator
    B = EDOperator(name[1],self).SO # create second operator
    h = self.get_hamiltonian() # Hamiltonian as matrix

    if T>0.: # finite-temperature Lehmann sum -- the exact reference used
        # to validate the finite-T MPS/METTS dynamical correlator (see
        # vevtk/mettsdynamicalcorrelator.py), the dynamical analog of
        # thermal_vev_ex's static Boltzmann sum in vevtk/thermalvev.py.
        if submode!="ED":
            raise NotImplementedError(
                "get_dynamical_correlator: T>0 is only implemented for "
                "submode='ED' (the exact Lehmann sum) so far, got "
                "submode=%r" % (submode,))
        if not is_hermitian(h):
            raise NotImplementedError(
                "get_dynamical_correlator: T>0 is not implemented for "
                "non-Hermitian Hamiltonians")
        emu,vs = self.get_diagonalized_hamiltonian()
        return dynamical_correlator_finite_T(h,A,B,T,emu=emu,vs=vs,**kwargs)

    if wf0 is None:  wf0 = self.get_gs_array() # compute ground state
    else:
        wf0 = wf0.v.copy()
    if not is_hermitian(h): # non Hermitian Hamiltonians
        if submode=="KPM": # non-Hermitian KPM
            from ..nonhermitian.kpm import dynamical_correlator_nhkpm_ed
            return dynamical_correlator_nhkpm_ed(self,name=name,**kwargs)
        from ..nonhermitian.dynamics import dynamical_correlator_non_hermitian
        return dynamical_correlator_non_hermitian(self,name=name,**kwargs)
#    print(wf0)
#    print(np.round(wf0,2))
    # for Hermitian Hamiltonians, continue
    if submode=="KPM":
        return dynamical_correlator_kpm(h,self.e0,wf0,A,B,chain=self,**kwargs)
    elif submode=="ED":
        emu,vs = self.get_diagonalized_hamiltonian()
        return dynamical_correlator_ED(h,A,B,emu=emu,vs=vs,**kwargs)
    elif submode=="EX":
        from .. import dcex
        return dcex.dynamical_correlator(self,name=name,**kwargs)
    elif submode=="INV":
      return dynamical_correlator_inv(h,wf0,self.e0,A,B,mode="full",**kwargs)
    elif submode=="CVM":
      return dynamical_correlator_inv(h,wf0,self.e0,A,B,mode="cv",**kwargs)
    elif submode=="ROOTN":
      return dynamical_correlator_rootn(h,self.e0,wf0,A,B,**kwargs)
    elif submode=="TD":
      from .. import timedependent
      return timedependent.dynamical_correlator(self,mode="ED",
              name=name,**kwargs)
    else:
        # A bare `raise` here (no active exception to re-raise) used to
        # crash with a confusing "RuntimeError: No active exception to
        # reraise" instead of naming the actual problem -- confirmed
        # directly: reachable via mode.py's own silent DMRG->ED fallback
        # whenever the requested submode (e.g. "TDZ", a DMRG/TDVP-only
        # complex-time-evolution method with no ED counterpart at all)
        # isn't one of the ones implemented above.
        raise NotImplementedError(
            "get_dynamical_correlator: submode=%r has no ED implementation"
            % (submode,))


def get_kpm_emax(chain,m,e0):
    """Return the KPM energy-scale estimate (top of the spectrum of the
    shifted Hamiltonian m=h-e0*Id), cached on the chain. This depends
    only on the Hamiltonian and its ground-state energy, not on the
    operators being correlated, so callers that sweep the same chain
    over several sites/operators (e.g. dynamicstk/spincorrelators.py)
    would otherwise repeat this ARPACK/dense eigenvalue solve on every
    single call."""
    if chain is not None and chain._kpm_emax_cache is not None:
        (cached_e0,cached_emax) = chain._kpm_emax_cache
        if cached_e0==e0: return cached_emax
    emax = -algebra.lowest_eigenvalues(-m,n=3)[0] # lowest energy
    if chain is not None: chain._kpm_emax_cache = (e0,emax)
    return emax


def dynamical_correlator_kpm(h,e0,wf0,A,B,chain=None,
        delta=1e-1,es=np.linspace(-1.,10,400)):
    A = np.conjugate(A.T)
    vi = B@wf0 # first wavefunction
    vj = A@wf0 # second wavefunction
    from scipy.sparse import identity
    m = -identity(h.shape[0])*e0+h # matrix to use
    emax = get_kpm_emax(chain,m,e0) # lowest energy (cached per chain)
    scale = np.max([np.abs(e0),np.abs(emax)])*3.0
    n = int(2*scale/delta) # number of polynomials
    (xs,ys) = kpm.dm_vivj_energy(m,vi,vj,scale=scale,
                                npol=n*4,ne=n*10,x=es)
    return xs,np.conjugate(ys)*scale/np.pi # return correlator




def dynamical_correlator_rootn(h,e0,wf0,A,B,delta=1e-1,
        es=np.linspace(-1.,10,400),N=8,nkry=20):
    """Compute the dynamical correlator with the root-N Krylov-space
    correction-vector method (Nocera & Alvarez, arXiv:2204.03165), against
    the exact ED Hamiltonian -- see algebra/rootn.py for the algorithm.
    N is the number of root-N steps (N=1 recovers the "conventional"
    Krylov-space correction-vector method of Nocera, PRE 2016); nkry is
    the Lanczos subspace dimension used at every step (the ED analogue of
    the bond dimension m in the paper's MPS/DMRG setting).

    Follows the same <GS|A(omega-H+E0+i*eta)^{-1}B|GS> convention (A used
    as-is, no dagger) as dynamical_correlator_ED/dynamical_correlator_inv
    above, so results are directly comparable to submode="ED"/"CVM"/"INV"
    -- unlike dynamical_correlator_kpm, which conjugate-transposes A for
    its own (different) vi/vj construction."""
    from ..algebra.rootn import rootn_correction_vector
    out = [rootn_correction_vector(h,wf0,e0,A,B,e,delta,N=N,nkry=nkry)
            for e in es]
    return es,np.array(out)


def check_dex_sensitivity(ex,dex,factor=3.):
    """Warn when the equal-weight `dex` cutoff of dynamical_correlator_ED
    is placed inside a cluster of levels rather than in a clean gap.

    The cutoff is a step function of the excitation energy, so the answer
    only depends on its exact value when some eigenvalue sits close to it.
    "Close" here is within `factor` on either side, i.e. in
    [dex/factor, dex*factor]: a level just below is included with full
    weight and one just above is dropped entirely, and an infinitesimal
    change of the swept parameter can move one across. `ex` is the array
    of excitation energies (ascending, already relative to the ground
    state)."""
    near = ex[(ex>=dex/factor) & (ex<=dex*factor)] # levels straddling it
    if len(near)==0: return # cutoff sits in a gap, result is robust
    import warnings
    warnings.warn(
        "dynamical_correlator_ED: %d eigenvalue(s) lie within a factor "
        "%g of the dex=%g cutoff (lowest such level at %g), so the equally-weighted "
        "initial manifold (nex=%d) depends on exactly where dex was "
        "placed and can jump discontinuously as a swept parameter moves a "
        "level across it. Choose dex relative to the splittings being "
        "resolved, or pass a small T to get_dynamical_correlator for the "
        "Boltzmann-weighted finite-temperature sum instead."
        %(len(near),factor,dex,near[0],len(ex[ex<dex])),
        RuntimeWarning,stacklevel=3)


def dynamical_correlator_ED(h,a0,b0,delta=2e-2,
        emu=None,vs = None,
        dex = 1e-5, # this is a tolerancy to consider something a GS
        es=np.linspace(-1.0,10.0,600)):
    """Compute a zero-temperature dynamical correlator.

    `dex` is a *hard, equal-weight* cutoff: every eigenstate with an
    excitation energy below it is taken as part of the initial manifold
    and given the same weight 1/nex (see dynamical_sum below). That is
    only physical when the manifold is genuinely degenerate. When it is
    split by an amount comparable to `dex` -- e.g. a Zeeman-field sweep
    whose splitting grows through the cutoff -- the result jumps
    discontinuously as states cross the threshold: a state at 0.99*dex
    counts fully, one at 1.01*dex not at all, and the sweep partially
    averages over the very splitting it is meant to resolve. `dex` has to
    be chosen relative to the splittings being resolved, and no single
    value works for both ends of such a sweep.

    For that regime use the finite-temperature route instead: pass a
    small `T` to get_dynamical_correlator, which dispatches to
    dynamical_correlator_finite_T and populates the manifold with exact
    Boltzmann weights, reproducing this function's degenerate-manifold
    average smoothly as T->0 rather than through a step. The equal-weight
    cutoff is kept as the T=0 default so existing callers' results are
    unchanged; check_dex_sensitivity below warns when the answer actually
    depends on where the cutoff was placed."""
    if emu is None or vs is None: # if not provided
        emu,vs = algebra.eigh(h) # compute them
    ex = emu-np.min(emu) # excitations
    check_dex_sensitivity(ex,dex) # warn if the cutoff is doing real work

    # crop to the needed states
    emax = np.max(es) # maximum energy
    vs = vs[:,emu<emax] # restrict
    emu = emu[emu<emax] # restrict
    # finnish cropping

    # compute the needed matrix elements
    U = np.array(vs) # matrix
    Uh = np.conjugate(np.transpose(U)) # dagger
    A = Uh@a0@U # get the matrix elements
    B = Uh@b0@U # get the matrix elements
    nex = len(ex[ex<dex]) # number of excited states
    out = dynamical_sum(emu,es,delta,A,B,nex=nex) # perform the summation
    return (es,-out.imag/(2*np.pi)) # return correlator

from numba import jit

@jit(nopython=True)
def dynamical_sum(es,ws,delta,A,B,nex=1):
    """Return the sum giving the dynamical correlator, i.e.
    sum_iw(ws+i*delta) - sum_iw(ws-i*delta). Both terms share the same
    matrix element A[i,j]*B[j,i], so they are accumulated together in a
    single pass over (i,j,iw) instead of two separate passes, each
    recomputing that matrix element from scratch."""
    es = es-np.min(es) # remove minimum energy (ground state)
    n = len(es) # number of eigenenergies
    out = np.zeros(len(ws),dtype=np.complex128) # initialize
    for i in range(nex): # loop over excited states
      for iw in range(len(ws)): # loop over frequencies
#        i = 0 # initial wavefunction (Ground state)
        # this will not work properly if there are degeneracies
        w = ws[iw]
        acc = 0.0+0.0j
        for j in range(n): # loop over eigenenergies
            tmp = A[i,j]*B[j,i] # relevant matrix element
            acc = acc + tmp*(1./(w+1j*delta+es[i]-es[j])
                            - 1./(w-1j*delta+es[i]-es[j]))
        out[iw] = out[iw] + acc
    return out/nex # return dynamical correlator


def dynamical_correlator_finite_T(h,a0,b0,T,delta=2e-2,
        emu=None,vs=None,
        es=np.linspace(-1.0,10.0,600)):
    """Finite-temperature dynamical correlator via the exact Lehmann sum
    (arXiv:2405.18484, Eq. before Sec. II.D):

        C_AB(omega) = (1/Z) sum_{n,m} exp(-beta*E_n) <n|A|m><m|B|n>
                      * [same +-i*delta kernel dynamical_sum uses at
                         omega = E_m - E_n]

    the finite-T generalization of dynamical_correlator_ED's T=0 sum
    (dynamical_sum's `nex`-averaged near-degenerate-ground-state sum)
    to a full thermal average over every eigenstate weighted by its
    exact Boltzmann factor. Unlike thermal_vev_ex (vevtk/thermalvev.py),
    which sums over a possibly-truncated set of excited states obtained
    from sparse/partial diagonalization and must therefore check the
    left-out Boltzmann weight is negligible, emu/vs here always come from
    a full dense diagonalization (get_diagonalized_hamiltonian ->
    algebra.eigh), so the Boltzmann weight is exact and no truncation
    safety check is needed -- the only limit is algebra.maxsize's
    existing hard cap on how large a Hilbert space can be fully
    diagonalized at all."""
    if emu is None or vs is None: # if not provided
        emu,vs = algebra.eigh(h) # compute them
    ex = emu-np.min(emu) # excitations, ascending

    beta = 1./T
    weights = np.exp(-beta*ex)
    weights = weights/weights.sum() # exact Boltzmann weight, full spectrum

    # Crop only the *final*-state (j) index space to the requested
    # frequency window, same approximation dynamical_correlator_ED already
    # makes -- the *initial*-state (i) sum below always ranges over the
    # FULL spectrum with its exact weight. Cropping i's index space too
    # (as an earlier version of this function did) would silently drop
    # any thermally-populated state above emax from the numerator while
    # `weights` (computed from the full spectrum, above) still counted it
    # in the normalization -- a real, silent under-count at any T large
    # enough for excited states above emax to still carry non-negligible
    # weight, confirmed by code review.
    emax = np.max(es) # maximum energy
    keep_j = emu<emax
    ex_j = ex[keep_j] # restrict

    Ufull = np.array(vs) # (dim, n) -- every eigenstate, for the i axis
    Uj = np.array(vs[:,keep_j]) # (dim, nj) -- cropped, for the j axis
    A = np.conjugate(np.transpose(Ufull))@a0@Uj # (n,nj): <i|A|j>
    B = np.conjugate(np.transpose(Uj))@b0@Ufull # (nj,n): <j|B|i>
    out = dynamical_sum_thermal(ex,ex_j,es,delta,A,B,weights) # perform the summation
    return (es,-out.imag/(2*np.pi)) # return correlator


@jit(nopython=True)
def dynamical_sum_thermal(ex_i,ex_j,ws,delta,A,B,weights):
    """Thermal generalization of dynamical_sum() above: sums over *every*
    initial state i (not just the first `nex` near-degenerate ones),
    weighting each by its own Boltzmann factor `weights[i]` (already
    normalized to sum to 1 over the full spectrum) instead of a uniform
    1/nex. A is (n_i,n_j) = <i|A|j>, B is (n_j,n_i) = <j|B|i>; ex_i/ex_j
    are already expressed relative to the same shared zero (the global
    ground-state energy, both derived from the same `ex` in the caller)
    -- deliberately *not* independently re-zeroed against their own
    separate minima here, since ex_i and ex_j only agree on their zero
    point because both were shifted the same way before being split."""
    ni = len(ex_i)
    nj = len(ex_j)
    out = np.zeros(len(ws),dtype=np.complex128) # initialize
    for i in range(ni): # loop over every initial state
      wi = weights[i]
      if wi==0.0: continue # skip states with no Boltzmann weight at all
      for iw in range(len(ws)): # loop over frequencies
        w = ws[iw]
        acc = 0.0+0.0j
        for j in range(nj): # loop over the cropped final-state space
            tmp = A[i,j]*B[j,i] # relevant matrix element
            acc = acc + tmp*(1./(w+1j*delta+ex_i[i]-ex_j[j])
                            - 1./(w-1j*delta+ex_i[i]-ex_j[j]))
        out[iw] = out[iw] + wi*acc
    return out # return dynamical correlator


def dynamical_correlator_inv(h0,wf0,e0,A,B,es=np.linspace(-1,10,600),
        delta=3e-2,mode="cv",**kwargs):
  """Calculate a correlation function AB in a frequency window"""
  ## default method
#  iden = np.identity(h0.shape[0],dtype=np.complex128) # identity
  from scipy.sparse import identity
  iden = identity(h0.shape[0],dtype=np.complex128) # matrix to use
  out = []
  for e in es: # loop over energies
      if mode=="full": # using exact inversion
          g1 = algebra.inv(iden*(e+e0+1j*delta)-h0)
          g2 = algebra.inv(iden*(e+e0-1j*delta)-h0)
          g = 1j*(g1-g2)/2.
          op = A@g@B # operator
          o = algebra.braket_wAw(wf0,op) # correlator
      elif mode=="cv": # correction vector algorithm
          o1 = solve_cv(h0,wf0,A,B,e+e0,delta=delta,**kwargs) # conjugate gradient
          o2 = solve_cv(h0,wf0,A,B,e+e0,delta=-delta,**kwargs) # conjugate gradient
          o = 1j*(o1 - o2)/2. # substract
      else: raise # not recognised
      out.append(o)
  return es,np.array(out)/np.pi # return result



def solve_cv(h0,wf0,si,sj,w,delta=0.0,rtol=1e-6):
    """Solve the dynamical correlator using conjugate gradient method"""
    ## This function may need some benchmarking
    from scipy.sparse import identity
    iden = identity(h0.shape[0],dtype=np.complex128) # matrix to use
    b = -delta*sj@wf0 # create the b vector
    A = (h0 - w*iden)@(h0-w*iden) + iden*delta*delta # define A matrix
    x,info = slg.cg(A,b,rtol=rtol) # solve the equation
    x = 1j*x + (h0 - w*iden)@x/delta # full correction vector
    x = si@x # apply second operator
    o = np.dot(np.conjugate(wf0),x) # compute the braket
    return o





