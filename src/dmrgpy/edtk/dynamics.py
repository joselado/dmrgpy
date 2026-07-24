from ..algebra import algebra
from .. import multioperator
import scipy.sparse.linalg as slg
from ..algebra import kpm
import numpy as np
from .edchain import EDOperator


#from numba import jit

is_hermitian = algebra.is_hermitian


def get_dynamical_correlator(self,name=None,submode="KPM",
        wf0=None,**kwargs):
    """
    Compute the dynamical correlator
    """
    if name is None: raise

    A = EDOperator(name[0],self).SO # create first operator
    B = EDOperator(name[1],self).SO # create second operator
    h = self.get_hamiltonian() # Hamiltonian as matrix

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
    else: raise


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


def dynamical_correlator_ED(h,a0,b0,delta=2e-2,
        emu=None,vs = None,
        dex = 1e-5, # this is a tolerancy to consider something a GS
        es=np.linspace(-1.0,10.0,600)):
    """Compute a dynamical correlator"""
    if emu is None or vs is None: # if not provided
        emu,vs = algebra.eigh(h) # compute them
    ex = emu-np.min(emu) # excitations

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





