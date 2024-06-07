from ..algebra import algebra
from .. import multioperator
import scipy.sparse.linalg as slg
from ..algebra import kpm
import numpy as np
#from numba import jit

is_hermitian = algebra.is_hermitian


def get_dynamical_correlator(self,name=None,submode="KPM",**kwargs):
    """
    Compute the dynamical correlator
    """
    if name is None: raise
    if type(name[0])==multioperator.MultiOperator and type(name[1])==multioperator.MultiOperator: # multioperator
      A = self.get_operator(name[0])
      B = self.get_operator(name[1])
    else:
      raise # this is no longer used
    h = self.get_operator(self.hamiltonian) # Hamiltonian in matrix form
    if not is_hermitian(h): # non Hermitian Hamiltonians
        from ..nonhermitian.dynamics import dynamical_correlator_non_hermitian
        return dynamical_correlator_non_hermitian(self,name=name,**kwargs)
    wf0 = self.get_gs_array() # compute ground state
    # for Hermitian Hamiltonians, continue
    if submode=="KPM":
      return dynamical_correlator_kpm(h,self.e0,wf0,A,B,**kwargs)
    elif submode=="ED":
      return dynamical_correlator_ED(h,A,B,**kwargs)
    elif submode=="INV":
      return dynamical_correlator_inv(h,wf0,self.e0,A,B,**kwargs)
    elif submode=="CVM":
      return dynamical_correlator_inv(h,wf0,self.e0,A,B,mode="cv",**kwargs)
    elif submode=="TD":
      from .. import timedependent
      return timedependent.dynamical_correlator(self,mode="ED",
              name=name,**kwargs)
    else: raise


def dynamical_correlator_kpm(h,e0,wf0,A,B,
        delta=1e-1,es=np.linspace(-1.,10,400)):
    A = np.conjugate(A.T)
    vi = B@wf0 # first wavefunction
    vj = A@wf0 # second wavefunction
    from scipy.sparse import identity
    m = -identity(h.shape[0])*e0+h # matrix to use
    emax = -algebra.lowest_eigenvalues(-m,n=3)[0] # lowest energy
    scale = np.max([np.abs(e0),np.abs(emax)])*3.0
    n = int(2*scale/delta) # number of polynomials
    (xs,ys) = kpm.dm_vivj_energy(m,vi,vj,scale=scale,
                                npol=n*4,ne=n*10,x=es)
    return xs,np.conjugate(ys)*scale/np.pi # return correlator




def dynamical_correlator_ED(h,a0,b0,delta=2e-2,
        es=np.linspace(-1.0,10.0,600)):
    """Compute a dynamical correlator"""
    raise # this is buggy
    emu,vs = algebra.eigh(h)
    U = np.array(vs) # matrix
    Uh = np.conjugate(np.transpose(U)) # Hermitian
    b0 = np.conjugate(b0.T)
    A = Uh@a0@U # get the matrix elements
    B = Uh@b0@U # get the matrix elements
    out = 0.0+es*0.0*1j # initialize
    out = dynamical_sum(emu,es+1j*delta,A,B,out) # perform the summation
    out -= dynamical_sum(emu,es-1j*delta,A,B,out) # perform the summation
    return (es,-out.imag/(2*np.pi)) # return correlator

#@jit(nopython=True)
def dynamical_sum(es,ws,A,B,out):
    """Return the sum giving the dynamical correlator"""
    out = out*0.0 # initialize
    es = es-np.min(es) # remove minimum
    n = len(es) # number of energies
    for iw in range(len(ws)): # loop over frequencies
        i = 0
        for j in range(n): # loop over energies
            tmp = A[i,j]*B[j,i]
            tmp *= 1./(ws[iw]+es[i] - es[j])
            out[iw] = out[iw] + tmp
    return out # return dynamical correlator


def dynamical_correlator_inv(h0,wf0,e0,A,B,es=np.linspace(-1,10,600),
        delta=3e-2,mode="cv"):
  """Calculate a correlation function AB in a frequency window"""
  ## default method
#  iden = np.identity(h0.shape[0],dtype=np.complex_) # identity
  from scipy.sparse import identity
  iden = identity(h0.shape[0],dtype=np.complex_) # matrix to use
  out = []
  for e in es: # loop over energies
      if mode=="full": # using exact inversion
        g1 = algebra.inv(iden*(e+e0+1j*delta)-h0)
        g2 = algebra.inv(iden*(e+e0-1j*delta)-h0)
        g = 1j*(g1-g2)/2./np.pi
        op = A@g@B # operator
        o = algebra.braket_wAw(wf0,op) # correlator
      elif mode=="cv": # correction vector algorithm
          o1 = solve_cv(h0,wf0,A,B,e+e0,delta=delta) # conjugate gradient
          o2 = solve_cv(h0,wf0,A,B,e+e0,delta=-delta) # conjugate gradient
          o = 1j*(o1 - o2)/2. # substract
      else: raise # not recognised
      out.append(o)
  return es,np.array(out)/np.pi # return result



def solve_cv(h0,wf0,si,sj,w,delta=0.0):
    """Solve the dynamical correlator using conjugate gradient method"""
    ## This function may need some benchmarking
    from scipy.sparse import identity
    iden = identity(h0.shape[0],dtype=np.complex_) # matrix to use
    b = -delta*sj@wf0 # create the b vector
    A = (h0 - w*iden)@(h0-w*iden) + iden*delta*delta # define A matrix
    x,info = slg.cg(A,b,tol=1e-10) # solve the equation
    x = 1j*x + (h0 - w*iden)@x/delta # full correction vector
    x = si@x # apply second operator
    o = np.dot(np.conjugate(wf0),x) # compute the braket
    return o





