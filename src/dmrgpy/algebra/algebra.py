from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import scipy.sparse.linalg as slg
import numpy as np


det = dlg.det

maxsize = 1000
tol = 1e-7 # tolerancy
maxiter = 1e6 # maximum number of iterations



def braket_wAw(w,A,wi=None):
  """
  Compute the braket of a wavefunction
  """
  if wi is None: wi = w
  if issparse(A): # sparse matrices
    return (np.conjugate(wi)@A@w) # modern way
  else: # matrices and arrays
      if len(w.shape)==1: return (np.conjugate(wi)@A@w) # modern way
      else: return (np.conjugate(wi)@A@w)[0,0] # modern way



def dagger(m):
  return np.conjugate(m.T)


def braket_ww(w,wi):
  """
  Compute the braket of two wavefunctions
  """
  w = matrix2vector(w) # convert to vector
  wi = matrix2vector(wi) # convert to vector
  return (np.conjugate(w)@wi) # modern way




def disentangle_manifold(wfs,A):
  """
  Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator
  """
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  evals,evecs = dlg.eigh(ma) # diagonalize

  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0j
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout



def get_representation(wfs,A):
    """
    Gets the matrix representation of a certain operator
    """
    n = len(wfs) # number of eigenfunctions
    ma = np.zeros((n,n),dtype=np.complex128) # representation of A
    sa = csc(A) # sparse matrix
    for i in range(n):
      vi = csc(np.conjugate(wfs[i])) # first wavefunction
      for j in range(n):
        vj = csc(wfs[j]).transpose() # second wavefunction
        data = (vi@sa@vj).todense()[0,0]
        ma[i,j] = data
    return ma





## routines for diagonalization ##

error = 1e-7



accelerate = False

def eigh(m):
    """Wrapper for linalg"""
    m = todense(m)
    return dlg.eigh(m)

def eigvalsh(m):
    """Wrapper for linalg"""
    m = todense(m)
    return dlg.eigvalsh(m)


def matrix2vector(v):
    """Transform a matrix into a vector"""
    if issparse(v): # sparse matrix
      v = v.todense() # convert to conventional matrix
    v = np.array(v) # convert to array
    if len(v.shape)==1: return v
    else: return v.reshape(v.shape[0]*v.shape[1])


def ground_state(h,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    try: 
      eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=int(maxiter))
    except:
      eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=int(maxiter),tol=tol)
    eig = np.sort(eig)
  else:
    if info: print("Full diagonalization")
    eig,eigvec = dlg.eigh(todense(h))
  return eig[0],eigvec.transpose()[0]


def todense(m):
    """Turn a matrix dense"""
    if issparse(m):
        if m.shape[0]>maxsize: raise
        else: return m.todense()
    else: return m


def lowest_eigenvalues(h,n=10):
  """Get a ground state"""
  info = False
  if h.shape[0]>maxsize: # for sparse use arpack
      eig,vs = lowest_states(h,n=n)
  else:
    if info: print("Full diagonalization")
    if ishermitian(h):
      eig = dlg.eigvalsh(h.todense())
    else:
      eig = dlg.eigvals(h.todense())
      eig,eigvec = sorteigen(eig,eig)
  return np.array(eig[0:n])


def lowest_states(h,n=10,**kwargs):
  """Get a ground state"""
  nmax = maxsize
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    if ishermitian(h): # Hermitian matrix
      eig,eigvec = slg.eigsh(h,k=n,which="SA",maxiter=int(maxiter),tol=tol)
      eig,eigvec = sorteigen(eig,eigvec.T)
      return (eig,eigvec)
    else: 
      eig,eigvec = slg.eigs(h,k=n,which="SR",maxiter=int(maxiter),tol=tol)
      eig,eigvec = sorteigen(eig,eigvec.T)
      return (eig,eigvec)
  else:
    if info: print("Full diagonalization")
    if ishermitian(h): # Hermitian matrix
      eig,vs = dlg.eigh(h.todense())
      return eig[0:n],vs.T[0:n] 
    else: # non Hermitian matrix
      eig,vs = dlg.eig(h.todense())
      eig,vs = sorteigen(eig,vs.T)
      return eig[0:n],vs[0:n]

lowest_eigenvectors = lowest_states

def sorteigen(eig,vs):
    """Return sorted eigenvalues and eigenvectors"""
    w = eig - np.min(eig.real) # smallest real part
    imweight = 1e-3 # how much to weight im when sorting 
    wz = w.real + imweight*w.imag
    vs = [y for (x,y) in sorted(zip(wz,vs),key=lambda x: x[0])]
    eig = [y.copy() for (x,y) in sorted(zip(wz,eig),key=lambda x: x[0])]
    return np.array(eig),vs



def ishermitian(m):
    """Check if a matrix is Hermitian"""
    d = m - np.conjugate(m.T)
    if np.max(np.abs(d))>1e-6: return False
    return True

def expm(m):
    """Compute exponential"""
    m = todense(m)
    return dlg.expm(m)


def inv(m):
    """Inverse"""
    m = todense(m)
    return dlg.inv(m)



def expm(m):
    m = todense(m)
    return dlg.expm(m) # exponential matrix

#def expm(m):
#    m = todense(m)
#    es,vs = dlg.eig(m)
#    d = np.zeros(m.shape,dtype=np.complex128)
#    for i in range(len(es)): d[i,i] = np.exp(es[i])
#    R = vs.T
#    Rh = np.conjugate(R.T)
#    U = Rh@d@R
#    return U



def ismatrix(m):
    return type(m)==np.ndarray or issparse(m) or type(m)==np.matrix





def smooth_gauge(w1,w2):
  """Perform a gauge rotation so that the second set of waves are smooth
  with respect to the first one"""
  m = uij(w1,w2) # matrix of wavefunctions
  U, s, V = np.linalg.svd(m, full_matrices=True) # sing val decomp
  R = np.conjugate(U@V).T # rotation matrix
  wnew = [w.copy()*0. for w in w2] # new WF
  wold = [w.copy() for w in w2] # old WF
  for ii in range(R.shape[0]):
    for jj in range(R.shape[0]):
      wnew[ii] = wnew[ii] + R[jj,ii]*wold[jj]
  return wnew



def uij(wf1,wf2):
  m = np.matrix(np.zeros((len(wf1),len(wf2)),dtype=np.complex128))
  for i in range(len(wf1)):
    for j in range(len(wf2)):
      m[i,j] = wf1[i].dot(wf2[j])
  return m


def trace(A):
    """Compute trace"""
    if issparse(A): return A.trace()
    else: return np.trace(A)



def is_hermitian(h,**kwargs):
    h = h - dagger(h) # difference
    return is_zero_matrix(h,**kwargs)


def is_zero_matrix(h,tol=1e-8):
    h = h@dagger(h)
    t = np.abs(trace(h)) # this should be a trace
    return t<tol



def applyinverse(A,b):
    """Apply A^-1 to b"""
    if A.shape[0]<30: return inv(A)@b
    else: return slg.spsolve(A,b)









