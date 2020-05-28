from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import scipy.sparse.linalg as slg
import numpy as np


maxsize = 3000



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
  ma = np.zeros((n,n),dtype=np.complex) # representation of A
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
    eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=100000)
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


def lowest_eigenvalues(h,n=10,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    eig,eigvec = slg.eigsh(h,k=n,which="SA",maxiter=100000)
    eig = np.sort(eig)
  else:
    if info: print("Full diagonalization")
    ishermitian(h)
    eig = dlg.eigvalsh(h.todense())
  return eig[0:n] # return eigenvalues



def ishermitian(m):
    d = m - np.conjugate(m.T)
    if np.max(np.abs(d))>1e-6: 
        print("Hamiltonian is not Hermitian")
        raise

def expm(m):
    m = todense(m)
    return dlg.expm(m)


def inv(m):
    """Inverse"""
    m = todense(m)
    return dlg.inv(m)



def lowest_eigenvectors(h,n=10,nmax=maxsize):
  """Get a ground state"""
  info = False
  if h.shape[0]>nmax:
    if info: print("Calling ARPACK")
    eig,eigvec = slg.eigsh(h,k=n,which="SA",maxiter=100000)
#    eigvec = [v for (e,v) in zip(eig,eigvec.T)]
  else:
    if info: print("Full diagonalization")
    eig,eigvec = dlg.eigh(h.todense())
    eigvec = eigvec.T # transpose
#  print(sorted(eig))
#  eigevec = [v for (e,v) in sorted(zip(eig,eigvec))]
  return eigvec[0:n] # return eigenvectors


def expm(m):
    m = todense(m)
    return dlg.expm(m) # exponential matrix

