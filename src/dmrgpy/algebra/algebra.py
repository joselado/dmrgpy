from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import numpy as np


def braket_wAw(w,A,wi=None):
  """
  Compute the braket of a wavefunction
  """
  if wi is None: wi = w
  if issparse(A): # sparse matrices
    return (np.conjugate(wi)@A@w) # modern way
  else: # matrices and arrays
    return (np.conjugate(wi)@A@w)[0,0] # modern way




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
      data = (vi*sa*vj).todense()[0,0]
      ma[i,j] = data
  return ma





## routines for diagonalization ##

error = 1e-7



accelerate = False

def eigh(m):
    """Wrapper for linalg"""
    from . import algebraf90
    if not accelerate: return dlg.eigh(m)
    # check if doing slices helps
    n = m.shape[0] # size of the matrix
    mo = m[0:n:2,1:n:2] # off diagonal is zero
#    if False: # assume block diagonal
    if np.max(np.abs(mo))<error: # assume block diagonal
        # detected block diagonal
        (es0,vs0) = eigh(m[0:n:2,0:n:2]) # recall
        (es1,vs1) = eigh(m[1:n:2,1:n:2]) # recall
        es = np.concatenate([es0,es1]) # concatenate array
        vs0 = algebraf90.todouble(vs0.T,0)
        vs1 = algebraf90.todouble(vs1.T,1)
        vs = np.concatenate([vs0,vs1])
        return (es,vs.T) # return the eigenvaleus and eigenvectors

    else:
      if np.max(np.abs(m.imag))<error: # assume real
          return dlg.eigh(m.real) # diagonalize real matrix
      else: return dlg.eigh(m) # diagonalize complex matrix


def eigvalsh(m):
    """Wrapper for linalg"""
    if not accelerate: return dlg.eigvalsh(m)
    # check if doing slices helps
    n = m.shape[0] # size of the matrix
    mo = m[0:n:2,1:n:2] # off diagonal is zero
#    if False: # assume block diagonal
    if np.max(np.abs(mo))<error: # assume block diagonal
        # detected block diagonal
        es0 = eigvalsh(m[0:n:2,0:n:2]) # recall
        es1 = eigvalsh(m[1:n:2,1:n:2]) # recall
        es = np.concatenate([es0,es1]) # concatenate array
        return es

    else:
      if np.max(np.abs(m.imag))<error: # assume real
          return dlg.eigvalsh(m.real) # diagonalize real matrix
      else: return dlg.eigvalsh(m) # diagonalize complex matrix



def matrix2vector(v):
    """Transform a matrix into a vector"""
    if issparse(v): # sparse matrix
      v = v.todense() # convert to conventional matrix
    v = np.array(v) # convert to array
    if len(v.shape)==1: return v
    else: return v.reshape(v.shape[0]*v.shape[1])


