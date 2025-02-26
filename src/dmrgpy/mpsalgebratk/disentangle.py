import numpy as np
import scipy.linalg as dlg


def disentangle_manifold(wfs,A,**kwargs):
  """
  Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator
  """
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  if A.is_hermitian(): evals,evecs = dlg.eigh(ma) # diagonalize
  else: evals,evecs = dlg.eig(ma) # diagonalize
  # if the Hamiltonian is non-Hermitian and the eigenvectors
  # are not perpendicular, like when there are degeneracies
  # this will not work ok
  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0j
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout



def get_representation(wfs,A,**kwargs):
    """
    Gets the matrix representation of a certain operator
    - plain: A
    - exponential exp(A)
    """
    n = len(wfs) # number of eigenfunctions
    ma = np.zeros((n,n),dtype=np.complex128) # representation of A
    for i in range(n):
        vi = wfs[i]
        for j in range(n):
            vj = wfs[j]
            data = vi.dot(A*vj)
#            if opmode=="exponential": 
#                data = vi.dot(vi.MBO.exponential(A,vj)) # compute exponential
            ma[i,j] = data
    return ma


