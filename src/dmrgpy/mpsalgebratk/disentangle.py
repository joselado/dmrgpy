import numpy as np
import scipy.linalg as dlg


def _is_hermitian(wfs,A):
  """Whether A is Hermitian, decided the way every other Hermiticity gate
  in dmrgpy decides it: by the chain's own is_hermitian (the canonical-form
  proof first, then a random-witness probe), falling back to the bare
  proof only when the states carry no chain. The bare proof alone is
  one-sided, "not proven" is not "not Hermitian": a same-site identity
  (1j*Sx0*Sy0 is exactly -Sz0/2), a term mixing a C-type and an A-type
  name, or any name off canonical.py's _PARITY table (every parafermionic
  operator) is Hermitian and unprovable (2026-09-24b audit, finding 1).
  An ED State's MBO is its EDchain, whose is_hermitian decides exactly on
  the operator's matrix; before it had one, every ED manifold raised
  AttributeError here (2026-09-24c audit, finding 18)."""
  mbo = getattr(wfs[0],"MBO",None)
  if mbo is not None and hasattr(mbo,"is_hermitian"): return mbo.is_hermitian(A)
  return A.is_hermitian()


def disentangle_manifold(wfs,A,**kwargs):
  """
  Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator
  """
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  if _is_hermitian(wfs,A):
      # the Hermitian part of ma: on a truncated manifold ma is Hermitian
      # only to the truncation level, and eigh would otherwise read one
      # triangle of it and ignore the other
      evals,evecs = dlg.eigh((ma + ma.conj().T)/2.)
  else: evals,evecs = dlg.eig(ma) # diagonalize
  # eig does not orthogonalize inside a degenerate eigenspace, so for a
  # Hermitian A with a degenerate spectrum on the manifold it returns a
  # basis that is non-orthonormal by an O(1), roundoff-dependent amount
  # (0.32 to 0.79 measured). That is why the branch above must see every
  # Hermitian A; this one is right only for a genuinely non-Hermitian A,
  # whose eigenvectors need not be orthogonal to begin with.
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


