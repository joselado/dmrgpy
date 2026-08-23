from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import scipy.sparse.linalg as slg
import numpy as np


det = dlg.det

maxsize = 2000
# ARPACK convergence tolerance for the ED path. 0 means machine precision,
# which is also scipy's own eigsh default; see _deflated_lowest_hermitian
# below for why anything looser is actively harmful here.
tol = 0.0 # ARPACK tolerance; 0 = machine precision


# Spare eigenpairs asked of ARPACK beyond what the caller wants, so the last
# requested Ritz pair is not the one being relied on. Used ONLY by the
# non-Hermitian branch of lowest_states, which cannot deflate (see there);
# it was measured to recover degenerate copies that a bare request loses.
# The Hermitian branch does not over-request: once deflation is doing the
# work, spare pairs are pure cost, since a round accepts only its own lowest
# level anyway. Measured on the ferromagnetic chain below, n=3 at dim 4096
# went 1.71 s -> 0.21 s on removing it, with all 18 checked (model, N, n)
# cases still exact -- and ED ground states run through lowest_states(n=3).
_extra_states = 10


def _ncv(k,n):
  """Krylov subspace size for an ARPACK call asking for `k` eigenpairs of an
  `n`-dimensional matrix.

  4*k, floored at scipy's own default of 20 (never go BELOW it: forcing
  ncv=16 for k=4 measured worse than leaving it alone) and capped at the
  matrix dimension."""
  return min(max(4*k,20),n)


def _arpack_k(n,dim):
  """How many eigenpairs to ask ARPACK for, in order to return `n`."""
  return min(n+_extra_states,dim-1)


def _matvec(h,x):
  """h@x as a flat 1D array, whatever h is (sparse, ndarray or matrix)."""
  return np.asarray(h@x).reshape(-1)


def _deflated_lowest_hermitian(h,n):
  """Lowest `n` eigenpairs of a Hermitian `h`, DEGENERATE LEVELS INCLUDED.

  A plain slg.eigsh(h,k=n) silently loses members of a degenerate
  multiplet: it stops once the Ritz pairs it holds are converged, which a
  single copy of a degenerate eigenvalue satisfies, so the remaining copies
  are never produced. Measured on a ferromagnetic spin-1/2 Heisenberg
  chain, whose ground multiplet is exactly (N+1)-fold by SU(2), at N=12
  (dim 4096, above `maxsize` and so on this path): asking for n=4 returned
  1 of the 4 ground-level copies, n=8 returned 3 of 8, n=16 returned 4 of
  13.

  This is the ED path -- the one used to decide whether a DMRG result is
  right -- so dropping levels HERE is worse than dropping them in the
  method under test. It already caused a wrong bug report: mode="ED"
  appeared to return "distinct levels only" while mode="DMRG" returned
  multiplet members, and DMRG was blamed for duplicating states (see
  docs/ed_vs_dmrg_degenerate_multiplets.md).

  Tuning does not fix it. Across (N=12,14) x (n=8,16), tol, ncv and
  over-requesting each repaired some cases and broke others -- ncv=8k gave
  8 of 8 at N=12/n=8 but 6 of 8 at N=14/n=8, and ncv=k+40 did the reverse.
  There is no setting of them that is uniformly right, because none of them
  gives ARPACK a REASON to produce a second copy of an eigenvalue it has
  already converged.

  Deflation does. Each round solves for the lowest level of

      P h P + sigma |V><V|,       P = 1 - |V><V|

  where V spans the eigenvectors found so far and sigma sits far above the
  region of interest. Any remaining copy of an already-found level is now
  the LOWEST eigenvalue of that operator, i.e. exactly the extremal case
  Lanczos is reliable at, so it has to come out. Three details are
  load-bearing, each measured:

    * accept only the lowest verified level per round (plus its exact
      degenerate partners), not every Ritz pair returned. Taking the higher
      ones on faith admits an excited state while a ground copy is still
      missing -- that alone was 6 of 8 instead of 8 of 8.
    * accept a pair only if it passes a real residual check against `h`.
      eigsh returns k Ritz pairs whether or not the trailing ones
      converged.
    * re-orthogonalize against everything already found before accepting.

  Verified elementwise against dense eigvalsh (max error ~5e-14, vectors
  orthonormal and eigen-residual <1e-13) at dim 4096 for n=1,4,8,16 on a
  ferromagnetic chain, an antiferromagnetic one, a random-field chain with
  no degeneracy at all, and an interacting spinless-fermion chain.

  The cost is one ARPACK call per level rather than one in total: n=8 takes
  ~0.8 s at dim 4096 and ~4.4 s at dim 16384. Accepted deliberately -- a
  fast wrong spectrum is worth less than a slow right one on the reference
  path.
  """
  dim = h.shape[0]
  es,vs = [],[] # levels and (orthonormal) eigenvectors found so far
  sigma = None # where deflated levels get parked; set after the first round
  while len(es)<n:
    m = len(es)
    k = min(n-m,dim-1-m) # no over-request: see _extra_states
    if k<1: break
    if m==0: A = h # nothing to deflate yet
    else:
      V = np.array(vs).T # (dim,m)
      Vh = np.conjugate(V.T)
      def mv(x,V=V,Vh=Vh,sg=sigma):
        x = np.asarray(x).reshape(-1)
        c = Vh@x
        y = _matvec(h,x-V@c) # h P x
        y = y - V@(Vh@y) # P h P x
        return y + sg*(V@c) # ... + sigma |V><V| x
      A = slg.LinearOperator((dim,dim),matvec=mv,dtype=np.complex128)
    e,vec = slg.eigsh(A,k=k,which="SA",maxiter=int(maxiter),tol=tol,
                      ncv=_ncv(k,dim))
    got,emin = 0,None
    for j in np.argsort(e.real):
      if sigma is not None and e[j].real>sigma-1.0: continue # a deflated copy
      w = vec[:,j].astype(np.complex128)
      if vs: # re-orthogonalize against everything already found
        V = np.array(vs).T
        w = w - V@(np.conjugate(V.T)@w)
      nw = np.linalg.norm(w)
      if nw<1e-6: continue # numerically inside the found space already
      w = w/nw
      hw = _matvec(h,w)
      ev = np.real(np.conjugate(w)@hw) # Rayleigh quotient
      if np.linalg.norm(hw-ev*w)>1e-7*(1.0+abs(ev)): continue # not converged
      if emin is None: emin = ev # the lowest level of this round
      elif ev>emin+1e-8*(1.0+abs(emin)): break # a higher one: next round
      es.append(ev) ; vs.append(w) ; got += 1
      if len(es)>=n: break
    if sigma is None and len(es)>0: sigma = es[0]+10.0*(1.0+abs(es[0]))
    if got==0: break # no progress this round, stop rather than spin
  if len(es)<n: # say so rather than hand back a short list silently
    import warnings
    warnings.warn("lowest_states: only %d of the %d requested levels "
            "converged; the returned list is short, not padded"%(len(es),n))
  order = np.argsort(np.array(es))
  return np.array([es[i] for i in order]),[vs[i] for i in order]


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
    o = dlg.eigh(m)
    return o

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
      eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=int(maxiter),
                             ncv=_ncv(10,h.shape[0]))
    except:
      eig,eigvec = slg.eigsh(h,k=10,which="SA",maxiter=int(maxiter),tol=tol,
                             ncv=_ncv(10,h.shape[0]))
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
      # Deflation, NOT a single eigsh call: the latter drops members of
      # degenerate levels. See _deflated_lowest_hermitian.
      return _deflated_lowest_hermitian(h,n)
    else: 
      # No deflation here: the projector trick above relies on eigenvectors
      # being orthogonal, which a non-Hermitian matrix's are not. This
      # branch keeps the over-request-and-truncate approach, so a degenerate
      # level can still lose copies -- non-Hermitian ED is not the path used
      # to validate DMRG, and a correct non-Hermitian version needs
      # biorthogonal (left/right) deflation rather than this one.
      k = _arpack_k(n,h.shape[0])
      eig,eigvec = slg.eigs(h,k=k,which="SR",maxiter=int(maxiter),tol=tol,
                            ncv=_ncv(k,h.shape[0]))
      eig,eigvec = sorteigen(eig,eigvec.T)
      return (eig[0:n],eigvec[0:n])
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


def biorthogonal_ground_state(h,**kwargs):
    """Return (e0,vr,vl), the ground state of a (possibly non-Hermitian)
    matrix h as a biorthogonal right/left pair: h@vr = e0*vr and
    dagger(h)@vl = conj(e0)*vl, normalized so that conj(vl)@vr = 1 (the
    convention needed by braket_ww-style expectation values). vr/vl pick
    the state with smallest real part of the eigenvalue, matching
    lowest_states' own non-Hermitian ordering."""
    er,vsr = lowest_states(h,n=1,**kwargs)
    el,vsl = lowest_states(dagger(h),n=1,**kwargs)
    vr = matrix2vector(vsr[0])
    vl = matrix2vector(vsl[0])
    norm = np.conjugate(vl)@vr
    if np.abs(norm)<1e-8:
        raise ValueError("Right and left ground states are (near) "
                "orthogonal, <vl|vr> = "+str(norm)+" - biorthogonal "
                "normalization failed")
    vl = vl/np.conjugate(norm)
    return er[0],vr,vl

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









