from __future__ import print_function
import numpy as np
from . import multioperator


def dynamical_correlator(self,es=np.linspace(-1.,10,100),
        delta=1e-1,name=None,N=8,nkry=20,**kwargs):
    """
    Compute the dynamical correlator using the root-N Krylov-space
    correction-vector method (Nocera & Alvarez, arXiv:2204.03165).

    This mirrors src/dmrgpy/algebra/rootn.py's ED implementation exactly
    (same recursion, same Lanczos-then-diagonalize-then-apply-fractional-
    power recipe), but builds the Krylov subspace out of *global* MPS
    vectors (obtained via repeated MPO application, i.e. Hmpo*v with
    truncation) instead of numpy vectors -- reusing the same backend-
    agnostic MPS algebra (self.toMPO, MPS +/scalar-*/.dot()) that cvm.py
    already uses, so this runs on itensor_version in (2,3,"python")
    exactly like CVM, with no backend-specific code at all.

    Why not a per-bond local sweep (which is how DMRG/TDVP usually
    implement "apply an operator to an MPS")? An earlier attempt tried
    exactly that in mpscpp3/chain_session.h (a TDVP-style single-MPS
    sweep, building each bond's effective Hamiltonian from that same
    MPS's own current tensors) and it produced wrong (sign-flipping,
    unstable) results when cross-checked against algebra/rootn.py's exact
    ED answer on a small chain. The root cause: unlike TDVP (a well-
    defined local Trotter step) or ground-state DMRG (repeated local
    energy minimization provably converges to the global minimum),
    reapplying a Lanczos-based function of H one bond at a time, using
    the *same* partially-updated MPS as both the source of the local
    environment and the thing being progressively overwritten, does not
    correctly realize "apply f(H) once" globally -- later bonds see a
    mix of already-updated and not-yet-updated neighbors. Getting a
    per-bond version right would need something like ITensor's own
    fitApplyMPO (mpoalgs.cc), which keeps the fixed input and the
    evolving ansatz as two separate objects with their own environments
    -- but that machinery is written for an explicit site-factorized MPO
    operator, not a Krylov-subspace-approximated nonlinear function of H,
    so it doesn't transfer directly either. The global-Krylov approach
    here sidesteps all of this by only ever using already-tested,
    already-correct whole-MPS primitives (MPO application, inner
    products, MPS addition) -- exactly the ones cvm.py's conjugate
    gradient already relies on -- at the cost of every Lanczos step
    costing a full (truncated) MPO application rather than a cheap local
    tensor contraction.
    """
    if delta<0.0: raise
    if type(name[0])!=multioperator.MultiOperator: raise
    A,B = name[0],name[1] # <GS|A(...)B|GS>, A used as-is (no dagger),
                          # matching cvm.py's / algebra.rootn's convention
    wf0 = self.get_gs() # ground state (also sets self.e0)
    Hmpo = self.toMPO(self.hamiltonian) # built once, shared across frequencies/steps
    out = []
    for e in es:
        o = rootn_correction_vector(self,A,B,e,delta,wf0=wf0,Hmpo=Hmpo,
                N=N,nkry=nkry)
        out.append(o)
    return es,np.array(out)


def rootn_correction_vector(self,A,B,omega,eta,N=8,nkry=20,
        wf0=None,Hmpo=None):
    """
    Correction vector <GS|A(omega-H+E0+i*eta)^{-1}B|GS>, built as N
    sequential fractional-power Lanczos steps (see module docstring).
    """
    if wf0 is None: wf0 = self.get_gs()
    if Hmpo is None: Hmpo = self.toMPO(self.hamiltonian)
    e0 = self.e0
    v = B*wf0 # phi = O_j|GS>, the p=0 seed
    for p in range(N):
        v = _apply_fractional_resolvent_mps(self,Hmpo,v,omega,e0,eta,N,nkry)
    G = wf0.dot(A*v)
    return -G.imag/np.pi


def _lanczos_basis_mps(self,Hmpo,v,nkry):
    """MPS analogue of algebra/rootn.py's lanczos_basis: build an
    orthonormal Lanczos basis of MPS vectors for Hmpo seeded by v, and
    the resulting real tridiagonal (alphas,betas). Every MPO application
    (Hmpo*q) truncates via the chain's ordinary self.maxm/self.cutoff,
    exactly like every other MPS operation in this codebase -- nkry plays
    the role of the paper's bond dimension m (see algebra/rootn.py's own
    module docstring)."""
    nrm = np.sqrt(v.dot(v).real)
    q = (1./nrm)*v
    Q = [q]
    alphas = []
    betas = []
    qprev = None
    beta = 0.0
    for it in range(nkry):
        w = Hmpo*Q[-1]
        if qprev is not None: w = w - beta*qprev
        alpha = Q[-1].dot(w).real # real for Hermitian H
        alphas.append(alpha)
        w = w - alpha*Q[-1]
        for qi in Q: w = w - qi.dot(w)*qi # full reorthogonalization
        if it==nkry-1: break
        beta = np.sqrt(abs(w.dot(w).real))
        if beta<1e-10: break # invariant subspace reached
        betas.append(beta)
        qprev = Q[-1]
        Q.append((1./beta)*w)
    return Q,np.array(alphas),np.array(betas)


def _apply_fractional_resolvent_mps(self,Hmpo,v,omega,e0,eta,N,nkry):
    """Apply (omega-H+e0+i*eta)^{1/N} to the MPS v, approximated within a
    Lanczos subspace of dimension nkry seeded by v itself -- same
    recipe/derivation as algebra/rootn.py's apply_fractional_resolvent,
    operating on MPS objects instead of numpy vectors."""
    nrm = np.sqrt(v.dot(v).real)
    Q,alphas,betas = _lanczos_basis_mps(self,Hmpo,v,nkry)
    m = len(alphas)
    T = np.diag(alphas)
    for i in range(m-1):
        T[i,i+1] = betas[i]
        T[i+1,i] = betas[i]
    es,V = np.linalg.eigh(T) # T real symmetric -> V real orthogonal
    coords = np.zeros(m,dtype=complex)
    coords[0] = nrm # v's coordinates in this basis: Q[0]=v/nrm by construction
    gamma = V.T@coords # V real, no conjugate needed
    f = 1./((omega-es+e0+1j*eta)**(1./N))
    gamma = f*gamma
    coords = V@gamma
    vnew = coords[0]*Q[0]
    for i in range(1,m): vnew = vnew + coords[i]*Q[i]
    return vnew
