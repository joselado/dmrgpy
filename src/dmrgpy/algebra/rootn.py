# Root-N Krylov-space correction-vector method
# Nocera & Alvarez, "Root-N Krylov-space correction-vectors for spectral
# functions with the density matrix renormalization group", arXiv:2204.03165
#
# Standard correction vector: x(w+i eta) = (w-H+E0+i eta)^{-1} O|GS>.
# The root-N idea builds this in N steps instead of one: at each step the
# *current* vector seeds a small Lanczos/Krylov subspace of H, and the
# fractional-power resolvent (w-H+E0+i eta)^{-1/N} is applied within that
# subspace, feeding the result forward as the seed for the next step. With
# N=1 this is exactly the "conventional" Krylov-space correction-vector
# method (Nocera, PRE 2016) that root-N is compared against in the paper.
#
# This module implements the method directly against a dense/sparse
# Hamiltonian (as used by the ED backend) rather than against an MPS: no
# multi-target state-averaging is needed here because there is no bond
# dimension to manage, so this is a numerically exact realization of the
# root-N recursion up to the Krylov-subspace truncation of dimension nkry
# at each step (nkry playing the role the paper's bond dimension m plays
# for the MPS/DMRG case).
import numpy as np


def lanczos_basis(H,v0,k):
    """Build an orthonormal Lanczos basis of size <=k for the Hermitian
    operator H (anything supporting H@v), starting from seed v0, and
    return (Q,T): Q's columns are the basis vectors, T is the resulting
    real symmetric tridiagonal projection of H onto that basis. Full
    reorthogonalization against every previous vector is used (not just
    the classic three-term recurrence) since k is always small (tens of
    vectors) here, so the O(k) extra work per step is negligible and it
    avoids the well known loss of orthogonality of plain Lanczos.
    If an invariant subspace is hit before k vectors are built (beta=0),
    the basis is simply truncated there -- this is the exact subspace
    already, not an error."""
    n = v0.shape[0]
    q = v0/np.linalg.norm(v0)
    Q = [q]
    alphas = []
    betas = []
    beta = 0.0
    q_prev = np.zeros(n,dtype=complex)
    for j in range(k):
        w = H@Q[j] - beta*q_prev
        alpha = np.vdot(Q[j],w).real # real for Hermitian H
        alphas.append(alpha)
        w = w - alpha*Q[j]
        for q_i in Q: # full reorthogonalization
            w = w - np.vdot(q_i,w)*q_i
        beta = np.linalg.norm(w)
        if j==k-1: break # T already has k rows/cols, no need for q_{k}
        if beta<1e-12: break # invariant subspace reached
        betas.append(beta)
        q_prev = Q[j]
        Q.append(w/beta)
    nb = len(Q)
    T = np.diag(alphas)
    for j in range(nb-1):
        T[j,j+1] = betas[j]
        T[j+1,j] = betas[j]
    Q = np.array(Q).T # (n,nb)
    return Q,T


def apply_fractional_resolvent(H,v,omega,e0,eta,N,nkry):
    """Apply (omega-H+e0+i*eta)^{1/N} to v, approximated within a Lanczos
    subspace of dimension nkry seeded by v itself."""
    Q,T = lanczos_basis(H,v,nkry)
    es,V = np.linalg.eigh(T) # T real symmetric -> V real orthogonal
    coords = np.conjugate(Q).T@v # v's coordinates in the Lanczos basis
    gamma = V.T@coords # coordinates in T's eigenbasis (V real, no conjugate needed)
    f = 1./((omega-es+e0+1j*eta)**(1./N))
    gamma = f*gamma
    coords = V@gamma
    return Q@coords


def rootn_correction_vector(H,wf0,e0,A,B,omega,eta,N=8,nkry=20):
    """Correction vector <GS|A (omega-H+E0+i*eta)^{-1} B|GS> computed with
    the root-N Krylov method: N sequential applications of the 1/N-power
    resolvent, each within its own nkry-dimensional Lanczos subspace
    re-seeded from the previous step's vector."""
    v = B@wf0 # phi = O_j|GS>, the p=0 seed
    for p in range(N):
        v = apply_fractional_resolvent(H,v,omega,e0,eta,N,nkry)
    G = np.vdot(wf0,A@v) # <GS|A|x>
    return -G.imag/np.pi
