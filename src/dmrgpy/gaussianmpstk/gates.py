"""Two-site gates on an MPS held as site arrays.

Each array has the shape (left link, site, right link), the layout of
mps_site_arrays/mps_from_site_arrays on every backend
(mpsalgebratk/hadamard.py). The orthogonality centre is moved onto the
bond of a gate before its SVD, so the truncation of every gate is the
optimal one for the whole state."""

import numpy as np


class SiteArrays:
    """An MPS as a list of arrays with a known orthogonality centre."""

    def __init__(self, arrays, center=0):
        self.A = [np.array(a, dtype=complex) for a in arrays]
        self.center = center
        self.discarded = 0.0  # summed discarded weight of all the SVDs

    @classmethod
    def product(cls, vectors):
        """A normalized product state from one local vector per site: every
        site array is then both left and right orthonormal."""
        arrays = []
        for v in vectors:
            v = np.asarray(v, dtype=complex)
            arrays.append((v/np.linalg.norm(v)).reshape(1, len(v), 1))
        return cls(arrays, center=0)

    def move_center(self, j):
        A = self.A
        while self.center < j:
            c = self.center
            Dl, d, Dr = A[c].shape
            Q, R = np.linalg.qr(A[c].reshape(Dl*d, Dr))
            A[c] = Q.reshape(Dl, d, Q.shape[1])
            A[c+1] = np.tensordot(R, A[c+1], axes=(1, 0))
            self.center += 1
        while self.center > j:
            c = self.center
            Dl, d, Dr = A[c].shape
            Q, R = np.linalg.qr(A[c].reshape(Dl, d*Dr).T)
            A[c] = Q.T.reshape(Q.shape[1], d, Dr)
            A[c-1] = np.tensordot(A[c-1], R.T, axes=(2, 0))
            self.center -= 1

    def apply(self, G, j, maxm, cutoff):
        """Apply the two-site gate G (rows the output, index d*s_j+s_{j+1})
        on sites j, j+1, keeping at most maxm states and dropping a
        discarded weight up to cutoff, relative to the weight of the
        state. The centre ends on site j+1."""
        self.move_center(j)
        A = self.A
        Dl, d, _ = A[j].shape
        _, d2, Dr = A[j+1].shape
        theta = np.tensordot(A[j], A[j+1], axes=(2, 0))  # (Dl,d,d2,Dr)
        G4 = np.asarray(G).reshape(d, d2, d, d2)
        theta = np.tensordot(G4, theta, axes=([2, 3], [1, 2]))  # (d,d2,Dl,Dr)
        theta = theta.transpose(2, 0, 1, 3).reshape(Dl*d, d2*Dr)
        U, S, Vh = np.linalg.svd(theta, full_matrices=False)
        w = S**2
        total = np.sum(w)
        if total == 0.0:
            raise ValueError("a gate annihilated the state")
        tail = np.cumsum(w[::-1])[::-1]/total  # tail[k]: weight from k on
        keep = len(S)
        for k in range(1, len(S)):
            if tail[k] <= cutoff:
                keep = k
                break
        keep = max(1, min(keep, int(maxm)))
        self.discarded += float(np.sum(w[keep:])/total)
        A[j] = U[:, :keep].reshape(Dl, d, keep)
        A[j+1] = (S[:keep, None]*Vh[:keep, :]).reshape(keep, d2, Dr)
        self.center = j+1

    def norm(self):
        return float(np.linalg.norm(self.A[self.center]))

    def normalize(self):
        nrm = self.norm()
        if not np.isfinite(nrm) or nrm == 0.0:
            raise ValueError("the state has zero norm")
        self.A[self.center] = self.A[self.center]/nrm
        return nrm

    def maxlinkdim(self):
        return max(a.shape[2] for a in self.A)
