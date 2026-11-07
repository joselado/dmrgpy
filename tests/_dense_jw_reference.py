"""Small, from-scratch, dense NumPy Jordan-Wigner reference, independent
of multioperatortk/jordanwigner.py and of the mixed-chain ED backend
(edtk/mixededchain.py).

Used as an ED ground-truth oracle for mixedchain.Mixed_Spin_Fermion_Chain
in test_mixed_chain.py and test_mixed_chain_ed.py, so the
Jordan-Wigner-string construction of the oracle is only maintained in one
place. Every site has local dimension 2 unless `dims` gives the local
dimension of each physical site (3 for a spin 1)."""
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = 0.5*np.array([[0, 1], [1, 0]], dtype=complex)
SY = 0.5*np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = 0.5*np.array([[1, 0], [0, -1]], dtype=complex)
_A = np.array([[0, 1], [0, 0]], dtype=complex)   # annihilation, basis (|0>,|1>)
_ADAG = _A.T.conj()
_Z = I2 - 2*(_ADAG @ _A)                          # fermion parity


def spin_one():
    """Dense Sx, Sy, Sz of a spin 1, basis m = 1, 0, -1"""
    r = 1/np.sqrt(2)
    sx = r*np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=complex)
    sy = r*np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]], dtype=complex)
    sz = np.diag([1., 0., -1.]).astype(complex)
    return sx, sy, sz


def kron_all(mats):
    out = np.array([[1.0]], dtype=complex)
    for m in mats: out = np.kron(out, m)
    return out


def _identities(nsites, dims):
    if dims is None: dims = [2]*nsites
    return [np.eye(d, dtype=complex) for d in dims]


def embed(mat, pos, nsites, dims=None):
    """Embed a local operator at physical site `pos` (0-based) of an
    nsites chain, identity elsewhere."""
    mats = _identities(nsites, dims)
    mats[pos] = mat
    return kron_all(mats)


def fermion_ops(fermionic_mask, dims=None):
    """Dense C/Cdag/N with an explicit Jordan-Wigner string, where
    fermionic_mask[i] is True at a fermionic-mode physical site and False
    at a spin site (whose parity contribution is the identity -- the same
    physics itensor_version=3's "F"-at-non-fermionic-site fallback
    implements). Entries at spin sites are None."""
    nsites = len(fermionic_mask)
    C, Cdag, N = [None]*nsites, [None]*nsites, [None]*nsites
    for p in range(nsites):
        if not fermionic_mask[p]: continue
        mats_c, mats_cdag = _identities(nsites, dims), _identities(nsites, dims)
        for k in range(p):
            if fermionic_mask[k]:
                mats_c[k] = _Z
                mats_cdag[k] = _Z
        mats_c[p] = _A
        mats_cdag[p] = _ADAG
        C[p] = kron_all(mats_c)
        Cdag[p] = kron_all(mats_cdag)
        N[p] = Cdag[p] @ C[p]
    return C, Cdag, N


def gs(H):
    """Ground-state energy and eigenvector of a dense Hermitian matrix"""
    es, vs = np.linalg.eigh(H)
    return es[0], vs[:, 0]
