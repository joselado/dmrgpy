"""Small, from-scratch, dense NumPy Jordan-Wigner reference (local
dimension 2 everywhere), independent of multioperatortk/jordanwigner.py.

Used as an ED ground-truth oracle for mixedchain.Mixed_Spin_Fermion_Chain
(no dedicated ED backend exists for a mixed spin/spinful-fermion Hilbert
space yet, see mixedchain.py's docstring), shared between
test_mixed_chain.py and examples/mixed_spin_fermion_chain/main.py so the
Jordan-Wigner-string construction is only maintained in one place."""
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = 0.5*np.array([[0, 1], [1, 0]], dtype=complex)
SY = 0.5*np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = 0.5*np.array([[1, 0], [0, -1]], dtype=complex)
_A = np.array([[0, 1], [0, 0]], dtype=complex)   # annihilation, basis (|0>,|1>)
_ADAG = _A.T.conj()
_Z = I2 - 2*(_ADAG @ _A)                          # fermion parity


def kron_all(mats):
    out = np.array([[1.0]], dtype=complex)
    for m in mats: out = np.kron(out, m)
    return out


def embed(mat, pos, nsites):
    """Embed a 2x2 operator at physical site `pos` (0-based) of an
    nsites chain, identity elsewhere."""
    mats = [I2]*nsites
    mats[pos] = mat
    return kron_all(mats)


def fermion_ops(fermionic_mask):
    """Dense C/Cdag/N with an explicit Jordan-Wigner string, for a chain
    of local dimension 2 where fermionic_mask[i] is True at a
    fermionic-mode physical site and False at a bosonic/spin site (whose
    parity contribution is the identity -- the same physics
    itensor_version=3's "F"-at-non-fermionic-site fallback implements)."""
    nsites = len(fermionic_mask)
    C, Cdag, N = [None]*nsites, [None]*nsites, [None]*nsites
    for p in range(nsites):
        mats_c, mats_cdag = [I2]*nsites, [I2]*nsites
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
