"""The circuit of Givens rotations of a fermionic Gaussian state, from its
correlation matrix (Fishman and White, PRB 92, 075132, arXiv:1504.07701).

The algorithm follows ITensorGaussianMPS.jl (Apache-2.0, commit 6761dd6,
src/gmps.jl: givens_rotations, isolate_subblock_eig and
correlation_matrix_to_gmps, number-conserving branch), rewritten here in
NumPy with the conventions derived and tested on their own:

  Lambda[i,j] = <c_i^dagger c_j>, a Hermitian matrix whose eigenvalues
  are 0 or 1 for a pure Gaussian state (a Slater determinant).

For each mode i we take the smallest block Lambda[i:i+B,i:i+B] that has an
eigenvalue within eigval_cutoff of 0 or 1, and rotate its eigenvector onto
mode i with B-1 nearest-neighbor Givens rotations. Mode i is then
decoupled from the rest up to eigval_cutoff, and we go on with mode i+1.
At the end V Lambda V^dagger is diagonal, with the occupations ns on the
diagonal, V being the product of all the rotations.

A rotation (n, c, s) is the 2x2 matrix g = [[c, s], [-conj(s), c]] acting
on modes (n, n+1), c real. Its many-body unitary U_g maps
c_n^dagger -> c c_n^dagger + s c_{n+1}^dagger and
c_{n+1}^dagger -> -conj(s) c_n^dagger + c c_{n+1}^dagger. With
V = g_K ... g_1 (g_1 found first), the state is
U_{g_1} U_{g_2} ... U_{g_K} |ns>, so the gates act on the product state
in the reverse of the order they were found."""

import numpy as np


def givens(a, b):
    """(c, s, r) with [[c, s], [-conj(s), c]] @ [a, b] = [r, 0], c real
    and non-negative."""
    ra, rb = abs(a), abs(b)
    if rb == 0.0:
        return 1.0, 0.0, a
    if ra == 0.0:
        return 0.0, 1.0, b
    rho = np.hypot(ra, rb)
    phase = a/ra
    return ra/rho, phase*np.conj(b)/rho, phase*rho


def givens_rotations(v):
    """Rotations (n, c, s), in the order they are applied to v, that bring
    the vector v onto its first entry: the product of all of them, last
    one leftmost, maps v to r e_0 with |r| = |v|."""
    v = np.array(v, dtype=complex)
    rotations = []
    for n in range(len(v)-2, -1, -1):
        c, s, r = givens(v[n], v[n+1])
        v[n], v[n+1] = r, 0.0
        rotations.append((n, c, s))
    return rotations


def rotate(L, n, c, s):
    """g L g^dagger for the rotation (n, c, s), in place."""
    g = np.array([[c, s], [-np.conj(s), c]])
    L[n:n+2, :] = g @ L[n:n+2, :]
    L[:, n:n+2] = L[:, n:n+2] @ g.conj().T
    return L


def _occupation_error(x):
    return min(abs(x), abs(1.0-x))


def correlation_matrix_to_circuit(L, eigval_cutoff=1e-8, maxblocksize=None,
                                  minblocksize=1):
    """(ns, rotations, error) for the correlation matrix L.

    ns: the occupations on the diagonal of V L V^dagger (close to 0 or 1).
    rotations: the list of (n, c, s), in the order they were found.
    error: the largest distance to 0 or 1 of an eigenvalue taken to
        decouple a mode, at most eigval_cutoff unless a block reached
        maxblocksize or the end of the chain first."""
    L = np.array(L, dtype=complex)
    N = L.shape[0]
    if maxblocksize is None:
        maxblocksize = N
    rotations = []
    error = 0.0
    for i in range(N):
        for B in range(minblocksize, maxblocksize+1):
            j = min(i+B, N)
            n, u = np.linalg.eigh(L[i:j, i:j])
            errs = [_occupation_error(x) for x in n]
            p = int(np.argmin(errs))
            if errs[p] <= eigval_cutoff or j == N:
                break
        error = max(error, errs[p])
        block = givens_rotations(u[:, p])
        for (k, c, s) in block:
            rotate(L, i+k, c, s)
            rotations.append((i+k, c, s))
    ns = np.real(np.diag(L))
    return ns, rotations, error


def gate_matrix(c, s):
    """The 4x4 matrix of U_g on two neighboring modes, rows the output,
    in the occupation basis |n_a n_b> ordered as 2*n_a+n_b.

    U_g = 1 + (c-1)(N_a+N_b-2 N_a N_b) + s Cdag_b C_a - conj(s) Cdag_a C_b,
    and on neighboring modes Cdag_b C_a = Adag_b A_a and
    Cdag_a C_b = Adag_a A_b, with no Jordan-Wigner sign."""
    G = np.zeros((4, 4), dtype=complex)
    G[0, 0] = 1.0
    G[3, 3] = 1.0
    G[2, 2] = c          # |10> -> c|10> + s|01>
    G[1, 2] = s
    G[1, 1] = c          # |01> -> -conj(s)|10> + c|01>
    G[2, 1] = -np.conj(s)
    return G


# Gaussian states with pairing (Bogoliubov-de Gennes ground states).
#
# The operators are interleaved, gamma = (c_0, c_0^dag, c_1, c_1^dag, ...),
# and the generalized correlation matrix is G[a,b] = <gamma_a^dag gamma_b>,
# so G[2i,2j] = <c_i^dag c_j> and G[2i,2j+1] = <c_i^dag c_j^dag>. Particle-
# hole symmetry reads X G X = 1 - conj(G), X the swap of 2i and 2i+1.
#
# A 2N x 2N unitary Q of the same symmetry (conj(Q) = X Q X) maps G to
# Q G Q^dagger, the correlation matrix of gamma'_b = sum_d conj(Q[b,d])
# gamma_d = U_Q gamma_b U_Q^dagger, so that a G' that is diagonal with
# occupations n gives |psi> = U_Q |n>, and U_{Q2 Q1} = U_{Q1} U_{Q2} as in
# the number-conserving case. Two kinds of local Q are used, on the four
# operators (c_k, c_k^dag, c_{k+1}, c_{k+1}^dag):
#   "hop"  (c, s): g = [[c, s], [-conj(s), c]] on (c_k, c_{k+1}) and
#          conj(g) on (c_k^dag, c_{k+1}^dag), the Givens rotation above;
#   "pair" (c, s): the same 2x2 matrix on (c_k, c_{k+1}^dag) and its
#          conjugate on (c_k^dag, c_{k+1}), a Bogoliubov rotation.
# The construction is the one of ITensorGaussianMPS.jl's parity branch
# (givens_rotations(::ConservesNfParity)), derived again here: the
# quasiparticle d = sum_k (a_k c_k + b_k c_k^dag) that annihilates the
# state obeys d^2 = sum_k a_k b_k = 0, which is what lets one Bogoliubov
# rotation finish the job after the Givens rotations.


def _local_q(kind, c, s):
    Q = np.zeros((4, 4), dtype=complex)
    if kind == "hop":
        Q[0, 0], Q[0, 2], Q[2, 0], Q[2, 2] = c, s, -np.conj(s), c
        Q[1, 1], Q[1, 3], Q[3, 1], Q[3, 3] = c, np.conj(s), -s, c
    elif kind == "pair":
        Q[0, 0], Q[0, 3], Q[3, 0], Q[3, 3] = c, s, -np.conj(s), c
        Q[1, 1], Q[1, 2], Q[2, 1], Q[2, 2] = c, np.conj(s), -s, c
    else:
        raise ValueError("unknown gate kind %r" % (kind,))
    return Q


def rotate_bdg(G, kind, k, c, s):
    """Q G Q^dagger for a local gate on modes k, k+1, in place."""
    Q = _local_q(kind, c, s)
    sl = slice(2*k, 2*k+4)
    G[sl, :] = Q @ G[sl, :]
    G[:, sl] = G[:, sl] @ Q.conj().T
    return G


def _apply_local_q(w, kind, k, c, s):
    sl = slice(2*k, 2*k+4)
    w[sl] = _local_q(kind, c, s) @ w[sl]


def bdg_rotations(w):
    """Gates (kind, k, c, s), in the order they are applied, that bring the
    coefficient vector w of an annihilator (interleaved, length 2B) onto
    c_0 or onto c_0^dag. The second case is a filled mode 0."""
    w = np.array(w, dtype=complex)
    B = len(w)//2
    gates = []

    def add(kind, k, c, s):
        _apply_local_q(w, kind, k, c, s)
        gates.append((kind, k, c, s))

    for k in range(B-2, -1, -1):  # the c coefficients onto mode 0
        c, s, _ = givens(w[2*k], w[2*k+2])
        add("hop", k, c, s)
    if abs(w[1]) > abs(w[0]):  # then w[0] = 0: d is a creation operator
        for k in range(B-2, -1, -1):  # the c^dag coefficients onto mode 0
            c, s, _ = givens(w[2*k+1], w[2*k+3])
            add("hop", k, c, np.conj(s))
        return gates
    for k in range(B-2, 0, -1):  # the c^dag coefficients onto mode 1
        c, s, _ = givens(w[2*k+1], w[2*k+3])
        add("hop", k, c, np.conj(s))
    if B >= 2:  # w[1] = 0 here, since w[0]*w[1] = 0
        c, s, _ = givens(w[0], w[3])
        add("pair", 0, c, s)
    return gates


def _swap_pairs(v):
    out = np.empty_like(v)
    out[0::2], out[1::2] = v[1::2], v[0::2]
    return out


def bdg_correlation_matrix_to_circuit(G, eigval_cutoff=1e-8, maxblocksize=None,
                                      minblocksize=1):
    """(ns, gates, error) for the generalized correlation matrix G, as
    correlation_matrix_to_circuit does for a number-conserving state;
    gates are (kind, n, c, s), kind "hop" or "pair"."""
    G = np.array(G, dtype=complex)
    N = G.shape[0]//2
    if maxblocksize is None:
        maxblocksize = N
    gates = []
    error = 0.0
    for i in range(N):
        for B in range(minblocksize, maxblocksize+1):
            j = min(i+B, N)
            n, u = np.linalg.eigh(G[2*i:2*j, 2*i:2*j])
            errs = [_occupation_error(x) for x in n]
            p = int(np.argmin(errs))
            if errs[p] <= eigval_cutoff or j == N:
                break
        error = max(error, errs[p])
        w = u[:, p]
        if n[p] > 0.5:  # the particle-hole partner, of occupation 1-n
            w = _swap_pairs(w.conj())
        for (kind, k, c, s) in bdg_rotations(w):
            rotate_bdg(G, kind, i+k, c, s)
            gates.append((kind, i+k, c, s))
    ns = np.real(np.diag(G))[0::2]
    return ns, gates, error


def pair_gate_matrix(c, s):
    """The 4x4 matrix of U_Q of a "pair" gate on neighboring modes a, b,
    in the basis 2*n_a+n_b: U c_a U^dag = c c_a + conj(s) c_b^dag. It acts
    on |00>, |11> only (|11> = c_a^dag c_b^dag |00>, with no sign on
    neighboring modes) and leaves |01>, |10> as they are."""
    G = np.eye(4, dtype=complex)
    G[0, 0], G[3, 0] = c, -np.conj(s)
    G[0, 3], G[3, 3] = s, c
    return G
