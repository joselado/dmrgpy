"""Fermionic Gaussian states (Slater determinants and BdG ground states) as
states of a chain.

A Slater determinant is fixed by its correlation matrix
Lambda[i,j] = <c_i^dagger c_j>, and a state with pairing by the
generalized one, G[a,b] = <gamma_a^dagger gamma_b> with the interleaved
gamma = (c_0, c_0^dagger, c_1, c_1^dagger, ...). On the MPS backends it is built by the
Fishman-White circuit (Fishman and White, PRB 92, 075132,
arXiv:1504.07701; the algorithm follows ITensorGaussianMPS.jl, Apache-2.0,
commit 6761dd6): nearest-neighbor Givens rotations found from Lambda
(gaussianmpstk/circuit.py), and with pairing also Bogoliubov rotations
of c_k and c_{k+1}^dagger, act on a product state of occupations, as
two-site gates on the site arrays (gaussianmpstk/gates.py), which go in
and out of the session through the same bridge as the Hadamard product
(mpsalgebratk/hadamard.py), so one kernel serves itensor_version=3,
"python" and "julia_live". mode="ED" builds the determinant directly, as
the product of sum_i Phi[i,k] Cdag_i over the occupied orbitals on the
vacuum, with no circuit, which makes it an independent reference; with
pairing it is the ground state of the parent Hamiltonian
(1/2) gamma^dagger (1 - 2 G^T) gamma, whose quasiparticles are +-1.

The modes are the chain's spinless sites in order. For
Spinful_Fermionic_Chain mode 2i is spin up and mode 2i+1 spin down of
orbital i, an n x n hopping matrix is taken for both spins, and an n x n
pairing matrix is the singlet pairing sum_ij Delta_ij c_iup^dag c_jdn^dag.

With pairing the state is the ground state of
H = sum_ij h_ij c_i^dag c_j + (1/2) sum_ij (Delta_ij c_i^dag c_j^dag + h.c.),
Delta antisymmetric, in the grand-canonical ensemble: the chemical
potential is in the diagonal of h, and nf is not used."""

import warnings

import numpy as np

from .gaussianmpstk import circuit
from .gaussianmpstk.gates import SiteArrays

# Above this summed discarded weight of the gate SVDs the state is not the
# Gaussian state any more, and gaussian_state() warns naming maxm.
DISCARDED_WARNING = 1e-6

# Distance to 0 or 1 above which an eigenvalue of the correlation matrix
# is not an occupation of a pure Gaussian state.
PURITY_TOLERANCE = 1e-6

_EMPTY = np.array([1.0, 0.0], dtype=complex)  # index 0: empty, 1: occupied
_OCCUPIED = np.array([0.0, 1.0], dtype=complex)


def _is_spinful(chain):
    from .fermionchain import Spinful_Fermionic_Chain
    return isinstance(chain, Spinful_Fermionic_Chain)


def _scale(m):
    return max(float(np.max(np.abs(m))) if m.size else 0.0, 1e-300)


def _square(m, name):
    m = np.array(m, dtype=complex)
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError("get_gaussian_state: %s must be a square matrix"
                         % name)
    return m


def _mode_hopping(chain, h):
    """The hopping matrix on the chain's modes."""
    h = _square(h, "h")
    nmodes = chain.ns
    if _is_spinful(chain) and 2*h.shape[0] == nmodes:
        h = np.kron(h, np.eye(2))
    if h.shape[0] != nmodes:
        raise ValueError("get_gaussian_state: h is %dx%d but the chain has "
                         "%d fermionic modes" % (h.shape[0], h.shape[1],
                                                 nmodes))
    if np.max(np.abs(h-h.conj().T)) > 1e-12*_scale(h):
        raise ValueError("get_gaussian_state: h must be Hermitian")
    return h


def _mode_pairing(chain, delta):
    """The antisymmetric pairing matrix on the chain's modes; on a spinful
    chain an n x n matrix is the singlet pairing Delta_ij c_iup^dag
    c_jdn^dag, and must be symmetric."""
    d = _square(delta, "delta")
    nmodes = chain.ns
    if _is_spinful(chain) and 2*d.shape[0] == nmodes:
        if np.max(np.abs(d-d.T)) > 1e-12*_scale(d):
            raise ValueError("get_gaussian_state: an n x n delta on a "
                             "spinful chain is a singlet pairing and must "
                             "be symmetric")
        D = np.zeros((nmodes, nmodes), dtype=complex)
        D[0::2, 1::2] = d
        D[1::2, 0::2] = -d.T
        return D
    if d.shape[0] != nmodes:
        raise ValueError("get_gaussian_state: delta is %dx%d but the chain "
                         "has %d fermionic modes" % (d.shape[0], d.shape[1],
                                                     nmodes))
    if np.max(np.abs(d+d.T)) > 1e-12*_scale(d):
        raise ValueError("get_gaussian_state: delta must be antisymmetric, "
                         "Delta_ij = -Delta_ji")
    return d


def bdg_correlation_matrix(chain, h, delta):
    """The generalized correlation matrix G[a,b] = <gamma_a^dag gamma_b>,
    gamma = (c_0, c_0^dag, c_1, ...), of the ground state of
    H = sum h_ij c_i^dag c_j + (1/2) sum (Delta_ij c_i^dag c_j^dag + h.c.):
    the quasiparticles of negative energy of the Bogoliubov-de Gennes
    matrix are filled."""
    h = _mode_hopping(chain, h)
    D = _mode_pairing(chain, delta)
    N = h.shape[0]
    M = np.zeros((2*N, 2*N), dtype=complex)
    M[0::2, 0::2] = h
    M[1::2, 1::2] = -h.T
    M[0::2, 1::2] = D
    M[1::2, 0::2] = -D.conj()
    E, W = np.linalg.eigh(M)
    tol = 1e-10*_scale(E)
    zero = np.abs(E) <= tol
    nzero = int(np.sum(zero))
    if nzero == 0:
        Wo = W[:, :N]
    else:
        warnings.warn("get_gaussian_state: the BdG spectrum has %d zero "
                      "mode(s) (|E| <= %.3e); the state is one of the "
                      "degenerate ground states" % (nzero//2, tol),
                      stacklevel=4)
        Wo = np.hstack([W[:, E < -tol], _zero_mode_annihilators(W[:, zero])])
    if Wo.shape[1] != N:
        raise RuntimeError("get_gaussian_state: %d filled quasiparticles "
                           "for %d modes" % (Wo.shape[1], N))
    G = Wo.conj() @ Wo.T
    _check_purity(G, "the BdG ground state")
    return G


def _zero_mode_annihilators(Z):
    """Half of the zero-energy subspace Z (columns), chosen orthogonal to
    its own particle-hole partners, so that it can be filled.

    eigh returns an arbitrary basis of Z, for instance a single Majorana
    mode, which is not a fermion that can be filled. The way we do this is
    by building a basis of Z invariant under the particle-hole conjugation
    P v = X conj(v) (Majorana-like vectors m_k, with real overlaps) and
    pairing them, s = (m_1 + i m_2)/sqrt(2), so that P s is orthogonal to
    every s."""
    two_m = Z.shape[1]
    if two_m % 2:
        raise ValueError("get_gaussian_state: an odd number (%d) of zero "
                         "modes of the BdG matrix; add a small term that "
                         "splits them" % two_m)
    P = np.empty_like(Z)  # P Z = X conj(Z), X swaps rows 2i and 2i+1
    P[0::2], P[1::2] = Z[1::2].conj(), Z[0::2].conj()
    C = np.hstack([Z+P, 1j*(Z-P)])  # P-invariant vectors
    K = np.real(C.conj().T @ C)
    lam, V = np.linalg.eigh(K)
    keep = lam > 1e-8*np.max(lam)
    Mj = (C @ V[:, keep])/np.sqrt(lam[keep])
    if Mj.shape[1] != two_m:
        raise RuntimeError("get_gaussian_state: the zero modes of the BdG "
                           "matrix could not be paired")
    return (Mj[:, 0::2]+1j*Mj[:, 1::2])/np.sqrt(2)


def _check_purity(m, what):
    n = np.linalg.eigvalsh(m)
    dist = np.max(np.minimum(np.abs(n), np.abs(1.0-n))) if n.size else 0.0
    if dist > PURITY_TOLERANCE:
        raise ValueError("get_gaussian_state: %s has a correlation matrix "
                         "with an eigenvalue %.3e away from 0 and 1, so it "
                         "is not a pure Gaussian state" % (what, dist))


def correlation_matrix(chain, h, nf):
    """The correlation matrix <c_i^dagger c_j> of the Slater determinant
    that fills the nf lowest orbitals of the hopping matrix h, where
    H = sum_ij h[i,j] c_i^dagger c_j. On a spinful chain an n x n matrix
    is taken for both spins (h kron 1)."""
    h = _mode_hopping(chain, h)
    nmodes = chain.ns
    if nf is None or int(nf) != nf or not 0 <= nf <= nmodes:
        raise ValueError("get_gaussian_state: nf must be an integer from 0 "
                         "to %d (got %r)" % (nmodes, nf))
    nf = int(nf)
    e, u = np.linalg.eigh(h)
    if 0 < nf < nmodes:
        gap = e[nf]-e[nf-1]
        if gap <= 1e-10*max(np.max(np.abs(e)), 1e-300):
            warnings.warn("get_gaussian_state: the level nf=%d is degenerate "
                          "with the next one (gap %.3e); the state is one "
                          "of the degenerate Slater determinants" % (nf, gap),
                          stacklevel=4)
    phi = u[:, :nf]
    return phi.conj() @ phi.T


def _check_correlation_matrix(chain, cm):
    """("normal", Lambda) or ("bdg", G) for a given correlation matrix,
    n x n or 2n x 2n; a generalized one with no anomalous part is taken as
    the normal one."""
    cm = np.array(cm, dtype=complex)
    nmodes = chain.ns
    if cm.shape not in ((nmodes, nmodes), (2*nmodes, 2*nmodes)):
        raise ValueError("get_gaussian_state: cm must be %dx%d, one row per "
                         "fermionic mode, or %dx%d with the interleaved "
                         "(c_i, c_i^dag) (got %s)" % (nmodes, nmodes,
                         2*nmodes, 2*nmodes, cm.shape))
    if np.max(np.abs(cm-cm.conj().T)) > 1e-12:
        raise ValueError("get_gaussian_state: cm must be Hermitian")
    _check_purity(cm, "cm")
    if cm.shape[0] == nmodes:
        return "normal", cm
    X = np.zeros((2*nmodes, 2*nmodes))
    X[np.arange(0, 2*nmodes, 2), np.arange(1, 2*nmodes, 2)] = 1.0
    X = X+X.T
    if np.max(np.abs(X @ cm @ X-(np.eye(2*nmodes)-cm.conj()))) > 1e-10:
        raise ValueError("get_gaussian_state: a 2n x 2n cm must have the "
                         "particle-hole symmetry of <gamma_a^dag gamma_b>, "
                         "gamma = (c_0, c_0^dag, c_1, ...)")
    if np.max(np.abs(cm[0::2, 1::2])) <= 1e-12:
        return "normal", cm[0::2, 0::2]
    return "bdg", cm


def _check_chain(chain, gutzwiller):
    from .fermionchain import Fermionic_Chain, Majorana_Chain
    if not isinstance(chain, Fermionic_Chain) or isinstance(chain,
                                                            Majorana_Chain):
        raise TypeError("get_gaussian_state needs a Fermionic_Chain or a "
                        "Spinful_Fermionic_Chain (spinless sites); the "
                        "native spinful chain is not supported")
    if getattr(chain, "conserved_sector", None):
        raise ValueError("get_gaussian_state: conserved-sector mode (%s) is "
                         "not supported; build the state on dense sites, "
                         "with the sector cleared by set_conserved_sector()"
                         % (chain.conserved_sector,))
    if gutzwiller and not _is_spinful(chain):
        raise ValueError("get_gaussian_state: gutzwiller=True needs a "
                         "Spinful_Fermionic_Chain")


def _bridge(chain):
    from .mpsalgebratk import hadamard
    if chain.itensor_version == "julia_live":
        return hadamard._JuliaStates(chain)
    session = getattr(chain, "_session", None)
    if session is None or not hasattr(session, "mps_from_site_arrays"):
        raise NotImplementedError("get_gaussian_state needs "
                "itensor_version=3, 'python' or 'julia_live' (got "
                "itensor_version=%r); use mode='ED' for a small chain"
                % (chain.itensor_version,))
    return hadamard._SessionStates(chain, session)


def _gaussian_state_dmrg(chain, kind, cm, eigval_cutoff, maxblocksize,
                         gutzwiller, maxm, cutoff):
    if kind == "normal":
        ns, rotations, err = circuit.correlation_matrix_to_circuit(
                cm, eigval_cutoff=eigval_cutoff, maxblocksize=maxblocksize)
        gates = [(circuit.gate_matrix(c, s), n) for (n, c, s) in rotations]
    else:
        ns, rotations, err = circuit.bdg_correlation_matrix_to_circuit(
                cm, eigval_cutoff=eigval_cutoff, maxblocksize=maxblocksize)
        gates = [(circuit.gate_matrix(c, s) if g == "hop" else
                  circuit.pair_gate_matrix(c, s), n)
                 for (g, n, c, s) in rotations]
    occ = np.rint(ns).astype(int)
    st = SiteArrays.product([_OCCUPIED if o == 1 else _EMPTY for o in occ])
    for (G, n) in reversed(gates):
        st.apply(G, n, maxm, cutoff)
    try:
        if gutzwiller:
            P = np.diag([1.0, 1.0, 1.0, 0.0]).astype(complex)  # 1-Nup Ndn
            for i in range(chain.ns//2):
                st.apply(P, 2*i, maxm, cutoff)
        st.normalize()
    except ValueError:
        raise ValueError("get_gaussian_state: the Gutzwiller projection "
                         "annihilates the state (no configuration without "
                         "double occupancy)") from None
    if st.discarded > DISCARDED_WARNING:
        warnings.warn("get_gaussian_state: the gates discarded a weight "
                      "%.3e at maxm=%d; raise maxm for the exact Gaussian "
                      "state" % (st.discarded, maxm), stacklevel=4)
    return _bridge(chain).to_mps(st.A)


def _parent_ground_state(chain, ed, G):
    """The ground state of (1/2) gamma^dag (1 - 2 G^T) gamma on ED, whose
    quasiparticle energies are -1 for the filled modes and +1 otherwise,
    so that it is unique, with a gap of 1."""
    from .edtk.edchain import State
    N = chain.ns
    gam = []
    gamd = []
    for i in range(N):
        gam += [chain.C[i], chain.Cdag[i]]
        gamd += [chain.Cdag[i], chain.C[i]]
    Mp = np.eye(2*N)-2.0*G.T
    Hp = 0
    for a in range(2*N):
        for b in range(2*N):
            if abs(Mp[a, b]) > 1e-14:
                Hp = Hp + 0.5*Mp[a, b]*gamd[a]*gam[b]
    M = ed.MO2matrix(Hp)
    if M.shape[0] <= 4096:
        M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
        v = np.linalg.eigh(M)[1][:, 0]
    else:
        from scipy.sparse.linalg import eigsh
        v = eigsh(M, k=1, which="SA")[1][:, 0]
    return State(np.asarray(v, dtype=complex), ed)


def _gaussian_state_ed(chain, kind, cm, gutzwiller):
    from .edtk.edchain import State
    if chain.hamiltonian is None:
        raise ValueError("get_gaussian_state: mode='ED' needs "
                         "set_hamiltonian() first, since an ED state lives "
                         "on the ED object of the Hamiltonian")
    ed = chain.get_ED_obj()
    if kind == "bdg":
        return _finish_ed(chain, ed, _parent_ground_state(chain, ed, cm),
                          gutzwiller)
    diag = np.real(ed.MO2matrix(sum(chain.N)).diagonal())
    vac = int(np.argmin(np.abs(diag)))
    if abs(diag[vac]) > 1e-9:
        raise RuntimeError("get_gaussian_state: no vacuum in the ED basis")
    v = np.zeros(len(diag), dtype=complex)
    v[vac] = 1.0
    wf = State(v, ed)
    n, u = np.linalg.eigh(cm)
    phi = u[:, n > 0.5].conj()  # <c_i^dag c_j> = sum_k phi*_ik phi_jk
    for k in range(phi.shape[1]):
        orbital = sum(phi[i, k]*chain.Cdag[i] for i in range(chain.ns)
                      if phi[i, k] != 0)
        wf = ed.applyoperator(orbital, wf)
    return _finish_ed(chain, ed, wf, gutzwiller)


def _finish_ed(chain, ed, wf, gutzwiller):
    from .edtk.edchain import State
    if gutzwiller:
        for i in range(chain.ns//2):
            wf = ed.applyoperator(1-chain.N[2*i]*chain.N[2*i+1], wf)
    nrm = np.sqrt(abs(wf.dot(wf)))
    if nrm == 0.0 or not np.isfinite(nrm):
        raise ValueError("get_gaussian_state: the Gutzwiller projection "
                         "annihilates the state (no configuration without "
                         "double occupancy)")
    return State(wf.v/nrm, ed)


def gaussian_state(chain, h=None, nf=None, cm=None, delta=None,
                   eigval_cutoff=1e-8, maxblocksize=None, gutzwiller=False,
                   maxm=None, cutoff=None, mode=None):
    """The normalized Gaussian state given by the hopping matrix h and the
    number of fermions nf (a Slater determinant), by h and the pairing
    matrix delta (the BdG ground state, no nf), or by its correlation
    matrix cm, n x n (<c_i^dag c_j>) or 2n x 2n (<gamma_a^dag gamma_b>,
    gamma = (c_0, c_0^dag, c_1, ...)), as a state of the chain.

    eigval_cutoff: how close to 0 or 1 an occupation of a block must be
        for its mode to be decoupled; it sets the accuracy of the circuit.
    maxblocksize: the largest block of the circuit (default: no limit).
        A smaller one gives fewer gates and a smaller bond dimension, and
        a less accurate state.
    gutzwiller: on a Spinful_Fermionic_Chain, project out every doubly
        occupied orbital, prod_i (1 - Nup_i Ndn_i), and normalize.
    maxm, cutoff: the truncation of every gate (default: the chain's).
    mode: "DMRG" (itensor_version=3, "python", "julia_live") or "ED"."""
    from .mode import resolve_mode
    mode = resolve_mode(chain, mode="DMRG" if mode is None else mode)
    _check_chain(chain, gutzwiller)
    if (h is None) == (cm is None):
        raise ValueError("get_gaussian_state: give either h (with nf or "
                         "delta) or cm")
    if cm is not None:
        if nf is not None or delta is not None:
            raise ValueError("get_gaussian_state: nf and delta go with h, "
                             "not cm")
        kind, cm = _check_correlation_matrix(chain, cm)
    elif delta is not None:
        if nf is not None:
            raise ValueError("get_gaussian_state: with delta the state is "
                             "grand canonical; put the chemical potential in "
                             "the diagonal of h, and drop nf")
        kind, cm = _check_correlation_matrix(
                chain, bdg_correlation_matrix(chain, h, delta))
    else:
        kind, cm = "normal", correlation_matrix(chain, h, nf)
    if mode == "ED":
        return _gaussian_state_ed(chain, kind, cm, gutzwiller)
    if maxm is None:
        maxm = chain.maxm
    if cutoff is None:
        cutoff = chain.cutoff
    return _gaussian_state_dmrg(chain, kind, cm, eigval_cutoff, maxblocksize,
                                gutzwiller, maxm, cutoff)
