"""The shells of an NRG read off the bond basis of a DMRG state, and the
hand-over from a DMRG head to an NRG tail (milestones N2 and N3 of
docs/nrg_mps_kondo_plan.md).

A matrix product state of a Wilson chain is an NRG in disguise
(Verstraete et al., cond-mat/0504305): its left isometries A_0 ... A_p map
the Fock space of sites 0..p onto chi_p states, the left basis of bond p,
and that basis plays the part of the kept states of the shell p. Saberi,
Weichselbaum and von Delft (arXiv:0804.0193) showed that a DMRG ground
state of a Wilson chain holds the NRG kept states of every shell when its
bond dimension is large enough. Here the chain is ysr.spin_chain_on_metal
(ysr.kondo_chain for ns = 1 has the same layout): spins on sites
0..ns-1, the Wilson chain on ns..ns+n-1, so Wilson site k is site
p = ns + k of the MPS, and the shell k of nrgspectral (the spins and the
Wilson sites 0..k) is the left block of bond p.

What is computed. Any operator O that acts on sites 0..p only is written
in the left basis by contracting the MPO the backend exports
(task.mpo_arrays on "python", Chain.mpo_site_arrays on v3, the arrays
TASK uses) with the left isometries of task.mixed_gauge, after a
compression that brings the MPO to bond dimension 1 right of p; each MPO
tensor right of p must then be a multiple of the identity, which is
checked. So the site conventions and the Jordan-Wigner strings are those of
the backend, with no copy of them here. The block Hamiltonian H_p is the
sum of the terms of H with every factor on sites <= p.

Labels. The tail of N3 needs the left basis graded by the charge and S^z
of the block. N_p and S^z_p are written in the left basis and diagonalized
together; their eigenvalues must be integers (half-integers for S^z) to
GRADE_TOL, or graded_basis raises. The left basis is graded when the state
has a definite total N and S^z and the cut at maxm does not split a
degenerate set of Schmidt values across labels. A plain solve of a doublet
or a triplet is an arbitrary member of the multiplet, so head_state solves
in the sector (Nf, 2Sz) and promotes the state to dense indices.

What a shell from DMRG is and is not. The block Hamiltonian in the left
basis is P H_p P, P the projector on the chi_p states. Its eigenvalues are
upper bounds of the lowest exact ones of H_p (Cauchy interlacing), and the
sum of its k lowest equals that of H_p exactly when the basis holds the k
lowest eigenstates. flow_diagram reports these levels, in units of the
scale omega_p of the shell, so a comparison with the exact NRG of
nrgspectral at the same shell measures the basis error. The left basis
holds at most the dimension of the right part of the chain, 4^m for m
sites after p.

The hand-over (handover_shell). At bond p the block Hamiltonian, the
labels, the two matrices of f_{p s} the next hopping needs and the tip
vertices X^dag_sigma of ysr.tunneling_operators become an
nrgspectral.StartShell; the tail is then nrgspectral.spectral_poles on the
Wilson sites after p. The anticommutator {X, X^dag} is projected as one
operator, so the sum rule of the result is <{X, X^dag}> on the state the
shell is cut from, and the part of X^dag|GS> outside the left basis is
missing from the poles: missing_weight = sum_rule - total.
"""

import copy

import numpy as np

from . import nrgspectral

GRADE_TOL = 1e-8 # the eigenvalues of N_p and 2 S^z_p must be integers to this
IDENTITY_TOL = 1e-10 # an MPO tensor right of the block is c*identity to this
BLOCK_TOL = 1e-4 # the block Hamiltonian couples two labels by at most this
# (first order in the error of the state, whose energies it moves at second
# order; the size dropped is kept as info["off_block"] of the shell)


def _layout(fc):
    """ns (the spin sites) and the sites of the Wilson chain of fc."""
    imp = list(getattr(fc, "impurity", []))
    wil = getattr(fc, "wilson", None)
    if wil is None:
        ch = getattr(fc, "channels", None)
        if not ch or len(ch) != 1:
            raise ValueError("nrgmps needs a chain of ysr.spin_chain_on_metal "
                             "or ysr.kondo_chain (one Wilson chain after the "
                             "spins)")
        wil = ch[0]
    wil = list(wil)
    if wil != list(range(len(imp), len(imp) + len(wil))):
        raise ValueError("nrgmps needs the Wilson chain right after the spins")
    if len(imp) != 1 and getattr(fc, "spin", None) is not None and \
            not hasattr(fc, "tip_spin"):
        raise ValueError("nrgmps needs spins 1/2 on single sites "
                         "(ysr.spin_chain_on_metal or a spin-1/2 kondo_chain)")
    return len(imp), wil


def head_state(fc, Nf=None, Sz=None):
    """A copy of fc solved by DMRG in the sector (Nf electrons, 2S^z = Sz)
    and promoted to dense indices, so that its left bases are graded.
    The default is half filling, Nf = number of sites (one electron on
    every spin site and on average one per Wilson site), and Sz = Nf mod 2,
    the smallest |2S^z|, which holds a member of the ground multiplet at
    B = 0. Sz="auto" solves every 2S^z from -(ns+1) to ns+1 of the parity
    of Nf and keeps the lowest (the default when the spins are in a
    field, fc.B_imp != 0). On itensor_version=3 the solve starts from the
    warm solve of ysr._seed_impurity_charge, since no term moves charge
    onto a spin site; every solved state is checked to hold one electron
    on every spin site. Returns (chain, E0)."""
    from . import ysr
    if fc.itensor_version not in ("python", 3):
        raise NotImplementedError("nrgmps: itensor_version=\"python\" or 3, "
                                  "got %r" % (fc.itensor_version,))
    ns, wil = _layout(fc)
    if Nf is None:
        Nf = ns + len(wil)
    if Sz is None:
        Sz = "auto" if getattr(fc, "B_imp", 0.0) != 0.0 else Nf % 2
    if isinstance(Sz, str):
        if Sz != "auto":
            raise ValueError("Sz must be an integer (2S^z) or \"auto\", got %r"
                             % (Sz,))
        best = None
        for s2 in range(-(ns + 1), ns + 2):
            if (Nf - s2) % 2:
                continue
            c, e = head_state(fc, Nf=Nf, Sz=s2)
            if best is None or e < best[1] - 1e-12*max(1.0, abs(e)):
                best = (c, e)
        return best
    if (Nf - Sz) % 2:
        raise ValueError("Nf and 2Sz must have the same parity, got %r and %r"
                         % (Nf, Sz))
    c = copy.deepcopy(fc)
    c.set_conserved_sector(Nf=int(Nf), Sz=int(Sz))
    if c.itensor_version == 3:
        ysr._seed_impurity_charge(c, ns)
    e0 = float(np.real(c.gs_energy()))
    ysr._check_impurity_charge(c, ns, "DMRG")
    c.promote_to_dense()
    c.head_sector = dict(Nf=int(Nf), Sz=int(Sz))
    return c, e0


class Export(object):
    """The ground state of a chain as left isometries (As[j], shape
    (chi_{j-1}, d, chi_j)) and a builder of MPO arrays for any
    MultiOperator (mpo(op)), on "python" and on 3. pad_bonds goes to
    task.mixed_gauge ("full" completes every left basis to the whole
    space of its block, for tests)."""

    def __init__(self, fc, pad_bonds=None):
        from . import taskdmrg
        from .pyitensor import task, backend
        if fc.itensor_version not in ("python", 3):
            raise NotImplementedError("nrgmps: itensor_version=\"python\" or "
                                      "3, got %r" % (fc.itensor_version,))
        taskdmrg.ground_state(fc)
        s = fc._session
        if fc.itensor_version == "python":
            from .pyitensor.chain import _mps_arrays_lpr
            Ms = _mps_arrays_lpr(s.wf0)
            self._mpo = lambda terms: task.mpo_arrays(s._mpo(terms))
        else:
            Ms = s.mps_site_arrays(s.gs_wavefunction())
            self._mpo = lambda terms: s.mpo_site_arrays(
                s.build_operator(terms))
        As, _, Cs = task.mixed_gauge([np.asarray(M) for M in Ms])
        As = [np.asarray(backend.to_host(A)) for A in As]
        Cs = [np.asarray(backend.to_host(C)) for C in Cs]
        As, self.real_residual = real_gauge(As, Cs)
        if pad_bonds is not None: # padded after the real gauge, so in reals
            As, _, _ = task.mixed_gauge(As, pad_bonds=pad_bonds)
            As = [np.ascontiguousarray(np.real(backend.to_host(A)))
                  for A in As]
        self.As = As
        self.L = len(self.As)

    def mpo(self, op):
        from .pyitensor import task, backend
        W = task.compress_mpo([np.asarray(w) for w in self._mpo(op.to_terms())])
        return [np.asarray(backend.to_host(w)) for w in W]

    def chi(self, p):
        return self.As[p].shape[2]

    def flipped(self, u):
        """The same Export with the state R|GS>, R the product of the local
        orthogonal matrix u on every site (spin_flip_matrix): the left
        isometries are u applied to the physical index of each tensor."""
        out = copy.copy(self)
        out.As = [np.ascontiguousarray(np.einsum("st,atb->asb", u, A))
                  for A in self.As]
        return out

    def overlap(self, other, p):
        """<a|b> between the left basis of bond p of self (bra) and of
        other (ket), a chi x chi' matrix."""
        Lm = np.ones((1, 1))
        for j in range(p + 1): # two pair contractions, chi^3 d (a single
            # three-operand einsum loops over all five indices, chi^4 d)
            X = np.tensordot(Lm, self.As[j].conj(), axes=(0, 0)) # b, s, c
            Lm = np.tensordot(X, other.As[j], axes=([0, 1], [0, 1])) # c, d
        return Lm

    def block_matrix(self, op, p, other=None):
        """<a|O|b> in the left basis of bond p (between sites p and p+1),
        O a MultiOperator with every factor on sites <= p. other: an Export
        whose left basis is the ket (default self)."""
        ket = self.As if other is None else other.As
        W = self.mpo(op)
        r = np.ones(1)
        for j in range(self.L - 1, p, -1): # right of the block: c * identity
            w = W[j]
            d = w.shape[1]
            T = np.einsum("asst->at", w)/d
            dev = w - T[:, None, None, :]*np.eye(d)[None, :, :, None]
            if np.max(np.abs(dev)) > IDENTITY_TOL*max(1.0, np.max(np.abs(w))):
                raise ValueError("block_matrix: the operator acts on site %d, "
                                 "right of the block 0..%d" % (j, p))
            r = T @ r
        Lm = np.ones((1, 1, 1))
        for j in range(p + 1):
            A = self.As[j]
            X = np.tensordot(Lm, A.conj(), axes=(0, 0)) # w, a', s, b
            X = np.tensordot(X, W[j], axes=([0, 2], [0, 1])) # a', b, t, w'
            Lm = np.tensordot(X, ket[j], axes=([0, 2], [0, 1])) # b, w', b'
        M = np.tensordot(Lm, r, axes=(1, 0))
        if np.iscomplexobj(M):
            if np.max(np.abs(M.imag)) > 1e-10*max(1.0, np.max(np.abs(M))):
                raise ValueError("block_matrix: a complex matrix in the left "
                                 "basis; nrgmps needs a real state")
            M = M.real
        return np.ascontiguousarray(M)


REAL_TOL = 1e-4 # the weighted left basis is closed under conjugation to this


def real_gauge(As, Cs=None, tol=REAL_TOL):
    """The left isometries As of a state that is real up to a global
    phase, in a real gauge: at every bond a real orthonormal basis close
    to the span of the left basis, and the matrix that relates the two
    moved into the next tensor; the last tensor loses the global phase.
    The "python" backend returns such a state in a complex gauge.

    The basis is the leading left singular vectors of [Re M, Im M], M the
    left basis of the bond written in the real gauge of the bonds before
    and weighted by the bond matrix of the state (Cs, the site-centre
    tensors of task.mixed_gauge, None for no weight). The weight matters:
    a direction with a Schmidt value of 1e-6 is fixed by the state only to
    about 1e-15/1e-6, and without it such directions fail the test below
    (residual 3.7e-4 on a 6-site chain whose state is real to 1e-15).

    Returns (real As, residual): the residual is the largest relative
    singular value beyond the bond dimension (zero for a weighted span
    closed under conjugation) and, last, the imaginary part left on the
    last tensor. It measures how far the state is from real, which for a
    solved state is the accuracy of the solve: the exact ground state of a
    real H in its sector is real, and a 7-site "python" solve at 20
    sweeps was real to 2.5e-7 with a residual of 1.1e-6. The real basis
    is the closest one to the weighted span, so it moves the left basis by
    about as much. Raises above tol (1e-4), which no converged state
    reaches."""
    if all(not np.iscomplexobj(A) or not np.any(A.imag) for A in As):
        return [np.ascontiguousarray(np.real(A)) for A in As], 0.0
    out, res = [], 0.0
    G = np.ones((1, 1), dtype=complex)
    L = len(As)
    for j, A in enumerate(As):
        X = np.tensordot(G, A, axes=(1, 0))
        a, d, b = X.shape
        M = X.reshape(a*d, b)
        if j == L - 1:
            k = np.argmax(np.abs(M))
            ph = M.flat[k]/abs(M.flat[k])
            M = M/ph
            res = max(res, float(np.max(np.abs(M.imag))/np.max(np.abs(M))))
            out.append(np.ascontiguousarray(M.real.reshape(a, d, b)))
            break
        Mw = M
        if Cs is not None: # the bond matrix: C_j = A_j Lambda_j
            A2 = A.reshape(-1, b)
            Lam = A2.conj().T @ np.asarray(Cs[j]).reshape(-1, Cs[j].shape[2])
            Mw = M @ Lam
        U, S, _ = np.linalg.svd(np.concatenate([Mw.real, Mw.imag], axis=1),
                                full_matrices=False)
        if len(S) > b and S[0] > 0:
            res = max(res, float(S[b]/S[0]))
        Q = U[:, :b]
        G = Q.T @ M
        out.append(np.ascontiguousarray(Q.reshape(a, d, b)))
    if res > tol:
        raise RuntimeError("real_gauge: the left basis is not closed under "
                           "complex conjugation (residual %.2e): the state "
                           "is not real up to a phase" % res)
    return out, res


def block_hamiltonian(fc, p):
    """The terms of fc.hamiltonian with every factor on sites <= p, as a
    MultiOperator."""
    h = fc.hamiltonian.copy()
    h.op = [term for term in fc.hamiltonian.op
            if all(f[1] <= p for f in term[1:])]
    h.i = len(h.op) - 1
    return h


def block_charges(fc, p):
    """N and S^z of the sites 0..p."""
    N = fc.Nup[0] + fc.Ndn[0]
    Sz = fc.Sz[0]
    for j in range(1, p + 1):
        N = N + fc.Nup[j] + fc.Ndn[j]
        Sz = Sz + fc.Sz[j]
    return N, Sz


def graded_basis(Nm, Szm, tol=GRADE_TOL):
    """Diagonalize the block charge N and S^z together. Returns (V,
    labels): V orthogonal, its columns the graded basis, labels the
    (Q, 2Sz) of each column. Raises when an eigenvalue is not an integer
    (half-integer for S^z) to tol, meaning the left basis mixes labels."""
    q, V = np.linalg.eigh(0.5*(Nm + Nm.T))
    qi = np.round(q)
    if np.max(np.abs(q - qi), initial=0.0) > tol:
        raise RuntimeError("graded_basis: the charge of the left basis is "
                           "not an integer to %.1e (off by %.2e): the basis "
                           "mixes charges; solve in a sector (head_state) or "
                           "raise maxm" % (tol, np.max(np.abs(q - qi))))
    cols, labels = [], []
    for Q in np.unique(qi):
        ix = np.where(qi == Q)[0]
        Vq = V[:, ix]
        s = Vq.T @ Szm @ Vq
        s2, U = np.linalg.eigh(s + s.T) # 2 S^z
        si = np.round(s2)
        if np.max(np.abs(s2 - si), initial=0.0) > tol:
            raise RuntimeError("graded_basis: 2S^z of the left basis is not "
                               "an integer to %.1e (off by %.2e) at Q = %d: "
                               "the basis mixes S^z; solve in a sector "
                               "(head_state) or raise maxm"
                               % (tol, np.max(np.abs(s2 - si)), Q))
        cols.append(Vq @ U)
        labels += [(int(Q), int(x)) for x in si]
    return np.concatenate(cols, axis=1), labels


def _scale(eps, t, k, Lambda):
    """The scale of shell k, as nrgspectral.chain_poles takes it."""
    if k == 0:
        return np.sqrt(Lambda)*abs(t[0]) if len(t) and t[0] != 0 else 1.0
    return abs(t[k - 1])


def wilson_parameters(fc):
    """(eps, t) of the Wilson chain of fc (fc.wilson_chain, stored by
    ysr.spin_chain_on_metal and ysr.kondo_chain)."""
    wc = getattr(fc, "wilson_chain", None)
    if wc is None:
        raise ValueError("nrgmps needs a chain of ysr.spin_chain_on_metal or "
                         "ysr.kondo_chain, which store fc.wilson_chain")
    return np.array(wc[0]), np.array(wc[1])


def flow_diagram(fc, bonds=None, nlev=40, Lambda=None, pad_bonds=None,
                 export=None):
    """The shells of the DMRG state of fc (head_state first, for a graded
    basis): for every Wilson site k in bonds (default all but the last),
    the lowest nlev eigenvalues of the block Hamiltonian in the left basis
    of bond p = ns + k, relative to its lowest one and in units of the
    scale omega_k of nrgspectral.chain_poles, with their labels (Q - Q0,
    2Sz), Q0 = ns + k + 1 being half filling of the block. Lambda (default
    fc's, from the ratio of the hoppings) sets omega_0. Returns a dict
    k -> dict(levels, labels, chi, omega, E0) with E0 the absolute lowest
    eigenvalue."""
    ns, wil = _layout(fc)
    eps, t = wilson_parameters(fc)
    if Lambda is None:
        Lambda = (abs(t[0]/t[2]) if len(t) > 2 else 4.0)
    ex = export if export is not None else Export(fc, pad_bonds=pad_bonds)
    if bonds is None:
        bonds = range(len(wil) - 1)
    out = {}
    for k in bonds:
        p = ns + k
        H = ex.block_matrix(block_hamiltonian(fc, p), p)
        Nop, Szop = block_charges(fc, p)
        V, labels = graded_basis(ex.block_matrix(Nop, p),
                                 ex.block_matrix(Szop, p))
        Hg = V.T @ H @ V
        E, lab = [], []
        for L in sorted(set(labels)):
            ix = [i for i, l in enumerate(labels) if l == L]
            e = np.linalg.eigvalsh(Hg[np.ix_(ix, ix)])
            E += list(e)
            lab += [(L[0] - (ns + k + 1), L[1])]*len(e)
        E = np.asarray(E)
        order = np.argsort(E, kind="stable")
        om = _scale(eps, t, k, Lambda)
        E0 = float(E[order[0]])
        out[k] = dict(levels=(E[order[:nlev]] - E0)/om,
                      labels=[lab[i] for i in order[:nlev]],
                      chi=ex.chi(p), omega=om, E0=E0)
    return out


JOIN_TOL = 1e-3 # a flipped direction joins the hand-over basis when at least this far outside it


def _joined_basis(fc, ex, p, Nop, Szop):
    """The hand-over basis of handover_shell(symmetric=True): the graded
    left basis of bond p of the state, and, label by label (Q, 2Sz), the
    graded left basis of its flipped partner R|GS> orthogonalized against
    it. A flipped direction is kept when its part outside the span is at
    least JOIN_TOL of its norm: the labels of each basis are integers only
    to the accuracy of the solve, and a nearly dependent direction would
    amplify that error by the inverse of its weight. The basis of (Q, -2Sz)
    is the exact image under R of that of (Q, 2Sz), and the one of 2Sz = 0
    is cut by the eigenvalues of an overlap matrix that commutes with R, so
    the joined basis is closed under R exactly, whatever was cut: a
    splitting of the partner spectra acts on a free moment as a field,
    which the flow makes relevant (with the cut alone, 60% of A off the
    exact route on 3 spins with the metal under the middle spin, without
    the flip average). Returns (bm, labels): bm(op) the matrix of a block
    operator in the joined basis, which is graded by construction."""
    exr = ex.flipped(spin_flip_matrix(fc, ex))
    O = ex.overlap(exr, p) # <a|R b>
    VA, la = graded_basis(ex.block_matrix(Nop, p), ex.block_matrix(Szop, p))
    # R|a> has the charge of |a> and the opposite 2Sz, so the graded basis
    # of the flipped partner is VA itself with the labels (Q, -2Sz)
    lr = [(q, -s) for (q, s) in la]
    n = VA.shape[0]
    def cols(labs, lab):
        return VA[:, [i for i, x in enumerate(labs) if x == lab]]
    blocks = {}
    for lab in sorted(set(la) | set(lr)):
        q, s = lab
        if s < 0:
            continue
        a, b = cols(la, lab), cols(lr, lab)
        if s > 0:
            # the state's own vectors, and the flipped ones orthogonalized
            # against them, kept when at least JOIN_TOL outside the span
            S = a.T @ O @ b
            w, U = np.linalg.eigh(np.eye(b.shape[1]) - S.T @ S)
            keep = w > JOIN_TOL**2
            Cc = U[:, keep]/np.sqrt(w[keep])
            top = np.hstack([a, -a @ S @ Cc])
            bot = np.hstack([np.zeros((n, a.shape[1])), b @ Cc])
            blocks[lab] = (top, bot)
            # the partner label is the image under R, exactly: R(A x + R A y)
            # = (-1)^Q A y + R A x, since R^2 = (-1)^N
            blocks[(q, -s)] = ((-1.0)**q*bot, top)
        else:
            # 2Sz = 0: R maps the label onto itself, and the overlap matrix of
            # [a, R a] commutes with R, so cutting it by its eigenvalues keeps
            # the span closed under R
            S = a.T @ O @ b
            g, Vg = np.linalg.eigh(np.block([[np.eye(a.shape[1]), S],
                                             [S.T, np.eye(b.shape[1])]]))
            keep = g > JOIN_TOL**2*np.max(g)
            C = Vg[:, keep]/np.sqrt(g[keep])
            blocks[lab] = (a @ C[:a.shape[1]], b @ C[a.shape[1]:])
    labels, top, bot = [], [], []
    for lab in sorted(blocks):
        t, b = blocks[lab]
        top.append(t)
        bot.append(b)
        labels += [lab]*t.shape[1]
    Cj = np.vstack([np.hstack(top), np.hstack(bot)])
    def bm(op):
        M = np.block([[ex.block_matrix(op, p), ex.block_matrix(op, p, exr)],
                      [exr.block_matrix(op, p, ex), exr.block_matrix(op, p)]])
        return Cj.T @ M @ Cj
    # the joined basis is orthonormal and graded, to the accuracy of the solve
    G = Cj.T @ np.block([[np.eye(n), O], [O.T, np.eye(n)]]) @ Cj
    Q = np.array([x[0] for x in labels], dtype=float)
    S2 = np.array([x[1] for x in labels], dtype=float)
    dev = max(np.max(np.abs(G - np.eye(len(labels)))),
              np.max(np.abs(bm(Nop) - np.diag(Q))),
              np.max(np.abs(2*bm(Szop) - np.diag(S2))))
    if dev > BLOCK_TOL:
        raise RuntimeError("handover_shell: the joined basis is orthonormal "
                           "and graded only to %.2e" % dev)
    return bm, labels


def spin_flip_matrix(fc, export):
    """The matrix u of R = exp(-i pi S^y_tot) on one site of fc, every site
    of ysr.spin_chain_on_metal being a native spinful site: R c_up R^dag =
    c_dn and R c_dn R^dag = -c_up. It is found on the backend's own basis
    and on-site signs, from the tensors of Cup[0] and Cdn[0] that export
    gives (site 0, where no Jordan-Wigner string enters): |0> and the
    doubly occupied state go to themselves and the two singly occupied
    ones to each other, with the signs that satisfy both relations. Real
    and orthogonal, and it keeps the parity, so the product of u on every
    site is R on the chain."""
    import itertools
    def local(op):
        W = [np.asarray(w) for w in export.mpo(op)]
        if any(w.shape[0] != 1 or w.shape[-1] != 1 for w in W):
            raise RuntimeError("spin_flip_matrix: the MPO of a one-site "
                               "operator is not of bond dimension 1")
        M = W[0][0, :, :, 0]
        for w in W[1:]:
            M = M*np.trace(w[0, :, :, 0])/w.shape[1]
        if np.max(np.abs(np.imag(M))) > 1e-12:
            raise RuntimeError("spin_flip_matrix: a complex on-site matrix")
        return np.real(M)
    Mu, Md = local(fc.Cup[0]), local(fc.Cdn[0])
    nu = np.round(np.diag(Mu.T @ Mu)).astype(int)
    nd = np.round(np.diag(Md.T @ Md)).astype(int)
    def state(a, b):
        ix = np.where((nu == a) & (nd == b))[0]
        if len(ix) != 1:
            raise RuntimeError("spin_flip_matrix: the site is not a spinful "
                               "fermion site of four states")
        return ix[0]
    e00, e10, e01, e11 = state(0, 0), state(1, 0), state(0, 1), state(1, 1)
    for s1, s2, s3 in itertools.product((1.0, -1.0), repeat=3):
        u = np.zeros((len(nu), len(nu)))
        u[e00, e00], u[e01, e10], u[e10, e01], u[e11, e11] = 1.0, s1, s2, s3
        if np.allclose(u @ Mu @ u.T, Md, atol=1e-12) and \
                np.allclose(u @ Md @ u.T, -Mu, atol=1e-12):
            return u
    raise RuntimeError("spin_flip_matrix: no sign choice gives R c_up R^dag = "
                       "c_dn and R c_dn R^dag = -c_up")


def handover_shell(fc, k0, U=0.0, imp=None, export=None, pad_bonds=None,
                   Lambda=None, symmetric=False):
    """The StartShell of nrgspectral at the Wilson site k0 of the DMRG
    state of fc: the block Hamiltonian, the labels, f_{p s} and the tip
    vertices X^dag_sigma (ysr.tunneling_operators with U and imp, the
    electron at the first Wilson site) in the left basis of bond
    p = ns + k0, graded and diagonalized. anti holds the projected
    anticommutator {X, X^dag} of each vertex, for the sum rule. Q in the
    labels is the charge of the whole block (spin sites included), whose
    parity is the Jordan-Wigner sign of the next site.

    symmetric=True (B = 0 only) hands over the span of the left bases of
    the state and of its partner R|GS> under the flip of all spins
    (R = exp(-i pi S^y_tot), spin_flip_matrix), orthonormalized through
    their overlap, so that the basis is closed under the flip also for a
    head with a moment (an odd number of electrons), whose own left basis
    is not."""
    from .ysr import tunneling_operators
    ns, wil = _layout(fc)
    if not 0 <= k0 < len(wil) - 1:
        raise ValueError("k0 must be a Wilson site before the last one, "
                         "0..%d, got %r" % (len(wil) - 2, k0))
    p = ns + k0
    eps, t = wilson_parameters(fc)
    if Lambda is None:
        Lambda = (abs(t[0]/t[2]) if len(t) > 2 else 4.0)
    ex = export if export is not None else Export(fc, pad_bonds=pad_bonds)
    Nop, Szop = block_charges(fc, p)
    if symmetric:
        if getattr(fc, "B_imp", 0.0) != 0.0:
            raise ValueError("handover_shell: symmetric=True needs B = 0, "
                             "where H commutes with the flip of all spins")
        bm, labels = _joined_basis(fc, ex, p, Nop, Szop)
        V = np.eye(len(labels))
    else:
        bm = lambda op: ex.block_matrix(op, p)
        V, labels = graded_basis(bm(Nop), bm(Szop))
    chi = len(labels)
    rot = lambda op: V.T @ bm(op) @ V
    H = rot(block_hamiltonian(fc, p))
    c = {1: rot(fc.Cup[p]), -1: rot(fc.Cdn[p])}
    where = {} if imp is None else {"imp": imp}
    Xs = tunneling_operators(fc, U=U, site=ns, **where)
    ops = {"up": rot(Xs[0]), "dn": rot(Xs[1])}
    anti = {}
    for name, X in zip(("up", "dn"), Xs):
        Xd = X.get_dagger()
        anti[name] = rot(X*Xd + Xd*X)
    shift = {"up": (1, 1), "dn": (1, -1)}
    om = _scale(eps, t, k0, Lambda)
    sh = nrgspectral._shell_from_dense(H, labels, c, ops, shift, om,
                                       anti=anti, block_tol=BLOCK_TOL,
                                       info=dict(k0=k0, p=p, chi=chi,
                                                 U=U, B=getattr(fc, "B_imp",
                                                                0.0)))
    return sh, ex


MIN_HANDOVER = 100 # a hand-over basis below this warns (README section 5)
N_MAX = 120 # the longest Wilson chain the tail is given, before its stop rule
OMEGA_MIN = 1e-30 # and the chain ends before its hoppings fall below this


def max_sites(Lambda, D=1.0):
    """The longest tail chain: N_MAX sites, or fewer, so that its scale
    D Lambda^(-n/2) stays above OMEGA_MIN. wilsonchain.wilson_chain loses
    its precision below a scale of about 1e-35 (at Lambda = 4 the
    120-site chain is broken from site 113 on, |eps|/t grows to 8e8)."""
    return int(min(N_MAX, 2.0*np.log(D/OMEGA_MIN)/np.log(Lambda)))


def tail_chain(fc, n_tail=None):
    """The Wilson chain the tail runs on: n_tail sites (default N_MAX) of
    the mesh of fc (fc.wilson_params), whose first sites must be the
    Wilson chain of fc itself. The default is max_sites(Lambda, D)."""
    from . import wilsonchain
    wp = getattr(fc, "wilson_params", None)
    eps_h, t_h = wilson_parameters(fc)
    if wp is None:
        raise ValueError("the head was built with wilson= by hand; pass the "
                         "tail chain as wilson_tail=(eps, t)")
    n = max_sites(wp["Lambda"], wp["D"]) if n_tail is None else int(n_tail)
    if n < len(eps_h):
        raise ValueError("n_tail=%d is shorter than the head (%d Wilson sites)"
                         % (n, len(eps_h)))
    eps, t = wilsonchain.wilson_chain(n, Lambda=wp["Lambda"],
                                      band=(-wp["D"], wp["D"]), z=wp["z"],
                                      mesh=wp["mesh"])
    return eps, t


def lost_wilson_bonds(fc):
    """The Wilson bonds k (between the Wilson sites k and k+1) whose
    channels the MPO of H on fc's backend has dropped: those whose link
    dimension in the MPO is below that of the first Wilson bond. The last
    bond is left out, since its link dimension is smaller at any |t| (the
    end of the MPO). On v3, ITensor's svdMPO cuts a squared singular value
    below an absolute 1e-13 at unit scale, so a hopping below about 3e-7 of
    the unit of H leaves the MPO and the deeper Wilson sites are nearly
    decoupled from the chain (measured on 1 spin at rho JK = 0.15: link
    dimension 7 down to |t| = 4.8e-7, 6 at 3.0e-7 to 3.5e-7, 3 from 1.5e-7
    on, at Lambda = 2.5, 4 and 6; README section 9). "python" truncates the
    return sweep of its builder at a relative 1e-14 and loses channels one
    bond earlier at Lambda = 4. The MPO is built again from fc.hamiltonian
    with the session's own builder. Returns a sorted list, empty when every
    Wilson bond is whole."""
    from .pyitensor import task
    ns, wil = _layout(fc)
    if len(wil) < 3:
        return []
    s = fc._session
    terms = fc.hamiltonian.to_terms()
    if fc.itensor_version == "python":
        W = task.mpo_arrays(s._mpo(terms))
    elif fc.itensor_version == 3:
        W = s.mpo_site_arrays(s.build_operator(terms))
    else:
        raise NotImplementedError("lost_wilson_bonds: itensor_version=\"python\" "
                                  "or 3, got %r" % (fc.itensor_version,))
    link = [np.shape(W[ns + k])[-1] for k in range(len(wil) - 2)]
    return [k for k in range(len(link)) if link[k] < link[0]]


def _check_prefix(eps, t, eps_h, t_h):
    if len(eps) < len(eps_h) or \
            np.max(np.abs(eps[:len(eps_h)] - eps_h), initial=0.0) > 1e-12 or \
            np.max(np.abs(t[:len(t_h)] - t_h), initial=0.0) > 1e-12:
        raise ValueError("the tail chain does not start with the Wilson chain "
                         "of the head")


class HeadTailSpectrum(object):
    """Result of head_tail_poles: omega, weight (the joined poles: the tail
    below w_join, TASK above), tail (the NRGSpectrum of the tail), task
    (the TASK poles of the head, the dict of head_task_poles, or None),
    w_join, omega_k0 (the scale of the hand-over shell), k0, chi (the
    size of the hand-over basis), E0 (of the head), sum_rule (sum over
    sigma of <{X, X^dag}> on the head state), missing_weight (sum_rule
    minus the weight of the tail poles at every energy: the part of
    X^dag|GS> outside the hand-over basis), overlap (the weights of both
    routes in the window w_join/window to w_join*window, by sign of w)
    and the head chain (chain)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def head_tail_poles(fc, m=4, k0=None, n_tail=None, keep=800, T=0.0, U=0.0,
                    nkry=90, w_join=None, stop_tol=1e-4, stop_extra=6,
                    wilson_tail=None, task=True, Nf=None, Sz=None,
                    window=None, solved=False, pad_bonds=None, symmetric="auto",
                    flip=None):
    """Poles of the STM spectrum A(w) of ysr.py on the chain fc (a
    ysr.spin_chain_on_metal, or a spin-1/2 ysr.kondo_chain at Delta = 0),
    with a DMRG head and an NRG tail (docs/nrg_mps_kondo_plan.md, N3
    and N4).

    The head is fc solved in its sector (head_state(fc, Nf, Sz), skipped
    with solved=True, when fc is already such a state). Its left basis at
    the Wilson site k0 (default the m-th site before the last one of fc)
    becomes the start shell of an NRG (handover_shell), which runs on the
    Wilson sites after k0 of the same mesh (tail_chain: up to n_tail
    sites, default N_MAX, or wilson_tail=(eps, t)) with keep states per
    shell. At T = 0 it stops stop_extra shells after the fixed point
    (stop_tol, nrgspectral.iterate); at T > 0 at the first shell below
    T/100. task=True adds the TASK poles of the head state
    (head_task_poles, nkry blocks of four) above w_join (default omega_k0,
    the scale of the hand-over shell) and keeps the tail below it. T > 0
    applies to the tail only: TASK is a T = 0 method, so T has to be far
    below w_join (a RuntimeWarning above w_join/10). pad_bonds goes to
    handover_shell ("full" completes the hand-over basis, for tests). It
    warns when the hand-over basis holds fewer than MIN_HANDOVER states, and
    when the MPO of H has dropped a Wilson bond inside the hand-over block
    (lost_wilson_bonds, README section 9). symmetric goes to handover_shell,
    a hand-over basis closed under the flip of all spins: "auto" (the
    default) takes it for a head with a moment (2Sz != 0) at B = 0, the
    case whose own left basis is not closed (README section 5), True asks
    for it and raises at B != 0, False is the plain left basis. flip goes
    to nrgspectral.spectral_poles (None: the average of partner spectra at
    B = 0, kept on with either basis).
    Returns a HeadTailSpectrum."""
    import warnings
    from .ysr import tunneling_operators
    from .ysrtk import zaverage
    if fc.itensor_version not in ("python", 3):
        raise NotImplementedError("head_tail_poles: itensor_version=\"python\""
                                  " or 3, got %r" % (fc.itensor_version,))
    if getattr(fc, "Delta_sc", 0.0) != 0.0:
        raise ValueError("head_tail_poles: a normal metal only (Delta = 0)")
    if symmetric not in (True, False, "auto"):
        raise ValueError("symmetric must be True, False or \"auto\", got %r"
                         % (symmetric,))
    if solved:
        c, e0 = fc, float(np.real(fc.gs_energy()))
    else:
        c, e0 = head_state(fc, Nf=Nf, Sz=Sz)
    if symmetric == "auto":
        s2 = (getattr(c, "head_sector", None) or {}).get("Sz")
        if s2 is None:
            Sz_tot = sum(c.Sz[i] for i in range(c.ns))
            s2 = int(round(2*float(np.real(c.vev(Sz_tot)))))
        symmetric = getattr(c, "B_imp", 0.0) == 0.0 and s2 != 0
    eps_h, t_h = wilson_parameters(c)
    nh = len(eps_h)
    if k0 is None:
        k0 = nh - 1 - int(m)
    if not 0 <= k0 < nh - 1:
        raise ValueError("k0 = %r: the hand-over needs a Wilson site of the "
                         "head before its last one (0..%d)" % (k0, nh - 2))
    if wilson_tail is None:
        eps, t = tail_chain(c, n_tail)
    else:
        eps, t = (np.asarray(x, dtype=float) for x in wilson_tail)
    _check_prefix(eps, t, eps_h, t_h)
    wp = getattr(c, "wilson_params", None) or {}
    Lambda = wp.get("Lambda", abs(t[0]/t[2]) if len(t) > 2 else 4.0)
    sh, ex = handover_shell(c, k0, U=U, Lambda=Lambda, pad_bonds=pad_bonds,
                            symmetric=symmetric)
    omega_k0 = sh.omega
    lost = [k for k in lost_wilson_bonds(c) if k < k0]
    if lost:
        warnings.warn("head_tail_poles: the MPO of H on this backend has dropped "
                      "channels of the Wilson bond %d (|t| = %.3g), inside the "
                      "block that hands over at the Wilson site %d, so that block "
                      "is not the one of the tail (lost_wilson_bonds; on v3 a "
                      "hopping below about 3e-7 of the unit of H is cut); hand "
                      "over at k0 <= %d" % (lost[0], abs(t_h[lost[0]]), k0, lost[0]),
                      RuntimeWarning, stacklevel=2)
    if sh.info["chi"] < MIN_HANDOVER:
        warnings.warn("head_tail_poles: the hand-over basis at the Wilson site "
                      "%d holds %d states, below the about %d the tail needs "
                      "(measured on one spin at rho J = 0.2: 35%% low with 16 "
                      "states, 3%% with 64); raise fc.maxm or m"
                      % (k0, sh.info["chi"], MIN_HANDOVER), RuntimeWarning,
                      stacklevel=2)
    if w_join is None:
        w_join = omega_k0
    if task and T > 0.1*w_join:
        warnings.warn("head_tail_poles: T = %.3g is not far below w_join = "
                      "%.3g, where the T = 0 TASK spectrum of the head takes "
                      "over" % (T, w_join), RuntimeWarning, stacklevel=2)
    r = nrgspectral.spectral_poles(
        sh, eps[k0+1:], t[k0:], keep=keep, T=T, flip=flip,
        stop_tol=stop_tol if T == 0.0 else None, stop_extra=stop_extra,
        stop_omega=T/100.0 if T > 0.0 else None)
    Xs = tunneling_operators(c, U=U, site=len(c.impurity))
    sum_rule = 0.0
    for X in Xs:
        Xd = X.get_dagger()
        sum_rule += float(np.real(c.vev(X*Xd + Xd*X)))
    if window is None:
        window = Lambda
    sel = np.abs(r.omega) < w_join
    om, wt = [r.omega[sel]], [r.weight[sel]]
    p = None
    if task:
        p = head_task_poles(c, U=U, nkry=nkry)
        selp = np.abs(p["omega"]) >= w_join
        om.append(p["omega"][selp])
        wt.append(p["weight"][selp])
    lo, hi = w_join/window, w_join*window
    ov = {}
    for name, (o, w) in (("tail", (r.omega, r.weight)),
                         ("task", (None, None) if p is None
                          else (p["omega"], p["weight"]))):
        if o is None:
            continue
        a = np.abs(o)
        for sgn, side in ((1, "positive"), (-1, "negative")):
            s = (sgn*o > 0) & (a >= lo) & (a <= hi)
            ov[(name, side)] = float(np.sum(w[s]))
        ov[(name, "below")] = float(np.sum(w[a < w_join]))
        ov[(name, "above")] = float(np.sum(w[a >= w_join]))
    return HeadTailSpectrum(omega=np.concatenate(om), weight=np.concatenate(wt),
                            tail=r, task=p, w_join=w_join, omega_k0=omega_k0,
                            k0=k0, chi=sh.info["chi"], E0=e0,
                            sum_rule=sum_rule,
                            missing_weight=r.sum_rule - r.total,
                            overlap=ov, window=(lo, hi), chain=c, T=T,
                            omega_N=r.omega_n[-1])


def stm_spectrum(es, ns=1, n_head=12, site=None, Js=0.5, JK=0.2, Lambda=4.0,
                 D=1.0, B=0.0, mesh="mean", nz=4, b=0.6, omega0=None, U=0.0,
                 m=4, keep=800, T=0.0, nkry=90, maxm=160, nsweeps=30,
                 cutoff=1e-12, itensor_version=3, continuum_norm=True,
                 task=True, tip_fermi=True, **kwargs):
    """A(w) at the energies es with a DMRG head and an NRG tail, z-averaged
    over nz meshes (z = k/nz): for each mesh, ysr.spin_chain_on_metal(ns,
    n_head, site, Js, JK, Lambda, D, B, z, mesh) solved by DMRG
    (itensor_version, maxm, nsweeps, cutoff), its head_tail_poles (m,
    keep, T, U, nkry, task and the extra keywords), the poles broadened by
    nrgspectral.curve with b and omega0 (at T = 0 the default is the
    scale of the last shell of the tail; at T > 0 with tip_fermi the
    result is dI/dV with the Fermi function of the tip), and the meshes
    averaged. With continuum_norm the result is multiplied by
    ysr.lambda_factor(Lambda) on mesh="mean" (1 on "zp"), as in
    nrgspectral.stm_spectrum. Returns (es, A, info): info["A_z"],
    info["results"] (the HeadTailSpectrum of every mesh, without its
    chain), info["z"] and info["dos_factor"]."""
    from .ysr import spin_chain_on_metal, lambda_factor
    from .ysrtk.zaverage import twists
    es = np.asarray(es, dtype=float)
    norm = lambda_factor(Lambda) if (continuum_norm and mesh == "mean") else 1.0
    A = np.zeros(len(es))
    info = dict(z=twists(nz), A_z=[], results=[], dos_factor=norm)
    for z in info["z"]:
        fc = spin_chain_on_metal(ns, n_head, site=site, Js=Js, JK=JK,
                                 Lambda=Lambda, D=D, B=B, z=z, mesh=mesh,
                                 itensor_version=itensor_version)
        fc.maxm, fc.nsweeps, fc.cutoff = maxm, nsweeps, cutoff
        r = head_tail_poles(fc, m=m, keep=keep, T=T, U=U, nkry=nkry,
                            task=task, **kwargs)
        r.chain = None
        info["results"].append(r)
    floor = nrgspectral.FLOOR*max(r.omega_N for r in info["results"])
    for r in info["results"]:
        Az = nrgspectral.curve(es, r.omega, r.weight, T=T, omega_N=r.omega_N,
                               b=b, omega0=omega0, tip_fermi=tip_fermi,
                               floor=floor)*norm
        A = A + Az
        info["A_z"].append(Az)
    return es, A/nz, info


def twin_chain(fc, z):
    """A chain with the model of fc on the mesh of twist z: a
    ysr.spin_chain_on_metal with fc's spins, Js, J_chain, B_imp, tip_spin,
    Wilson length and mesh (fc.wilson_params), and fc's backend, maxm,
    nsweeps, cutoff and noise. A spin-1/2 kondo_chain becomes the
    spin_chain_on_metal of one spin, the same Hamiltonian."""
    from .ysr import spin_chain_on_metal
    ns, wil = _layout(fc)
    wp = getattr(fc, "wilson_params", None)
    if wp is None:
        raise ValueError("the chain was built with wilson= by hand, so its "
                         "twists cannot be rebuilt; use nz=1")
    c = spin_chain_on_metal(ns, len(wil), site=getattr(fc, "tip_spin", ns - 1),
                            Js=getattr(fc, "Js", 0.0), JK=fc.J_chain,
                            Lambda=wp["Lambda"], D=wp["D"],
                            B=getattr(fc, "B_imp", 0.0), z=z, mesh=wp["mesh"],
                            itensor_version=fc.itensor_version)
    for key in ("maxm", "nsweeps", "cutoff", "noise"):
        if hasattr(fc, key):
            setattr(c, key, getattr(fc, key))
    return c


def get_stm_spectrum(fc, es, U=0.0, head="dmrg", nz=1, b=None, omega0=None,
                     T=0.0, keep=800, m=4, k0=None, n_tail=None, nkry=90,
                     stop_tol=1e-4, stop_extra=6, tip_fermi=True,
                     return_poles=False, symmetric="auto"):
    """ysr.get_stm_spectrum(fc, es, submode="NRG"): A(w) at the energies es,
    in the units of the TASK route (the chain's own normalization;
    multiply by fc.dos_factor for the continuum), averaged over nz twists
    of the mesh, z = fc's twist + k/nz (twin_chain builds the others from
    the parameters fc stores, with fc's solver settings).

    head="dmrg" (the default) is head_tail_poles(fc, m, k0, n_tail, keep,
    T, U, nkry, stop_tol, stop_extra, symmetric) for every twist: the
    DMRG head (fc's own maxm, nsweeps and cutoff), the NRG tail and the
    TASK poles of the
    head above the scale of the hand-over. head="exact" skips DMRG: the
    spins and the first Wilson site are the exact start shell
    (nrgspectral.spin_cluster_shell with fc's Js, J_chain, B_imp and
    tip_spin) and the NRG runs on the Wilson chain continued to n_tail
    sites, with the same stop rule. The poles are broadened by
    nrgspectral.curve with b (default 0.8 for one mesh at Lambda = 4, 0.4
    for nz of 4 or more, 0.6 in between) and omega0 (default: the last
    shell's scale at T = 0), all meshes read with one floor. At T > 0
    with tip_fermi (the default) the result is dI/dV with the Fermi
    function of a normal tip at T (nrgspectral.thermal_didv), else A(w).
    With return_poles the result also holds the list of the
    HeadTailSpectrum or NRGSpectrum of every twist, as (es, A, results)."""
    from .ysrtk.zaverage import twists
    if head not in ("dmrg", "exact"):
        raise ValueError("head must be \"dmrg\" or \"exact\", got %r"
                         % (head,))
    if getattr(fc, "Delta_sc", 0.0) != 0.0:
        raise ValueError("submode=\"NRG\": a normal metal only (Delta = 0)")
    if b is None:
        b = 0.8 if nz == 1 else (0.6 if nz < 4 else 0.4)
    ns, wil = _layout(fc)
    es = np.asarray(es, dtype=float)
    z0 = (getattr(fc, "wilson_params", None) or {}).get("z", 0.0)
    results, scales = [], []
    for k, dz in enumerate(twists(nz)):
        c = fc if k == 0 else twin_chain(fc, (z0 + dz) % 1.0)
        if head == "exact":
            eps, t = tail_chain(c, n_tail)
            wp = c.wilson_params
            om_start = np.sqrt(wp["Lambda"])*abs(t[0]) if len(t) else wp["D"]
            start = nrgspectral.spin_cluster_shell(
                ns=ns, site=getattr(c, "tip_spin", ns - 1),
                Js=getattr(c, "Js", 0.0), JK=c.J_chain,
                B=getattr(c, "B_imp", 0.0), U=U, eps0=eps[0], omega=om_start)
            r = nrgspectral.spectral_poles(
                start, eps[1:], t, keep=keep, T=T,
                stop_tol=stop_tol if T == 0.0 else None, stop_extra=stop_extra,
                stop_omega=T/100.0 if T > 0.0 else None)
            scales.append(r.omega_n[-1])
        else:
            r = head_tail_poles(c, m=m, k0=k0, n_tail=n_tail, keep=keep, T=T,
                                U=U, nkry=nkry, stop_tol=stop_tol,
                                stop_extra=stop_extra, symmetric=symmetric)
            scales.append(r.omega_N)
        results.append(r)
    floor = nrgspectral.FLOOR*max(scales)
    A = np.zeros(len(es))
    for r, wN in zip(results, scales):
        A = A + nrgspectral.curve(es, r.omega, r.weight, T=T, omega_N=wN, b=b,
                                  omega0=omega0, tip_fermi=tip_fermi,
                                  floor=floor)
    A = A/len(results)
    if return_poles:
        return es, A, results
    return es, A


def _session_arrays(c, ops):
    """The ground state, H and the MPO of every MultiOperator of ops as
    site arrays, on "python" and on 3, as taskdmrg exports them."""
    from . import taskdmrg
    from .pyitensor import task
    taskdmrg.ground_state(c)
    s = c._session
    if c.itensor_version == "python":
        from .pyitensor.chain import _mps_arrays_lpr
        Ms = _mps_arrays_lpr(s.wf0)
        Ws = task.mpo_arrays(s.H)
        Wo = [task.mpo_arrays(s._mpo(op.to_terms())) for op in ops]
    else:
        Ms = s.mps_site_arrays(s.gs_wavefunction())
        Ws = s.mpo_site_arrays(s.hamiltonian_mpo())
        Wo = [s.mpo_site_arrays(s.build_operator(op.to_terms()))
              for op in ops]
    return Ms, Ws, Wo


def head_task_poles(c, U=0.0, nkry=90):
    """TASK poles of the STM spectrum of the head state c (a solved chain,
    one of definite parity), from one block Krylov space
    (pyitensor.task.block_krylov, nkry blocks) started from P X^dag_sigma
    |GS> and P X_sigma|GS> of both tip spins, so that four starts share one
    basis and every application of H acts on a block of four. Returns a
    dict with "omega" (the electron poles at +w, the hole poles at -w),
    "weight", "missing_weight" (the sum rule <{X, X^dag}> minus the
    weights: the part of X^dag|GS> outside the one-site tangent space and
    the Krylov space) and "dimension" (of the Krylov basis). Only the
    column <alpha|O|GS> is computed, not the operator between Ritz states
    that task.projected_spectrum also builds."""
    from .ysr import tunneling_operators
    from .pyitensor import task, backend
    Xd = tunneling_operators(c, U=U, site=len(c.impurity))
    X = [x.get_dagger() for x in Xd]
    Ms, Ws, Wo = _session_arrays(c, Xd + X)
    ts = task.TangentSpace(Ms, Ws)
    starts = [ts.project_operator(task.compress_mpo(W)) for W in Wo]
    store, n, Hq = task.block_krylov(ts, starts, int(nkry))
    ws, V = np.linalg.eigh(Hq)
    Vh = V.conj().T
    om, wt = [], []
    for k, b in enumerate(starts):
        amp = Vh @ backend.to_host(backend.row_overlaps(store, ts.flat(b),
                                                        n))[:n]
        sign = 1.0 if k < len(Xd) else -1.0 # electron at +w, hole at -w
        om.append(sign*ws)
        wt.append(np.abs(amp)**2)
    om, wt = np.concatenate(om), np.concatenate(wt)
    total = 0.0
    for x in Xd:
        xd = x.get_dagger()
        total += float(np.real(c.vev(x*xd + xd*x)))
    return dict(omega=om, weight=wt, missing_weight=total - float(np.sum(wt)),
                dimension=n)
