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
BLOCK_TOL = 1e-6 # the block Hamiltonian couples two labels by at most this


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
    the smallest |2S^z|. On itensor_version=3 the solve starts from the
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
        Sz = Nf % 2
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
        As, _, _ = task.mixed_gauge([np.asarray(M) for M in Ms],
                                    pad_bonds=pad_bonds)
        As = [np.asarray(backend.to_host(A)) for A in As]
        self.As, self.real_residual = real_gauge(As)
        self.L = len(self.As)

    def mpo(self, op):
        from .pyitensor import task, backend
        W = task.compress_mpo([np.asarray(w) for w in self._mpo(op.to_terms())])
        return [np.asarray(backend.to_host(w)) for w in W]

    def chi(self, p):
        return self.As[p].shape[2]

    def block_matrix(self, op, p):
        """<a|O|b> in the left basis of bond p (between sites p and p+1),
        O a MultiOperator with every factor on sites <= p."""
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
            Lm = np.tensordot(X, A, axes=([0, 2], [0, 1])) # b, w', b'
        M = np.tensordot(Lm, r, axes=(1, 0))
        if np.iscomplexobj(M):
            if np.max(np.abs(M.imag)) > 1e-10*max(1.0, np.max(np.abs(M))):
                raise ValueError("block_matrix: a complex matrix in the left "
                                 "basis; nrgmps needs a real state")
            M = M.real
        return np.ascontiguousarray(M)


REAL_TOL = 1e-8 # the left basis is closed under complex conjugation to this


def real_gauge(As, tol=REAL_TOL):
    """The left isometries As of a state that is real up to a global
    phase, in a real gauge: at every bond a real orthonormal basis of the
    span of the left basis (the left singular vectors of [Re M, Im M], M
    the isometry of that bond in the gauge of the bonds before), and the
    unitary that relates the two moved into the next tensor; the last
    tensor loses the global phase. The "python" backend returns such a
    state in a complex gauge. Returns (real As, residual), the residual
    the largest relative singular value beyond the bond dimension (zero
    for a span closed under conjugation) and, last, the imaginary part
    left on the last tensor. Raises above tol."""
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
        U, S, _ = np.linalg.svd(np.concatenate([M.real, M.imag], axis=1),
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


def handover_shell(fc, k0, U=0.0, imp=None, export=None, pad_bonds=None,
                   Lambda=None):
    """The StartShell of nrgspectral at the Wilson site k0 of the DMRG
    state of fc: the block Hamiltonian, the labels, f_{p s} and the tip
    vertices X^dag_sigma (ysr.tunneling_operators with U and imp, the
    electron at the first Wilson site) in the left basis of bond
    p = ns + k0, graded and diagonalized. anti holds the projected
    anticommutator {X, X^dag} of each vertex, for the sum rule. Q in the
    labels is the charge of the whole block (spin sites included), whose
    parity is the Jordan-Wigner sign of the next site."""
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
    V, labels = graded_basis(ex.block_matrix(Nop, p), ex.block_matrix(Szop, p))
    rot = lambda op: V.T @ ex.block_matrix(op, p) @ V
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
                                       info=dict(k0=k0, p=p, chi=ex.chi(p),
                                                 U=U))
    return sh, ex


N_MAX = 120 # the longest Wilson chain the tail is given, before its stop rule


def tail_chain(fc, n_tail=None):
    """The Wilson chain the tail runs on: n_tail sites (default N_MAX) of
    the mesh of fc (fc.wilson_params), whose first sites must be the
    Wilson chain of fc itself."""
    from . import wilsonchain
    wp = getattr(fc, "wilson_params", None)
    eps_h, t_h = wilson_parameters(fc)
    if wp is None:
        raise ValueError("the head was built with wilson= by hand; pass the "
                         "tail chain as wilson_tail=(eps, t)")
    n = N_MAX if n_tail is None else int(n_tail)
    if n < len(eps_h):
        raise ValueError("n_tail=%d is shorter than the head (%d Wilson sites)"
                         % (n, len(eps_h)))
    eps, t = wilsonchain.wilson_chain(n, Lambda=wp["Lambda"],
                                      band=(-wp["D"], wp["D"]), z=wp["z"],
                                      mesh=wp["mesh"])
    return eps, t


def _check_prefix(eps, t, eps_h, t_h):
    if len(eps) < len(eps_h) or \
            np.max(np.abs(eps[:len(eps_h)] - eps_h), initial=0.0) > 1e-12 or \
            np.max(np.abs(t[:len(t_h)] - t_h), initial=0.0) > 1e-12:
        raise ValueError("the tail chain does not start with the Wilson chain "
                         "of the head")


class HeadTailSpectrum(object):
    """Result of head_tail_poles: omega, weight (the joined poles: the tail
    below w_join, TASK above), tail (the NRGSpectrum of the tail), task
    (the TASK poles of the head, a dict of zaverage.stm_poles, or None),
    w_join, omega_k0 (the scale of the hand-over shell), k0, chi (the
    size of the hand-over basis), E0 (of the head), sum_rule (sum over
    sigma of <{X, X^dag}> on the head state), missing_weight (sum_rule
    minus the weight of the tail poles at every energy: the part of
    X^dag|GS> outside the hand-over basis), overlap (the weights of both
    routes in the window w_join/window to w_join*window, by sign of w)
    and the head chain (chain)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def head_tail_poles(fc, m=2, k0=None, n_tail=None, keep=400, T=0.0, U=0.0,
                    nkry=300, w_join=None, stop_tol=1e-4, stop_extra=6,
                    wilson_tail=None, task=True, Nf=None, Sz=None,
                    window=None, solved=False):
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
    (zaverage.stm_poles with nkry) above w_join (default omega_k0, the
    scale of the hand-over shell) and keeps the tail below it. T > 0
    applies to the tail only: TASK is a T = 0 method, so T has to be far
    below w_join (a RuntimeWarning above w_join/10). Returns a
    HeadTailSpectrum."""
    import warnings
    from .ysr import tunneling_operators
    from .ysrtk import zaverage
    if fc.itensor_version not in ("python", 3):
        raise NotImplementedError("head_tail_poles: itensor_version=\"python\""
                                  " or 3, got %r" % (fc.itensor_version,))
    if getattr(fc, "Delta_sc", 0.0) != 0.0:
        raise ValueError("head_tail_poles: a normal metal only (Delta = 0)")
    if solved:
        c, e0 = fc, float(np.real(fc.gs_energy()))
    else:
        c, e0 = head_state(fc, Nf=Nf, Sz=Sz)
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
    sh, ex = handover_shell(c, k0, U=U, Lambda=Lambda)
    omega_k0 = sh.omega
    if w_join is None:
        w_join = omega_k0
    if task and T > 0.1*w_join:
        warnings.warn("head_tail_poles: T = %.3g is not far below w_join = "
                      "%.3g, where the T = 0 TASK spectrum of the head takes "
                      "over" % (T, w_join), RuntimeWarning, stacklevel=2)
    r = nrgspectral.spectral_poles(
        sh, eps[k0+1:], t[k0:], keep=keep, T=T,
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
        p = zaverage.stm_poles(c, U=U, nkry=nkry, parity=False)
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
                            missing_weight=sum_rule - r.total,
                            overlap=ov, window=(lo, hi), chain=c, T=T,
                            omega_N=r.omega_n[-1])
