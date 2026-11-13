"""Spectral functions of a spin cluster on a Wilson chain by the numerical
renormalization group (NRG), with the complete basis of discarded states
and the full density matrix: the STM spectrum dI/dV of ysr.py at T = 0 and
T > 0. Milestones N0 and N1 of docs/nrg_mps_kondo_plan.md.

The model. ns spins 1/2 with a Heisenberg exchange Js, one of them (site)
Kondo coupled by JK to the first site of a Wilson chain,

    H = Js sum_i S_i.S_{i+1} + JK S_site.s_0 + B sum_i S^z_i
        + sum_n eps_n N_n + sum_n t_n (f^dag_n f_{n+1} + h.c.),

the Hamiltonian of ysr.spin_chain_on_metal, and for ns = 1 that of
ysr.kondo_chain at Delta = 0 (the spin is a bare doublet here, the chains
of ysr hold it as a penalized fermion site, whose low states are the
same). JK is the coupling on the chain, as ysr takes it. The operators are
the tunneling vertices of ysr.tunneling_operators,

    X^dag_sigma = sum_s' f^dag_{0 s'} [U delta + sigma.S_site]_{s' sigma}/2,

and the result is

    A(w) = sum_sigma sum_ab (p_a + p_b) |<a|X^dag_sigma|b>|^2
                                        delta(w - E_a + E_b),

the spectral function of the anticommutator, which at T = 0 is
sum_sigma [A^e_sigma(w) + A^h_sigma(-w)] of ysr.py, the dI/dV of a normal
tip in units of 2 pi T0^2 per unit of tip density of states. At T > 0 it
is the spectrum of the sample at T; the thermal broadening by the Fermi
function of the tip is not included.

The iteration. The start shell is the spins with the Wilson site 0,
diagonalized exactly in blocks of the bath charge Q and 2 S^z_tot. Site
n+1 is added as in nrgthermo (its block Hamiltonian, with the
Jordan-Wigner sign (-1)^Q of the old state on f_{n+1}); every shell keeps
all its eigenvalues and eigenvectors, the lowest `keep` states are kept
(the cut moved up past near degenerate partners), and the rest are
discarded. An operator on the first sites is carried as O x 1 into the
next product basis, with no sign. Together, the discarded states of every
shell and all states of the last one are a complete basis of the Fock
space of the whole chain (Anders and Schiller, cond-mat/0604517), in
which each state is approximated as an eigenstate with the energy of its
shell.

The density matrix. With the complete basis the density matrix is the
full one of Weichselbaum and von Delft (cond-mat/0607497),

    rho = sum_n sum_{s in D_n, e} 4^(N-n) exp(-E_s^n/T)/Z |s e; n><s e; n|,

e running over the 4^(N-n) states of the sites after shell n, and
Z = sum_n 4^(N-n) sum_{s in D_n} exp(-E_s^n/T). A backward pass gives its
reduced density matrix R^n on every shell: diagonal on the discarded
states, and on the kept ones the trace over the next site of R^{n+1}. At
T = 0 rho is the equal mixture of the ground manifold of the last shell,
the states within deg_tol*omega_N of its ground energy, which on a
symmetric band can hold several values of Q.

The Lehmann sum. Since B = sum_n (P_n B P_n + P_n B Q_n + Q_n B P_n) for
the projectors P_n on the discarded states of shell n and Q_n on its kept
ones, every pair (a, b) of states of shell n that are not both kept
contributes one pole at E_a - E_b with the weight

    X^dag_ab [ (X^dag R^n)_ab + (R^n X^dag)_ab ],

the electron part and the hole part of the anticommutator. The weights of
the kept-discarded pairs carry the off-diagonal elements of R^n, so a
single weight can be negative. The sum of all weights is
sum_sigma Tr[rho {X_sigma, X^dag_sigma}] exactly, at any keep and T
(Peters, Pruschke and Anders, cond-mat/0607494): the result holds the sum
rule computed on the start shell from rho (sum_rule) next to the total
weight, and their difference is roundoff.

The poles of all shells are binned on a logarithmic mesh of |w|
(bins_per_decade), the poles below w_min into one bin at w = 0, and
broadened with the log-Gaussian kernel of ysrtk.zaverage.broaden_poles.
The twist average (nz meshes z = k/nz) is stm_spectrum.

Long chains. At B = 0 the spectra of the blocks (Q, 2Sz) and (Q, -2Sz)
are set to their mean at every shell, as in nrgthermo, since roundoff
would otherwise split the two members of a free moment and act as a
field once omega_N falls below it.

What was checked (tests/test_nrgspectral.py): with no truncation the
poles are those of a dense diagonalization of the same model, at T = 0
and T > 0, one spin and several; with a decoupled tail (t = 0 after the
first sites) the truncated result is the untruncated one of the short
chain, which tests the backward pass and the kept-discarded terms; the
sum rule holds to roundoff at small keep; the poles of the tip vertex on
kondo_chain agree with its ED dynamical correlator; at J = 0 the spectrum
is (U^2/2 + 3/8) times the density of states per spin of the first site.
"""

import numpy as np

from . import wilsonchain
from .nrgthermo import _NLOC, _S2LOC, _CLOC, _components, _block_hamiltonian

LOG4 = np.log(4.0)
DEG_TOL = 1e-4 # ground manifold at T = 0: within DEG_TOL*omega_N (module docstring)
MAX_START_SPINS = 10 # the start shell is dense, 4*2^ns states


class StartShell(object):
    """The first shell of the iteration, every eigenstate of it, by label
    lab = (Q, 2Sz).

    E: {lab: energies, ascending}. F: {s: {(lower, upper): matrix}} with
    <lower|f_{last s}|upper> between all states, f of the last site of the
    shell, s = +1 (up) or -1 (dn), upper = lower + (1, s). ops: {name:
    {(out, in): matrix}}, the creation-like operators whose spectra are
    computed, between all states, and shift: {name: (dQ, d2Sz)}, with
    out = in + shift. anti: {name: {lab: matrix}}, the anticommutator
    {O, O^dag} of each operator in the same basis, for the sum rule (None
    when it cannot be computed exactly on the shell). omega is the energy
    scale of the shell. The parity (-1)^Q must be the fermion parity of
    the shell, in the ordering in which the next sites come after it."""

    def __init__(self, E, F, ops, shift, anti=None, omega=1.0, info=None):
        self.E, self.F, self.ops, self.shift = E, F, ops, shift
        self.anti, self.omega = anti, float(omega)
        self.info = {} if info is None else info


def _local_site():
    """c_up, c_dn on |0>, |up>, |dn>, |up dn> (nrgthermo's basis)."""
    c = {}
    for s in (1, -1):
        m = np.zeros((4, 4))
        for (a1, a2), v in _CLOC[s].items():
            m[a1, a2] = v
        c[s] = m
    return c


def _heisenberg_cluster(ns, Js, B):
    """H of ns spins 1/2 (Js S_i.S_{i+1} + B sum S^z_i) in the product
    basis (index 0 = up), its S^z_i and S^+_i, and 2S^z of every basis
    state."""
    sz = np.diag([0.5, -0.5])
    sp = np.array([[0.0, 1.0], [0.0, 0.0]]) # S^+|dn> = |up>

    def spin(o, i):
        return np.kron(np.kron(np.eye(2**i), o), np.eye(2**(ns - i - 1)))

    Sz = [spin(sz, i) for i in range(ns)]
    Sp = [spin(sp, i) for i in range(ns)]
    H = np.zeros((2**ns, 2**ns))
    for i in range(ns - 1):
        H = H + Js*(Sz[i] @ Sz[i+1] + 0.5*(Sp[i] @ Sp[i+1].T
                                            + Sp[i].T @ Sp[i+1]))
    for i in range(ns):
        H = H + B*Sz[i]
    bits = (np.arange(2**ns)[:, None] >> np.arange(ns)[::-1][None, :]) & 1
    s2 = np.sum(1 - 2*bits, axis=1)
    return H, Sz, Sp, [int(x) for x in s2]


def spin_cluster_shell(ns=1, site=None, Js=0.5, JK=0.2, B=0.0, U=0.0,
                       eps0=0.0, omega=1.0, operators="tip", project=None):
    """The exact start shell: ns spins 1/2 (index 0 = up) and the first
    Wilson site, with the Hamiltonian of the module docstring restricted
    to them (eps0 the onsite energy of the Wilson site). site (default
    ns-1, the end spin, as in ysr.spin_chain_on_metal) is the spin under
    the tip and the one Kondo coupled. operators="tip" gives the two
    vertices "up" and "dn" of ysr.tunneling_operators (X^dag_sigma, with
    U); "ldos" gives f^dag_{0 sigma}, the local density of states of the
    first Wilson site. omega is the scale of the shell, for its
    truncation.

    project=k keeps only the k lowest states of the bare spin chain (the
    shortcut of the plan's section on where exact NRG stops: the cut
    states then take no part in the Kondo flow); S_site is projected onto
    them. None keeps all 2^ns."""
    ns = int(ns)
    if not 1 <= ns <= MAX_START_SPINS:
        raise ValueError("spin_cluster_shell: 1 <= ns <= %d, got %r"
                         % (MAX_START_SPINS, ns))
    if site is None:
        site = ns - 1
    if site != int(site) or not 0 <= site < ns:
        raise ValueError("site must be one of the spins 0..%d, got %r"
                         % (ns - 1, site))
    site = int(site)
    Hc, Sz, Sp, s2 = _heisenberg_cluster(ns, Js, B)
    Szs, Sps = Sz[site], Sp[site]
    if project is not None:
        Hc, Szs, Sps, s2 = _project_cluster(Hc, Szs, Sps, s2, int(project))
    return cluster_shell(Hc, s2, Szs, Sps, JK=JK, U=U, eps0=eps0,
                         omega=omega, operators=operators,
                         info=dict(ns=ns, site=site, Js=Js, JK=JK, B=B, U=U,
                                   eps0=eps0, operators=operators,
                                   project=project))


def _project_cluster(Hc, Sz, Sp, s2, k):
    """The k lowest eigenstates of the cluster Hamiltonian Hc, found in
    blocks of 2S^z (s2 the 2S^z of each basis state), and Sz, Sp and Hc
    written in them. The cut is not moved past degenerate partners: k has
    to close a multiplet, which is checked."""
    s2 = np.asarray(s2)
    E, vecs, lab = [], [], []
    for m in sorted(set(s2.tolist())):
        ix = np.where(s2 == m)[0]
        e, v = np.linalg.eigh(Hc[np.ix_(ix, ix)])
        for j in range(len(e)):
            u = np.zeros(len(s2))
            u[ix] = v[:, j]
            E.append(e[j])
            vecs.append(u)
            lab.append(int(m))
    order = np.argsort(E, kind="stable")
    E = np.asarray(E)[order]
    if not 1 <= k <= len(E):
        raise ValueError("project must be between 1 and %d" % len(E))
    if k < len(E) and E[k] - E[k-1] < 1e-10*max(1.0, abs(E[k])):
        raise ValueError("project=%d splits a degenerate multiplet of the "
                         "spin chain" % k)
    P = np.array([vecs[j] for j in order[:k]]).T
    return (np.diag(E[:k]), P.T @ Sz @ P, P.T @ Sp @ P,
            [lab[j] for j in order[:k]])


def cluster_shell(Hc, s2, Sz, Sp, JK=0.2, U=0.0, eps0=0.0, omega=1.0,
                  operators="tip", info=None):
    """The start shell of any spin cluster: Hc its Hamiltonian in a basis
    of S^z eigenstates (s2 their 2S^z), Sz and Sp the S^z and S^+ of the
    spin under the tip, which is also the one coupled by JK S.s_0 to the
    first Wilson site (onsite energy eps0). operators, omega: see
    spin_cluster_shell."""
    if operators not in ("tip", "ldos"):
        raise ValueError("operators must be \"tip\" or \"ldos\", got %r"
                         % (operators,))
    m = Hc.shape[0]
    I4 = np.eye(4)
    cl = _local_site()
    c = {s: np.kron(np.eye(m), cl[s]) for s in (1, -1)}
    nloc = c[1].T @ c[1] + c[-1].T @ c[-1]
    szl = 0.5*(c[1].T @ c[1] - c[-1].T @ c[-1])
    spl = c[1].T @ c[-1] # s^+ = f^dag_up f_dn
    S_z, S_p = np.kron(Sz, I4), np.kron(Sp, I4)
    S_m = S_p.T
    H = np.kron(Hc, I4) + eps0*nloc
    H = H + JK*(S_z @ szl + 0.5*(S_p @ spl.T + S_m @ spl))
    one = np.eye(4*m)
    cd = {s: c[s].T for s in (1, -1)}
    if operators == "tip":
        ops = {"up": 0.5*(cd[1] @ (U*one + S_z) + cd[-1] @ S_p),
               "dn": 0.5*(cd[1] @ S_m + cd[-1] @ (U*one - S_z))}
    else:
        ops = {"up": cd[1], "dn": cd[-1]}
    shift = {"up": (1, 1), "dn": (1, -1)}
    labels = [(_NLOC[a], int(s2[k]) + _S2LOC[a])
              for k in range(m) for a in range(4)]
    return _shell_from_dense(H, labels, c, ops, shift, omega, info=info)


def _shell_from_dense(H, labels, c, ops, shift, omega, info=None, anti=None,
                      block_tol=1e-12):
    """Block-diagonalize a dense shell Hamiltonian H by the labels of its
    basis and write f (c[s], the last site) and the operators in the
    eigenbasis. anti ({name: dense matrix}) is the anticommutator of each
    operator in the same basis; by default it is computed from the
    matrices, which is exact when the basis is the whole space of the
    shell. H may couple two labels by at most block_tol relative to its
    largest element, and that part is dropped (info["off_block"])."""
    idx = {}
    for k, lab in enumerate(labels):
        idx.setdefault(lab, []).append(k)
    idx = {lab: np.array(v) for lab, v in idx.items()}
    lab_of = np.empty(len(labels), dtype=object)
    for k, lab in enumerate(labels):
        lab_of[k] = lab
    same = np.array([[a == b for b in labels] for a in labels])
    off = float(np.max(np.abs(H[~same]), initial=0.0))
    if off > block_tol*max(1.0, np.max(np.abs(H))):
        raise RuntimeError("the start shell Hamiltonian is not block "
                           "diagonal in (Q, 2Sz): %.2e off the blocks" % off)
    info = {} if info is None else dict(info)
    info["off_block"] = off
    E, V = {}, {}
    for lab, ix in idx.items():
        e, v = np.linalg.eigh(H[np.ix_(ix, ix)])
        E[lab], V[lab] = e, v

    def blocks(O, dlab):
        out = {}
        for li, ixi in idx.items():
            lo = (li[0] + dlab[0], li[1] + dlab[1])
            if lo not in idx:
                continue
            M = O[np.ix_(idx[lo], ixi)]
            if np.any(M):
                out[(lo, li)] = V[lo].T @ M @ V[li]
        return out

    F = {}
    for s in (1, -1): # <lower|f_s|upper>, upper = lower + (1, s)
        Fs = {}
        for (lo, li), M in blocks(c[s], (-1, -s)).items():
            Fs[(lo, li)] = M
        F[s] = Fs
    opsb = {name: blocks(O, shift[name]) for name, O in ops.items()}
    antib = {}
    for name, O in ops.items():
        A = O @ O.T + O.T @ O if anti is None else anti[name]
        antib[name] = {lab: V[lab].T @ A[np.ix_(ix, ix)] @ V[lab]
                       for lab, ix in idx.items()}
    return StartShell(E, F, opsb, dict(shift), anti=antib, omega=omega,
                      info=info)


class _Shell(object):
    """One shell of the forward pass: e (energies relative to the shell
    ground state, all states), e0 (absolute ground energy), nk (kept
    states per label), U and comps (the eigenvectors in the product basis
    of the kept states of the previous shell and one site, None for the
    start shell), ops (every operator between all states), omega."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def _truncate(e, keep, omega):
    """Kept states per label: the lowest keep, the cut moved up past near
    degenerate partners (nrgthermo's rule)."""
    allE = np.concatenate([v for v in e.values()])
    if keep is None or len(allE) <= keep:
        return {lab: len(v) for lab, v in e.items()}
    ecut = np.sort(allE)[keep - 1]
    ecut = ecut + 1e-8*max(ecut, omega)
    return {lab: int(np.sum(v <= ecut)) for lab, v in e.items()}


def _carry(blocks, nk, comps, U, shift):
    """An operator on the first sites, O x 1 in the product basis of the
    next shell, between all its eigenstates. blocks are the matrices of
    the previous shell (all states), of which the kept ones are used."""
    out = {}
    pos = {lab: {(l, a): (o, k) for l, a, o, k in cl}
           for lab, (cl, d) in comps.items()}
    for li, (cl, d) in comps.items():
        lo = (li[0] + shift[0], li[1] + shift[1])
        if lo not in comps:
            continue
        M = np.zeros((comps[lo][1], d))
        hit = False
        for l, a, o, k in cl:
            lp = (l[0] + shift[0], l[1] + shift[1])
            if (lp, a) not in pos[lo] or (lp, l) not in blocks:
                continue
            o2, k2 = pos[lo][(lp, a)]
            M[o2:o2+k2, o:o+k] = blocks[(lp, l)][:k2, :k]
            hit = True
        if hit:
            out[(lo, li)] = U[lo].T @ M @ U[li]
    return out


def _site_f(comps, U, nk):
    """<lower|f_{N s}|upper> of the site just added, between kept states,
    with the Jordan-Wigner sign (-1)^Q of the old state."""
    F = {1: {}, -1: {}}
    for s in (1, -1):
        for lab, k in nk.items():
            lab2 = (lab[0] - 1, lab[1] - s)
            if k == 0 or nk.get(lab2, 0) == 0:
                continue
            cl, d = comps[lab]
            cl2, d2 = comps[lab2]
            pos2 = {(l, a): (o, kk) for l, a, o, kk in cl2}
            M = np.zeros((d2, d))
            for (a1, a2), cv in _CLOC[s].items():
                for l, a, o, kk in cl:
                    if a == a2 and (l, a1) in pos2:
                        o1 = pos2[(l, a1)][0]
                        sign = -1.0 if l[0] % 2 else 1.0
                        M[o1:o1+kk, o:o+kk] += cv*sign*np.eye(kk)
            F[s][(lab2, lab)] = U[lab2][:, :nk[lab2]].T @ M @ U[lab][:, :k]
    return F


def _flip_average(eig):
    """B = 0: one spectrum for (Q, 2Sz) and (Q, -2Sz), their mean."""
    for (q, s2), (e, U) in list(eig.items()):
        p = (q, -s2)
        if s2 > 0 and p in eig and len(eig[p][0]) == len(e):
            avg = 0.5*(e + eig[p][0])
            eig[(q, s2)] = (avg, U)
            eig[p] = (avg.copy(), eig[p][1])


def _levels(sh, nlev):
    e = np.sort(np.concatenate([v for v in sh.e.values()]))
    return e[:nlev]/sh.omega


def iterate(start, eps, t, keep=400, flip=False, stop_tol=None, stop_nlev=16,
            stop_extra=4, stop_omega=None):
    """The forward pass: from the StartShell start, add one site per
    (eps[k], t[k]), t[k] the hopping from the last site of the previous
    shell to the new one. Returns the list of _Shell, start first.

    The stop rule (both off by default, so that every site is added):
    with stop_tol, the lowest stop_nlev levels of shell N in units of its
    scale omega_N are compared with those of shell N-2; once two
    consecutive shells agree with their partners to stop_tol (the fixed
    point, which has period two in N), stop_extra more shells are added
    and the pass stops. With stop_omega, it stops at the first shell with
    omega_N <= stop_omega. The last shell records why it stopped in
    stop ("fixed point", "omega", or "chain" when eps ran out)."""
    if len(eps) != len(t):
        raise ValueError("iterate needs one hopping per added site")
    if keep is not None and (isinstance(keep, bool) or not
                             isinstance(keep, (int, np.integer)) or keep < 1):
        raise ValueError("keep must be None (keep all) or a positive integer,"
                         " got %r" % (keep,))
    e0 = min(v[0] for v in start.E.values() if len(v))
    erel = {lab: np.asarray(v) - e0 for lab, v in start.E.items()}
    sh = _Shell(e=erel, e0=e0, de=0.0, nk=_truncate(erel, keep, start.omega),
                U=None, comps=None, ops=start.ops, omega=start.omega)
    shells = [sh]
    F = start.F
    conv, countdown, why = [], None, "chain"
    for ek, tk in zip(eps, t):
        nk = sh.nk
        kE = {lab: sh.e[lab][:k] for lab, k in nk.items() if k > 0}
        kF = {s: {key: M[:nk[key[0]], :nk[key[1]]]
                  for key, M in F[s].items()
                  if key[0] in kE and key[1] in kE}
              for s in (1, -1)}
        comps = _components(kE)
        old = {"E": kE, "F": kF}
        eig = {}
        for new, (cl, d) in comps.items():
            H = _block_hamiltonian(old, cl, d, ek, "hop", tk)
            eig[new] = np.linalg.eigh(H)
        if flip:
            _flip_average(eig)
        e0 = min(e[0] for e, _ in eig.values())
        erel = {lab: e - e0 for lab, (e, _) in eig.items()}
        omega = abs(tk) if tk != 0.0 else sh.omega # a decoupled site
        nk_new = _truncate(erel, keep, omega)
        U = {lab: V for lab, (_, V) in eig.items()}
        ops = {}
        for name, blocks in sh.ops.items():
            kept = {key: M for key, M in blocks.items()
                    if nk.get(key[0], 0) > 0 and nk.get(key[1], 0) > 0}
            ops[name] = _carry(kept, nk, comps, U, start.shift[name])
        F = _site_f(comps, U, nk_new)
        sh = _Shell(e=erel, e0=sh.e0 + e0, de=e0, nk=nk_new, U=U, comps=comps,
                    ops=ops, omega=omega)
        shells.append(sh)
        if stop_omega is not None and omega <= stop_omega:
            why = "omega"
            break
        if stop_tol is not None and countdown is None and len(shells) >= 3:
            a = _levels(shells[-1], stop_nlev)
            b = _levels(shells[-3], stop_nlev)
            m = min(len(a), len(b))
            conv.append(m > 0 and float(np.max(np.abs(a[:m] - b[:m])))
                        < stop_tol)
            if len(conv) >= 2 and conv[-1] and conv[-2]:
                countdown = int(stop_extra)
        if countdown is not None:
            if countdown <= 0:
                why = "fixed point"
                break
            countdown -= 1
    shells[-1].stop = why
    return shells


def _logsumexp(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def density_matrices(shells, T=0.0, deg_tol=DEG_TOL):
    """The reduced density matrix R^n of the full density matrix on every
    shell, {lab: matrix} over all states of the shell (kept first), and
    the weight of the discarded states of each shell, sum_s R^n_ss over
    D_n (its total over the shells is 1)."""
    N = len(shells) - 1
    RD = []
    if T > 0.0:
        # the ground energy of shell n relative to that of the last shell,
        # a sum of the small increments of the later shells, so that it
        # keeps a relative precision at low T (the absolute energies would
        # carry 1e-16 of the energy of the first shell into exp(-E/T))
        de = np.array([getattr(sh, "de", 0.0) for sh in shells])
        rel = -np.concatenate([np.cumsum(de[::-1])[::-1][1:], [0.0]])
        logs = []
        for n, sh in enumerate(shells):
            d = {}
            for lab, e in sh.e.items():
                k = sh.nk[lab] if n < N else 0
                d[lab] = (N - n)*LOG4 - (e[k:] + rel[n])/T
            RD.append(d)
            logs += [v for v in d.values() if len(v)]
        logZ = _logsumexp(np.concatenate(logs))
        RD = [{lab: np.exp(v - logZ) for lab, v in d.items()} for d in RD]
    else:
        last = shells[N]
        tol = deg_tol*last.omega
        g = sum(int(np.sum(e <= tol)) for e in last.e.values())
        for n, sh in enumerate(shells):
            if n < N:
                RD.append({lab: np.zeros(len(e) - sh.nk[lab])
                           for lab, e in sh.e.items()})
            else:
                RD.append({lab: (e <= tol)/float(g) for lab, e in sh.e.items()})
    R = [None]*(N + 1)
    R[N] = {lab: np.diag(v) for lab, v in RD[N].items()}
    for n in range(N - 1, -1, -1):
        sh, nxt = shells[n], shells[n + 1]
        RK = {lab: np.zeros((k, k)) for lab, k in sh.nk.items()}
        for lab, (cl, d) in nxt.comps.items():
            Y = nxt.U[lab] @ R[n + 1][lab] @ nxt.U[lab].T
            for l, a, o, k in cl:
                RK[l] += Y[o:o+k, o:o+k]
        R[n] = {}
        for lab, e in sh.e.items():
            k = sh.nk[lab]
            M = np.zeros((len(e), len(e)))
            M[:k, :k] = RK[lab]
            M[k:, k:] = np.diag(RD[n][lab])
            R[n][lab] = 0.5*(M + M.T)
    wD = np.array([sum(float(np.sum(v)) for v in d.values()) for d in RD])
    return R, wD


class _Collector(object):
    """Poles of all shells: raw (bins_per_decade None) or summed into
    logarithmic bins of |w| between w_min and w_max, the poles below
    w_min into one bin at w = 0."""

    def __init__(self, bins_per_decade, w_min, w_max):
        self.raw = bins_per_decade is None
        self.w, self.W = [], []
        if not self.raw:
            self.lw0 = np.log10(w_min)
            self.nper = float(bins_per_decade)
            self.nb = int(np.ceil((np.log10(w_max) - self.lw0)*self.nper)) + 1
            self.pos = np.zeros(self.nb)
            self.neg = np.zeros(self.nb)
            self.zero = 0.0
            self.w_min = w_min

    def add(self, w, W):
        sel = W != 0.0
        w, W = w[sel], W[sel]
        if self.raw:
            self.w.append(w)
            self.W.append(W)
            return
        a = np.abs(w)
        z = a < self.w_min
        self.zero += float(np.sum(W[z]))
        w, W, a = w[~z], W[~z], a[~z]
        i = np.clip(((np.log10(a) - self.lw0)*self.nper).astype(int), 0,
                    self.nb - 1)
        p = w > 0
        self.pos += np.bincount(i[p], weights=W[p], minlength=self.nb)
        self.neg += np.bincount(i[~p], weights=W[~p], minlength=self.nb)

    def poles(self):
        if self.raw:
            if not self.w:
                return np.zeros(0), np.zeros(0)
            return np.concatenate(self.w), np.concatenate(self.W)
        c = 10.0**(self.lw0 + (np.arange(self.nb) + 0.5)/self.nper)
        om = np.concatenate([-c[::-1], [0.0], c])
        wt = np.concatenate([self.neg[::-1], [self.zero], self.pos])
        sel = wt != 0.0
        return om[sel], wt[sel]


class NRGSpectrum(object):
    """Result of spectral_poles: omega and weight (the poles, binned or
    raw), per operator in by_operator, sum_rule (sum over the operators of
    Tr[rho {O, O^dag}] on the start shell, None if the start shell has no
    exact anticommutator) and total (the sum of the weights), T, the
    shells' scales omega_n, ground energies E0 and kept counts nkept,
    trace (Tr rho), the weight of the discarded states per shell
    (discarded_weight) and the spectrum of the last shell (spectrum,
    relative to its ground energy, with its labels)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def spectral_poles(start, eps, t, keep=400, T=0.0, flip=None,
                   bins_per_decade=200, w_min=None, deg_tol=DEG_TOL,
                   shells=None, stop_tol=None, stop_nlev=16, stop_extra=4,
                   stop_omega=None):
    """Poles of A(w) of every operator of start (a StartShell) on the
    chain that continues it with the onsite energies eps and hoppings t
    (t[k] links the last site of the previous shell to the site of
    eps[k]), by the complete-basis NRG of the module docstring, keeping
    keep states per shell (None keeps all), at temperature T (0 is the
    ground manifold of the last shell). flip (default: when the start
    shell is symmetric under the flip of all spins, which holds at B = 0)
    sets the spectra of (Q, 2Sz) and (Q, -2Sz) to their mean on every
    shell. bins_per_decade=None returns every pole as it is (for short
    chains and tests); otherwise the poles are binned, those with |w|
    below w_min (default 1e-4 times the smaller of omega_N and T>0) into
    one bin at w = 0. shells, the output of iterate, skips the forward
    pass. stop_tol, stop_nlev, stop_extra and stop_omega are the stop rule
    of iterate (off by default: every site of eps is added). Returns an
    NRGSpectrum; its stop says why the forward pass ended."""
    eps, t = np.asarray(eps, dtype=float), np.asarray(t, dtype=float)
    if T < 0.0:
        raise ValueError("T must be >= 0, got %r" % (T,))
    if flip is None:
        flip = _flip_symmetric(start)
    if shells is None:
        shells = iterate(start, eps, t, keep=keep, flip=flip,
                         stop_tol=stop_tol, stop_nlev=stop_nlev,
                         stop_extra=stop_extra, stop_omega=stop_omega)
    N = len(shells) - 1
    R, wD = density_matrices(shells, T=T, deg_tol=deg_tol)
    trace = sum(float(np.trace(M)) for M in R[0].values())
    spread = max(float(max(np.max(e) for e in sh.e.values())) for sh in shells)
    if w_min is None:
        w_min = 1e-4*(min(shells[-1].omega, T) if T > 0 else shells[-1].omega)
    by_op, total = {}, 0.0
    allc = _Collector(bins_per_decade, w_min, 2.0*spread + w_min)
    for name in start.ops:
        col = _Collector(bins_per_decade, w_min, 2.0*spread + w_min)
        for n, sh in enumerate(shells):
            last = n == N
            for (lo, li), M in sh.ops[name].items():
                W = M*(M @ R[n][li] + R[n][lo] @ M)
                if not last:
                    W[:sh.nk[lo], :sh.nk[li]] = 0.0
                w = sh.e[lo][:, None] - sh.e[li][None, :]
                col.add(w.ravel(), W.ravel())
                allc.add(w.ravel(), W.ravel())
        om, wt = col.poles()
        by_op[name] = (om, wt)
        total += float(np.sum(wt))
    om, wt = allc.poles()
    sum_rule = None
    if start.anti is not None:
        sum_rule = 0.0
        for name, blocks in start.anti.items():
            for lab, M in blocks.items():
                sum_rule += float(np.sum(R[0][lab]*M))
    last = shells[-1]
    labs = sorted(last.e)
    return NRGSpectrum(omega=om, weight=wt, by_operator=by_op,
                       sum_rule=sum_rule, total=total, T=T, trace=trace,
                       omega_n=np.array([sh.omega for sh in shells]),
                       E0=np.array([sh.e0 for sh in shells]),
                       nkept=np.array([sum(sh.nk.values()) for sh in shells]),
                       discarded_weight=wD,
                       spectrum=np.concatenate([last.e[l] for l in labs]),
                       labels=[l for l in labs for _ in last.e[l]],
                       keep=keep, flip=flip,
                       stop=getattr(shells[-1], "stop", "chain"),
                       nshells=len(shells))


def _flip_symmetric(start, tol=1e-13):
    """Whether the spectra of (Q, 2Sz) and (Q, -2Sz) may be averaged: when
    the start shell records its field (info["B"], set by
    spin_cluster_shell and nrgmps.handover_shell), exactly when B == 0;
    otherwise when every block of the start shell has the spectrum of its
    partner to tol relative to its largest energy."""
    if "B" in start.info:
        return start.info["B"] == 0.0
    scale = max(1.0, max(float(np.max(np.abs(e))) for e in start.E.values()))
    for (q, s2), e in start.E.items():
        p = start.E.get((q, -s2))
        if p is None or len(p) != len(e) or np.max(np.abs(p - e)) > tol*scale:
            return False
    return True


def chain_poles(ns=1, n=20, site=None, Js=0.5, JK=0.2, Lambda=4.0, D=1.0,
                B=0.0, z=0.0, mesh="mean", wilson=None, U=0.0, keep=400,
                T=0.0, operators="tip", project=None, **kwargs):
    """Poles of the STM spectrum of ysr.spin_chain_on_metal(ns, n, site,
    Js, JK, Lambda, D, B, z, mesh, wilson) (ysr.kondo_chain(n, J=JK, B=B)
    for ns = 1) by spectral_poles, from the exact start shell of the spins
    and the first Wilson site. wilson=(eps, t) replaces the Wilson chain
    (n values and n-1). The scale of the start shell is sqrt(Lambda)|t_0|,
    as in nrgthermo, and that of shell k is |t_{k-1}|. Extra keywords go
    to spectral_poles; project goes to spin_cluster_shell. Returns an
    NRGSpectrum, with eps and t."""
    if wilson is None:
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D), z=z,
                                          mesh=mesh)
    else:
        eps, t = (np.asarray(x, dtype=float) for x in wilson)
        if len(eps) != n or len(t) != n - 1:
            raise ValueError("wilson=(eps, t) needs %d energies and %d "
                             "hoppings" % (n, n - 1))
    omega = np.sqrt(Lambda)*abs(t[0]) if n > 1 and t[0] != 0 else D
    start = spin_cluster_shell(ns=ns, site=site, Js=Js, JK=JK, B=B, U=U,
                               eps0=eps[0], omega=omega, operators=operators,
                               project=project)
    r = spectral_poles(start, eps[1:], t, keep=keep, T=T, **kwargs)
    r.eps, r.t = eps, t
    return r


TINY = 1e-150 # w = 0 is the mean of the limits from both sides (TINY^2 > 0)


def broaden(es, omega, weight, b=0.6, omega0=None):
    """ysrtk.zaverage.broaden_poles at gap 0, except at es == 0, where
    broaden_poles gives 0 and this gives the mean of the limits from
    w > 0 and w < 0."""
    from .ysrtk.zaverage import broaden_poles
    es = np.asarray(es, dtype=float)
    A = broaden_poles(es, omega, weight, b=b, omega0=omega0)[0]
    z = es == 0.0
    if np.any(z):
        side = broaden_poles(np.array([TINY, -TINY]), omega, weight, b=b,
                             omega0=omega0)[0]
        A[z] = 0.5*(side[0] + side[1])
    return A


# the window and the grid of the thermal convolution, in units of T
THERMAL_WINDOW = 40.0 # -f' is below 4e-18 of its peak beyond 40 T
THERMAL_POINTS = 16 # points per T
THERMAL_FAR = 60.0 # |V| beyond this, in units of T: dI/dV = A_b(V)


def thermal_didv(es, omega, weight, T, b=0.6, omega0=None):
    """dI/dV(V) = int dw A_b(w) [-f'(w - V)] at the biases es: the poles
    broadened by broaden (b, omega0, default T/4) and convolved with the
    derivative of the Fermi function of a normal tip at the temperature
    T, which is what STM measures at T > 0. Within |w| < T the NRG at
    T > 0 has the right weight but not its shape (the module docstring),
    and the convolution averages over that window. For |V| above
    THERMAL_FAR*T, where A_b changes on the scale b|V|, the convolution is
    left out (its correction is of order (T/V)^2)."""
    es = np.asarray(es, dtype=float)
    if T <= 0.0:
        raise ValueError("thermal_didv needs T > 0")
    if omega0 is None:
        omega0 = T/4.0
    out = broaden(es, omega, weight, b=b, omega0=omega0)
    m = int(THERMAL_WINDOW*THERMAL_POINTS)
    x = (np.arange(-m, m) + 0.5)*(T/THERMAL_POINTS) # never w - V = 0
    k = 1.0/(4.0*T*np.cosh(x/(2.0*T))**2)
    dx = x[1] - x[0]
    for i in np.where(np.abs(es) < THERMAL_FAR*T)[0]:
        out[i] = np.sum(broaden(es[i] + x, omega, weight, b=b,
                                omega0=omega0)*k)*dx
    return out


def stm_spectrum(es, ns=1, n=40, nz=4, b=0.6, omega0=None,
                 continuum_norm=True, Lambda=4.0, mesh="mean", T=0.0,
                 tip_fermi=True, **kwargs):
    """A(w) at the energies es, z-averaged over nz meshes (z = k/nz): the
    poles of chain_poles for each mesh, broadened by broaden (the kernel
    of ysrtk.zaverage.broaden_poles) with width b and crossover omega0
    (default: the scale omega_N of the last shell at T = 0), and
    averaged; at T = 0, |w| below FLOOR times the largest omega_N is read
    there (curve). At T > 0 with tip_fermi (the default) the result is dI/dV,
    the poles convolved with -f' of a normal tip at T (thermal_didv,
    omega0 default T/4); with tip_fermi=False it is A(w) itself, broadened
    with omega0 default T, whose shape within |w| < omega0 is that of the
    kernel and not of the spectrum. At T > 0 the chain stops at the first
    shell below T/100 unless stop_omega is given. With continuum_norm
    the result is multiplied by ysr.lambda_factor(Lambda) on mesh="mean"
    (1 on "zp"), so that it stands for the continuum, as in
    zaverage_stm_spectrum. Extra keywords go to chain_poles. Returns
    (es, A, info) with info["A_z"], info["results"] (the NRGSpectrum of
    each mesh), info["z"] and info["dos_factor"]."""
    from .ysrtk.zaverage import twists
    from .ysr import lambda_factor
    es = np.asarray(es, dtype=float)
    norm = lambda_factor(Lambda) if (continuum_norm and mesh == "mean") else 1.0
    if T > 0.0:
        kwargs.setdefault("stop_omega", T/100.0)
    A = np.zeros(len(es))
    info = dict(z=twists(nz), A_z=[], dos_factor=norm)
    info["results"] = [chain_poles(ns=ns, n=n, z=z, Lambda=Lambda, mesh=mesh,
                                   T=T, **kwargs) for z in info["z"]]
    floor = FLOOR*max(r.omega_n[-1] for r in info["results"])
    for r in info["results"]:
        Az = curve(es, r.omega, r.weight, T=T, omega_N=r.omega_n[-1], b=b,
                   omega0=omega0, tip_fermi=tip_fermi, floor=floor)*norm
        A = A + Az
        info["A_z"].append(Az)
    return es, A/nz, info


FLOOR = 30.0 # at T = 0, |w| below FLOOR*omega_N is read at FLOOR*omega_N


def curve(es, omega, weight, T=0.0, omega_N=None, b=0.6, omega0=None,
          tip_fermi=True, floor=None):
    """The broadened spectrum of the poles. At T = 0: broaden with omega0
    (default omega_N, the scale of the last shell), and every |w| below
    floor (default FLOOR*omega_N) read at sign(w)*floor (w = 0 as the
    mean of both sides), since below the last shell the curve is that of
    the kernel, while past the fixed point (the stop rule of iterate) the
    spectrum no longer changes. A twist average has to pass one floor
    for all its meshes, or each mesh is read at the same place of its own
    comb. At T > 0: thermal_didv with omega0 (default T/4) when
    tip_fermi, else broaden with omega0 (default T)."""
    es = np.asarray(es, dtype=float)
    if T > 0.0 and tip_fermi:
        return thermal_didv(es, omega, weight, T, b=b, omega0=omega0)
    if omega0 is None:
        omega0 = omega_N if T == 0.0 else T
    if T > 0.0 or omega_N is None:
        return broaden(es, omega, weight, b=b, omega0=omega0)
    lo = FLOOR*omega_N if floor is None else float(floor)
    low = np.abs(es) < lo
    ev = np.where(low, np.sign(es)*lo, es)
    A = broaden(ev, omega, weight, b=b, omega0=omega0)
    z = es == 0.0
    if np.any(z):
        side = broaden(np.array([lo, -lo]), omega, weight, b=b, omega0=omega0)
        A[z] = 0.5*(side[0] + side[1])
    return A
