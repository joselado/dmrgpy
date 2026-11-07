"""Thermodynamics of the Kondo impurity on a Wilson chain by the numerical
renormalization group (NRG), and the Kondo temperature in Wilson's
definition.

The model is the spin-1/2 Kondo impurity of ysr.kondo_chain at Delta = 0,
on exactly the same Wilson chain (wilsonchain.wilson_chain, mesh "mean" or
"zp"),

    H = J S.s_0 + B S^z_imp + sum_n eps_n N_n
        + sum_n t_n (f^dag_{n s} f_{n+1 s} + h.c.),

with the impurity a bare doublet (two states, no charge; kondo_chain holds
it as a fermion site with a penalty, whose lowest states are these). J is
the coupling on the chain, as kondo_chain takes it: J*ysr.lambda_factor
for a continuum J on mesh="mean", the continuum J bare on mesh="zp".

Two susceptibilities, two Kondo temperatures. The local one,
chi_loc = -d<S^z_imp>/dB with B on the impurity only, is what the DMRG
route of ysr measures at T = 0 from 2(E0 - E(B))/B^2. Wilson's impurity
susceptibility chi_imp is the response of the total moment to a uniform
field, with the impurity minus without it,

    T chi_imp(T) = [<(S^z_tot)^2> - <S^z_tot>^2]_T - (same, free chain),

(Wilson, RMP 47, 773; Bulla, Costi and Pruschke, RMP 80, 395, section
III.A), which is what the Kondo temperature of NRG papers is defined from:
chi_imp(T -> 0) = w/(4 T_K), w = 0.4128, and then T_K chi_imp(T_K) = 0.0701
on the universal curve. The two agree in the scaling limit and differ at
O(rho J). chi_imp cannot be read from a T = 0 ground state: S^z_tot
commutes with H, so a uniform field moves every level linearly and the
ground-state energy of a finite chain has no curvature. It lives in the
spectrum of the shells at temperatures of order their level spacing,
which is what this module computes.

The iteration. Shell N is the impurity with sites 0..N. Its eigenstates
are labelled by the bath charge Q and 2 S^z_tot, so every block is
diagonalized on its own. Site N+1 is added in the product basis
|r_N> x |alpha>, alpha in {0, up, dn, up dn}, the hopping uses the matrix
of f_{N s} carried in the eigenbasis of shell N, with the Jordan-Wigner
sign (-1)^Q of the old state, and S^z_imp is carried the same way. After
diagonalization the lowest `keep` states are kept, the cut moved up past
any (near) degenerate partner so that no SU(2) multiplet is split. The
energies are kept relative to the ground state of each shell (the
running offset is the absolute ground energy). The thermodynamics of
shell N uses all its states at T_N = omega_N/beta_bar, with
omega_N = |t_{N-1}| the energy scale of the last bond (the RMP's
Lambda^-(N-1)/2 up to the constant of the chain) and beta_bar a keyword.
The free chain is exact free fermions of the same N+1 sites:
Var(S^z) = (1/2) sum_k f(e_k)(1-f(e_k)) and S = 2 sum_k s(f(e_k)) over the
single-particle levels e_k, since S^z = (N_up - N_dn)/2 and the two spins
are independent. The levels e_k are grown site by site in the eigenbasis
of the shorter chain, dropping those above 1e4 omega_N (frozen and empty
at every later shell), so that they keep a relative precision: a dense
diagonalization of the whole chain has an absolute error of about
1e-16 D, which is larger than the signal once omega_N is below about
1e-13 D. S_imp is the entropy with the impurity minus without it, ln 2
for a free moment and 0 once it is screened.

Long chains. With (Q, 2Sz) labels and no SU(2), nothing protects the
degeneracy of a free moment: roundoff of order 1e-16 D from the first
shells would split the Sz = +-1/2 partners and act as a field on the
impurity once T_N falls below it. At B = 0 the spectra of the blocks
(Q, 2Sz) and (Q, -2Sz) are set equal at every shell (their mean), which
is the spin-flip symmetry of H, and removes that splitting.

Even and odd shells. Far below T_K the shell spectrum is the strong
coupling fixed point, a free chain one site shorter than the reference,
so T chi_imp of one shell carries an offset of order 1e-2 (6.8e-3 at
Lambda = 4, beta_bar = 0.75; 1.4e-3 at Lambda = 2, beta_bar = 1.5) whose sign
alternates with N, on top of the physical chi_imp(0) T. kondo_temperature
fits a + b T_N on the even and on the odd shells separately and takes the
mean slope; the alternating offset drops out (the two intercepts sum to
2e-7 at Lambda = 4), and so does an alternating term linear in T (the
even-odd average of the RMP is the same statement). The crossing 0.0701 is read on the even-odd averaged curve of
the RMP.

What was checked (tests/test_nrgthermo.py): with no truncation (3 bath
sites, 128 states) the ground energy equals the ED ground energy of
kondo_chain to 3e-15 and the full spectrum a dense ED of the same model to
4e-15, with a field and a twist and on both meshes; at J = 0,
T chi_imp = 1/4 and S_imp = ln 2 at every shell (1e-15 untruncated, 2e-7
and 1.5e-5 at keep=300). At Lambda = 4, rho J = 0.1, chi_loc(T -> 0) is
8605.4 against the DMRG 8604 of kondo_chain(27) on the mean mesh and
9946 against 9938 on zp. At Lambda = 4 keep = 400 is enough (keep 400,
800 and 1600 give T_K to 1e-5); beta_bar = 0.5, 0.75 and 1 give the same
T_K to 1.5e-3, and the slopes of the even and of the odd shells differ by
15% at beta_bar = 0.75, which is why their mean is taken. At Lambda = 2
beta_bar = 1.5 is needed: keep = 800 and 1600 give 1.13920e-5 and
1.13904e-5 there (keep 3200 equals keep 1600 to 1e-6), while keep = 800
at beta_bar = 0.75 is 5.6% off and warns. Long chains: at rho J = 0.1 the
T_K of n = 40, 50, 60 and 80 agree to 3e-6 (T_N down to 3e-24 D), and at
J = 0 T chi_imp stays 1/4 to 7e-8 and S_imp ln 2 to 4e-6 down to n = 60
(keep 400).

The crossing T chi_imp = 0.0701 is a finite-Lambda quantity here: it sits
at 1.111 T_K at Lambda = 4 (z-averaged or not), 1.073 T_K at Lambda = 3 and
1.034 T_K at Lambda = 2, linear in Lambda with the limit 1 at Lambda = 1.
The slope is the Lambda-independent T_K (the RMP's 4 T_K chi_imp(0) =
0.413) and the default. The crossing is definition="crossing", because it
is how NRG Ljubljana reads "Wilson's T_K" (Zitko and Pruschke,
arXiv:0810.4765: chi_imp(0) = 0.103/T_K, extracted with
T_K chi_imp(T_K) = 0.07), and so the convention of the NRG tables we
compare with, which are at Lambda = 2.
"""

import numpy as np

from . import wilsonchain

W_WILSON = 0.4128 # Wilson's number, chi_imp(0) = w/(4 T_K)
CROSSING = 0.0701 # T_K chi_imp(T_K) on the universal curve

# local basis of a spinful site: |0>, |up>, |dn>, |up dn> = c^dag_up c^dag_dn |0>
_NLOC = (0, 1, 1, 2)
_S2LOC = (0, 1, -1, 0)
# <alpha'|c_s|alpha> as {(alpha', alpha): value}, s = +1 (up), -1 (dn)
_CLOC = {+1: {(0, 1): 1.0, (2, 3): 1.0},
         -1: {(0, 2): 1.0, (1, 3): -1.0}}


class NRGResult(object):
    """Shell by shell output of kondo_nrg. Arrays over the shells N:

    N, T (the temperature T_N), omega (the scale omega_N), E0 (absolute
    ground energy of the shell), Tchi_imp, S_imp, Tchi_tot, Tchi_free and,
    when B != 0, Sz_imp (<S^z_imp>_T) and chi_loc (-<S^z_imp>_T/B). nkept
    is the number of states kept after each shell. spectrum holds the
    absolute energies of every state of the last shell, labels their
    (Q, 2 S^z_tot)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def _impurity_shell(B):
    """The bare impurity doublet as a shell: labels (Q=0, 2Sz=+-1)."""
    E = {(0, 1): np.array([0.5*B]), (0, -1): np.array([-0.5*B])}
    Sz = {(0, 1): np.array([[0.5]]), (0, -1): np.array([[-0.5]])}
    Sp = {((0, 1), (0, -1)): np.array([[1.0]])} # S^+ from 2Sz=-1 to +1
    return {"E": E, "Sz": Sz, "F": None, "Sp": Sp}


def _components(E):
    """For every label of the next shell, the list of (old label, alpha,
    offset, size) of its product basis, and its dimension."""
    comps = {}
    for lab in sorted(E):
        for a in range(4):
            new = (lab[0] + _NLOC[a], lab[1] + _S2LOC[a])
            comps.setdefault(new, []).append([lab, a])
    out = {}
    for new, cl in comps.items():
        d = 0
        for c in cl:
            k = len(E[c[0]])
            c += [d, k]
            d += k
        out[new] = (cl, d)
    return out


def _block_hamiltonian(shell, cl, d, eps, coupling, t):
    """Hamiltonian of one block of the next shell in its product basis."""
    E = shell["E"]
    H = np.zeros((d, d))
    where = {}
    for lab, a, o, k in cl:
        H[o:o+k, o:o+k] += np.diag(E[lab] + eps*_NLOC[a])
        where[(lab, a)] = (o, k)
    if coupling == "kondo": # J S.s_0 with the impurity alone as old shell
        for lab, a, o, k in cl:
            if a in (1, 2): # S^z s^z_0
                H[o:o+k, o:o+k] += t*np.diag(np.diag(shell["Sz"][lab]))*0.5*_S2LOC[a]
        for (lp, lm), Sp in shell["Sp"].items(): # (S^+ s^-_0 + h.c.)/2
            if (lp, 2) in where and (lm, 1) in where:
                o1, k1 = where[(lp, 2)]
                o2, k2 = where[(lm, 1)]
                H[o1:o1+k1, o2:o2+k2] += 0.5*t*Sp
                H[o2:o2+k2, o1:o1+k1] += 0.5*t*Sp.T
        return H
    for s in (1, -1): # t (f^dag_{N s} f_{N+1 s} + h.c.)
        F = shell["F"][s]
        for (a1, a2), c in _CLOC[s].items():
            for lab2, a, o2, k2 in cl:
                if a != a2:
                    continue
                lab1 = (lab2[0] + 1, lab2[1] + s)
                if (lab1, a1) not in where or (lab2, lab1) not in F:
                    continue
                o1, k1 = where[(lab1, a1)]
                sign = -1.0 if lab2[0] % 2 else 1.0
                blk = t*c*sign*F[(lab2, lab1)].T # <r1|f^dag|r2>
                H[o1:o1+k1, o2:o2+k2] += blk
                H[o2:o2+k2, o1:o1+k1] += blk.T
    return H


def _step(shell, eps, t, coupling, keep, omega, need_sz, flip_symmetric):
    """Add one site: diagonalize every block, return the thermodynamic
    data of the new shell (all states) and the truncated shell."""
    comps = _components(shell["E"])
    eig = {}
    for new, (cl, d) in comps.items():
        H = _block_hamiltonian(shell, cl, d, eps, coupling, t)
        e, U = np.linalg.eigh(H)
        eig[new] = (e, U)
    if flip_symmetric: # B = 0: the blocks (Q, 2Sz) and (Q, -2Sz) have one spectrum
        for (q, s2), (e, U) in list(eig.items()):
            p = (q, -s2)
            if s2 > 0 and p in eig and len(eig[p][0]) == len(e):
                avg = 0.5*(e + eig[p][0])
                eig[(q, s2)] = (avg, U)
                eig[p] = (avg.copy(), eig[p][1])
    e0 = min(e[0] for e, U in eig.values())
    allE = np.concatenate([e for e, U in eig.values()]) - e0
    if keep is None or len(allE) <= keep:
        ecut = np.inf
    else:
        srt = np.sort(allE)
        ecut = srt[keep-1]
        ecut = ecut + 1e-8*max(ecut, omega) # never split a multiplet
    kept = {lab: int(np.sum(e - e0 <= ecut)) for lab, (e, U) in eig.items()}
    kept = {lab: k for lab, k in kept.items() if k > 0}
    # S^z_imp in the new eigenbasis: all diagonal elements, kept block
    szdiag, Sz = {}, {}
    if need_sz:
        for lab, (e, U) in eig.items():
            cl, d = comps[lab]
            M = np.zeros((d, d))
            for l, a, o, k in cl:
                M[o:o+k, o:o+k] = shell["Sz"][l]
            X = M @ U
            szdiag[lab] = np.sum(U*X, axis=0)
            if lab in kept:
                k = kept[lab]
                Sz[lab] = U[:, :k].T @ X[:, :k]
    # f_{N+1 s} between kept states of the new shell
    F = {1: {}, -1: {}}
    for s in (1, -1):
        for lab, k in kept.items():
            lab2 = (lab[0] - 1, lab[1] - s) # f lowers Q by 1, 2Sz by s
            if lab2 not in kept:
                continue
            cl, d = comps[lab]
            cl2, d2 = comps[lab2]
            pos2 = {(l, a): (o, kk) for l, a, o, kk in cl2}
            M = np.zeros((d2, d))
            for (a1, a2), c in _CLOC[s].items():
                for l, a, o, kk in cl:
                    if a == a2 and (l, a1) in pos2:
                        o1 = pos2[(l, a1)][0]
                        sign = -1.0 if l[0] % 2 else 1.0
                        M[o1:o1+kk, o:o+kk] += c*sign*np.eye(kk)
            U, U2 = eig[lab][1][:, :k], eig[lab2][1][:, :kept[lab2]]
            F[s][(lab2, lab)] = U2.T @ M @ U
    new = {"E": {lab: eig[lab][0][:k] - e0 for lab, k in kept.items()},
           "Sz": Sz, "F": F}
    full = {lab: (eig[lab][0] - e0, szdiag.get(lab)) for lab in eig}
    return new, full, e0


FREE_CUT = 1e4 # single-particle levels kept below FREE_CUT*omega_N


def _free_chain_levels(eps, t):
    """Single-particle levels of the free chains 0..N, N = 0..n-1, with
    relative precision: the chain is grown site by site in the eigenbasis
    of the shorter one, as the NRG does, dropping the levels above
    FREE_CUT*omega_N, which are frozen and thermally empty at every later
    shell. A dense eigvalsh of the whole chain has an absolute error of
    order 1e-16 D, larger than the signal below omega_N of about 1e-13 D."""
    e, u = np.array([eps[0]]), np.array([1.0]) # levels, amplitude on site N
    out = [e]
    for N in range(1, len(eps)):
        m = len(e)
        h = np.diag(np.append(e, eps[N]))
        h[:m, m] = t[N-1]*u
        h[m, :m] = t[N-1]*u
        w, V = np.linalg.eigh(h)
        out.append(w)
        sel = np.abs(w) <= FREE_CUT*abs(t[N-1])
        e, u = w[sel], V[m, sel]
    return out


def _free_thermo(levels, beta):
    """Var(S^z) and entropy of free spinful fermions with single-particle
    levels `levels` at inverse temperature beta, grand canonical at mu = 0."""
    x = beta*np.abs(levels)
    ex = np.exp(-x)
    f = ex/(1.0 + ex)
    var = 0.5*np.sum(f*(1.0 - f))
    S = 2.0*np.sum(np.log1p(ex) + x*f)
    return var, S


def kondo_nrg(n, J=0.2, Lambda=4.0, D=1.0, z=0.0, mesh="mean", B=0.0,
              keep=800, beta_bar=0.75):
    """Iterative diagonalization of the spin-1/2 Kondo impurity on the
    Wilson chain of n sites (the chain of ysr.kondo_chain with the same n,
    Lambda, D, z and mesh, at Delta = 0), keeping keep states per shell
    (None keeps all). J is the coupling on the chain and B the field on the
    impurity, B S^z_imp. Thermodynamics at T_N = omega_N/beta_bar, omega_N =
    |t_{N-1}| (omega_0 = sqrt(Lambda)|t_0|). Returns an NRGResult."""
    if n < 1:
        raise ValueError("the Wilson chain needs at least one site")
    if beta_bar <= 0.0:
        raise ValueError("beta_bar must be positive")
    if keep is not None and (isinstance(keep, bool) or not
                             isinstance(keep, (int, np.integer)) or keep < 1):
        raise ValueError("keep must be None (keep all) or a positive integer, "
                         "got %r" % (keep,))
    eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D), z=z,
                                      mesh=mesh)
    need_sz = B != 0.0
    free_levels = _free_chain_levels(eps, t)
    shell = _impurity_shell(B)
    offset = 0.0
    rows = {k: [] for k in ("N", "T", "omega", "E0", "Tchi_tot", "Tchi_free",
                            "S_tot", "S_free", "Sz_imp", "nkept")}
    full = None
    for N in range(n):
        if N == 0:
            omega = np.sqrt(Lambda)*abs(t[0]) if n > 1 else D
            shell, full, e0 = _step(shell, eps[0], J, "kondo", keep, omega,
                                    need_sz, not need_sz)
        else:
            omega = abs(t[N-1])
            shell, full, e0 = _step(shell, eps[N], t[N-1], "hop", keep, omega,
                                    need_sz, not need_sz)
        offset += e0
        beta = beta_bar/omega
        Z, m1, m2, eav, szi = 0.0, 0.0, 0.0, 0.0, 0.0
        for lab, (e, szd) in full.items():
            w = np.exp(-beta*e)
            sz = 0.5*lab[1]
            Z += np.sum(w)
            m1 += sz*np.sum(w)
            m2 += sz*sz*np.sum(w)
            eav += np.sum(w*e)
            if need_sz:
                szi += np.sum(w*szd)
        m1, m2, eav, szi = m1/Z, m2/Z, eav/Z, szi/Z
        vfree, sfree = _free_thermo(free_levels[N], beta)
        rows["N"].append(N)
        rows["T"].append(omega/beta_bar)
        rows["omega"].append(omega)
        rows["E0"].append(offset)
        rows["Tchi_tot"].append(m2 - m1*m1)
        rows["Tchi_free"].append(vfree)
        rows["S_tot"].append(beta*eav + np.log(Z))
        rows["S_free"].append(sfree)
        rows["Sz_imp"].append(szi)
        rows["nkept"].append(sum(len(e) for e in shell["E"].values()))
    r = {k: np.array(v) for k, v in rows.items()}
    r["Tchi_imp"] = r["Tchi_tot"] - r["Tchi_free"]
    r["S_imp"] = r["S_tot"] - r["S_free"]
    r["chi_loc"] = -r["Sz_imp"]/B if need_sz else None
    if not need_sz:
        r["Sz_imp"] = None
    labs = sorted(full)
    r["spectrum"] = np.concatenate([full[l][0] for l in labs]) + offset
    r["labels"] = [l for l in labs for _ in full[l][0]]
    return NRGResult(J=J, Lambda=Lambda, D=D, z=z, mesh=mesh, B=B, keep=keep,
                     beta_bar=beta_bar, **r)


def even_odd_average(T, O):
    """The even-odd average of the RMP (section III.A): at shell N, half
    of O_N plus O_{N-1} and O_{N+1} interpolated linearly in T to T_N.
    Returns T and the averaged O on the inner shells."""
    T, O = np.asarray(T), np.asarray(O)
    Ta = T[1:-1]
    Oi = O[:-2] + (O[2:] - O[:-2])*(Ta - T[:-2])/(T[2:] - T[:-2])
    return Ta, 0.5*(O[1:-1] + Oi)


def _crossing(T, y, level):
    """The temperature where y(T) first rises through level
    from below, scanning from low T, interpolated in ln T."""
    for i in range(len(T)-1, 0, -1): # T decreases with the index
        lo, hi = y[i], y[i-1]
        if lo < level <= hi:
            f = (level - lo)/(hi - lo)
            return float(np.exp(np.log(T[i]) + f*(np.log(T[i-1])
                                                  - np.log(T[i]))))
    return None


def kondo_temperature(res, definition="wilson", window=None):
    """Kondo temperature of an NRGResult.

    definition="wilson" (the default): chi_imp(0) from the shells with T_N
    in window (default: below 0.03 times the 0.0701 crossing), as the mean
    of the slopes of T chi_imp = a + b T_N fitted on the even and on the
    odd shells, and T_K = w/(4 chi_imp(0)). This T_K does not depend on
    Lambda: at rho J = 0.1, z = 0, T_K/T_K(local) is 1.1176 at Lambda = 2
    (zp) and 1.1174 at Lambda = 4 (mean). The mean intercept (a_even + a_odd)/2 is zero in a converged
    run, since the even-odd offset alternates; a warning says when it is
    above 1% of T chi_imp at the top of the window, which is what a
    truncation with too small keep (or beta_bar) leaves.

    definition="crossing": T_K is the temperature where the even-odd
    averaged T chi_imp crosses 0.0701, the prescription of NRG Ljubljana,
    which writes 0.07 (the difference moves T_K by well under 1%). This
    one depends on Lambda: it is 1.034 T_K(wilson) at Lambda = 2 and 1.111
    at Lambda = 4, and it goes to T_K(wilson) as Lambda -> 1.

    "wilson" and "crossing" are meant for a run with B = 0, and warn
    otherwise (a field on the impurity of 0.09 T_K moves T_K by 0.16%).

    definition="local" (needs B != 0): chi_loc(0) as the mean of chi_loc
    over the window (default: below 0.01 T_K estimated from the last
    shell), T_K = w/(4 chi_loc(0)).

    Returns a dict with T_K and chi0 (= w/(4 T_K) for "crossing"); for
    "wilson" and "crossing" also T_cross, slope_even, slope_odd,
    intercept_even, intercept_odd, chi0_averaged (the slope of the
    even-odd averaged curve), nshells and window; for "local" spread
    (max - min of chi_loc over the window, relative) and nshells."""
    import warnings
    T = res.T
    if definition == "local":
        if res.chi_loc is None:
            raise ValueError("definition=\"local\" needs a run with B != 0")
        if window is None:
            tk0 = W_WILSON/(4.0*res.chi_loc[-1])
            window = (0.0, 0.01*tk0)
        sel = (T >= window[0]) & (T <= window[1])
        if np.sum(sel) < 2:
            raise ValueError("fewer than two shells in the window; take a "
                             "longer chain")
        c = res.chi_loc[sel]
        chi0 = float(np.mean(c))
        return {"T_K": W_WILSON/(4.0*chi0), "chi0": chi0,
                "spread": float((c.max() - c.min())/chi0),
                "nshells": int(np.sum(sel))}
    if definition not in ("wilson", "crossing"):
        raise ValueError("definition must be \"wilson\", \"crossing\" or "
                         "\"local\", got %r" % (definition,))
    if res.B != 0.0:
        warnings.warn("kondo_temperature: definition=%r on a run with a field "
                      "B = %g on the impurity; Wilson's T_K needs B = 0"
                      % (definition, res.B))
    Ta, ya = even_odd_average(T, res.Tchi_imp)
    tc = _crossing(Ta, ya, CROSSING)
    if window is None:
        if tc is None:
            raise ValueError("T chi_imp never crosses 0.0701: the chain does "
                             "not reach below T_K")
        window = (0.0, 0.03*tc)
    slopes, icpts = [], []
    nsh, tmax = 0, 0.0
    for p in (0, 1):
        sel = (T >= window[0]) & (T <= window[1]) & (res.N % 2 == p)
        if np.sum(sel) < 3:
            raise ValueError("fewer than three shells of one parity in the "
                             "window; take a longer chain")
        b, a = np.polyfit(T[sel], res.Tchi_imp[sel], 1)
        slopes.append(b)
        icpts.append(a)
        nsh += int(np.sum(sel))
        tmax = max(tmax, float(np.max(T[sel])))
    chi0 = 0.5*(slopes[0] + slopes[1])
    sel = (Ta >= window[0]) & (Ta <= window[1])
    chi_avg = np.polyfit(Ta[sel], ya[sel], 1)[0] # the averaged curve, a + b T
    a_mean = 0.5*(icpts[0] + icpts[1])
    if abs(a_mean) > 0.01*abs(chi0)*tmax:
        warnings.warn("kondo_temperature: the mean intercept of T chi_imp is "
                      "%.2e, %.1f%% of T chi_imp at the top of the window; "
                      "the run is truncated, raise keep (or beta_bar)"
                      % (a_mean, 100*abs(a_mean)/(abs(chi0)*tmax)))
    out = {"T_K": W_WILSON/(4.0*chi0), "chi0": float(chi0), "T_cross": tc,
           "slope_even": float(slopes[0]), "slope_odd": float(slopes[1]),
           "intercept_even": float(icpts[0]), "intercept_odd": float(icpts[1]),
           "chi0_averaged": float(chi_avg), "nshells": nsh,
           "window": tuple(window)}
    if definition == "crossing":
        if tc is None:
            raise ValueError("T chi_imp never crosses 0.0701: the chain does "
                             "not reach below T_K")
        out["T_K"], out["chi0"] = tc, W_WILSON/(4.0*tc)
    return out
