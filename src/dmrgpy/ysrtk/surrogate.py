"""Surrogate BCS baths: several Kondo spins on a superconductor, each
coupled to a few discrete BCS levels.

This is the surrogate model of Baran, Frost and Paaske, arXiv:2307.11646.
A flat BCS band of half-width D and gap Delta, seen from one point, has
the Nambu hybridization Gamma g(w_n) (i w_n + Delta tau_x) on the
Matsubara axis, with

    g(w) = (1/pi) int_{-D}^{D} dxi / (xi^2 + Delta^2 + w^2)
         = (2/pi) arctan(D/sqrt(Delta^2+w^2)) / sqrt(Delta^2+w^2).

Every BCS level xi with the same gap gives one term 1/(xi^2+Delta^2+w^2),
so a few levels can fit it:

    g~(w) = sum_l gamma_l / (xi_l^2 + Delta^2 + w^2),

L levels in pairs +-xi_l with equal weights, plus one level at xi=0 when
L is odd. gamma_l and xi_l minimize sum_j |g(w_j) - g~(w_j)|^2 on 1000
points log-spaced in [1e-3 Delta, wc]. This is exactly the fit of the
paper (done there with Mathematica), and the names gamma, xi are those of
the two public codes of the authors, both MIT licensed:
github.com/virgilb91/surrogate_models (DMRG, ITensors.jl, with the fitted
parameters) and github.com/emilfrost/fermpy (exact diagonalization).
fit_surrogate gives back their parameters at D = wc = 10 Delta to 1e-7 for
L = 1 to 8 (tests/test_ysr_surrogate.py holds L = 1 to 4).

The Kondo exchange. A spin couples by J S.s to the local orbital psi of
the substrate, whose Green function is pi rho g(w) (i w + Delta tau_x)
with rho = 1/(2D). The surrogate reproduces it with
psi = sum_l sqrt(w_l) a_l, w_l = pi rho gamma_l, and psi is normalized
only if W = sum_l w_l = pi rho sum_l gamma_l is one. With f_1 = psi/sqrt(W)
the normalized orbital, J S.s_psi = J W S.s_{f_1}, so the coupling on the
surrogate is J*W. Said differently, the surrogate is a band of width
D~ = pi sum_l gamma_l / 2 (the effective bandwidth of the paper) and
rho~ J W = rho J. At D = wc = 10 Delta, W is 0.15, 0.40, 0.66, 0.83 for
L = 1 to 4 and 0.99 at L = 8; it falls when D >> Delta (see below).

All levels share one Delta, which is the identity in level space, so the
star of L levels tridiagonalizes into a chain of L sites with uniform
onsite pairing, the same argument as for the Wilson chain of ysr.py. The
spin couples to the first site only. With symmetric levels the onsite
energies vanish.

The spins are spin 1/2, each a singly occupied native spinful site with
the charge penalty of ysr.kondo_chain, so that Sz sectors are bath parity
sectors: with N spins, 2Sz = N + N_bath mod 2. The sites are ordered in
blocks [spin i, f_1 .. f_L of bath i], so a hopping t_b between the first
sites of neighbouring baths and a direct exchange J12 S_i.S_{i+1} both have
a range of L+1 sites. A potential scattering V at each spin breaks
particle-hole symmetry, which the hybridization of two YSR states needs.
t_b is a phenomenological stand-in for the substrate between the atoms,
which in a real metal decays as sin(k_F r)/(k_F r) exp(-r/xi_0); here it
is a parameter.

What a few levels cannot give is the Kondo scale. T_K comes from the
logarithm between Delta and D, and the fit weights w < Delta, so at
D >> Delta the surrogate is a band of width D~ = W D << D (W = 0.83, 0.25,
0.03 at L = 4 for D/Delta = 10, 100, 1000) and the logarithm between D~
and D is missing. Measured against the Wilson chain of ysr.py (D = 1):
at rho J = 0.15 and 0.2 no L from 1 to 4 reaches the screened side at any
Delta from 1e-4 to 1, where the Wilson chain has Delta_c = 2.0e-3 and
1.4e-2 (Lambda = 4); at rho J = 0.3, L = 4 puts Delta_c at 0.073 against
0.091 (Lambda = 2), and at rho J = 0.45 at 0.33 against 0.34 (Lambda =
1.5). For Delta >= 0.1 D, against Lambda = 1.5, L = 4 gives E_odd - E_even
within 0.06 Delta at Delta = 0.1 D and within 0.01 Delta above 0.5 D, for
rho J = 0.15 to 0.45. So the
surrogate is a model for Delta not much smaller than D, or for a spin
whose exchange is given rather than derived from a band; the classical
YSR energy, which has no logarithm, is right to 0.004 to 0.05 Delta at
L = 4 for any D/Delta (tests/test_ysr_surrogate.py).

Sector solves need care. H conserves the charge of each spin site, so
a DMRG solve keeps the impurity charges of its start; sector_energies
starts from a state with one electron on each (_warm_sector_solve) and
checks it. On itensor_version=3 the odd bath sector of a chain converges
slowly, since one quasiparticle has to spread over all the baths: four
spins at L = 1 needed 80 sweeps for 1e-12. Near a degeneracy of the
ground state the plain solve that get_stm_spectrum uses can mix in a
state of the other parity, which shows up in dI/dV (1e-2 of the peak
height at a gap of 1.7e-3 on a dimer); get_stm_spectrum(parity=True)
starts from the sector ground state instead (ysrtk/paritystate.py).
"""

import fractions
import functools

import numpy as np

from .. import fermionchain


def g_function(w, Delta=1.0, D=10.0):
    """The g-function of a flat BCS band at the Matsubara frequencies w."""
    s = np.sqrt(Delta*Delta + np.asarray(w, dtype=float)**2)
    return 2.0/np.pi*np.arctan(D/s)/s


def g_surrogate(w, gamma, xi, Delta=1.0):
    """The g-function of the levels xi with weights gamma."""
    w = np.asarray(w, dtype=float)
    gamma, xi = np.asarray(gamma, dtype=float), np.asarray(xi, dtype=float)
    den = xi[None, :]**2 + Delta*Delta + w[:, None]**2
    return np.sum(gamma[None, :]/den, axis=1)


def _unpack(p, L):
    """Levels from the log-parameters: K pair weights, K pair energies and
    the weight of the zero level when L is odd. The order is that of the
    authors' files: the zero level first, then +xi, -xi for each pair."""
    K = L//2
    x = np.exp(p)
    gamma, xi = [], []
    if L % 2 == 1:
        gamma.append(x[2*K])
        xi.append(0.0)
    for k in range(K):
        gamma += [x[k], x[k]]
        xi += [x[K+k], -x[K+k]]
    return np.array(gamma), np.array(xi)


@functools.lru_cache(maxsize=256)
def _fit_unit_gap(L, D, wc, wmin, npoints, nstart, seed):
    from scipy.optimize import least_squares
    w = np.logspace(np.log10(wmin), np.log10(wc), npoints)
    y = g_function(w, 1.0, D)
    K = L//2
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(nstart):
        xi0 = np.sort(np.exp(rng.uniform(np.log(0.2), np.log(max(D, 2.0)), K)))
        ga0 = np.exp(rng.uniform(-2.0, 1.0, K))
        p0 = np.concatenate([np.log(ga0), np.log(xi0)]
                            + ([[np.log(0.5)]] if L % 2 == 1 else []))
        r = least_squares(lambda p: g_surrogate(w, *_unpack(p, L)) - y, p0,
                          xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    gamma, xi = _unpack(best.x, L)
    return tuple(gamma), tuple(xi), 2.0*best.cost


def fit_surrogate(L, Delta=1.0, D=10.0, wc=None, wmin=1e-3, npoints=1000,
                  nstart=40, seed=0):
    """Fit L surrogate levels to the g-function of a flat BCS band of
    half-width D and gap Delta, as in arXiv:2307.11646.

    wc is the upper end of the Matsubara grid (default D; the paper uses
    wc = D = 10 Delta) and wmin its lower end, both in units of Delta for
    wmin and of energy for wc. The fit is a least squares in the
    logarithms of gamma and xi from nstart random starts with a fixed
    seed, so it is deterministic. The cost function is the one of the
    paper, sum |g - g~|^2, dominated by w < Delta where g is largest, so
    the tail up to D, which holds the bandwidth, is fitted loosely and W
    falls below one when D >> Delta (0.83, 0.25, 0.03 at L=4 for D/Delta
    = 10, 100, 1000). Returns a dict with gamma and xi (arrays
    of L values, energies), chi2 (in units of 1/Delta^2), W = pi rho
    sum(gamma), the weight of the local orbital in the surrogate (see the
    module docstring), and D_tilde = pi sum(gamma)/2, the effective
    bandwidth."""
    L = int(L)
    if L < 1:
        raise ValueError("the surrogate needs at least one level")
    if Delta <= 0:
        raise ValueError("the surrogate bath needs a gap, Delta > 0")
    if wc is None:
        wc = D
    g, x, chi2 = _fit_unit_gap(L, float(D)/Delta, float(wc)/Delta,
                               float(wmin), int(npoints), int(nstart),
                               int(seed))
    gamma, xi = np.array(g)*Delta, np.array(x)*Delta
    return dict(gamma=gamma, xi=xi, chi2=chi2, L=L, Delta=Delta, D=D, wc=wc,
                W=np.pi*np.sum(gamma)/(2.0*D),
                D_tilde=np.pi*np.sum(gamma)/2.0)


def star_to_chain(xi, gamma):
    """Lanczos map of the star of levels xi with weights gamma, started
    from the normalized orbital sqrt(gamma/sum(gamma)): onsite energies
    eps (L values) and hoppings t (L-1 values), with full
    reorthogonalization. A spectrum symmetric under xi -> -xi with equal
    weights has eps = 0 exactly, and it is set so."""
    xi, gamma = np.asarray(xi, dtype=float), np.asarray(gamma, dtype=float)
    n = len(xi)
    V = np.zeros((n, n))
    V[:, 0] = np.sqrt(gamma/np.sum(gamma))
    eps, t = np.zeros(n), np.zeros(n-1)
    for k in range(n):
        u = xi*V[:, k]
        eps[k] = np.dot(V[:, k], u)
        if k == n-1:
            break
        u = u - eps[k]*V[:, k]
        if k > 0:
            u = u - t[k-1]*V[:, k-1]
        for _ in range(2):
            u = u - V[:, :k+1] @ (V[:, :k+1].T @ u)
        t[k] = np.linalg.norm(u)
        V[:, k+1] = u/t[k]
    order = np.argsort(xi)
    if np.allclose(xi[order], -xi[order][::-1]) and \
       np.allclose(gamma[order], gamma[order][::-1]):
        eps[:] = 0.0
    return eps, t


def bath_hamiltonian(N, eps, t, t_b=0.0):
    """Single-particle matrix of N surrogate chains, blocks of L sites,
    joined by t_b between their first sites."""
    L = len(eps)
    h1 = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    h = np.kron(np.eye(N), h1)
    for i in range(N-1):
        h[i*L, (i+1)*L] = h[(i+1)*L, i*L] = t_b
    return h


def surrogate_kondo_chain(N=1, J=0.2, Delta=0.1, L=None, D=1.0, t_b=0.0,
                          J12=0.0, B=0.0, V=0.0, wc=None,
                          Upen=None, fit=None, **kwargs):
    """N spin-1/2 Kondo impurities on a superconductor, each exchange
    coupled to its own surrogate BCS bath of L levels, as a
    Spinful_Fermionic_Chain_Native.

    H = sum_i [J W S_i.s_{f_1^i} + V W n_{f_1^i} + B S_i^z + Upen (N_i - 1)^2]
        + sum_i H_bath^i + t_b sum_i sum_s (f^dag_{1,s}^i f_{1,s}^{i+1} + h.c.)
        + J12 sum_i S_i.S_{i+1},
    H_bath = sum_n t_n (f^dag_n f_{n+1} + h.c.)
             + Delta sum_n (f^dag_{n up} f^dag_{n dn} + h.c.).

    J is the continuum exchange (rho J = J/(2D)); the chain gets J*W, see
    the module docstring. J may be a list of N values. V is a potential
    scattering at the spin, mapped like J; at V = 0 the model is
    particle-hole symmetric and the YSR states of two spins do not split
    under t_b (Schmid et al., arXiv:2203.15011, eq. A8, in the
    zero-bandwidth model: the splitting is t_b |u^2 - v^2|, with u and v
    the electron and hole amplitudes of the YSR state on the coupled
    site); at V != 0 it is linear in t_b. L is the number of levels per
    bath (default 2). fit (default fit_surrogate(L, Delta, D, wc)) can be
    passed to reuse a fit; it must have been made at this Delta, D, L (if
    given) and wc (if given), or ValueError is raised. Spin i
    is on site i(L+1) and its bath on the L sites after it. The chain
    stores fc.impurities, fc.bath_sites, fc.surrogate (the fit) and the
    attributes ysr.py reads (impurity, spin, channels, Delta_sc, qp_gap),
    so for N = 1 ysr.subgap_energies works unchanged, and
    ysr.get_stm_spectrum works at any spin given imp and site (or through
    get_stm_spectrum below)."""
    N = int(N)
    if N < 1:
        raise ValueError("at least one spin")
    if fit is None:
        fit = fit_surrogate(2 if L is None else L, Delta=Delta, D=D, wc=wc)
    else:
        for k, v in (("Delta", Delta), ("D", D), ("L", L), ("wc", wc)):
            if v is not None and not np.isclose(fit[k], v, rtol=1e-12,
                                                atol=0.0):
                raise ValueError("the fit was made at %s=%g, not at the "
                                 "%s=%g given here" % (k, fit[k], k, v))
    L = len(fit["xi"])
    eps, t = star_to_chain(fit["xi"], fit["gamma"])
    Js = np.broadcast_to(np.asarray(J, dtype=float), (N,))
    Jc = Js*fit["W"]
    if Upen is None:
        Upen = 2.0*(D + np.max(np.abs(Jc)) + abs(B) + abs(J12) + abs(t_b)
                    + abs(V))
    fc = fermionchain.Spinful_Fermionic_Chain_Native(N*(L+1), **kwargs)
    imps = [i*(L+1) for i in range(N)]
    baths = [list(range(i*(L+1)+1, (i+1)*(L+1))) for i in range(N)]
    def dot(a, b):
        return fc.Sx[a]*fc.Sx[b] + fc.Sy[a]*fc.Sy[b] + fc.Sz[a]*fc.Sz[b]
    def hop(a, b, v):
        return _hop(fc, a, b, v)
    h = 0
    for i in range(N):
        s, f = imps[i], baths[i]
        h = h + Jc[i]*dot(s, f[0]) + B*fc.Sz[s]
        if V != 0.0:
            h = h + V*fit["W"]*(fc.Nup[f[0]] + fc.Ndn[f[0]])
        n0 = fc.Nup[s] + fc.Ndn[s] - 1
        h = h + Upen*n0*n0
        for k in range(L):
            if eps[k] != 0.0:
                h = h + eps[k]*(fc.Nup[f[k]] + fc.Ndn[f[k]])
            h = h + Delta*(fc.Cdagup[f[k]]*fc.Cdagdn[f[k]]
                           + fc.Cdn[f[k]]*fc.Cup[f[k]])
        for k in range(L-1):
            h = h + hop(f[k], f[k+1], t[k])
        if i < N-1:
            if t_b != 0.0:
                h = h + hop(f[0], baths[i+1][0], t_b)
            if J12 != 0.0:
                h = h + J12*dot(s, imps[i+1])
    fc.set_hamiltonian(h)
    fc.impurities, fc.bath_sites = imps, baths
    fc.impurity, fc.channels = [imps[0]], [baths[0]]
    fc.spin = fractions.Fraction(1, 2)
    fc.Delta_sc, fc.Dz, fc.B_imp = Delta, 0.0, B
    fc.surrogate, fc.J_chain, fc.Upen = fit, Jc, Upen
    hb = bath_hamiltonian(N, eps, t, t_b)
    nb = len(hb)
    HB = np.block([[hb, Delta*np.eye(nb)], [Delta*np.eye(nb), -hb]])
    fc.qp_gap = float(np.min(np.abs(np.linalg.eigvalsh(HB))))
    fc.reaches_gap = False # no bound on E_odd - E_even in a small bath
    return fc


def _hop(fc, a, b, v):
    return v*(fc.Cdagup[a]*fc.Cup[b] + fc.Cdagdn[a]*fc.Cdn[b]
              + fc.Cdagup[b]*fc.Cup[a] + fc.Cdagdn[b]*fc.Cdn[a])


def _impurity_charge_error(fc, mode):
    """sum_i <(N_i - 1)^2> over the spin sites: zero in every physical
    state."""
    q = 0
    for s in fc.impurities:
        n0 = fc.Nup[s] + fc.Ndn[s] - 1
        q = q + n0*n0
    return float(np.real(fc.vev(q, mode=mode)))


def _warm_sector_solve(c):
    """Ground state of a sector by DMRG, started where every spin site
    holds one electron.

    H conserves the charge of each spin site exactly, so a DMRG solve
    keeps the impurity charges of its start. The sector start of
    itensor_version=3 is a sum of at most 16 random product states, so
    with N spins it often has no component with all of them singly
    occupied, and the solve stops Upen or more above the answer (measured
    at N = 4). So we first solve with a hopping 0.1 Upen between each spin
    site and its bath, under which the penalty drives the charges to one,
    and then sweep the true H from that state, where the wrong-charge part
    is Upen higher and drops out of the local solves. The first solve gets
    the full nsweeps: with a quarter of them the odd sector of four spins
    on v3 ended 1.2e-3 above ED instead of 4.7e-4."""
    h = c.hamiltonian
    w = 0
    for s, f in zip(c.impurities, c.bath_sites):
        w = w + _hop(c, s, f[0], 0.1*c.Upen)
    c.set_hamiltonian(h + w)
    c.gs_energy()
    wf = c.get_gs()
    c.set_hamiltonian(h)
    c.set_initial_wf_guess(wf)
    return c.gs_energy()


def sector_energies(fc, mode="DMRG", sectors=None, nstates=1,
                    observables=None):
    """Lowest nstates energies in each Sz sector (2*Sz units) of a
    surrogate_kondo_chain, from sector solves on copies of fc.

    With N spins of 1/2 the bath parity of a sector is (2Sz - N) mod 2. At
    B = 0 every multiplet has a member at 2Sz = 0 or 1, so sectors (0, 1)
    hold the whole spectrum, and 2Sz = 2 picks out the states with S >= 1
    (the triplets of a dimer). The default sectors are (0, 1) at B = 0,
    and every 2Sz from -(N+1) to N+1 at B != 0, where the lowest member of
    a multiplet sits at the largest |Sz| of one sign. Returns "sectors" (2Sz -> array of
    energies), "bath_parity" (2Sz -> "even" or "odd"), "E_even" and "E_odd"
    (lowest energy with 2Sz even and odd, the total fermion parity),
    "E_bath_even", "E_bath_odd", "E0" and "Sz_ground" (the sector of E0;
    within 1e-10 max(1, |E0|) of it, the smallest |2Sz|, as in
    ysr.subgap_energies). observables (name ->
    operator) are measured in the lowest state of each sector, as
    "vevs"[name][2Sz]. Raises at Delta = 0, where an Sz sector does not fix
    the particle number.

    On DMRG every sector ground state is solved from a start with one
    electron on each spin site (_warm_sector_solve), and a state that ends
    with another charge there raises RuntimeError. Excited states
    (nstates > 1) start from random states, and an energy more than Upen/2
    above the lowest one of its sector, the signature of a wrong impurity
    charge, raises too; on itensor_version="python" they are refused."""
    import copy
    if getattr(fc, "Delta_sc", 0.0) == 0.0:
        raise ValueError("sector_energies needs a gap, Delta > 0")
    N = len(fc.impurities)
    if sectors is None:
        if getattr(fc, "B_imp", 0.0) == 0.0:
            sectors = (0, 1)
        else:
            sectors = range(-(N + 1), N + 2)
    sectors = sorted(int(s) for s in sectors)
    ens, vevs = {}, {}
    for sz in sectors:
        c = copy.deepcopy(fc)
        c.set_conserved_sector(Sz=sz)
        dmrg = c.get_mode(mode=mode) == "DMRG"
        if dmrg and nstates > 1 and c.itensor_version == "python":
            raise NotImplementedError("nstates > 1 on itensor_version="
                "'python': its excited-state solve keeps the charge of its "
                "random start on the impurity sites, which H conserves "
                "exactly, and stops at Upen above the true state; use "
                "itensor_version=3 or mode='ED'")
        if dmrg:
            e0 = _warm_sector_solve(c)
            q = _impurity_charge_error(c, mode)
            if q > 1e-6:
                raise RuntimeError("sector 2Sz=%d: the DMRG state has "
                    "sum_i <(N_i-1)^2> = %.2e on the spin sites; raise "
                    "nsweeps or maxm" % (sz, q))
        else:
            e0 = c.gs_energy(mode=mode)
        if nstates == 1:
            es = np.array([e0])
        else:
            es = np.sort(np.real(c.get_excited(n=nstates, mode=mode)))
            if dmrg and np.max(es) > es[0] + 0.5*c.Upen:
                raise RuntimeError("sector 2Sz=%d: an excited state is "
                    "%.3g above the lowest one, more than Upen/2, so it "
                    "holds a wrong charge on a spin site" %
                    (sz, np.max(es) - es[0]))
        ens[sz] = es
        if observables:
            for name, op in observables.items():
                vevs.setdefault(name, {})[sz] = np.real(c.vev(op, mode=mode))
    par = {s: ("even" if (s - N) % 2 == 0 else "odd") for s in sectors}
    def lowest(cond):
        v = [ens[s][0] for s in sectors if cond(s)]
        return min(v) if v else None
    out = dict(sectors=ens, bath_parity=par)
    out["E_even"] = lowest(lambda s: s % 2 == 0)
    out["E_odd"] = lowest(lambda s: s % 2 == 1)
    out["E_bath_even"] = lowest(lambda s: par[s] == "even")
    out["E_bath_odd"] = lowest(lambda s: par[s] == "odd")
    out["E0"] = min(e[0] for e in ens.values())
    tie = 1e-10*max(1.0, abs(out["E0"])) # a degenerate multiplet
    out["Sz_ground"] = min((s for s in sectors
                            if ens[s][0] <= out["E0"] + tie),
                           key=lambda s: (abs(s), s))
    if observables:
        out["vevs"] = vevs
    return out


def get_stm_spectrum(fc, es, spin=0, parity=False, **kwargs):
    """dI/dV(eV) with the tip over spin number spin of a
    surrogate_kondo_chain: ysr.get_stm_spectrum with the composite
    tunneling operator of that spin and the first site of its bath.

    With parity=True the correlator starts from the ground state solved
    inside its parity sector, paritystate.parity_ground_state(fc) (the
    sectors of sector_energies, warm-started and charge-checked), promoted
    to dense indices; fc itself is not solved then. Use it near a
    degeneracy of the ground state with the other parity, where the plain
    solve can stop at a mixture of the two and the mixture shows up in
    dI/dV. With mode="ED" the keyword changes nothing, since the ED ground
    state has a definite parity. On a dimer with J12 = t_b = 0.01 (L=3,
    J=0.671, Delta=0.1, maxm=64, nsweeps=20) the parity route is off ED by
    2.5e-6 of the peak on itensor_version="python" and 4.5e-9 on v3; the
    plain route on "python" is 5.8e-7 off there.

    With the default submode="TASK" the Krylov size nkry (default 100) has
    to grow with the number of spins. At N = 4 (L = 3, maxm 128) nkry=100
    put the four YSR peaks on the sector energies; at N = 8 it gave six
    peaks for eight states, some at |eV| differing by up to 0.0125 Delta
    between the two bias signs, while nkry=300 gave all eight, at the same
    |eV| on both signs and at both spins measured (nothing above 300 was
    tried). A
    peak that sits at different |eV| on the two sides is the sign of a
    Krylov space that is too small."""
    from .. import ysr
    if parity:
        from .paritystate import parity_ground_state
        fc = parity_ground_state(fc, mode=kwargs.get("mode", "DMRG"))[0]
    return ysr.get_stm_spectrum(fc, es, imp=fc.impurities[spin],
                                site=fc.bath_sites[spin][0], **kwargs)
