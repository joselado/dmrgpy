"""A magnetic impurity on a superconducting substrate: the Kondo impurity
coupled to the Wilson chain of a BCS band, its subgap (Yu-Shiba-Rusinov)
states, and the STM spectrum dI/dV.

The substrate is a BCS band of half-width D with s-wave gap Delta,

    H_band = sum_k eps_k c^dag_k c_k
             + Delta sum_k (c^dag_{k up} c^dag_{-k dn} + h.c.),

the model of Zitko, Bodensiek and Pruschke, arXiv:1008.1734. The pairing
is diagonal in energy and the same for both spins, and the Wilson map is
one real orthogonal matrix O applied to both spins, so
sum_q a^dag_{q up} a^dag_{q dn} = sum_n f^dag_{n up} f^dag_{n dn} since
O O^T = 1. The Wilson chain of the BCS band is the normal Wilson chain of
wilsonchain.wilson_chain plus the uniform onsite pairing on every site
(tests/test_ysr.py checks it on the Bogoliubov-de Gennes Green function).

The mesh. On the default mean-energy mesh the chain has the low-energy
density of states rho/A_Lambda (lambda_factor), so the continuum J is
J*A_Lambda on the chain. kondo_chain(mesh="zp") takes the
representative-energy mesh of Zitko and Pruschke (arXiv:0810.4765, see
wilsonchain.py), whose density of states is rho at low energy and at the
band edge, and there the continuum J goes on the chain bare
(tests/test_ysr_mesh.py checks both on the classical YSR energy and on
the z-averaged density of states of the BdG chain).

The impurity is a spin 1/2, J S.s_0 with s_0 = f_0^dag (sigma/2) f_0,
held as a native spinful site with one electron on it: it has no hopping,
so its charge is conserved and a penalty Upen (N_imp - 1)^2 removes the
empty and doubly occupied states exactly. Total Sz is conserved too, and
since the impurity carries one electron, Sz = 0 (in 2*Sz units) holds the
states of even total parity (the screened singlet, whose bath holds the
bound quasiparticle and so an odd number of electrons) and Sz = 1 the odd
ones (the free doublet over the BCS vacuum, an even bath). So the energy
of the subgap state is a difference of two sector ground states, with no
excited-state solve.

A spin S > 1/2 with uniaxial anisotropy, J S.s_0 + D_z (S^z)^2 + B S^z,
is held as 2S such sites (sites 0 to 2S-1, the bath from site 2S on), each
with one electron, bound by a ferromagnetic Hund coupling
-J_H sum_{a<b} s_a.s_b = -(J_H/2) (S_imp^2 - 3S/2), S_imp = sum_a s_a.
Every physical term is a function of the vector S_imp, so it commutes
with S_imp^2 and the spin-S multiplet is an exactly decoupled block: J_H
only has to put it lowest, which it does once J_H S, the gap to the
multiplet S-1, exceeds the physical terms (tests/test_ysr_spin.py checks
<S_imp^2> = S(S+1) to roundoff on ED and DMRG). The impurity carries 2S
electrons, so the parity of the bath is (2 Sz - 2S) mod 2 in 2*Sz units:
the free spin over the BCS vacuum (bath even) has 2 Sz = 2S mod 2 and the
screened state, spin S-1/2 in one channel, the other parity. With D_z the
lowest state of a parity can sit at any |Sz| up to S+1/2, so
subgap_energies scans the sectors. This is the model of arXiv:1008.1734,
whose exchange is J S_imp.s with the same J for every S. No term changes
the charge of an impurity site, so a two-site solve keeps the impurity
charges of its starting state, and on itensor_version=3 a sector solve
starts from a short solve with a weak impurity-bath hopping
(_seed_impurity_charge); every solved state is checked to hold one
electron per impurity site.

The STM spectrum is the composite-operator form of
docs/kondo_higher_order_plan.md: at lowest order in the tip coupling the
tip electron enters through

    X_sigma^dag = sum_s' f^dag_{site,s'} [U delta + sigma.S]_{s' sigma}/2

(f at the substrate site under the tip, S the impurity spin, U the ratio
of direct to exchange tunneling), and with a normal tip at T=0

    dI/dV(V) ~ A(eV),  A(w) = sum_sigma [A^e_sigma(w) + A^h_sigma(-w)],
    A^e_sigma(w) = sum_n |<n|X_sigma^dag|GS>|^2 delta(w - (E_n - E_GS)),

and A^h the same with X_sigma. A superconducting tip is a convolution of
A(w) with the tip BCS density of states (tip_convolution).
"""

import copy
import fractions

import numpy as np

from . import fermionchain
from . import wilsonchain


def lambda_factor(Lambda):
    """Wilson's A_Lambda = (1+1/Lambda) ln(Lambda) / (2 (1-1/Lambda)).

    The logarithmic mesh has a smaller low-energy density of states than
    the continuum, by this factor: a coupling v on the chain acts on the
    impurity as v/A_Lambda in the continuum (the classical YSR energy of
    tests/test_ysr.py checks it to 5e-4 of Delta). So J*A_Lambda on the
    chain stands for J in the continuum. This is a property of the
    mean-energy mesh (mesh="mean" of wilsonchain.star_levels) only: on the
    representative-energy mesh (mesh="zp") the low-energy density of
    states is rho itself, and the continuum J goes on the chain bare
    (kondo_chain stores the factor that applies as fc.dos_factor)."""
    L = float(Lambda)
    return 0.5*(1.0 + 1.0/L)/(1.0 - 1.0/L)*np.log(L)


def two_spin(S):
    """2S as an integer, from "1/2", "3/2", 1, 1.5, ... (S >= 1/2)."""
    s2 = 2*fractions.Fraction(str(S))
    if s2.denominator != 1 or s2 < 1:
        raise ValueError("the impurity spin must be 1/2, 1, 3/2, ..., not %r"
                         % (S,))
    return int(s2)


def impurity_spin(fc, imp=None):
    """(S^x, S^y, S^z) of the impurity, the sum of the spins of the sites
    imp (an int or a list, default fc.impurity, else site 0)."""
    if imp is None:
        imp = getattr(fc, "impurity", 0)
    imps = list(np.atleast_1d(imp))
    out = []
    for Sa in (fc.Sx, fc.Sy, fc.Sz):
        s = Sa[int(imps[0])]
        for a in imps[1:]:
            s = s + Sa[int(a)]
        out.append(s)
    return tuple(out)


def kondo_chain(n, J=0.2, Delta=0.0, Lambda=2.0, D=1.0, B=0.0, z=0.0,
                Upen=None, S="1/2", Dz=0.0, JH=None, mesh="mean", **kwargs):
    """Kondo impurity of spin S, coupled by J S.s_0 to the first site of
    the Wilson chain of a BCS band, as a Spinful_Fermionic_Chain_Native.
    The impurity takes sites 0 to 2S-1 and the chain sites 2S to 2S+n-1.

    H = J S.s_0 + Dz (S^z)^2 + B S^z - JH sum_{a<b} s_a.s_b
        + Upen sum_a (N_a - 1)^2
        + sum_n eps_n N_n + sum_n t_n (f^dag_n f_{n+1} + h.c.)
        + Delta sum_n (f^dag_{n up} f^dag_{n dn} + h.c.)

    with S = sum_a s_a the impurity spin. J is the coupling on the chain,
    taken as given on both meshes. mesh is the Wilson mesh of the band,
    "mean" (default) or "zp" (wilsonchain.star_levels). On "mean" the
    chain has the low-energy density of states rho/A_Lambda, so pass
    J*lambda_factor(Lambda) to stand for J in the continuum; on "zp" it
    has rho, and the continuum J goes on the chain bare. The chain stores
    that factor as fc.dos_factor (lambda_factor(Lambda) or 1), the mesh as
    fc.mesh, the coupling as fc.J_chain and the last hopping |t_{n-1}| as
    fc.t_last (D for n=1); ysrtk.zaverage reads both. The band
    is flat on [-D, D], so rho = 1/(2D) and rho J is J/(2D). B is the field
    on the impurity only and Dz the uniaxial anisotropy, both in energy
    units (Dz (S^z)^2 is a constant for S=1/2). S is "1/2" (default), "1",
    "3/2", "2", ... or a number. Upen (default 2(D+|J|+|B|+2S|Dz|)) only
    has to push the empty and the doubly occupied impurity sites above the
    singly occupied ones; JH (default the same, for S > 1/2) puts the spin-S
    multiplet of the 2S sites lowest, since JH S is its gap to S-1 and no
    physical term mixes the two. Extra keywords go to the chain. The chain
    stores fc.impurity = [0..2S-1], fc.spin = S (a Fraction), fc.channels =
    [[2S..2S+n-1]] and fc.Delta_sc."""
    if n < 1:
        raise ValueError("the Wilson chain needs at least one site")
    ni = two_spin(S)
    eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D), z=z,
                                      mesh=mesh)
    if Upen is None:
        Upen = 2.0*(D + abs(J) + abs(B))
        if Dz != 0.0:
            Upen = Upen + 2.0*ni*abs(Dz)
    if JH is None:
        JH = Upen
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n+ni, **kwargs)
    Sx, Sy, Sz = impurity_spin(fc, list(range(ni)))
    h = J*(Sx*fc.Sx[ni] + Sy*fc.Sy[ni] + Sz*fc.Sz[ni])
    h = h + B*Sz
    if Dz != 0.0:
        h = h + Dz*Sz*Sz
    for a in range(ni): # Hund coupling between the impurity sites
        for b in range(a+1, ni):
            h = h - JH*(fc.Sx[a]*fc.Sx[b] + fc.Sy[a]*fc.Sy[b]
                        + fc.Sz[a]*fc.Sz[b])
    for a in range(ni): # one electron on every impurity site
        na = fc.Nup[a] + fc.Ndn[a] - 1
        h = h + Upen*na*na
    for k in range(n):
        f = k + ni
        h = h + eps[k]*(fc.Nup[f] + fc.Ndn[f])
        if Delta != 0.0:
            h = h + Delta*(fc.Cdagup[f]*fc.Cdagdn[f] + fc.Cdn[f]*fc.Cup[f])
    for k in range(n-1):
        a, b = k + ni, k + ni + 1
        h = h + t[k]*(fc.Cdagup[a]*fc.Cup[b] + fc.Cdagdn[a]*fc.Cdn[b]
                      + fc.Cdagup[b]*fc.Cup[a] + fc.Cdagdn[b]*fc.Cdn[a])
    fc.set_hamiltonian(h)
    fc.impurity = list(range(ni))
    fc.spin = fractions.Fraction(ni, 2)
    fc.channels = [list(range(ni, n+ni))]
    fc.Delta_sc = Delta
    fc.mesh = mesh
    fc.dos_factor = lambda_factor(Lambda) if mesh == "mean" else 1.0
    fc.J_chain = J
    fc.t_last = abs(t[-1]) if n > 1 else D
    fc.Dz = Dz
    fc.B_imp = B
    fc.qp_gap = bath_quasiparticle_gap(eps, t, Delta)
    fc.reaches_gap = n > 1 and Delta != 0 and t[-1] < 0.1*abs(Delta)
    return fc


def bath_quasiparticle_gap(eps, t, Delta):
    """Lowest Bogoliubov energy of the Wilson chain alone. It is Delta
    once the chain reaches below Delta, and larger on a short chain, whose
    lowest level sits at about D Lambda^-n/2: then the continuum of dI/dV
    starts above Delta and the chain is too short for the gap."""
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
    return float(np.min(np.abs(np.linalg.eigvalsh(H))))


def subgap_energies(fc, mode="DMRG", nrestart=1, sectors=None,
                    observables=None):
    """Lowest energies of the two parity sectors of a kondo_chain, from
    Sz-sector solves on copies of fc (the chain itself is left untouched).

    For S=1/2 these are E_even (Sz=0, the screened singlet) and E_odd (Sz=1
    in 2*Sz units, the free doublet). For any S, E_even and E_odd are the
    lowest energies with 2Sz even and odd, the total fermion parity, and
    E_free and E_screened the lowest energies with the bath even (2Sz = 2S
    mod 2, the free spin S over the BCS vacuum) and odd (the spin screened
    by a bound quasiparticle, to S-1/2 in one channel). So for S=1/2
    E_screened = E_even, and for integer S E_screened = E_odd.
    E_ysr = E_free - E_screened > 0 means a screened ground state, and
    |E_ysr| is the energy of the YSR excitation, the subgap pole of dI/dV
    (for S=1/2, E_ysr = E_odd - E_even).

    sectors are the 2Sz values solved. By default at Dz = 0 and B = 0 the
    lowest |2Sz| of each parity, since there the spin is isotropic and every
    multiplet has a member there; otherwise every 2Sz of each parity with
    |2Sz| <= 2S+1, both signs if B != 0, since with anisotropy the lowest
    state of a parity can sit at any |Sz| up to S+1/2. "sectors" in the
    result maps each 2Sz to its energy, and "Sz_ground" is the 2Sz of the
    lowest one, the smallest |2Sz| within 1e-10 of it on a degenerate
    multiplet (for S=1/2 "ground_state" is "singlet" or "doublet", for
    S > 1/2 "screened" or "free"). observables, a dict of name ->
    operator, are measured in the ground state of each sector and returned
    as "vevs"[name][2Sz].

    DMRG energies are upper bounds, so with nrestart > 1 each sector is
    solved nrestart times from fresh random starts and the lowest energy
    is kept; "spread" is the largest difference between the restarts of
    one sector, a measure of how far a single solve can be trusted. A
    sector solve can stop in a metastable state, for instance the
    singlet plus a free quasiparticle in the odd sector. On a chain that
    reaches well below Delta (fc.reaches_gap, the last hopping below
    Delta/10) the two parities differ by at most one quasiparticle far from
    the impurity, which costs fc.qp_gap = Delta, so a larger |E_ysr|
    flags it, with a warning. On a shorter chain there is no such bound:
    the lowest bath level is then strongly coupled to the impurity, and
    |E_ysr| can exceed Delta in the exact answer. Sectors of one parity are
    split by Dz and B, so their differences are not bounded by Delta."""
    import warnings
    if getattr(fc, "Delta_sc", 0.0) == 0.0:
        raise ValueError("subgap_energies needs a gap, Delta > 0: at "
            "Delta = 0 the particle number is conserved too, so an Sz "
            "sector alone does not fix the state, and a solve stays at the "
            "particle number it started from")
    ni = int(2*getattr(fc, "spin", fractions.Fraction(1, 2)))
    if sectors is None:
        if getattr(fc, "Dz", 0.0) == 0.0 and getattr(fc, "B_imp", 0.0) == 0.0:
            sectors = [0, 1]
        else:
            sectors = list(range(0, ni + 2))
            if getattr(fc, "B_imp", 0.0) != 0.0:
                sectors = [-s for s in sectors[::-1] if s > 0] + sectors
    sectors = sorted(int(s) for s in sectors)
    if not any((s - ni) % 2 == 0 for s in sectors) or \
       not any((s - ni) % 2 == 1 for s in sectors):
        raise ValueError("sectors must contain both parities of 2Sz")
    out, spread, ens, vevs = {}, 0.0, {}, {}
    for sz in sectors:
        es, cs = [], []
        for _ in range(nrestart if mode == "DMRG" else 1):
            c = copy.deepcopy(fc)
            c.set_conserved_sector(Sz=sz)
            if ni > 1 and mode == "DMRG" and c.itensor_version == 3:
                _seed_impurity_charge(c, ni)
            es.append(c.gs_energy(mode=mode))
            if ni > 1:
                _check_impurity_charge(c, ni, mode)
            cs.append(c)
        ens[sz] = min(es)
        spread = max(spread, max(es) - min(es))
        if observables:
            c = cs[int(np.argmin(es))]
            for name, op in observables.items():
                vevs.setdefault(name, {})[sz] = c.vev(op, mode=mode)
    def lowest(par):
        return min(e for s, e in ens.items() if s % 2 == par)
    out["E_even"], out["E_odd"] = lowest(0), lowest(1)
    out["E_free"], out["E_screened"] = lowest(ni % 2), lowest((ni + 1) % 2)
    out["E_ysr"] = out["E_free"] - out["E_screened"]
    out["spread"] = spread
    if ni == 1:
        out["ground_state"] = "singlet" if out["E_ysr"] > 0 else "doublet"
    else:
        out["ground_state"] = "screened" if out["E_ysr"] > 0 else "free"
    out["sectors"] = ens
    e0 = min(ens.values())
    tie = 1e-10*max(1.0, abs(e0)) # a degenerate multiplet: smallest |2Sz|
    out["Sz_ground"] = min((s for s, e in ens.items() if e <= e0 + tie),
                           key=lambda s: (abs(s), s))
    if observables:
        out["vevs"] = vevs
    gap = getattr(fc, "qp_gap", 0.0)
    if getattr(fc, "reaches_gap", False) and abs(out["E_ysr"]) > gap*(1 + 1e-2):
        warnings.warn("subgap_energies: |E_free - E_screened| = %.3e is "
            "larger than the lowest quasiparticle energy of the bath, %.3e, "
            "so one sector solve is not its ground state; raise maxm, "
            "nsweeps or nrestart" % (abs(out["E_ysr"]), gap))
    return out


def _default_site(fc, site):
    if site is not None:
        return site
    ch = getattr(fc, "channels", None)
    return ch[0][0] if ch else 1


def _impurity_hopping(fc, ni):
    """Hopping between every impurity site and the first bath site."""
    b = ni
    T = 0
    for a in range(ni):
        T = T + (fc.Cdagup[a]*fc.Cup[b] + fc.Cdagdn[a]*fc.Cdn[b]
                 + fc.Cdagup[b]*fc.Cup[a] + fc.Cdagdn[b]*fc.Cdn[a])
    return T


def _seed_impurity_charge(c, ni, t=0.1):
    """Starting state for a v3 sector solve of an impurity of 2S = ni
    sites. No term of H changes the charge of an impurity site, not even
    a factor of one across a bond, so a two-site solve keeps the impurity
    charges of its start, and the in-sector start of itensor_version=3 is
    a product state that can leave an impurity site empty (a spin 3/2 then
    lands at S(S+1) = 2 and an energy Upen too high). We solve first with a
    weak hopping t between the impurity and the bath, which connects the
    charges, at a small bond dimension, and start the real solve there."""
    h = c.hamiltonian
    maxm, nsweeps = c.maxm, c.nsweeps
    c.set_hamiltonian(h + t*_impurity_hopping(c, ni))
    c.maxm, c.nsweeps = min(maxm, 16), min(nsweeps, 4)
    c.gs_energy()
    wf = c.get_gs()
    c.set_hamiltonian(h)
    c.maxm, c.nsweeps = maxm, nsweeps
    c.set_initial_wf_guess(wf)


def _check_impurity_charge(c, ni, mode):
    """Every impurity site must hold one electron in a solved state."""
    w = sum(np.real(c.vev((c.Nup[a] + c.Ndn[a] - 1)*(c.Nup[a] + c.Ndn[a] - 1),
                          mode=mode)) for a in range(ni))
    if w > 1e-6:
        raise RuntimeError("subgap_energies: the solved state has weight "
            "%.3e outside one electron per impurity site, so it is not in "
            "the spin-S multiplet; raise nsweeps or Upen" % w)


def tunneling_operators(fc, imp=None, site=None, U=0.0):
    """The two tunneling vertices X_up^dag and X_dn^dag of the module
    docstring, for the impurity spin on imp (a site or a list of sites
    whose spins add up, default fc.impurity, else site 0) and the substrate
    site site (default the first site of the Wilson chain, else site 1)."""
    S = impurity_spin(fc, imp)
    site = _default_site(fc, site)
    Sp, Sm = S[0] + 1j*S[1], S[0] - 1j*S[1]
    cu, cd = fc.Cdagup[site], fc.Cdagdn[site]
    Xup = 0.5*(cu*(U + S[2]) + cd*Sp)
    Xdn = 0.5*(cu*Sm + cd*(U - S[2]))
    return [Xup, Xdn]


def get_stm_spectrum(fc, es, delta=None, imp=None, site=None, U=0.0,
                     submode="TASK", mode="DMRG", **kwargs):
    """dI/dV(eV) at the energies es, with a normal tip at T=0, in units of
    2 pi T0^2 per unit of tip density of states: A(w) of the module
    docstring, each pole broadened into a Lorentzian of half-width delta.

    The electron part A^e(w) and the hole part A^h(-w) are two dynamical
    correlators, (X, X^dag) and (X^dag, X), summed over the tip spin.
    submode and the remaining keywords go to get_dynamical_correlator. The
    ground state of fc is the one the correlator uses: on the free-doublet
    side at B=0 any member of the doublet gives the same spin-summed
    spectrum, since the sum over sigma is a scalar under spin rotation.
    With Dz != 0 the same holds within a degenerate pair +-Sz, since the
    spin-summed correlator conserves Sz and the two members differ in Sz.
    imp and site are as in tunneling_operators: for S > 1/2 the impurity
    spin is the total spin of the 2S impurity sites."""
    es = np.asarray(es, dtype=float)
    if delta is None:
        delta = 0.02*(np.max(es) - np.min(es))
    Jc = getattr(fc, "J_chain", None) # a number, or one per spin (surrogate)
    if (submode == "TASK" and mode != "ED" and getattr(fc, "mode", None) != "ED"
            and Jc is not None and np.all(np.asarray(Jc) == 0)
            and kwargs.get("pad_bonds") != "full"):
        import warnings
        warnings.warn("get_stm_spectrum: J = 0 on the chain, so the ground "
            "state is a product of the impurity and the bath, and the "
            "one-site tangent space of TASK misses the part of X|GS> that "
            "flips the impurity spin and adds an electron to the bath at "
            "once: the spectrum is incomplete, down to zero. Pass "
            "pad_bonds=\"full\" or use mode=\"ED\" on a short chain, or the "
            "BdG density of states of zaverage.bdg_chain_dos, for J = 0",
            RuntimeWarning)
    out = np.zeros(len(es))
    for X in tunneling_operators(fc, imp=imp, site=site, U=U):
        Xd = X.get_dagger() # the operator that removes the electron
        _, ce = fc.get_dynamical_correlator(mode=mode, submode=submode,
                    name=(Xd, X), es=es, delta=delta, **kwargs)
        _, ch = fc.get_dynamical_correlator(mode=mode, submode=submode,
                    name=(X, Xd), es=-es[::-1], delta=delta, **kwargs)
        out = out + np.real(ce) + np.real(ch)[::-1]
    return es, out


def bcs_dos(e, Delta, gamma=0.0):
    """Dynes density of states Re[(e+i gamma)/sqrt((e+i gamma)^2-Delta^2)],
    normalized to 1 far from the gap. At gamma=0 it is the BCS density
    |e|/sqrt(e^2-Delta^2) for |e| > |Delta| and exactly zero inside the
    gap; the edge |e| = |Delta| itself, an integrable singularity, is given
    zero too."""
    if gamma == 0.0:
        e = np.asarray(e, dtype=float)
        a = np.abs(e)
        g = abs(float(Delta))
        if g == 0.0:
            return np.ones(np.shape(e))[()]
        out = np.zeros(np.shape(e))
        o = a > g
        out[o] = a[o]/np.sqrt((a[o] - g)*(a[o] + g))
        return out[()] if out.ndim == 0 else out
    w = np.asarray(e, dtype=complex) + 1j*gamma
    r = w/np.sqrt(w*w - Delta*Delta)
    return np.abs(np.real(r))


def tip_convolution(es, A, Delta_tip, gamma=0.0, Vs=None):
    """dI/dV(V) for a superconducting tip with gap Delta_tip, from the
    sample spectrum A(w) on the uniform grid es, at T=0:

        I(V) = integral dw rho_tip(w - eV) A(w) [theta(eV - w) - theta(-w)],

    and dI/dV by finite differences. Vs (default es inside the grid by
    Delta_tip) are the biases. With Delta_tip = 0 this gives back A(eV).
    The grid must reach past max|V| + Delta_tip for the integral to be
    complete; A outside es is taken as zero."""
    es = np.asarray(es, dtype=float)
    A = np.asarray(A, dtype=float)
    dw = es[1] - es[0]
    if not np.allclose(np.diff(es), dw):
        raise ValueError("tip_convolution needs a uniform grid es")
    if Vs is None:
        Vs = es[np.abs(es) < np.max(np.abs(es)) - Delta_tip]
    Vs = np.asarray(Vs, dtype=float)
    def current(V):
        rho = bcs_dos(es - V, Delta_tip, gamma) if Delta_tip > 0 else 1.0
        occ = (es < V).astype(float) - (es < 0.0).astype(float)
        return np.sum(rho*A*occ)*dw
    h = dw
    I1 = np.array([current(V + h/2) for V in Vs])
    I0 = np.array([current(V - h/2) for V in Vs])
    return Vs, (I1 - I0)/h
