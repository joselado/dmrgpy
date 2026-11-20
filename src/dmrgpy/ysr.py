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

A spin chain on a superconducting island (island_chain) is the same
encoding with a different bath: ns spins 1/2 on singly occupied native
sites, coupled by Js S_i.S_{i+1}, the last one Kondo coupled by JK to the
first of ni island sites with the pairing Delta put by hand, and the last
island site hopping by t_c onto a normal Wilson chain. The metal is
gapless, so the two parities are split by a finite-size level of the
Wilson chain, not by a subgap state, and island_energies reports that
splitting, E_parity = E_odd - E_even, against its value at JK=0; the
sector solves, the charge seeding and the parity state are the ones
above, with ns penalized sites.

A spin chain on a normal metal (spin_chain_on_metal) is that encoding
with no island: ns spins 1/2 and the Wilson chain of the metal, Kondo
coupled to any one spin of the chain, the one under the tip. The Wilson
chain comes after all the spins also when that spin is in the middle, so
the Kondo term is then long-range (docs/development.md 4.8i says why this
order and not one that alternates the spins with the Wilson sites). With
a list of spins every spin of the list has a Wilson chain of its own,
independent baths, placed by order= (_chain_on_metal), and
anderson_spin_chain makes the spins of the list Anderson levels that
hop onto their metal instead of Kondo spins, with the level's own
d^dag as the tip vertex.

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
    imp (an int or a list). The default is fc.tip_spin when the chain has
    one (island_chain: the end spin, the one in contact with the island),
    else fc.impurity, else site 0."""
    if imp is None:
        imp = getattr(fc, "tip_spin", None)
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
    fc.mesh, the coupling as fc.J_chain, the last hopping |t_{n-1}| as
    fc.t_last (D for n=1) and the chain as fc.wilson_chain = (eps, t)
    with fc.wilson_params (Lambda, D, z, mesh); ysrtk.zaverage reads the
    first two. The band
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
    fc.wilson_chain = (np.array(eps), np.array(t))
    fc.wilson_params = dict(Lambda=Lambda, D=D, z=z, mesh=mesh)
    fc.Dz = Dz
    fc.B_imp = B
    fc.qp_gap = bath_quasiparticle_gap(eps, t, Delta)
    fc.reaches_gap = n > 1 and Delta != 0 and t[-1] < 0.1*abs(Delta)
    return fc


def bath_quasiparticle_gap(eps, t, Delta):
    """Lowest Bogoliubov energy of the Wilson chain alone. It is Delta
    once the chain reaches below Delta, and larger on a short chain, whose
    lowest level sits at about D Lambda^-n/2: then the continuum of dI/dV
    starts above Delta and the chain is too short for the gap. Delta is a
    number (uniform pairing) or one value per site (island_chain, whose
    pairing sits on the island only)."""
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    P = np.diag(np.broadcast_to(np.asarray(Delta, dtype=float), (n,)))
    H = np.block([[h, P], [P, -h]])
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
    if on_metal_leads(fc):
        raise ValueError("subgap_energies is for one gapped bath with the "
            "impurity on the first sites (kondo_chain): on a chain of "
            "spin_chain_on_metal the spins sit between the Wilson chains, "
            "and with several leads each lead keeps its own fermion "
            "parity, so a 2Sz sector alone does not fix the state")
    if hasattr(fc, "island"):
        raise ValueError("subgap_energies is for a gapped bath "
            "(kondo_chain): on an island_chain the bath is a metal, the two "
            "parities are split by a finite-size level of the Wilson chain "
            "and E_free - E_screened is not a subgap energy; use "
            "island_energies")
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
    out = {}
    ens, spread, vevs = _sector_ground_energies(fc, sectors, mode, nrestart,
                                                observables, ni, ni > 1)
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


def _sector_ground_energies(fc, sectors, mode, nrestart, observables, ni,
                            check_charge, seed_always=False):
    """Lowest energy of every 2Sz sector in sectors, each from nrestart
    solves (one on ED) on copies of fc, the loop of subgap_energies and
    island_energies. ni is the number of singly occupied penalized sites
    (sites 0..ni-1): on itensor_version=3 a DMRG sector solve with ni > 1
    first takes the warm start of _seed_impurity_charge, and with
    check_charge every solved state is checked to hold one electron per
    penalized site. With seed_always (island_chain) every DMRG sector
    solve, on every backend and for any ni, takes the warm start: on
    "python" a 16-site island chain (ns=4) stopped with a spin site empty or
    doubly occupied without it. Returns (energies, spread, vevs)."""
    spread, ens, vevs = 0.0, {}, {}
    for sz in sectors:
        es, cs = [], []
        for _ in range(nrestart if mode == "DMRG" else 1):
            c = copy.deepcopy(fc)
            c.set_conserved_sector(Sz=sz)
            if mode == "DMRG" and (seed_always or
                    (ni > 1 and c.itensor_version == 3)):
                _seed_impurity_charge(c, ni)
            es.append(c.gs_energy(mode=mode))
            if check_charge:
                _check_impurity_charge(c, ni, mode)
            cs.append(c)
        ens[sz] = min(es)
        spread = max(spread, max(es) - min(es))
        if observables:
            c = cs[int(np.argmin(es))]
            for name, op in observables.items():
                vevs.setdefault(name, {})[sz] = c.vev(op, mode=mode)
    return ens, spread, vevs


def _hop(fc, a, b):
    """The four-term spinful hopping between sites a and b."""
    return (fc.Cdagup[a]*fc.Cup[b] + fc.Cdagdn[a]*fc.Cdn[b]
            + fc.Cdagup[b]*fc.Cup[a] + fc.Cdagdn[b]*fc.Cdn[a])


def island_chain(ns, ni, n, Js=0.5, JK=0.5, Delta=0.1, t_island=0.5,
                 eps_island=0.0, t_c=0.5, Lambda=2.0, D=1.0, B=0.0, z=0.0,
                 Upen=None, mesh="mean", wilson=None, **kwargs):
    """A Heisenberg chain of ns spins 1/2 whose last spin is Kondo coupled
    to a superconducting island of ni sites, which leaks into a normal
    metal, the Wilson chain of n sites. One linear chain of ns+ni+n native
    spinful sites (Spinful_Fermionic_Chain_Native): spins 0..ns-1, island
    ns..ns+ni-1, Wilson chain ns+ni..ns+ni+n-1. With a0 = ns, a1 = ns+ni-1
    and w0 = ns+ni,

        H = Js sum_i S_i.S_{i+1} + JK S_{ns-1}.s_{a0} + B sum_i S^z_i
            + Upen sum_i (N_i - 1)^2                      (spin sites)
            + eps_island sum_k N_k + t_island sum_k (f^dag_k f_{k+1} + h.c.)
            + Delta sum_k (f^dag_{k up} f^dag_{k dn} + h.c.)   (island)
            + t_c (f^dag_{a1} f_{w0} + h.c.)
            + sum_m eps_m N_m + sum_m t_m (f^dag_m f_{m+1} + h.c.)  (metal)

    with s_a0 = f^dag_a0 (sigma/2) f_a0 and every spin a native site with
    one electron, as in kondo_chain. The pairing is put by hand on the
    island only, so the island gap is imposed, not self-consistent, and
    the metal softens it: there is no true subgap state on this chain,
    every "subgap" feature is a resonance whose width is the leak into the
    metal, and on a finite Wilson chain it is a few discrete poles.

    Js > 0 and JK > 0 are antiferromagnetic. B is a field on the spins
    only. The Wilson chain is wilsonchain.wilson_chain(n, Lambda, (-D, D),
    z, mesh), or wilson=(eps, t) by hand (n values and n-1). Couplings are
    taken as given. JK couples to the island, not to the Wilson chain, so
    no lambda_factor applies to it. t_c is a hopping into the metal, with
    hybridization Gamma = pi rho t_c^2 (rho = 1/(2D) on the flat band):
    on mesh="mean" the chain has the low-energy density of states
    rho/A_Lambda, so t_c*sqrt(lambda_factor(Lambda)) on the chain stands
    for t_c in the continuum (one square root); on mesh="zp" t_c goes bare.
    Checked on the occupation of a level at e_d << Gamma hopping onto a
    60-site chain, averaged over 8 twists: t_c sqrt(A_Lambda) gives 0.99
    of the continuum Gamma at Lambda = 2 and 4, bare t_c gives 1/A_Lambda
    of it (0.95 and 0.86).
    Upen (default 2(max(D, |t_island|, |eps_island|, |Delta|, |t_c|) + |JK|
    + |Js| + |B|)) only has to exceed the magnetic binding of a spin site,
    since no term moves charge on or off it; energy differences inside the
    one-electron block do not depend on it.

    Checked (tests/test_ysr_island.py, on d22): ns=1, ni=1 with the island
    taken as the first Wilson site gives kondo_chain at Delta=0 to
    roundoff; at JK=0 the energy is the Heisenberg chain plus the
    Bogoliubov-de Gennes ground state of island plus metal, and the island
    density of states is bdg_poles(*fc.bdg_params); sector energies and
    the STM spectrum agree between ED, "python" and 3. On 16 sites (ns=4,
    ni=2, n=10, maxm=64) "python" and 3 agree to 2e-7 in the sector
    energies; this is agreement between the backends, not accuracy: at
    JK=0 both are about 2e-6 above the exact energy at maxm=64
    (truncation). The STM spectrum at delta=0.03 to 1.5e-3 of its peak
    at nkry=300; the default nkry=100 of TASK was 1.5e-2 of the peak away
    from nkry=300 there, on both backends, so pass a larger nkry for
    sharp lines on a long chain.

    Stored on the chain: fc.impurity = [0..ns-1] (the penalized sites),
    fc.spin = ns/2, which is the kondo_chain convention "2*fc.spin
    penalized sites" read by the sector solvers, NOT the total spin of the
    chain; fc.tip_spin = ns-1, the default spin of the tunneling vertex;
    fc.island, fc.wilson and fc.channels = [island, wilson] (the vertex
    takes its electron at channels[0][0] = a0); fc.Delta_sc, fc.Dz = 0,
    fc.B_imp, fc.J_chain = JK, fc.Js, fc.t_island, fc.eps_island, fc.t_c,
    fc.Upen, fc.mesh, fc.dos_factor, fc.t_last; fc.bdg_params = (eps,
    t, Delta) of island plus metal, the single-particle problem at JK=0;
    fc.qp_gap, its lowest Bogoliubov level (the cost of one quasiparticle,
    not Delta; at eps=0 it alternates with the parity of ni+n, and is
    much smaller when ni+n is odd, since then a zero mode is gapped only
    by its weight on the island); fc.island_gap, the same for
    the isolated island (at least |Delta|); fc.reaches_gap = False."""
    ns, ni, n = int(ns), int(ni), int(n)
    if ns < 1 or ni < 1 or n < 1:
        raise ValueError("island_chain needs ns >= 1 spins, ni >= 1 island "
                         "sites and n >= 1 Wilson sites")
    if wilson is None:
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D),
                                          z=z, mesh=mesh)
    else:
        eps, t = (np.asarray(x, dtype=float) for x in wilson)
        if len(eps) != n or len(t) != n - 1:
            raise ValueError("wilson=(eps, t) needs %d energies and %d "
                             "hoppings" % (n, n - 1))
    if Upen is None:
        scale = max(D, abs(t_island), abs(eps_island), abs(Delta), abs(t_c))
        Upen = 2.0*(scale + abs(JK) + abs(Js) + abs(B))
    a0, a1, w0 = ns, ns + ni - 1, ns + ni
    fc = fermionchain.Spinful_Fermionic_Chain_Native(ns + ni + n, **kwargs)
    h = 0
    for i in range(ns - 1): # Heisenberg chain
        h = h + Js*(fc.Sx[i]*fc.Sx[i+1] + fc.Sy[i]*fc.Sy[i+1]
                    + fc.Sz[i]*fc.Sz[i+1])
    e = ns - 1 # the end spin, Kondo coupled to the first island site
    h = h + JK*(fc.Sx[e]*fc.Sx[a0] + fc.Sy[e]*fc.Sy[a0] + fc.Sz[e]*fc.Sz[a0])
    for i in range(ns):
        if B != 0.0:
            h = h + B*fc.Sz[i]
        ni_ = fc.Nup[i] + fc.Ndn[i] - 1 # one electron on every spin site
        h = h + Upen*ni_*ni_
    for k in range(a0, a1 + 1): # the island
        if eps_island != 0.0:
            h = h + eps_island*(fc.Nup[k] + fc.Ndn[k])
        if Delta != 0.0:
            h = h + Delta*(fc.Cdagup[k]*fc.Cdagdn[k] + fc.Cdn[k]*fc.Cup[k])
        if k < a1:
            h = h + t_island*_hop(fc, k, k + 1)
    h = h + t_c*_hop(fc, a1, w0) # island to metal
    for m in range(n): # the Wilson chain, no pairing
        h = h + eps[m]*(fc.Nup[w0+m] + fc.Ndn[w0+m])
        if m < n - 1:
            h = h + t[m]*_hop(fc, w0 + m, w0 + m + 1)
    fc.set_hamiltonian(h)
    fc.impurity = list(range(ns))
    fc.spin = fractions.Fraction(ns, 2)
    fc.tip_spin = ns - 1
    fc.island = list(range(a0, a1 + 1))
    fc.wilson = list(range(w0, w0 + n))
    fc.channels = [fc.island, fc.wilson]
    fc.Delta_sc = Delta
    fc.Dz = 0.0
    fc.B_imp = B
    fc.J_chain = JK
    fc.Js, fc.t_island, fc.eps_island, fc.t_c = Js, t_island, eps_island, t_c
    fc.Upen = Upen
    fc.mesh = mesh
    fc.dos_factor = lambda_factor(Lambda) if mesh == "mean" else 1.0
    fc.t_last = abs(t[-1]) if n > 1 else abs(t_c)
    eps_f = np.concatenate([np.full(ni, float(eps_island)), eps])
    t_f = np.concatenate([np.full(ni - 1, float(t_island)), [float(t_c)], t])
    Delta_f = np.concatenate([np.full(ni, float(Delta)), np.zeros(n)])
    fc.bdg_params = (eps_f, t_f, Delta_f)
    fc.qp_gap = bath_quasiparticle_gap(eps_f, t_f, Delta_f)
    fc.island_gap = bath_quasiparticle_gap(np.full(ni, float(eps_island)),
                                           np.full(ni - 1, float(t_island)),
                                           Delta)
    fc.reaches_gap = False
    return fc


def island_energies(fc, mode="DMRG", nrestart=1, sectors=None,
                    observables=None):
    """Lowest energies of the 2Sz sectors of an island_chain, from sector
    solves on copies of fc (the chain itself is left untouched), and the
    parity splitting E_parity = E_odd - E_even.

    Every spin site holds one electron, so 2Sz fixes the total fermion
    parity (2Sz mod 2, "even" or "odd") and the parity of island plus
    metal, (2Sz - ns) mod 2. The metal is gapless, so the two parities
    differ by the cost of one quasiparticle, the lowest Bogoliubov level of
    island plus metal, not Delta. At eps=0 it alternates with the parity
    of ni+n, and is much smaller when ni+n is odd (for ns=4, ni=2 at
    Lambda=2: 1.0e-1 at n=6, 7.5e-3 at n=7, 5.3e-2 at n=8). E_parity is a
    finite-size number that goes to zero with n, and fc.qp_gap is its value
    at JK=0 (+qp_gap for ns even, -qp_gap for ns odd). The ratio
    E_parity/qp_gap is the finite-size spectrum where the fixed point is
    read (end spin screened or free); "Sz_ground" and "parity_ground" are
    the robust statements. This is not a YSR energy.

    sectors (default [0, 1] at B=0, where the chain is isotropic and every
    multiplet has a member at the lowest |2Sz| of its parity; every 2Sz
    with |2Sz| <= ns+1, both signs, at B != 0) are the 2Sz solved. Every
    DMRG sector solve, on "python" and on 3, starts from the warm stage of
    _seed_impurity_charge (a weak hopping of every spin onto the first
    island site), and every solved state is checked to hold one electron
    per spin site. nrestart,
    "spread" and observables ("vevs"[name][2Sz]) are those of
    subgap_energies. Returns "sectors" (2Sz -> energy), "E_even", "E_odd",
    "E_parity", "Sz_ground" (the smallest |2Sz| within 1e-10 of the
    lowest), "parity_ground", "spread", "qp_gap" and "vevs". On "python"
    the restarts are identical (a seeded start and a seeded generator), so
    "spread" is 0 by construction there; it measures the start dependence
    on itensor_version=3 only. A warning
    says when the ground state is even and E_parity > 1.05 qp_gap: adding
    the lowest quasiparticle to the even state costs about qp_gap, so a
    larger splitting usually means the odd sector solve is not its ground
    state. The bound is heuristic once JK > 0: it holds only up to JK
    times the weight of that mode on a0, small for a mode at the far end
    of the Wilson chain and not small for n <= 2. An odd ground state can
    sit below the even one by more than qp_gap (the Kondo binding of the
    end spin, ratio -1.06 on ED at ns=2, ni=1, n=2, JK=1), so no warning
    is given on that side. At small Delta the island pairing is the only
    term that changes the particle number, and a sector solve can stall at
    the particle number it started from, which can get the ground-state
    parity wrong: compare with ED on a short chain. On "python" the
    restarts are identical, so raise nsweeps (at Delta=0.02 on 7 sites
    20 sweeps gave the wrong parity and 40 the right one); on 3,
    nrestart > 1 measures the start dependence."""
    import warnings
    if not hasattr(fc, "island"):
        raise ValueError("island_energies needs an island_chain; for a "
                         "kondo_chain use subgap_energies")
    if fc.Delta_sc == 0.0:
        raise ValueError("island_energies needs Delta != 0: at Delta = 0 "
            "the particle number is conserved too, and an Sz sector alone "
            "does not fix the state")
    ns = len(fc.impurity)
    if sectors is None:
        if fc.B_imp == 0.0:
            sectors = [0, 1]
        else:
            sectors = list(range(-(ns + 1), ns + 2))
    sectors = sorted(set(int(s) for s in sectors))
    if not any(s % 2 == 0 for s in sectors) or \
       not any(s % 2 == 1 for s in sectors):
        raise ValueError("sectors must contain both parities of 2Sz")
    ens, spread, vevs = _sector_ground_energies(fc, sectors, mode, nrestart,
                                                observables, ns, True, True)
    out = dict(sectors=ens, spread=spread, qp_gap=fc.qp_gap)
    out["E_even"] = min(e for s, e in ens.items() if s % 2 == 0)
    out["E_odd"] = min(e for s, e in ens.items() if s % 2 == 1)
    out["E_parity"] = out["E_odd"] - out["E_even"]
    e0 = min(ens.values())
    tie = 1e-10*max(1.0, abs(e0))
    out["Sz_ground"] = min((s for s, e in ens.items() if e <= e0 + tie),
                           key=lambda s: (abs(s), s))
    out["parity_ground"] = "even" if out["Sz_ground"] % 2 == 0 else "odd"
    if observables:
        out["vevs"] = vevs
    if out["E_parity"] > 1.05*fc.qp_gap:
        warnings.warn("island_energies: E_odd - E_even = %.3e is larger "
            "than the cost of one quasiparticle, %.3e, so the odd sector "
            "solve is likely not its ground state; raise maxm, nsweeps or "
            "nrestart" % (out["E_parity"], fc.qp_gap), RuntimeWarning)
    return out


def island_spin_correlations(fc, mode="DMRG"):
    """The (ns+1) x (ns+1) matrix <S_a.S_b> over the ns spins and the
    first island site a0 (last row and column), in the ground state of fc.
    <S_i.S_{i+1}> is -3/4 for an isolated singlet pair, and the last
    column measures how much the island electrons screen each spin; the
    diagonal is 3/4 on a spin site and 3/4 of the single occupancy on a0.
    Every entry is a scalar under spin rotation, so any member of a
    degenerate multiplet gives the same matrix."""
    sites = list(fc.impurity) + [fc.island[0]]
    m = len(sites)
    C = np.zeros((m, m))
    for a in range(m):
        for b in range(a, m):
            i, j = sites[a], sites[b]
            op = fc.Sx[i]*fc.Sx[j] + fc.Sy[i]*fc.Sy[j] + fc.Sz[i]*fc.Sz[j]
            C[a, b] = C[b, a] = np.real(fc.vev(op, mode=mode))
    return C


def spin_chain_on_metal(ns, n, site=None, Js=0.5, JK=0.5, Lambda=2.0, D=1.0,
                        B=0.0, z=0.0, Upen=None, mesh="mean", wilson=None,
                        order=None, Delta=0.0, **kwargs):
    """A Heisenberg chain of ns spins 1/2 on a normal metal: the spin
    site (any of 0..ns-1, default the end spin ns-1) is Kondo coupled by
    JK to the first site of a Wilson chain of n sites, the metal under
    that spin. With site a list of spins, every spin of the list has a
    Wilson chain of its own, see "Several spins on the metal" below.
    ns+n native spinful sites (Spinful_Fermionic_Chain_Native):
    spins 0..ns-1, Wilson chain ns..ns+n-1, for every site. With w0 = ns,

        H = Js sum_i S_i.S_{i+1} + JK S_site.s_{w0} + B sum_i S^z_i
            + Upen sum_i (N_i - 1)^2                      (spin sites)
            + sum_m eps_m N_m + sum_m t_m (f^dag_m f_{m+1} + h.c.)  (metal)

    with s_w0 = f^dag_w0 (sigma/2) f_w0 and every spin a native site with
    one electron, as in kondo_chain and island_chain. With site = ns-1
    this is one linear chain (island_chain with one island site and no
    pairing, to roundoff). With a spin in the middle the geometry is a T
    junction, and the Wilson chain still comes after all the spins, so
    the Kondo term jumps over the spins site+1..ns-1. Keep this order: a
    cut between two of those spins carries the entropy of the spin bond
    plus that of the spins with the metal (0.5 to 0.8 on 7 and 15 spins
    at Js=0.5, JK=1), while every bond of the Wilson chain carries 1.2
    to 1.6. Alternating the spins with the Wilson sites, which makes
    every term of range 2, puts a Wilson bond on every spin bond: on 15
    spins and 21 Wilson sites with the metal on spin 7 it was 2.7e-3
    above the best energy at maxm=64, against 5e-5 for this order
    (docs/development.md 4.8i).

    Js > 0 and JK > 0 are antiferromagnetic. B is a field on the spins
    only. The Wilson chain is wilsonchain.wilson_chain(n, Lambda, (-D, D),
    z, mesh), or wilson=(eps, t) by hand (n values and n-1). JK is the
    coupling on the chain, as J in kondo_chain: on mesh="mean" pass
    JK*lambda_factor(Lambda) to stand for JK in the continuum, on
    mesh="zp" pass it bare (fc.dos_factor is the factor that applies).
    Upen (default 2(D + |JK| + |Js| + |B|)) only has to exceed the
    magnetic binding of a spin site, since no term moves charge on or off
    it. Extra keywords go to the chain.

    The total spin of the ground state depends on site. On 3 and on 7
    spins with n odd it is a singlet with the metal on an even spin (0,
    2, ...), which the metal screens, and a triplet with the metal on an
    odd spin, which points against the total spin of the chain, so its
    exchange with the metal acts as ferromagnetic (S_tot^2 = 0 and 2,
    measured). With ns+n odd it is a doublet. The spin-summed dI/dV is a
    scalar under spin rotation, so every state of a multiplet gives the
    same one.

    The metal has no gap, so the sector tools (subgap_energies,
    island_energies, parity_ground_state) raise on this chain, and a
    plain solve gives the ground state. get_stm_spectrum puts the tip on
    the spin site and takes the electron at the first Wilson site.

    Stored on the chain: fc.impurity = [0..ns-1] (the penalized sites),
    fc.spin = ns/2 (the convention of island_chain, not the total spin),
    fc.tip_spin = site, fc.wilson and fc.channels = [wilson], fc.J_chain
    = JK, fc.Js, fc.B_imp, fc.Dz = 0, fc.Delta_sc = 0, fc.Upen, fc.mesh,
    fc.dos_factor, fc.t_last (the last hopping, the energy resolution of
    the chain), fc.wilson_chain = (eps, t), fc.wilson_params (Lambda, D,
    z and mesh, None for wilson= by hand), fc.qp_gap (the lowest level
    of the Wilson chain alone, zero for n odd on a symmetric band) and
    fc.reaches_gap = False.

    Several spins on the metal. With site a list (or tuple, or array) of
    spins, each of them has its own Wilson chain of n sites, the same
    chain for all, and JK is a number or one value per spin of the list.
    The baths are independent: this is the model of spins far apart on
    the scale of the Fermi wavelength, where the metal under one spin
    does not see the other (on a common metal the number of channels is
    the rank of the coupling matrix instead, Eickhoff and Anders,
    arXiv:2008.09034). order sets where the chains go: "sides" (the
    default for several spins) puts the chains of the left half of the
    chain on its left end, reversed, and those of the right half on its
    right end, "comb" puts every chain right after its spin, and
    "linear" all the chains after all the spins, one after the other;
    see _chain_on_metal for the measurement behind the default.
    The positions are then no longer the labels: fc.spin_sites[i] is the
    site of spin i. The chain stores fc.contacts (the spins of the list,
    in order), fc.impurities (their sites) and fc.bath_sites (their
    Wilson chains, first site first), the names of
    ysrtk.surrogate.surrogate_kondo_chain, fc.channels = fc.bath_sites,
    fc.tip_spin = the site of the first spin of the list, fc.J_chain
    (one value per chain when there are several), fc.order and
    fc.anderson = []. get_stm_spectrum(fc, es, imp=fc.impurities[k])
    puts the tip on spin fc.contacts[k] and takes the electron at the
    first site of its own chain; a tip on a spin with no chain raises.
    submode="NRG" takes one chain only. A list of one spin with
    order="linear" is the chain above, up to the order of the terms.

    Leads that differ, and superconducting leads. With site a list, D,
    Lambda and Delta are a number or one value per spin of the list, as
    JK, so every Wilson chain can be a different metal: wilson_chain(n,
    Lambda[k], (-D[k], D[k]), z, mesh) plus the uniform onsite pairing
    Delta[k] sum_m (f^dag_{m up} f^dag_{m dn} + h.c.) of kondo_chain on
    every site of chain k. site=[0, ns-1] (order "sides") is a spin chain
    between two semi-infinite leads, the left one reversed on the left end
    of the MPS and the right one on the right end. The leads couple to the
    spins by exchange only, which moves no charge in or out of a lead, so
    a phase of one pairing can be gauged away (f -> exp(i phi/2) f on that
    lead): the energy does not depend on it (Delta[k] -> -Delta[k] leaves
    the spectrum unchanged), and there is no Josephson coupling through
    the spins. Delta is real, and a complex one raises. Each lead keeps
    its own fermion parity. With JK on the mean mesh, pass JK[k]*lambda_factor(
    Lambda[k]). The chain stores fc.lead_Delta, fc.lead_D, fc.lead_Lambda
    (one value per spin of the list), fc.lead_chains (the (eps, t) of
    every chain) and fc.lead_qp_gap (the lowest Bogoliubov level of every
    chain alone); fc.Delta_sc is the largest |Delta[k]|, fc.qp_gap the
    smallest gap, fc.t_last the smallest last hopping, fc.dos_factor that
    of the leads when they share Lambda (None otherwise), and
    fc.wilson_chain and fc.wilson_params those of the leads when all are
    the same metal (None otherwise). With several leads the sector tools
    (subgap_energies, ysrtk.paritystate.parity_ground_state,
    ysrtk.surrogate.sector_energies) raise (see on_metal_leads), and with
    a nonzero Delta the STM spectrum takes get_stm_spectrum(mode="ED",
    submode="CVM") only: the DMRG routes raise, also through
    zaverage.stm_poles and surrogate.get_stm_spectrum (see
    refuse_dmrg_on_paired_leads for the measurement). One spin with one
    lead, spin_chain_on_metal(1, n, site=[0]), is kondo_chain term for
    term and is not refused. D, Lambda and Delta per lead need
    site as a list, and wilson= by hand takes D and Lambda as numbers."""
    if np.ndim(site) == 1:
        return _chain_on_metal(ns, n, site, order, Js=Js, B=B, Lambda=Lambda,
                               D=D, z=z, Upen=Upen, mesh=mesh, wilson=wilson,
                               JK=JK, kwargs=kwargs, Delta=Delta)
    if np.any(np.asarray(Delta) != 0.0) or np.ndim(D) > 0 or \
            np.ndim(Lambda) > 0:
        raise ValueError("Delta, and D or Lambda per lead, need site as a "
                         "list: pass site=[%s] for one spin on the metal, "
                         "or site=[0, ns-1] for a chain between two leads"
                         % ("ns-1" if site is None else site))
    if order is not None:
        raise ValueError("order= places several Wilson chains; with one "
                         "spin on the metal (site an int) the chain always "
                         "comes after the spins. Pass site as a list for "
                         "the general layout")
    ns, n = int(ns), int(n)
    if ns < 1 or n < 1:
        raise ValueError("spin_chain_on_metal needs ns >= 1 spins and "
                         "n >= 1 Wilson sites")
    if site is None:
        site = ns - 1
    if site != int(site) or not 0 <= site < ns:
        raise ValueError("site must be one of the spins 0..%d, got %r"
                         % (ns - 1, site))
    site = int(site)
    if wilson is None:
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D),
                                          z=z, mesh=mesh)
    else:
        eps, t = (np.asarray(x, dtype=float) for x in wilson)
        if len(eps) != n or len(t) != n - 1:
            raise ValueError("wilson=(eps, t) needs %d energies and %d "
                             "hoppings" % (n, n - 1))
    if Upen is None:
        Upen = 2.0*(D + abs(JK) + abs(Js) + abs(B))
    w0 = ns
    fc = fermionchain.Spinful_Fermionic_Chain_Native(ns + n, **kwargs)
    h = 0
    for i in range(ns - 1): # Heisenberg chain
        h = h + Js*(fc.Sx[i]*fc.Sx[i+1] + fc.Sy[i]*fc.Sy[i+1]
                    + fc.Sz[i]*fc.Sz[i+1])
    # the Kondo term, long-range when site is not the end spin
    h = h + JK*(fc.Sx[site]*fc.Sx[w0] + fc.Sy[site]*fc.Sy[w0]
                + fc.Sz[site]*fc.Sz[w0])
    for i in range(ns):
        if B != 0.0:
            h = h + B*fc.Sz[i]
        ni_ = fc.Nup[i] + fc.Ndn[i] - 1 # one electron on every spin site
        h = h + Upen*ni_*ni_
    for m in range(n): # the Wilson chain
        if eps[m] != 0.0:
            h = h + eps[m]*(fc.Nup[w0+m] + fc.Ndn[w0+m])
        if m < n - 1:
            h = h + t[m]*_hop(fc, w0 + m, w0 + m + 1)
    fc.set_hamiltonian(h)
    fc.impurity = list(range(ns))
    fc.spin = fractions.Fraction(ns, 2)
    fc.tip_spin = site
    fc.wilson = list(range(w0, w0 + n))
    fc.channels = [fc.wilson]
    fc.Delta_sc = 0.0
    fc.Dz = 0.0
    fc.B_imp = B
    fc.J_chain = JK
    fc.Js = Js
    fc.Upen = Upen
    fc.mesh = mesh
    fc.dos_factor = lambda_factor(Lambda) if mesh == "mean" else 1.0
    fc.t_last = abs(t[-1]) if n > 1 else D
    fc.wilson_chain = (np.array(eps), np.array(t))
    fc.wilson_params = (None if wilson is not None
                        else dict(Lambda=Lambda, D=D, z=z, mesh=mesh))
    fc.qp_gap = bath_quasiparticle_gap(eps, t, 0.0)
    fc.reaches_gap = False
    return fc


def anderson_spin_chain(ns, n, site=None, Js=0.5, ed=None, U=1.0, tD=0.2,
                        Lambda=2.0, D=1.0, B=0.0, z=0.0, Upen=None,
                        mesh="mean", wilson=None, order=None, **kwargs):
    """A Heisenberg chain of ns spins 1/2 on a normal metal in which the
    spin site (default the end spin ns-1) is an Anderson level: an
    orbital d with energy ed, repulsion U and hopping tD onto the first
    site of a Wilson chain of n sites, the metal under it,

        H = Js sum_i S_i.S_{i+1} + B sum_i S^z_i + Upen sum_{i != d} (N_i - 1)^2
            + ed N_d + U n_{d up} n_{d dn} + tD sum_s (d^dag_s f_{0 s} + h.c.)
            + sum_m eps_m N_m + sum_m t_m (f^dag_m f_{m+1} + h.c.),

    with S_d = d^dag (sigma/2) d in the Heisenberg term, so the exchange
    with the neighbours acts on the singly occupied level and vanishes on
    the empty and the doubly occupied one. The other spins are charge
    frozen, as in spin_chain_on_metal. ed defaults to -U/2, the particle
    hole symmetric point. The hybridization is Gamma = pi rho tD^2 with
    rho = 1/(2D); on mesh="mean" pass tD*sqrt(lambda_factor(Lambda)) to
    stand for tD in the continuum (the factor of island_chain's t_c), on
    mesh="zp" pass it bare. At U >> Gamma the level is a spin 1/2 Kondo
    coupled to the metal with J = 2 tD^2 (1/|ed| + 1/(ed+U)), 8 tD^2/U at
    ed = -U/2 (Schrieffer and Wolff), and the chain becomes
    spin_chain_on_metal with that JK.

    site may be a list: every spin of it is then an Anderson level with
    its own Wilson chain, and ed, U and tD are numbers or one value per
    spin of the list, and so are D and Lambda, the bandwidth and the mesh
    of each chain, as in spin_chain_on_metal (no pairing here); order is
    as in spin_chain_on_metal. The STM vertex
    of a tip on an Anderson site is d^dag_s (plus q f^dag_{0 s} for the
    direct path into the metal, get_stm_spectrum(q=)), so dI/dV is the
    spectral function of the level summed over spin, the ideal STM limit
    of Jacob, arXiv:1805.10624. On mesh="mean" pass
    q*sqrt(lambda_factor(Lambda)) to stand for q in the continuum, as for
    tD (a bare q acts as q/sqrt(A_Lambda), 7% low at Lambda = 4), and the
    spectral function of the level needs no fc.dos_factor. Stored on the chain as in
    spin_chain_on_metal with several spins (fc.spin_sites, fc.contacts,
    fc.impurities, fc.bath_sites, fc.channels, fc.tip_spin), plus
    fc.anderson (the sites of the levels), fc.ed, fc.U_d and fc.tD (one
    value per level), fc.penalized (the charge-frozen spins) and
    fc.J_chain = None. Extra keywords go to the chain."""
    return _chain_on_metal(ns, n, ns - 1 if site is None else site, order,
                           Js=Js, B=B, Lambda=Lambda, D=D, z=z, Upen=Upen,
                           mesh=mesh, wilson=wilson, kwargs=kwargs,
                           anderson=dict(ed=ed, U=U, tD=tD))


METAL_ORDERS = ("comb", "linear", "sides")
MISSING_RAISE = 0.1 # TASK on several Kondo baths raises above this missing
# share of the sum rule (Fable's choice of 2026-11-18, between the 0.02 of
# expand_bonds=8, expand_rounds=3 and the 0.62 to 0.95 of no expansion)
METAL_ORDER = "sides" # the default order for several Wilson chains


def _metal_contacts(ns, site):
    """The spins on the metal, as a list of labels in the order given."""
    out = []
    for s in np.atleast_1d(site):
        if s != int(s) or not 0 <= s < ns:
            raise ValueError("site must be one of the spins 0..%d, got %r"
                             % (ns - 1, s))
        out.append(int(s))
    if len(out) == 0:
        raise ValueError("site is an empty list: at least one spin must be "
                         "on the metal")
    if len(set(out)) != len(out):
        raise ValueError("site lists a spin twice: %r" % (list(site),))
    return out


def _metal_layout(ns, n, contacts, order):
    """The site of every spin and the sites of the Wilson chain of every
    contact (first Wilson site first), for the orders of METAL_ORDERS."""
    if order == "linear":
        spins = list(range(ns))
        baths = [list(range(ns + k*n, ns + (k + 1)*n))
                 for k in range(len(contacts))]
        return spins, baths
    if order == "comb":
        spins, own, p = [], {}, 0
        for i in range(ns):
            spins.append(p)
            p += 1
            if i in contacts:
                own[i] = list(range(p, p + n))
                p += n
        return spins, [own[c] for c in contacts]
    if order == "sides":
        left = sorted(c for c in contacts if 2*c < ns - 1)
        right = sorted((c for c in contacts if 2*c > ns - 1), reverse=True)
        for c in contacts: # the middle spin of an odd chain
            if 2*c == ns - 1:
                (left if len(left) < len(right) else right).append(c)
        nl, own = len(left), {}
        for k, c in enumerate(left): # depth k from the spins, reversed
            own[c] = [(nl - k)*n - 1 - m for m in range(n)]
        spins = [nl*n + i for i in range(ns)]
        r0 = nl*n + ns
        for k, c in enumerate(right):
            own[c] = [r0 + k*n + m for m in range(n)]
        return spins, [own[c] for c in contacts]
    raise ValueError("order must be one of %s, got %r" % (METAL_ORDERS, order))


def _per_contact(name, x, nc):
    x = np.asarray(x, dtype=float)
    if x.ndim > 1 or (x.ndim == 1 and len(x) != nc):
        raise ValueError("%s must be a number or one value per spin on the "
                         "metal (%d), got %r" % (name, nc, x))
    return np.broadcast_to(x, (nc,)).copy()


def _chain_on_metal(ns, n, site, order, Js, B, Lambda, D, z, Upen, mesh,
                    wilson, kwargs, JK=None, anderson=None, Delta=0.0):
    """The builder behind spin_chain_on_metal with a list of spins and
    behind anderson_spin_chain: ns spins 1/2, the spins of site on the
    metal, each with its own Wilson chain, Kondo coupled by JK or, with
    anderson = dict(ed, U, tD), an Anderson level hopping onto it.

    Why "sides" for several chains. A cut of the MPS carries the sum of
    the entropies of what it separates. With the chains one after the
    other ("linear"), a cut inside one chain also separates every spin
    from the chains still to come, so its entropy grows with the number
    of chains. "sides" splits them between the two ends, and "comb" puts
    every chain right after its spin, where a cut inside a chain carries
    that chain and one bond of the spin chain, but every Heisenberg bond
    that jumps a chain costs bond dimension (2 to 4 times for one chain,
    docs/development.md 4.8i). Measured on 5 spins with 13 Wilson sites
    per chain (docs/new_methods/chains_on_metal/README.md section 2):
    at JK = 0.3 sides is the lowest at maxm 64 and 128 with 2, 3 and 4
    chains, and its largest bond entropy the smallest (1.79 against 2.35
    for comb with 4 chains); comb leads only at JK = 1 with 3 and 4
    chains. At JK = 0, where the exact energy is the Heisenberg chain
    plus free chains, comb stopped 0.177 (3 chains) and 0.36 (4 chains)
    above it at maxm 128, since two-site DMRG does not build a spin bond
    that jumps an uncoupled chain, so an explicit order="comb" with a
    contact of zero coupling warns."""
    ns, n = int(ns), int(n)
    if ns < 1 or n < 1:
        raise ValueError("a chain on a metal needs ns >= 1 spins and n >= 1 "
                         "Wilson sites")
    contacts = _metal_contacts(ns, site)
    nc = len(contacts)
    if order is None:
        order = METAL_ORDER if nc > 1 else "linear"
    if order not in METAL_ORDERS:
        raise ValueError("order must be one of %s, got %r"
                         % (METAL_ORDERS, order))
    Dl = _per_contact("D", D, nc)
    Ll = _per_contact("Lambda", Lambda, nc)
    if np.iscomplexobj(Delta):
        raise ValueError("Delta must be real: with exchange coupling only, "
            "no term moves charge in or out of a lead, so the phase of its "
            "pairing is a gauge")
    Pl = _per_contact("Delta", Delta, nc)
    if wilson is None:
        made = {} # one chain per distinct metal
        for k in range(nc):
            key = (Dl[k], Ll[k])
            if key not in made:
                made[key] = wilsonchain.wilson_chain(n, Lambda=Ll[k],
                                                     band=(-Dl[k], Dl[k]),
                                                     z=z, mesh=mesh)
        chains = [made[(Dl[k], Ll[k])] for k in range(nc)]
    else:
        if np.ndim(D) > 0 or np.ndim(Lambda) > 0:
            raise ValueError("wilson=(eps, t) by hand is one chain for "
                             "every lead, so D and Lambda must be numbers")
        eps, t = (np.asarray(x, dtype=float) for x in wilson)
        if len(eps) != n or len(t) != n - 1:
            raise ValueError("wilson=(eps, t) needs %d energies and %d "
                             "hoppings" % (n, n - 1))
        chains = [(eps, t)]*nc
    spins, baths = _metal_layout(ns, n, contacts, order)
    level = [spins[c] for c in contacts]
    if anderson is None:
        Jc = _per_contact("JK", JK, nc)
        andersons, scale = [], np.max(np.abs(Jc))
    else:
        Ud = _per_contact("U", anderson["U"], nc)
        ed = (-Ud/2.0 if anderson["ed"] is None
              else _per_contact("ed", anderson["ed"], nc))
        tD = _per_contact("tD", anderson["tD"], nc)
        andersons = list(level)
        scale = max(np.max(np.abs(ed)), np.max(np.abs(Ud)),
                    np.max(np.abs(tD)))
    weak = Jc if anderson is None else tD
    if order == "comb" and nc > 1 and np.any(weak == 0.0):
        import warnings
        warnings.warn("order=\"comb\" with a spin of zero coupling to its "
                      "metal: a Heisenberg bond that jumps an uncoupled "
                      "Wilson chain is not built by two-site DMRG, and the "
                      "solve can stop far above the ground state (0.36 on "
                      "4 chains at JK = 0). comb was measured to lead only "
                      "at strong coupling; use order=\"sides\" (the "
                      "default)", RuntimeWarning, stacklevel=3)
    penalized = [p for p in spins if p not in andersons]
    if Upen is None:
        Upen = 2.0*(np.max(Dl) + scale + abs(Js) + abs(B))
    fc = fermionchain.Spinful_Fermionic_Chain_Native(ns + n*nc, **kwargs)
    def dot(a, b):
        return fc.Sx[a]*fc.Sx[b] + fc.Sy[a]*fc.Sy[b] + fc.Sz[a]*fc.Sz[b]
    h = 0
    for i in range(ns - 1): # Heisenberg chain
        h = h + Js*dot(spins[i], spins[i+1])
    for k in range(nc): # the coupling of each spin to its own metal
        p, w0 = level[k], baths[k][0]
        if anderson is None:
            h = h + Jc[k]*dot(p, w0)
        else:
            if ed[k] != 0.0:
                h = h + ed[k]*(fc.Nup[p] + fc.Ndn[p])
            if Ud[k] != 0.0:
                h = h + Ud[k]*fc.Nup[p]*fc.Ndn[p]
            h = h + tD[k]*_hop(fc, p, w0)
    for p in spins:
        if B != 0.0:
            h = h + B*fc.Sz[p]
    for p in penalized: # one electron on every charge-frozen spin
        q = fc.Nup[p] + fc.Ndn[p] - 1
        h = h + Upen*q*q
    for k, w in enumerate(baths): # the Wilson chains
        eps, t = chains[k]
        for m in range(n):
            if eps[m] != 0.0:
                h = h + eps[m]*(fc.Nup[w[m]] + fc.Ndn[w[m]])
            if Pl[k] != 0.0: # onsite pairing, as in kondo_chain
                h = h + Pl[k]*(fc.Cdagup[w[m]]*fc.Cdagdn[w[m]]
                               + fc.Cdn[w[m]]*fc.Cup[w[m]])
            if m < n - 1:
                h = h + t[m]*_hop(fc, w[m], w[m+1])
    fc.set_hamiltonian(h)
    fc.impurity = list(spins)
    fc.spin_sites = list(spins)
    fc.penalized = penalized
    fc.contacts = list(contacts)
    fc.impurities = list(level)
    fc.bath_sites = [list(w) for w in baths]
    fc.channels = fc.bath_sites
    if nc == 1:
        fc.wilson = list(baths[0])
    fc.order = order
    fc.anderson = andersons
    fc.spin = fractions.Fraction(ns, 2)
    fc.tip_spin = level[0]
    fc.lead_Delta, fc.lead_D, fc.lead_Lambda = Pl, Dl, Ll
    fc.lead_chains = [(np.array(e), np.array(x)) for e, x in chains]
    fc.lead_qp_gap = np.array([bath_quasiparticle_gap(e, x, Pl[k])
                               for k, (e, x) in enumerate(chains)])
    fc.Delta_sc = float(np.max(np.abs(Pl)))
    fc.Dz = 0.0
    fc.B_imp = B
    if anderson is None:
        fc.J_chain = float(Jc[0]) if nc == 1 else Jc
    else:
        fc.J_chain = None
        fc.ed, fc.U_d, fc.tD = ed, Ud, tD
    fc.Js = Js
    fc.Upen = Upen
    fc.mesh = mesh
    same = bool(np.all(Dl == Dl[0]) and np.all(Ll == Ll[0]))
    one_mesh = mesh != "mean" or bool(np.all(Ll == Ll[0]))
    fc.dos_factor = None # the factor depends on Lambda and the mesh only
    if one_mesh:
        fc.dos_factor = lambda_factor(Ll[0]) if mesh == "mean" else 1.0
    if same: # one metal for every lead: the attributes of one chain
        eps, t = chains[0]
        fc.wilson_chain = (np.array(eps), np.array(t))
        fc.wilson_params = (None if wilson is not None else
                            dict(Lambda=float(Ll[0]), D=float(Dl[0]), z=z,
                                 mesh=mesh))
    else: # no single value; the per-lead ones are fc.lead_*
        fc.wilson_chain = fc.wilson_params = None
    fc.t_last = float(min(abs(x[-1]) if n > 1 else Dl[k]
                          for k, (e, x) in enumerate(chains)))
    fc.qp_gap = float(np.min(fc.lead_qp_gap))
    fc.reaches_gap = bool(np.any(Pl != 0.0) and all(
        n > 1 and abs(x[-1]) < 0.1*abs(Pl[k])
        for k, (e, x) in enumerate(chains) if Pl[k] != 0.0))
    return fc


def on_metal_leads(fc):
    """True for a chain of _chain_on_metal (spin_chain_on_metal with a
    list of spins, anderson_spin_chain) other than one Kondo spin with one
    Wilson chain after it, spin_chain_on_metal(1, n, site=[0]), which is
    kondo_chain term for term. The sector tools of one impurity seed and
    check the charge of the sites 0..2S-1 and solve a total 2Sz sector,
    which takes the wrong sites when a chain sits before or between the
    spins and does not fix the parity of each lead when there are
    several. One spin with one lead is measured against ED with pairing
    (subgap_energies, parity_ground_state and the TASK STM spectrum on v3
    and "python"); several spins on one paired lead are not, and are
    refused as the int site path refuses Delta."""
    if getattr(fc, "order", None) not in METAL_ORDERS:
        return False
    return not (len(getattr(fc, "bath_sites", None) or ()) == 1
                and not getattr(fc, "anderson", None)
                and list(getattr(fc, "spin_sites", None) or ()) == [0])


def refuse_dmrg_on_paired_leads(fc, mode="DMRG", where="get_stm_spectrum"):
    """Raise when the STM spectrum of a chain of spin_chain_on_metal with
    pairing (on_metal_leads and a nonzero fc.Delta_sc) is asked of a DMRG
    route. Measured on 2 spins with 3 Wilson sites per lead at
    saturation: sides with expand_bonds 3e-4 to 6e-4 of the peak off ED
    (1e-8 at Delta = 0), linear tip-first plain TASK 1.6e-3 on "python"
    and 0.89 on v3, with the ground state exact and of definite parity.
    Called before any solve by ysr.get_stm_spectrum, zaverage.stm_poles
    and surrogate.get_stm_spectrum."""
    if (on_metal_leads(fc) and getattr(fc, "Delta_sc", 0.0) != 0.0
            and fc.get_mode(mode=mode) != "ED"):
        raise ValueError("%s: the DMRG routes are not verified on a chain "
            "of spin_chain_on_metal with pairing other than one spin on one "
            "lead (TASK "
            "was measured 0.89 of the peak off ED on v3 on 8 sites, with "
            "an exact ground state); use ysr.get_stm_spectrum with "
            "mode=\"ED\", submode=\"CVM\" on a short chain" % where)


def _tip_position(fc, imp):
    """The site of the tip spin, as impurity_spin resolves imp, or None
    when it is several sites (a spin S > 1/2 of kondo_chain)."""
    if imp is None:
        imp = getattr(fc, "tip_spin", None)
    if imp is None:
        imp = getattr(fc, "impurity", 0)
    imps = list(np.atleast_1d(imp))
    return int(imps[0]) if len(imps) == 1 else None


def _default_site(fc, site, imp=None):
    """The substrate site under the tip: site when given; with several
    Wilson chains (fc.bath_sites), the first site of the chain of the tip
    spin, which must have one; else the first site of fc.channels[0]."""
    if site is not None:
        return site
    baths = getattr(fc, "bath_sites", None)
    if baths is not None and len(baths) > 1:
        p = _tip_position(fc, imp)
        if p is None:
            raise ValueError("a tip on several sites at once has no single "
                             "chain under it on a chain with several Wilson "
                             "chains; pass site= by hand")
        imps = [int(a) for a in fc.impurities]
        if p not in imps:
            raise ValueError("the tip spin, site %r, has no Wilson chain of "
                             "its own: with several chains the electron is "
                             "taken at the first site of the chain of the "
                             "tip spin, so the tip must be on one of the "
                             "sites %s (or pass site= by hand)" % (p, imps))
        return baths[imps.index(p)][0]
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
    sites (and for every DMRG sector solve of an island_chain, ni = ns,
    where site ni is the first island site). No term of H changes the
    charge of an impurity site, not even a factor of one across a bond, so
    a two-site solve keeps the impurity charges of its start, and the
    in-sector start of itensor_version=3 is a product state that can
    leave an impurity site empty (a spin 3/2 then
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
    """Every penalized site 0..ni-1 (the impurity sites of kondo_chain,
    the spins of island_chain) must hold one electron in a solved state.
    ni may also be the list of the penalized sites (the spins of
    anderson_spin_chain other than the Anderson levels)."""
    sites = range(ni) if np.ndim(ni) == 0 else ni
    w = sum(np.real(c.vev((c.Nup[a] + c.Ndn[a] - 1)*(c.Nup[a] + c.Ndn[a] - 1),
                          mode=mode)) for a in sites)
    if w > 1e-6:
        raise RuntimeError("the solved state has weight %.3e outside one "
            "electron per penalized spin site (an impurity site of "
            "kondo_chain or a spin of island_chain), so a spin site is "
            "empty or doubly occupied; raise nsweeps or Upen" % w)


def tunneling_operators(fc, imp=None, site=None, U=0.0, q=0.0):
    """The two tunneling vertices X_up^dag and X_dn^dag of the module
    docstring, for the impurity spin on imp (a site or a list of sites
    whose spins add up, default fc.tip_spin, else fc.impurity, else site
    0, see impurity_spin) and the substrate site site (default
    fc.channels[0][0], the first site of the Wilson chain of kondo_chain
    or the first island site of island_chain, else site 1; with several
    Wilson chains, fc.bath_sites, the first site of the chain of the tip
    spin). On an island_chain the electron is taken at that island site
    for every spin, since it is the only fermionic site a spin touches.

    A tip on an Anderson level (a site of fc.anderson, anderson_spin_chain)
    puts the electron on the level itself, X_s^dag = d^dag_s +
    q f^dag_{site s}, with q the direct tunneling into the metal relative
    to the level (0 by default, the ideal STM limit of the level's
    spectral function; on mesh="mean", q*sqrt(lambda_factor(Lambda))
    stands for q in the continuum); U belongs to the exchange vertex of a
    Kondo spin and raises there, as q does on a Kondo spin. On a chain of
    _chain_on_metal imp must be a spin (fc.spin_sites), since the
    positions are not the labels in its orders."""
    p = _tip_position(fc, imp)
    spins = getattr(fc, "spin_sites", None)
    if spins is not None:
        for a in np.atleast_1d(fc.tip_spin if imp is None else imp):
            if int(a) not in spins:
                raise ValueError("imp=%r: site %d is not a spin of this "
                                 "chain, whose spins are on the sites %s "
                                 "(fc.spin_sites)" % (imp, int(a), spins))
    if p is not None and p in getattr(fc, "anderson", ()):
        if U != 0.0:
            raise ValueError("U is the potential scattering of the exchange "
                             "vertex of a Kondo spin; a tip on an Anderson "
                             "level takes q=, the direct tunneling into the "
                             "metal")
        X = [fc.Cdagup[p], fc.Cdagdn[p]]
        if q != 0.0:
            w0 = _default_site(fc, site, imp)
            X = [X[0] + q*fc.Cdagup[w0], X[1] + q*fc.Cdagdn[w0]]
        return X
    if q != 0.0:
        raise ValueError("q= is the direct tunneling into the metal of a tip "
                         "on an Anderson level (anderson_spin_chain); on a "
                         "Kondo spin the direct path is U=")
    S = impurity_spin(fc, imp)
    site = _default_site(fc, site, imp)
    Sp, Sm = S[0] + 1j*S[1], S[0] - 1j*S[1]
    cu, cd = fc.Cdagup[site], fc.Cdagdn[site]
    Xup = 0.5*(cu*(U + S[2]) + cd*Sp)
    Xdn = 0.5*(cu*Sm + cd*(U - S[2]))
    return [Xup, Xdn]


def get_stm_spectrum(fc, es, delta=None, imp=None, site=None, U=0.0, q=0.0,
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
    spin is the total spin of the 2S impurity sites. With several Wilson
    chains (spin_chain_on_metal with a list of spins) imp picks the spin
    under the tip and the electron goes to the first site of its own
    chain. On an Anderson level (anderson_spin_chain) the vertex is
    d^dag_s + q f^dag_{0 s}, and dI/dV at q=0 is the spectral function of
    the level summed over spin.

    submode="NRG" (spin_chain_on_metal with one spin on the metal,
    anderson_spin_chain with one level and the linear order, or a
    spin-1/2 kondo_chain at Delta = 0; an Anderson level takes
    head="dmrg" only) is the numerical renormalization group of nrgmps
    instead,
    on the mesh of fc: a DMRG head and an NRG tail (head="dmrg", the
    default), or the exact start shell of the spins and the first Wilson
    site (head="exact", up to nrgspectral.MAX_START_SPINS spins), with
    the poles broadened on a logarithmic mesh rather than by delta; T > 0
    is allowed there. The tip is on fc.tip_spin and the first Wilson
    site, so imp and site must be None. See nrgmps.get_stm_spectrum."""
    es = np.asarray(es, dtype=float)
    if submode == "NRG":
        if imp is not None or site is not None:
            raise ValueError("get_stm_spectrum: submode=\"NRG\" puts the tip "
                             "on fc.tip_spin and the first Wilson site; imp "
                             "and site must be None")
        if mode != "DMRG":
            raise ValueError("get_stm_spectrum: submode=\"NRG\" takes no "
                             "mode=; use head=\"exact\" for the exact start "
                             "shell")
        if delta is not None:
            raise ValueError("get_stm_spectrum: submode=\"NRG\" broadens on "
                             "a logarithmic mesh, with b=, not delta=")
        if len(getattr(fc, "bath_sites", None) or ()) > 1:
            raise ValueError("get_stm_spectrum: submode=\"NRG\" runs one "
                             "Wilson chain, and this chain has %d (one per "
                             "spin on the metal); use submode=\"TASK\""
                             % len(fc.bath_sites))
        from . import nrgmps
        return nrgmps.get_stm_spectrum(fc, es, U=U, q=q, **kwargs)
    refuse_dmrg_on_paired_leads(fc, mode=mode)
    if delta is None:
        delta = 0.02*(np.max(es) - np.min(es))
    Jc = getattr(fc, "J_chain", None) # a number, or one per spin (surrogate)
    # the Wilson chains of _chain_on_metal (not the surrogate BCS baths,
    # which pairing and t_b do not leave with a charge of their own)
    nbath = len(getattr(fc, "contacts", None) or ())
    on_task = (submode == "TASK" and mode != "ED"
               and getattr(fc, "mode", None) != "ED"
               and kwargs.get("pad_bonds") != "full")
    kondo_baths = nbath > 1 and not getattr(fc, "anderson", None)
    plain = kwargs.get("expand_bonds") is None
    if on_task and kondo_baths and plain:
        import warnings
        p = _tip_position(fc, imp)
        if not (getattr(fc, "order", None) == "linear"
                and p == fc.impurities[0]):
            raise ValueError("get_stm_spectrum: plain TASK on a chain with "
                "%d Wilson chains misses most of the spectrum unless the "
                "tip spin is the first of site= in order=\"linear\": each "
                "chain keeps its own charge, and the one-site tangent space "
                "of any other layout cannot hold X^dag|GS> (62%% to 95%% of "
                "its weight missing with 15 and 21 Wilson sites in sides, "
                "and 3e-2 to 1e-1 of the peak against ED on short chains "
                "for the other tips). Build the chain again with "
                "site=[tip, others...] and order=\"linear\" (reorder a "
                "per-spin JK, D, Lambda and Delta with it), or pass "
                "expand_bonds=8, "
                "expand_rounds=3 and check the result against a stronger "
                "expansion or submode=\"CVM\" with "
                "fc.cvm_solver=\"variational\"" % nbath)
        warnings.warn("get_stm_spectrum: TASK with the tip spin first in "
            "order=\"linear\" was measured on two chains of up to 15 "
            "Wilson sites (nothing of X^dag|GS> missing, 1.7e-4 of the peak "
            "against ED on 3 sites); check it against a stronger expansion "
            "or submode=\"CVM\" with fc.cvm_solver=\"variational\". With "
            "three chains or more the linear order also costs about 1e-2 "
            "in the energy at maxm 128 (docs/new_methods/chains_on_metal)",
            RuntimeWarning)
    elif on_task and nbath > 1 and plain: # Anderson levels: not measured
        import warnings
        warnings.warn("get_stm_spectrum: this chain has %d Wilson chains, "
            "one per level on the metal; on several Kondo baths the "
            "one-site tangent space of TASK misses part of X^dag|GS>, and "
            "levels were not measured. Pass expand_bonds=8, "
            "expand_rounds=3 and check the result against a stronger "
            "expansion or against submode=\"CVM\" with "
            "fc.cvm_solver=\"variational\"" % nbath, RuntimeWarning)
    if (on_task and nbath <= 1
            and Jc is not None and np.all(np.asarray(Jc) == 0)
            and not (kwargs.get("expand_bonds") is not None
                     and kwargs.get("expand_rounds", 1) >= 2)):
        import warnings
        warnings.warn("get_stm_spectrum: J = 0 on the chain, so the ground "
            "state is a product of the impurity and the bath, and the "
            "one-site tangent space of TASK misses the part of X|GS> that "
            "flips the impurity spin and adds an electron to the bath at "
            "once: the spectrum is incomplete, down to zero. Pass "
            "expand_bonds=4 with expand_rounds=2, or pad_bonds=\"full\" or "
            "mode=\"ED\" on a short chain, or use the BdG density of states "
            "of zaverage.bdg_chain_dos, for J = 0",
            RuntimeWarning)
    # with several Kondo baths the poles come back with their missing
    # weight, and a TASK space that lost a tenth of the sum rule raises
    check = on_task and kondo_baths
    extra = dict(return_poles=True, variance=False) if check else {}
    out, lost, rule = np.zeros(len(es)), 0.0, 0.0
    for X in tunneling_operators(fc, imp=imp, site=site, U=U, q=q):
        Xd = X.get_dagger() # the operator that removes the electron
        for name, x in (((Xd, X), es), ((X, Xd), -es[::-1])):
            r = fc.get_dynamical_correlator(mode=mode, submode=submode,
                        name=name, es=x, delta=delta, **kwargs, **extra)
            c = np.real(r[1])
            out = out + (c if x is es else c[::-1])
            if check:
                info = r[2]
                m = abs(info["missing_weight"])
                lost += m
                rule += abs(info["missing_weight"] + info["elastic"]
                            + complex(np.sum(info["weight"])))
    if check and lost > MISSING_RAISE*rule:
        raise ValueError("get_stm_spectrum: the TASK space misses %.3g of "
                         "the sum rule %.3g of the tip vertices on this "
                         "chain with several Wilson chains, above %g of it; "
                         "raise expand_bonds or fc.maxm, or use the linear "
                         "order with the tip spin first"
                         % (lost, rule, MISSING_RAISE))
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
