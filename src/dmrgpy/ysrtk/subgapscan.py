"""Scans of the subgap energies of a magnetic impurity on a
superconductor (ysr.subgap_energies) over the gap Delta or the uniaxial
anisotropy Dz, and the bisection of the singlet-doublet transition
Delta_c, with every point after the first warm-started from the solved
states of its neighbour.

The Hamiltonian of ysr.kondo_chain is linear in both parameters,

    H(Delta, Dz) = H0 + Delta P + Dz Q,
    P = sum_n (f^dag_{n up} f^dag_{n dn} + f_{n dn} f_{n up}),
    Q = (S^z_imp)^2,

with n over the sites of the Wilson chain, as long as Upen and JH stay
fixed. Their defaults in kondo_chain depend on Dz (Upen = 2(D+|J|+|B|) +
2*2S|Dz| once Dz != 0, JH = Upen), so a scan resolves them once, with
the largest |Dz| of the scan, and passes them to every point. On the
spin-S multiplet the penalty is zero and the Hund term is the constant
-(JH/2)(S(S+1) - 3S/2), so every energy difference, E_ysr included, is
independent of them; the absolute sector energies of a Dz scan differ
from those of kondo_chain's defaults at the same Dz by that constant.

Each 2Sz sector is one chain (a deepcopy of the first kondo_chain with
the sector set), kept for the whole scan. The first point is solved as
subgap_energies solves it: a fresh start with nsweeps, after the seed
stage of ysr._seed_impurity_charge on itensor_version=3 for S > 1/2.
Every later point sets H(x) on the same chain and sweeps nsweeps_warm
times from the state of the nearest point already solved
(set_initial_wf_guess, which works on every backend; on "python" a
changed Hamiltonian drops the session's own state, so restart=False
alone is a cold start there). The states of neighbouring points are
close, and a warm state already holds one electron on every impurity
site, which no term of H changes, so the seed stage is not needed after
the first point. ysr._check_impurity_charge still checks every solved
state for S > 1/2. The points are visited in the order given, each from
the one before, so a scan should be monotonic and closely spaced.

The number of warm sweeps. Both backends take the bond dimension of the
state already present as the floor of the sweep schedule
(Chain::gs_energy in mpscpp3/chain_session.h and pyitensor/chain.py's
gs_energy, floor_dim into make_sweeps_ramped). A cold start spends
bond_ramp_fraction = 0.5 of its sweeps at a bond dimension growing from
10, which costs little; a warm start sweeps at maxm from the first sweep.
So nsweeps warm sweeps cost about 1.5 times a cold solve (measured on
26 sites), and the default is nsweeps_warm = min(nsweeps, max(8,
nsweeps//3)), the same on every backend (warm_sweeps). "python" mixes no
noise into a warm start (pyitensor passes noisy_start only on a fresh
start). v3 does: a state already at maxm takes make_sweeps, with the
chain's noise (1e-7) on the first half of the warm sweeps and none on
the second, so a smaller nsweeps_warm shortens the noisy stage too. A
converged metastable state is the same risk at 8 sweeps or at 20, and
cold_every is the guard in both cases.

A warm start carries along the state it starts from, a metastable one
too, for instance the singlet plus a free quasiparticle in the odd
sector, which a cold solve escapes only by chance. Two guards: the
warning of subgap_energies, raised at every point where the chain
reaches below the gap and |E_ysr| exceeds the bath's quasiparticle gap,
and cold_every=k, which solves every k-th point once more with
ysr.subgap_energies itself and returns the largest warm-minus-cold
energy as "warm_cold", with a warning once a warm sector is more than
1e-4 |Delta| above the cold one or the two disagree on the ground state.

Checked (tests/test_ysr_scan.py). On 4 sites at J=0.6, Lambda=4, maxm=64,
nsweeps=15, a Delta scan (0.05, 0.1, 0.2) for S=1/2 and for S=1 at
Dz=0.05, and a Dz scan (-0.1 to 0.1) for S=1 at Delta=0.1, give every
sector energy of ED to 3e-13 on v3 and 4e-11 on "python" at the default
nsweeps_warm=8, and the cold subgap_energies of the same point to
4e-11; the bisection of Delta_c on 3 sites visits the same points as
on ED. On 26 sites (v3, rho J=0.1,
maxm=128, nsweeps=20) a warm scan over five gaps around Delta_c gives
every sector energy of the cold protocol to 1e-12 with nsweeps_warm=20,
8, 6 and 3, and find_transition gives Delta_c = 5.2853e-5 at z=0 as the
cold bisection of docs/development.md 4.8f. Measured on d22 at load 5
to 13, a warm point costs 0.25 (nsweeps_warm=3) to 0.5 (8) of a cold one
for S=1/2, and 0.3 (5) to 0.5 (8) for S=1 at maxm=192. With
nsweeps_warm = nsweeps it costs 1.5 times a cold point. The floor of 8
is set by "python": on the 4-site S=1 scan at Dz=0.05, 4 warm sweeps
left 4e-5 and a wrong Sz_ground, 6 left 7e-9 and 8 left 4e-11, while v3
is at 3e-13 from 4 sweeps on. A large jump between two points at a
small nsweeps_warm is not caught by the quasiparticle-gap warning, only
by cold_every: on 14 sites at nsweeps_warm=2, Delta from 0.1 back to
0.001 gave E_ysr/Delta = 0.94 to 1.00 instead of 0.86, and even a
descending scan in steps of sqrt(10) left 1.4e-3 at Delta = 0.0032.
"""

import copy
import time
import warnings

import numpy as np

from .. import wilsonchain
from .. import ysr


WARM_FLOOR = 8 # default floor of the warm sweeps, set by "python"


def warm_sweeps(nsweeps, nsweeps_warm=None):
    """The sweeps of a warm point: nsweeps_warm if given, else the larger
    of WARM_FLOOR and nsweeps//3, never more than nsweeps."""
    if nsweeps_warm is not None:
        if int(nsweeps_warm) < 1:
            raise ValueError("nsweeps_warm must be a positive integer")
        return int(nsweeps_warm)
    return min(int(nsweeps), max(WARM_FLOOR, int(nsweeps)//3))


def resolve_penalties(J, D=1.0, B=0.0, S="1/2", Dzs=(0.0,), Upen=None,
                      JH=None):
    """(Upen, JH) of a scan: the kondo_chain defaults at the largest |Dz|
    of the scan, unless given."""
    ni = ysr.two_spin(S)
    if Upen is None:
        Upen = 2.0*(D + abs(J) + abs(B))
        dzmax = max(abs(float(d)) for d in Dzs)
        if dzmax != 0.0:
            Upen = Upen + 2.0*ni*dzmax
    if JH is None:
        JH = Upen
    return float(Upen), float(JH)


def default_sectors(S, Dzs, B):
    """The 2Sz sectors subgap_energies solves, made the same for every
    point: the lowest |2Sz| of each parity when Dz = 0 and B = 0 at every
    point, else every 2Sz up to 2S+1 (both signs if B != 0)."""
    ni = ysr.two_spin(S)
    if all(float(d) == 0.0 for d in Dzs) and B == 0.0:
        return [0, 1]
    sectors = list(range(0, ni + 2))
    if B != 0.0:
        sectors = [-s for s in sectors[::-1] if s > 0] + sectors
    return sectors


def summarize(ens, ni, qp_gap=0.0, reaches_gap=False, where=""):
    """The keys of ysr.subgap_energies from the sector energies ens
    (2Sz -> energy), with its metastable-sector warning."""
    def lowest(par):
        return min(e for s, e in ens.items() if s % 2 == par)
    out = {"E_even": lowest(0), "E_odd": lowest(1),
           "E_free": lowest(ni % 2), "E_screened": lowest((ni + 1) % 2)}
    out["E_ysr"] = out["E_free"] - out["E_screened"]
    if ni == 1:
        out["ground_state"] = "singlet" if out["E_ysr"] > 0 else "doublet"
    else:
        out["ground_state"] = "screened" if out["E_ysr"] > 0 else "free"
    out["sectors"] = dict(ens)
    e0 = min(ens.values())
    tie = 1e-10*max(1.0, abs(e0)) # a degenerate multiplet: smallest |2Sz|
    out["Sz_ground"] = min((s for s, e in ens.items() if e <= e0 + tie),
                           key=lambda s: (abs(s), s))
    if reaches_gap and abs(out["E_ysr"]) > qp_gap*(1 + 1e-2):
        warnings.warn("subgap_scan%s: |E_free - E_screened| = %.3e is "
            "larger than the lowest quasiparticle energy of the bath, %.3e, "
            "so one sector solve is not its ground state; raise maxm, "
            "nsweeps or nsweeps_warm, or check with cold_every"
            % (where, abs(out["E_ysr"]), qp_gap))
    return out


class _SectorChains(object):
    """The sector chains of one scan, H(Delta, Dz) and the solve of one
    point, cold or warm."""

    def __init__(self, n, J, Lambda, D, S, B, z, mesh, Upen, JH, maxm,
                 nsweeps, nsweeps_warm, sectors, mode, chain_kwargs, Delta0,
                 Dz0):
        if n < 1:
            raise ValueError("the Wilson chain needs at least one site")
        self.n, self.ni, self.mode = n, ysr.two_spin(S), mode
        self.nsweeps = int(nsweeps)
        self.nsweeps_warm = warm_sweeps(nsweeps, nsweeps_warm)
        self.maxm = int(maxm)
        self.kw = dict(J=J, Lambda=Lambda, D=D, B=B, z=z, mesh=mesh,
                       Upen=Upen, S=S, JH=JH)
        self.chain_kwargs = dict(chain_kwargs)
        # the mesh of the chain, for the quasiparticle gap of every point
        self.eps, self.t = wilsonchain.wilson_chain(n, Lambda=Lambda,
                                                    band=(-D, D), z=z,
                                                    mesh=mesh)
        fc = self.new_chain(Delta0, Dz0)
        # H0: the same chain with no pairing and no anisotropy
        h0 = ysr.kondo_chain(n, Delta=0.0, Dz=0.0, **self.kw,
                             **self.chain_kwargs).hamiltonian
        P = 0
        for f in fc.channels[0]:
            P = P + fc.Cdagup[f]*fc.Cdagdn[f] + fc.Cdn[f]*fc.Cup[f]
        Sz = ysr.impurity_spin(fc, fc.impurity)[2]
        self.h0, self.P, self.Q = h0, P, Sz*Sz
        self.sectors = sorted(int(s) for s in sectors)
        if not any((s - self.ni) % 2 == 0 for s in self.sectors) or \
           not any((s - self.ni) % 2 == 1 for s in self.sectors):
            raise ValueError("sectors must contain both parities of 2Sz")
        self.chains = {}
        for sz in self.sectors:
            c = copy.deepcopy(fc)
            c.set_conserved_sector(Sz=sz)
            self.chains[sz] = c
        self.first = (float(Delta0), float(Dz0)) # H of the chains now

    def new_chain(self, Delta, Dz):
        fc = ysr.kondo_chain(self.n, Delta=Delta, Dz=Dz, **self.kw,
                             **self.chain_kwargs)
        fc.maxm, fc.nsweeps = self.maxm, self.nsweeps
        return fc

    def hamiltonian(self, Delta, Dz):
        h = self.h0 + Delta*self.P
        if Dz != 0.0:
            h = h + Dz*self.Q
        return h

    def gap(self, Delta):
        """(qp_gap, reaches_gap) of the bath at this Delta, as in
        kondo_chain."""
        qp = ysr.bath_quasiparticle_gap(self.eps, self.t, Delta)
        reaches = self.n > 1 and Delta != 0 and \
            self.t[-1] < 0.1*abs(Delta)
        return qp, bool(reaches)

    def solve(self, Delta, Dz, start=None, nsweeps=None):
        """Sector energies at (Delta, Dz). start=None is the cold first
        solve of the chains as built; otherwise a dict 2Sz -> state to
        sweep from, nsweeps times (default nsweeps_warm). Returns
        (energies, states)."""
        if Delta == 0.0:
            raise ValueError("subgap_scan needs a gap, Delta != 0 at every "
                             "point (see ysr.subgap_energies)")
        dmrg = self.mode == "DMRG"
        ens, states = {}, {}
        for sz in self.sectors:
            c = self.chains[sz]
            if start is None:
                if (float(Delta), float(Dz)) != self.first:
                    raise RuntimeError("a cold solve needs the chains as "
                                       "built")
                c.nsweeps = self.nsweeps
                if self.ni > 1 and dmrg and c.itensor_version == 3:
                    ysr._seed_impurity_charge(c, self.ni)
            else:
                c.set_hamiltonian(self.hamiltonian(Delta, Dz), restart=False)
                if dmrg:
                    c.set_initial_wf_guess(start[sz])
                c.nsweeps = self.nsweeps_warm if nsweeps is None \
                    else int(nsweeps)
            # c.nsweeps stays at the value of this solve: the stored state
            # answers only for the parameters it was solved with, so
            # resetting it here made the vev and get_gs below solve again
            ens[sz] = c.gs_energy(mode=self.mode)
            if self.ni > 1:
                ysr._check_impurity_charge(c, self.ni, self.mode)
            if dmrg:
                states[sz] = c.get_gs()
        return ens, states

    def point(self, Delta, Dz, start=None, nsweeps=None):
        t0 = time.time()
        ens, states = self.solve(Delta, Dz, start=start, nsweeps=nsweeps)
        qp, reaches = self.gap(Delta)
        out = summarize(ens, self.ni, qp, reaches,
                        where=" (Delta = %.6g, Dz = %.6g)" % (Delta, Dz))
        out.update(Delta=float(Delta), Dz=float(Dz), qp_gap=qp,
                   reaches_gap=reaches, warm=start is not None,
                   time=time.time() - t0)
        return out, states

    def cold(self, Delta, Dz):
        """The reference: ysr.subgap_energies on a fresh kondo_chain."""
        return ysr.subgap_energies(self.new_chain(Delta, Dz), mode=self.mode,
                                   sectors=self.sectors)


def _compare(warm, cold, Delta, where):
    """warm - cold per sector, with a warning when the warm solve is
    above the cold one by more than 1e-4 |Delta| (E_ysr/Delta off by
    1e-4), or when the two disagree on the ground state: on its parity,
    or on its 2Sz where the cold sectors of the two differ by more than
    that."""
    tol = 1e-4*abs(Delta)
    d = {s: warm["sectors"][s] - cold["sectors"][s] for s in warm["sectors"]}
    s = max(d, key=lambda k: d[k])
    if d[s] > tol:
        warnings.warn("subgap_scan%s: the warm solve of sector 2Sz = %d is "
            "%.3e above a cold one, so the warm start carried a metastable "
            "or unconverged state along; raise nsweeps_warm or solve this "
            "point cold" % (where, s, d[s]))
    sw, sc = warm["Sz_ground"], cold["Sz_ground"]
    split = abs(cold["sectors"][sw] - cold["sectors"][sc])
    if warm["ground_state"] != cold["ground_state"] or \
       (sw != sc and split > tol):
        warnings.warn("subgap_scan%s: the warm solve gives the ground state "
            "%s at 2Sz = %d, a cold one %s at 2Sz = %d; raise nsweeps_warm "
            "or solve this point cold" % (where, warm["ground_state"], sw,
                                          cold["ground_state"], sc))
    return d


def subgap_scan(n, Deltas=None, Dzs=None, J=0.2, Lambda=4.0, D=1.0,
                S="1/2", Dz=None, B=0.0, z=0.0, mesh="mean", Delta=None,
                Upen=None, JH=None, maxm=64, nsweeps=20, nsweeps_warm=None,
                sectors=None, cold_every=None, mode="DMRG", **chain_kwargs):
    """ysr.subgap_energies at every point of a scan over the gap (Deltas,
    at fixed Dz) or over the anisotropy (Dzs, at fixed Delta), the points
    after the first warm-started from the previous one.

    Exactly one of Deltas and Dzs is a sequence. A Deltas scan takes Dz
    (default 0) and refuses Delta; a Dzs scan needs Delta and refuses Dz.
    The chain is ysr.kondo_chain(n, J=J, Lambda=Lambda, D=D, S=S, B=B,
    z=z, mesh=mesh, Upen=Upen, JH=JH, **chain_kwargs), so J is the
    coupling on the chain (J*ysr.lambda_factor(Lambda) for J in the
    continuum on the "mean" mesh). Upen and JH, if not given, are
    kondo_chain's defaults at the largest |Dz| of the scan, fixed for the
    whole scan so that H is linear in the scanned parameter (see the
    module docstring). The first point is solved with nsweeps from a
    fresh start, every later one with nsweeps_warm sweeps from the states
    of the point before, in the order given, so the points should be
    monotonic and closely spaced. The default nsweeps_warm is the larger
    of 8 and nsweeps//3, at most nsweeps (warm_sweeps). A smaller
    nsweeps_warm is checked only by cold_every. sectors are the 2Sz
    values solved at every point; the default is subgap_energies's, made
    the same for all points (every 2Sz up to 2S+1 once any Dz of the scan
    is nonzero).

    Returns a dict with the arrays "Delta", "Dz", "E_even", "E_odd",
    "E_free", "E_screened", "E_ysr", "qp_gap", "time" (seconds per
    point), the lists "sectors" (2Sz -> energy at every point),
    "ground_state", "Sz_ground" and "points" (the per-point dicts, with
    the keys of subgap_energies), and "Upen", "JH", "nsweeps_warm" (the
    value used). With cold_every=k every k-th point (k, 2k, ...) is
    solved once more with ysr.subgap_energies on a fresh chain: "cold"
    maps the index to that result and "warm_cold" is the largest
    |E_warm - E_cold| over those points and sectors (None without
    cold_every). The metastable-sector warning of subgap_energies is
    raised at every point."""
    if (Deltas is None) == (Dzs is None):
        raise ValueError("give exactly one of Deltas and Dzs")
    if Deltas is not None:
        if Delta is not None:
            raise ValueError("a Deltas scan takes no Delta")
        Deltas = [float(x) for x in np.atleast_1d(Deltas)]
        Dzs = [float(0.0 if Dz is None else Dz)]*len(Deltas)
    else:
        if Delta is None:
            raise ValueError("a Dz scan needs Delta")
        if Dz is not None:
            raise ValueError("a Dzs scan takes no Dz")
        Dzs = [float(x) for x in np.atleast_1d(Dzs)]
        Deltas = [float(Delta)]*len(Dzs)
    if len(Deltas) == 0:
        raise ValueError("the scan has no points")
    if any(x == 0.0 for x in Deltas):
        raise ValueError("subgap_scan needs a gap, Delta != 0 at every "
                         "point (see ysr.subgap_energies)")
    if cold_every is not None and int(cold_every) < 1:
        raise ValueError("cold_every must be a positive integer")
    Upen, JH = resolve_penalties(J, D, B, S, Dzs, Upen, JH)
    if sectors is None:
        sectors = default_sectors(S, Dzs, B)
    sc = _SectorChains(n, J, Lambda, D, S, B, z, mesh, Upen, JH, maxm,
                       nsweeps, nsweeps_warm, sectors, mode, chain_kwargs,
                       Deltas[0], Dzs[0])
    points, cold, diffs, states = [], {}, [], None
    for i, (x, dz) in enumerate(zip(Deltas, Dzs)):
        p, states = sc.point(x, dz, start=states if i > 0 else None)
        points.append(p)
        if cold_every and i > 0 and i % int(cold_every) == 0:
            cold[i] = sc.cold(x, dz)
            d = _compare(p, cold[i], x, " (Delta = %.6g, Dz = %.6g)" % (x, dz))
            diffs.extend(abs(v) for v in d.values())
    out = {"Delta": np.array(Deltas), "Dz": np.array(Dzs)}
    for k in ("E_even", "E_odd", "E_free", "E_screened", "E_ysr", "qp_gap",
              "time"):
        out[k] = np.array([p[k] for p in points])
    for k in ("sectors", "ground_state", "Sz_ground"):
        out[k] = [p[k] for p in points]
    out.update(points=points, cold=cold, Upen=Upen, JH=JH,
               nsweeps_warm=sc.nsweeps_warm,
               warm_cold=max(diffs) if diffs else None)
    return out


def find_transition(n, J, Delta_lo, Delta_hi, ratio=1.04, Lambda=4.0,
                    D=1.0, S="1/2", Dz=0.0, B=0.0, z=0.0, mesh="mean",
                    Upen=None, JH=None, maxm=64, nsweeps=20, nsweeps_warm=None,
                    sectors=None, nexpand=8, check_cold=False, mode="DMRG",
                    **chain_kwargs):
    """Delta_c, the gap where E_ysr = E_free - E_screened changes sign, by
    bisection in ln Delta with warm starts.

    E_ysr > 0 (screened) is taken below Delta_c and E_ysr < 0 (free)
    above. If [Delta_lo, Delta_hi] does not bracket the sign change, the
    bracket is moved by factors of 1.5, at most nexpand times. Then it is
    halved in ln Delta down to Delta_hi/Delta_lo <= ratio, and Delta_c is
    the zero of x = E_ysr/Delta interpolated linearly in ln Delta between
    the bracket ends (the protocol of docs/development.md, 4.8f). Delta_lo
    is solved cold, every other point with nsweeps_warm sweeps (default
    warm_sweeps(nsweeps), the larger of 8 and nsweeps//3) from the states
    of the nearer bracket end. The other keywords are those of
    subgap_scan. check_cold=True solves the two final bracket ends once
    more with ysr.subgap_energies and returns the largest |E_warm -
    E_cold| as "warm_cold", and Delta_c from the cold energies as
    "Delta_c_cold".

    Returns a dict with "Delta_c", the bracket "lo", "hi" and its "x_lo",
    "x_hi", the arrays "Delta" and "x" of every point visited (sorted in
    Delta), "points" (the per-point dicts of subgap_scan, in the same
    order), "monotonic" (x decreasing along the visited points),
    "nsweeps_warm" (the value used) and "warm_cold"."""
    lo, hi = float(Delta_lo), float(Delta_hi)
    if not 0 < lo < hi:
        raise ValueError("find_transition needs 0 < Delta_lo < Delta_hi")
    Upen, JH = resolve_penalties(J, D, B, S, [Dz], Upen, JH)
    if sectors is None:
        sectors = default_sectors(S, [Dz], B)
    sc = _SectorChains(n, J, Lambda, D, S, B, z, mesh, Upen, JH, maxm,
                       nsweeps, nsweeps_warm, sectors, mode, chain_kwargs,
                       lo, Dz)
    pts, states = {}, {}

    def x(Delta):
        if Delta in pts:
            return pts[Delta]["E_ysr"]/Delta
        if not pts:
            p, st = sc.point(Delta, Dz)
        else:
            near = min(states, key=lambda d: abs(np.log(d/Delta)))
            p, st = sc.point(Delta, Dz, start=states[near])
        pts[Delta], states[Delta] = p, st
        return p["E_ysr"]/Delta

    def keep(*ds): # only the bracket ends are warm starts later
        for d in list(states):
            if d not in ds:
                del states[d]

    x(lo)
    for _ in range(int(nexpand) + 1):
        if x(lo) <= 0:
            lo, hi = lo/1.5, lo
        elif x(hi) >= 0:
            lo, hi = hi, hi*1.5
        else:
            break
        keep(lo, hi)
    else:
        raise RuntimeError("find_transition: no sign change of E_ysr found "
            "between %.4g and %.4g" % (min(pts), max(pts)))
    keep(lo, hi)
    while hi/lo > ratio:
        m = float(np.sqrt(lo*hi))
        if x(m) > 0:
            lo = m
        else:
            hi = m
        keep(lo, hi)

    def interpolate(xl, xh):
        a, b = np.log(lo), np.log(hi)
        return float(np.exp(a + (b - a)*xl/(xl - xh)))
    out = {"Delta_c": interpolate(x(lo), x(hi)), "lo": lo, "hi": hi,
           "x_lo": x(lo), "x_hi": x(hi), "Upen": Upen, "JH": JH,
           "nsweeps_warm": sc.nsweeps_warm,
           "warm_cold": None}
    ds = sorted(pts)
    out["Delta"] = np.array(ds)
    out["x"] = np.array([pts[d]["E_ysr"]/d for d in ds])
    out["points"] = [pts[d] for d in ds]
    out["monotonic"] = bool(np.all(np.diff(out["x"]) < 0))
    if check_cold:
        diffs, xc = [], []
        for d in (lo, hi):
            c = sc.cold(d, Dz)
            if pts[d]["warm"]:
                diffs.extend(abs(v) for v in _compare(pts[d], c, d,
                             " (Delta = %.6g)" % d).values())
            xc.append(c["E_ysr"]/d)
        out["warm_cold"] = max(diffs) if diffs else 0.0
        out["Delta_c_cold"] = interpolate(*xc)
    return out
