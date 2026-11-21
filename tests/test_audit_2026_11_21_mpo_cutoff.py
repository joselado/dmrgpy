"""Regressions of docs/audit_2026_11_21_hole_hunt.md, findings 18 and 19,
the per-bond cutoff of the MPO build (mo_terms.h bond_cutoff_factor on v2
and v3, pyitensor/mpobuilder.py _bond_cutoff_factor on "python").

Finding 18. The factor was m^2, m the strongest coefficient crossing the
weakest bond over the largest one, so it was 1 whenever every bond is
crossed by an O(1) term, and a weak multi-site term on such a bond was cut
as if the per-bond cutoff did not exist: a weak link of 1e-7 in a 6-site
ring, the same link next to a strong Sz0*Sz5, an Ising chain plus
1e-7*(XX+YY), and a geometric chain jumped by a strong S_0.S_11 read 0.000
of the weak part on v2, v3 and "python". m now also takes the weakest
multi-site coefficient below WEAK_TERM_TAU = 1e-5 of the largest one,
whatever else crosses its bonds.

Finding 19. On native spinful sites every fermionic factor carries its own
Jordan-Wigner string from site 0, and the factor read the span of a term
from its first listed site, so every hopping and every spin term crossed
every bond to its left: a Wilson chain reversed on the left end of the MPS
(the "sides" layout of the YSR chain on a metal) lost its hoppings of 1e-8
and below. The span now starts at the leftmost site whose factors are not
strings cancelling in pairs; a parity-odd term keeps its real string.

The probe is linearity: with X = big + s*small, [vev(X) - vev(big)] /
(s vev(small)) is exactly 1 on any fixed state, which ED confirms, and a
dropped weak part reads 0. Each MPS backend measures on one fixed state,
its own ground state of a uniform chain, so no eigensolver enters the
ratio. On "python" the factor is also compared with the factor before this
fix, kept verbatim below: it was wrong on every case above, and on the
reference chains of the 2026-11-14 per-bond cutoff it is unchanged, which
makes their MPOs the same to the bit. On v2 and v3 the factor is not bound
to Python, so that check is a comparison of the site arrays across the
rebuild, not a test here."""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, ysr
from dmrgpy.pyitensor import mpobuilder

from test_audit_2026_09_25_scale import _backend, _chain, _heisenberg, _hubbard


def _old_bond_cutoff_factor(ampo):
    """mpobuilder._bond_cutoff_factor before this fix, verbatim."""
    merged = {}
    for t in ampo.terms:
        key = tuple(t.ops)
        merged[key] = merged.get(key, 0.0) + t.coef
    cmax = max((abs(c) for c in merged.values()), default=0.0)
    if not cmax > 0.0:
        return 1.0
    bmax = {}
    for ops, c in merged.items():
        c = abs(c)
        if c <= 1e-14*cmax or len(ops) < 2:
            continue
        sites = [i for _, i in ops]
        for b in range(min(sites), max(sites)):
            bmax[b] = max(bmax.get(b, 0.0), c)
    m = min([1.0] + [v/cmax for v in bmax.values()])
    return m*m


def _factors(chain, op):
    """(old, new) per-bond factor of op, on a "python" chain"""
    ampo = chain._session._ampo(op.to_terms())
    return _old_bond_cutoff_factor(ampo), mpobuilder._bond_cutoff_factor(ampo)


def _ratio(chain, X, big, small, s, **kw):
    """[vev(X) - vev(big)]/(s vev(small)): exactly 1 when X = big + s*small"""
    vx, vb, vs = (np.real(chain.vev(op, **kw)) for op in (X, big, small))
    assert abs(vs) > 0.05  # the probe reads something
    return (vx - vb)/(s*vs)


def _dot(c, i, j):
    return c.Sx[i]*c.Sx[j] + c.Sy[i]*c.Sy[j] + c.Sz[i]*c.Sz[j]


def _xy(c, i, j):
    return c.Sx[i]*c.Sx[j] + c.Sy[i]*c.Sy[j]


def _heis(c, bonds):
    """sum of J S_i.S_j over bonds = [(i, j, J)]"""
    h = 0
    for (i, j, J) in bonds:
        h = h + J*_dot(c, i, j)
    return h


# ------------------------------------------------- finding 18, six spins

N = 6
OBC = [(i, i+1) for i in range(N-1)]
WEAK = (2, 3)
S = 1e-7
KINDS = ["ring", "end_to_end", "ising_xy"]


def _case(c, kind, s):
    """(X, big, small), X = big + s*small, as in the record's repro"""
    strong = [(i, j, 1.0) for (i, j) in OBC if (i, j) != WEAK]
    if kind == "ring":  # the closing bond (0,5) crosses the weak one
        big = _heis(c, strong + [(0, N-1, 1.0)])
        small = _dot(c, *WEAK)
    elif kind == "end_to_end":  # one strong Sz0*Sz5 crosses every bond
        big = _heis(c, strong) + 1.0*c.Sz[0]*c.Sz[N-1]
        small = _dot(c, *WEAK)
    else:  # ZZ at 1 and XY at s on every bond
        big, small = 0, 0
        for (i, j) in OBC:
            big = big + c.Sz[i]*c.Sz[j]
            small = small + _xy(c, i, j)
    return big + s*small, big, small


@pytest.mark.parametrize("kind", KINDS)
def test_a_weak_term_sets_the_factor_whatever_crosses_its_bonds(kind):
    """The factor before the fix was 1 on all three, since every bond is
    crossed by an O(1) term; it is s^2 now, and 1 again with s above the
    window, where nothing was lost."""
    c = _chain("python")
    old, new = _factors(c, _case(c, kind, S)[0])
    assert old == 1.0
    assert new == pytest.approx(S**2, rel=1e-9)
    old, new = _factors(c, _case(c, kind, 2e-5)[0])
    assert old == 1.0 and new == 1.0


SPIN = [_backend("ED"), _backend("python"), _backend(3), _backend(2)]


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("version", SPIN)
def test_a_weak_term_crossed_by_order_one_terms_keeps_its_channels(version, kind):
    """It read 0.000000 on v2, v3 and "python" at s = 1e-7, and 1.000000
    on ED; the state is the ground state of the uniform ring."""
    np.random.seed(1)
    c = _chain(version)
    c.set_hamiltonian(_heis(c, [(i, j, 1.0) for (i, j) in OBC] + [(0, N-1, 1.0)]))
    kw = dict(mode="ED") if version == "ED" else dict(wf=c.get_gs())
    X, big, small = _case(c, kind, S)
    assert _ratio(c, X, big, small, S, **kw) == pytest.approx(1.0, abs=1e-2)


# ---------------------------------------- finding 18, a jumped geometric chain

L = 12
DEEP = [7, 8, 9]  # t = 1e-6, 1e-7, 1e-8


def _t(k):
    return 10.0**(-(k - 1))


def _geometric(c, J, skip=None):
    """S_0.S_1 and S_10.S_11 at 1, t_k S_k.S_k+1 for k = 1..9 with
    t_k = 10^-(k-1), and J S_0.S_11 jumping the whole chain"""
    op = _dot(c, 0, 1) + _dot(c, L-2, L-1)
    for k in range(1, L-2):
        if k != skip:
            op = op + _t(k)*_dot(c, k, k+1)
    if J != 0:
        op = op + J*_dot(c, 0, L-1)
    return op


def test_the_factor_of_a_jumped_chain_follows_its_deepest_bond():
    """With the jump every bond held an O(1) term and the factor was 1;
    without it both rules read the deepest bond, t_9^2."""
    c = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
    old, new = _factors(c, _geometric(c, 1.0))
    assert old == 1.0
    assert new == pytest.approx(_t(9)**2, rel=1e-9)
    old, new = _factors(c, _geometric(c, 0.0))
    assert old == new == pytest.approx(_t(9)**2, rel=1e-9)


@pytest.mark.parametrize("version", SPIN)
def test_a_chain_jumped_by_a_strong_coupling_keeps_its_deep_bonds(version):
    """t = 1e-7 and 1e-8 read 0.0000 on "python" and v3 with J = 1, and
    1.0000 with J = 0 and on ED."""
    np.random.seed(4)
    c = spinchain.Spin_Chain(["S=1/2"]*L,
            itensor_version=("python" if version == "ED" else version))
    c.maxm = 32
    c.nsweeps = 6
    c.set_hamiltonian(_heis(c, [(i, i+1, 1.0) for i in range(L-1)]))
    kw = dict(mode="ED") if version == "ED" else dict(wf=c.get_gs())
    full = _geometric(c, 1.0)
    for k in DEEP:
        r = _ratio(c, full, _geometric(c, 1.0, skip=k), _dot(c, k, k+1),
                   _t(k), **kw)
        assert r == pytest.approx(1.0, abs=1e-2), (k, _t(k))


# ------------------------------------------- finding 19, a chain on a metal

NS = 2
TW = [0.5, 1e-4, 1e-8, 1e-9, 1e-10]  # the Wilson hoppings of the record


def _metal(order, version, t=TW):
    """Two spins, the first on a Wilson chain of len(t)+1 native sites"""
    np.random.seed(11)
    n = len(t) + 1
    fc = ysr.spin_chain_on_metal(NS, n, site=[0], order=order,
                                 wilson=([0.0]*n, list(t)),
                                 itensor_version=version)
    fc.maxm = 32
    fc.nsweeps = 6
    return fc


def test_the_span_starts_after_the_cancelling_strings():
    """With every hopping above the window only the span acts. Reversed
    on the left end ("sides") the chain was read as crossed by every
    fermionic term through the strings, and its factor sat far above the
    deepest hopping; now it is the factor of the same chain after the
    spins ("linear"), which the fix leaves as it was."""
    t = [0.5, 0.1, 1e-2, 1e-3]
    lin = _metal("linear", "python", t)
    sides = _metal("sides", "python", t)
    old_lin, new_lin = _factors(lin, lin.hamiltonian)
    old_sides, new_sides = _factors(sides, sides.hamiltonian)
    assert new_lin == old_lin
    assert new_sides == pytest.approx(new_lin, rel=1e-12)
    assert old_sides > 100*new_sides


def test_a_parity_odd_term_keeps_its_string():
    """A string of an odd number of F on a site is a real part of the
    operator, so a sum of a_i Cup_i keeps the span it had: its factor is
    the weakest a_i squared, as before, and not 1."""
    fn = fermionchain.Spinful_Fermionic_Chain_Native(4, itensor_version="python")
    a = [1.0, 1.0, 1.0, 1e-3]
    op = sum(a[i]*fn.Cup[i] for i in range(4))
    old, new = _factors(fn, op)
    assert new == old
    assert new == pytest.approx(1e-6, rel=1e-9)


METAL = [_backend("python"), _backend(3), _backend(2)]


@pytest.mark.parametrize("order", ["sides", "comb"])
@pytest.mark.parametrize("version", METAL)
def test_a_wilson_chain_off_the_right_end_keeps_its_deep_hoppings(version, order):
    """"sides" (finding 19) read 0, 0, 0 on "python" and 0.25, 0, 0 on v3
    for the hoppings of 1e-8, 1e-9 and 1e-10; "comb", the second spin
    after the chain, so that its Js jumps it (finding 18), read 0 for all
    three on both. The state is the ground state with every hopping at 1;
    the anchor is linearity, since ED on 8 native sites is 65536 states."""
    fc = _metal(order, version)
    w = list(fc.bath_sites[0])
    t = [float(x) for x in fc.lead_chains[0][1]]
    assert t == TW
    H = fc.hamiltonian
    hops = [ysr._hop(fc, w[m], w[m+1]) for m in range(len(t))]
    Hst = H
    for m in range(len(t)):
        Hst = Hst + (1.0 - t[m])*hops[m]
    fc.set_hamiltonian(Hst)
    wf = fc.get_gs()
    for m in range(len(t)):
        r = _ratio(fc, H, H - t[m]*hops[m], hops[m], t[m], wf=wf)
        assert r == pytest.approx(1.0, abs=1e-2), (m, t[m])


# ------------------------------------------------ the reference chains

def _reference(name, version="python"):
    """(chain, operator) for the chains whose MPOs the 2026-11-14 per-bond
    cutoff was checked to leave unchanged. The test reads them on
    "python"; with version=3 they are the chains whose v3 site arrays
    (Chain.mpo_site_arrays of build_operator) are compared across a
    rebuild of the extension."""
    def spins(n):
        return spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    if name == "heisenberg":
        c = spins(12)
        return c, _heisenberg(c, 12)
    if name == "xx":
        c = spins(12)
        return c, 0.5*sum(_xy(c, i, i+1) for i in range(11))
    if name == "dimerized":
        c = spins(12)
        return c, _heis(c, [(i, i+1, 1.0 if i % 2 == 0 else 0.5)
                            for i in range(11)])
    if name == "j1j2":
        c = spins(12)
        return c, _heisenberg(c, 12) + _heis(c, [(i, i+2, 0.5)
                                                 for i in range(10)])
    if name == "inverse_square":
        c = spins(16)
        return c, _heis(c, [(i, j, 1.0/(j - i)**2)
                            for i in range(16) for j in range(i+1, 16)])
    if name == "hubbard":
        fs = fermionchain.Spinful_Fermionic_Chain(4, itensor_version=version)
        return fs, _hubbard(fs, 4)
    if name == "hubbard_native":
        fn = fermionchain.Spinful_Fermionic_Chain_Native(6, itensor_version=version)
        h = sum(ysr._hop(fn, i, i+1) for i in range(5))
        return fn, h + 2.0*sum(fn.Nup[i]*fn.Ndn[i] for i in range(6))
    if name == "nrg_head":
        fc = ysr.spin_chain_on_metal(1, 12, Js=0.01,
                                     JK=0.3*ysr.lambda_factor(4.0),
                                     Lambda=4.0, itensor_version=version)
        return fc, fc.hamiltonian
    raise ValueError(name)


REFERENCE = ["heisenberg", "xx", "dimerized", "j1j2", "inverse_square",
             "hubbard", "hubbard_native", "nrg_head"]


@pytest.mark.parametrize("name", REFERENCE)
def test_the_reference_chains_build_as_before(name):
    """The build reads the operator through the factor only (the automaton
    is exact and the factor sets the cutoff of its one truncating sweep),
    so an equal factor is the same MPO to the bit."""
    c, op = _reference(name)
    old, new = _factors(c, op)
    assert new == old
