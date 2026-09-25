"""Regression tests for the `scale` cluster of the 2026-09-25 hole hunt,
items clean-threshold and small-units.

multioperator.clean_threshold was an absolute 1e-8: every term whose
|coefficient| was at or below it was dropped at every consumption point
(to_terms, MO2matrix, write) and in the canonical form, on every backend,
ED included. So 1e-8*Sz0 had no terms, vev and every correlator of it were
exactly 0, gs_energy() of a Hamiltonian written at that scale was 0, the
canonical form proved 1e-9*(Sz0+1j*Sx0) Hermitian and 1e-9*Sz0 zero, and
meanfield.py's own absolute 1e-10 emptied the mean-field Hamiltonian below
that scale. The drop is relative now: at consumption, to the largest
coefficient of the list; in the canonical form, to the larger of that and
the magnitudes summed into each signature, which is what keeps the rounding
dust of H - H^dagger provably zero.

The anchors are the same calculation at scale 1 and linearity, so no golden
number enters.

small-units: once the terms reached the backend, v2/v3 still got a small-unit
Hamiltonian wrong through two absolute thresholds inside ITensor. toMPO's
svdMPO drops the XX+YY channels of a Heisenberg bond once their squared
singular values sum under 1e-13 (the Neel energy, -1.25 against -2.4936 at
s <= 4e-7 on v3 and s <= 1e-7 on v2), and davidson() replaces every Krylov
direction whose residual is below 1e-10 by a random vector, which costs
accuracy continuously from about s=1e-6 down even on an exact MPO. Both are
fixed in mo_terms.h and chain_session.h by handing ITensor the operator
multiplied by the power of two that brings its largest coefficient into
[1,2), and dividing back; a largest coefficient of 1 or more takes the old
path unchanged. The tests below pin the two mechanisms separately (the MPO
alone through vev on a fixed state, the solver alone through a chain whose
MPO has one channel per bond) and then every consumer of the local solver.
"python"'s own Lanczos stops on an absolute value and breakdown test, not
fixed here, and its strict xfail at the end pins where that still fails.
"""

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, meanfield, multioperator, spinchain


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("ED"), _backend("python"), _backend(3), _backend(2)]


def _heisenberg(sc, n):
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h


def _hubbard(fs, n, phi=0.3, U=2.0):
    """Hubbard chain with a complex, site-dependent hopping phase, so that
    H - H^dagger cancels between conjugated coefficients."""
    h = 0
    for i in range(n-1):
        t = np.exp(1j*phi*(i+1))
        h = h + t*fs.Cdagup[i]*fs.Cup[i+1] + np.conj(t)*fs.Cdagup[i+1]*fs.Cup[i]
        h = h + t*fs.Cdagdn[i]*fs.Cdn[i+1] + np.conj(t)*fs.Cdagdn[i+1]*fs.Cdn[i]
    for i in range(n): h = h + U*fs.Nup[i]*fs.Ndn[i]
    return h


def _chain(version, n=6):
    sc = spinchain.Spin_Chain(["S=1/2"]*n,
            itensor_version=("python" if version == "ED" else version))
    sc.maxm = 30; sc.nsweeps = 10
    return sc


def _mode(version):
    return "ED" if version == "ED" else "DMRG"


# ------------------------------------------------------ the drop itself

@pytest.mark.parametrize("eps", [1e-8, 1e-9, 1e-12])
def test_an_operator_in_small_units_keeps_its_terms(eps):
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    A = eps*sc.Sz[0]
    assert A.to_terms() == [(complex(eps), [("Sz", 1)])]
    assert not A.is_zero()
    assert len(A.simplify().op) == 1
    assert not (eps*(sc.Sz[0] + 1j*sc.Sx[0])).is_hermitian()


def test_exact_zeros_are_still_dropped():
    """The 0*identity() placeholder of "h = 0; h = h + term" goes whatever
    else the list holds, and an all-zero list has no terms."""
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    assert len((0 + 1e-9*sc.Sz[0]).to_terms()) == 1
    assert multioperator.zero().to_terms() == []
    assert (0*sc.Sz[0]).to_terms() == []


def test_the_drop_is_relative_to_the_largest_coefficient():
    """A term twelve orders below the largest one is dropped, nine orders
    below is kept; an absolute 1e-8 kept the first and a relative 1e-8
    would drop the second."""
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    assert len((1e9*sc.Sz[0] + sc.Sz[1]).to_terms()) == 2
    assert len((1e9*sc.Sz[0] + sc.Sz[1]).simplify().op) == 2
    assert len((1e9*sc.Sz[0] + 1e-9*sc.Sz[1]).to_terms()) == 1
    assert len((1e9*sc.Sz[0] + 1e-9*sc.Sz[1]).simplify().op) == 1


# ------------------------------------------- the traps a relative rule has

def test_rounding_dust_of_a_cancellation_is_zero():
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    assert (0.1*sc.Sz[0] + 0.2*sc.Sz[0] - 0.3*sc.Sz[0]).is_zero()


@pytest.mark.parametrize("s", [1.0, 1e-6, 1e-9, 1e-12])
def test_ordinary_hamiltonians_are_proven_hermitian_at_every_scale(s):
    sc = spinchain.Spin_Chain(["S=1/2"]*6)
    fs = fermionchain.Spinful_Fermionic_Chain(4)
    for H, extra in ((_heisenberg(sc, 6), 0.1j*sc.Sz[0]),
                     (_hubbard(fs, 4), 0.1*fs.Cdagup[0]*fs.Cup[1])):
        X = s*H
        assert (X - X.get_dagger()).is_zero()
        assert X.is_hermitian()
        assert len(X.simplify().op) == len(H.simplify().op)
        assert not (s*(H + extra)).is_hermitian()


def test_a_long_accumulation_on_one_signature_is_proven_hermitian():
    """2000 copies of each of two terms minus their adjoints leave about
    2e-12 of the largest coefficient, above a floor on that alone, but
    1e-15 of what was summed into the signature."""
    sc = spinchain.Spin_Chain(["S=1/2"]*4)
    np.random.seed(7)
    cs = np.random.random(2000)
    H = multioperator.msum([c*sc.Sz[0]*sc.Sz[1] for c in cs])
    H = H + multioperator.msum([c*sc.Sx[1]*sc.Sx[2] for c in cs[::-1]])
    for s in (1.0, 1e-9):
        assert (s*H - (s*H).get_dagger()).is_zero()
        assert (s*H).is_hermitian()
        assert len((s*H).simplify().op) == 2


# ------------------------------------------------- through the public API

@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("eps", [1e-8, 1e-9, 1e-12])
def test_vev_is_linear_in_the_scale_of_the_operator(version, eps):
    sc = _chain(version)
    sc.set_hamiltonian(_heisenberg(sc, 6) + 0.3*sc.Sz[0])
    ref = sc.vev(sc.Sz[0], mode=_mode(version))
    assert abs(ref) > 0.1
    assert sc.vev(eps*sc.Sz[0], mode=_mode(version))/eps == \
        pytest.approx(ref, rel=1e-10)


@pytest.mark.parametrize("version", [_backend("ED"), _backend("python")])
def test_correlator_is_quadratic_in_the_scale_of_the_operators(version):
    sc = _chain(version)
    sc.set_hamiltonian(_heisenberg(sc, 6) + 0.3*sc.Sz[0])
    es = np.linspace(-0.5, 5.0, 200)
    eps = 1e-9
    kw = dict(mode=_mode(version), delta=0.2, es=es)
    x, y1 = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), **kw)
    x, y = sc.get_dynamical_correlator(name=(eps*sc.Sz[0], eps*sc.Sz[3]), **kw)
    y1, y = np.asarray(y1), np.asarray(y)
    assert np.max(np.abs(y1)) > 0.1
    assert np.max(np.abs(y/eps**2 - y1)) < 1e-6*np.max(np.abs(y1))


def _scaled_energy(version, s):
    sc = _chain(version)
    sc.set_hamiltonian(s*_heisenberg(sc, 6))
    return sc.gs_energy(mode=_mode(version))/s


@pytest.mark.parametrize("version,s,tol", [
    pytest.param("ED", 1e-8, 1e-10, id="ED-1e-08"),
    pytest.param("ED", 1e-9, 1e-10, id="ED-1e-09"),
    pytest.param("ED", 1e-12, 1e-10, id="ED-1e-12"),
    pytest.param("python", 1e-8, 1e-6, id="python-1e-08")])
def test_ground_state_energy_is_scale_covariant(version, s, tol):
    """It was exactly 0 at s=1e-8 on every backend. ED is exact at any
    scale; "python" is pinned where its own solver still converges."""
    assert _scaled_energy(version, s) == \
        pytest.approx(_scaled_energy(version, 1.0), abs=tol)


# ------------------------------------------------------------ small-units

CPP = [_backend(3), _backend(2)]


def _one_param(version, s, tol):
    """A (version, s, tol) case, skipped where the backend is not built."""
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = ("v%d" % version if version in (2, 3) else version) + "-%.0e" % s
    return pytest.param(version, s, tol, id=ident, marks=marks)


# gs_energy()/s against ED over the scales an ordinary model is written in.
# On v2/v3 it was 7.8e-10 to 1.1e-9 off at s=1e-6 (davidson's randomization)
# and the Neel -1.25 from s=4e-7 (v3) or 1e-7 (v2) down (svdMPO's
# truncation); it is now at the 1e-15 level everywhere. "python" is the
# control: its MPO is scale-free and its own Lanczos, not fixed here, is
# still 2.5e-9 off at s=1e-8, so it is held at 1e-6.
@pytest.mark.parametrize("version,s,tol",
    [_one_param(v, s, 1e-11) for v in (3, 2)
     for s in (1.0, 1e-2, 1e-4, 1e-6, 6e-7, 3e-7, 1e-7, 1e-8)]
    + [_one_param("python", s, 1e-6) for s in (1.0, 1e-4, 1e-6, 1e-8)])
def test_ground_state_energy_is_scale_invariant_on_every_mps_backend(version, s, tol):
    assert _scaled_energy(version, s) == \
        pytest.approx(_scaled_energy("ED", 1.0), abs=tol)


@pytest.mark.parametrize("version", CPP)
@pytest.mark.parametrize("s", [1e-10, 1e-12])
def test_small_units_ground_state_on_the_compiled_backends(version, s):
    """Below every threshold the solver used to meet: -1.06 and -0.90 on v3,
    -1.17 and -0.71 on v2, at these two scales in one measured run."""
    assert _scaled_energy(version, s) == \
        pytest.approx(_scaled_energy("ED", 1.0), abs=1e-11)


@pytest.mark.parametrize("version", CPP)
def test_the_mpo_of_a_small_unit_operator_keeps_every_channel(version):
    """The MPO alone, with no eigensolver: vev(s*X)/s on one fixed state.
    svdMPO's absolute cutoff read the Heisenberg H as 1/3 or 2/3 of itself
    and its XX+YY part as 1/2, from s=6e-7 on v3 and 3e-7 on v2."""
    sc = _chain(version)
    h = _heisenberg(sc, 6)
    xy = 0
    for i in range(5): xy = xy + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    sc.set_hamiltonian(h)
    for X in (h, xy):
        ref = sc.vev(X)
        assert abs(ref) > 1.0
        for s in (6e-7, 3e-7, 1e-7, 1e-10, 1e-12):
            assert sc.vev(s*X)/s == pytest.approx(ref, rel=1e-10)


@pytest.mark.parametrize("version", CPP)
@pytest.mark.parametrize("s", [1e-8, 1e-10])
def test_the_local_solver_is_scale_invariant(version, s):
    """The solver alone: the transverse-field Ising chain has one channel
    per bond, which svdMPO never truncates, and was 1.8e-5 to 4.4e-5 off
    at s=1e-8 and 0.24 to 0.64 at s=1e-10, over the runs measured, through
    davidson()'s absolute 1e-10 alone."""
    def ising(c):
        return (sum(c.Sz[i]*c.Sz[i+1] for i in range(5)) + 0.7*sum(c.Sx))
    sc = _chain(version)
    sc.set_hamiltonian(s*ising(sc))
    ref = _chain(version)
    ref.set_hamiltonian(ising(ref))
    assert sc.gs_energy()/s == \
        pytest.approx(ref.gs_energy(mode="ED"), abs=1e-10)


@pytest.mark.parametrize("version", CPP)
def test_excited_states_in_small_units_match_ed(version):
    """The overlap penalty is an energy and has to follow the operator the
    solver sees: the four lowest levels at s=1e-8 were the Neel -1.25 and
    -0.75 against -2.4936 and -2.0020."""
    s = 1e-8
    sc = _chain(version)
    sc.set_hamiltonian(s*_heisenberg(sc, 6))
    es = np.real(np.asarray(sc.get_excited(n=4)))/s
    ed = np.real(np.asarray(sc.get_excited(n=4, mode="ED")))/s
    assert np.max(np.abs(es - ed)) < 1e-8


@pytest.mark.parametrize("version", CPP)
def test_kpm_correlator_in_small_units_matches_ed(version):
    """KPM reads both band edges, the upper one from the -H solve, and sums
    a shift into H: at s=1e-8 it was 1.1 (v2) and 2.8 (v3) of the ED peak
    off, against 4.1e-4 at s=1, which is what it is at s=1e-8 now."""
    s = 1e-8
    sc = _chain(version)
    sc.set_hamiltonian(s*_heisenberg(sc, 6))
    kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=s*np.linspace(-0.5, 4.0, 46),
              delta=0.2*s)
    x, y = sc.get_dynamical_correlator(**kw)
    x, yed = sc.get_dynamical_correlator(mode="ED", **kw)
    y, yed = np.asarray(y), np.asarray(yed)
    assert np.max(np.abs(y - yed)) < 1e-2*np.max(np.abs(yed))


@pytest.mark.skipif(not cppext.available(3), reason="needs v3")
def test_generalized_ground_state_in_small_units():
    """H|psi> = lambda*A|psi> with A = 1 + 0.2*Sz0 at s=1e-8, against the
    dense generalized eigenproblem; it was -1.389 against -2.5748."""
    import scipy.linalg as sla
    s = 1e-8
    sc = _chain(3)
    h = _heisenberg(sc, 6)
    sc.set_hamiltonian(s*h)
    A = 1.0 + 0.2*sc.Sz[0]
    dense = lambda M: np.asarray(M.todense() if hasattr(M, "todense") else M)
    ed = sc.get_ED_obj()
    lref = np.min(sla.eigh(dense(ed.get_operator(h)), dense(ed.get_operator(A)),
                           eigvals_only=True))
    assert np.real(sc.gs_energy_generalized(A))/s == pytest.approx(lref, abs=1e-8)


# Where the same property still fails: "python"'s _lanczos_ground_state
# stops on an absolute value/breakdown test (2.8e-5 off at s=1e-10), which
# the fix of pyitensor/dmrg.py has to turn into a plain test. Strict.
@pytest.mark.xfail(strict=True, reason="small-units, python Lanczos, not fixed")
@pytest.mark.parametrize("version,s", [pytest.param("python", 1e-10, id="python-1e-10")])
def test_small_units_ground_state_still_open(version, s):
    assert _scaled_energy(version, s) == \
        pytest.approx(_scaled_energy(version, 1.0), abs=1e-6)


@pytest.mark.parametrize("s", [1e-9, 1e-11])
def test_meanfield_does_not_depend_on_units(s, capsys):
    def run(scale):
        sc = spinchain.Spin_Chain(["S=1/2"]*4)
        sc.set_hamiltonian(scale*(_heisenberg(sc, 4) + 0.3*sum(sc.Sz)))
        out = meanfield.spinchain_meanfield(sc, p=0.0, m0=[[0., 0., 0.5]]*4,
                                            maxite=100, mode="ED")
        mz = np.array(out.get_magnetization(mode="ED")).T[:, 2].real
        return out.gs_energy(mode="ED")/scale, mz
    e1, m1 = run(1.0)
    es, ms = run(s)
    assert es == pytest.approx(e1, abs=1e-8)
    assert np.max(np.abs(ms - m1)) < 1e-8


# ---------------------------------------------- small-units, repair pass
#
# KPM rescales H with a single-term shift*Id AutoMPO, and svdMPO skips a
# coefficient below an absolute 1e-14 (isZero(coef,1e-14)), so once the band
# centre 0.62*s fell below it the spectrum lost its position: 1.2 of the ED
# peak off at s=1.5e-14 against 4.1e-4 at 1.8e-14, on both backends. The
# shift (and CVM's z*Id) now goes through the same unit scale as the terms.

@pytest.mark.parametrize("version", CPP)
@pytest.mark.parametrize("s", [1.5e-14, 1e-20])
def test_kpm_shift_survives_below_the_absolute_coefficient_floor(version, s):
    sc = _chain(version)
    sc.set_hamiltonian(s*_heisenberg(sc, 6))
    kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=s*np.linspace(-0.5, 4.0, 46),
              delta=0.2*s)
    x, y = sc.get_dynamical_correlator(**kw)
    x, yed = sc.get_dynamical_correlator(mode="ED", **kw)
    y, yed = np.asarray(y), np.asarray(yed)
    assert np.max(np.abs(y - yed)) < 1e-2*np.max(np.abs(yed))


# Where one scale per operator does not reach: svdMPO truncates bond by bond,
# each bond against its own largest weight and the absolute 1e-13, so a bond
# whose strongest crossing term sits far below the operator's largest
# coefficient loses channels at any units. An O(1) energy offset (sent as
# ('Id', site)) or one-site field next to exchange at 1e-7, and a weak link
# J'=1e-7 in a J=1 chain, all read 1/3 of the exchange on v3 and v2, before
# the fix and after it; "python"'s builder drops them outright (0 and 0.09),
# through the relative cutoff of pyitensor/mpobuilder.py's return sweep. It
# is a property of the bond rather than of the units, and a per-bond cutoff
# is what would turn these into plain tests. Strict.

def _mixed(sc, kind, s):
    """(operator, operator without its small part, small part alone)"""
    h = _heisenberg(sc, 6)
    if kind == "offset":
        return s*h + 1.0, None, h
    if kind == "field":
        return s*h + sc.Sz[0], sc.Sz[0], h
    b23 = sc.Sx[2]*sc.Sx[3] + sc.Sy[2]*sc.Sy[3] + sc.Sz[2]*sc.Sz[3]
    return h - b23 + s*b23, h - b23, b23


@pytest.mark.xfail(strict=True, reason="small-units repair: bond-local truncation, not fixed")
@pytest.mark.parametrize("kind", ["offset", "field", "weak-link"])
@pytest.mark.parametrize("version", [_backend(3), _backend(2), _backend("python")])
def test_a_bond_far_below_the_largest_coefficient_still_open(version, kind):
    s = 1e-7
    sc = _chain(version)
    sc.set_hamiltonian(_heisenberg(sc, 6))
    X, big, small = _mixed(sc, kind, s)
    rest = 1.0 if big is None else np.real(sc.vev(big))
    assert (np.real(sc.vev(X)) - rest)/s == \
        pytest.approx(np.real(sc.vev(small)), rel=1e-6)


@pytest.mark.xfail(strict=True, reason="small-units repair: bond-local truncation, not fixed")
@pytest.mark.parametrize("kind", ["offset", "field"])
@pytest.mark.parametrize("version", CPP)
def test_small_units_next_to_an_order_one_term_still_open(version, kind):
    """gs_energy of s*H + 1 and s*H + Sz0 at s=1e-7: (E-1)/s was -1.25 and
    (E+0.5)/s -1.25 (v3) or -1.00 (v2), against -2.4936 and -2.0926."""
    s = 1e-7
    sc = _chain(version)
    X, big, small = _mixed(sc, kind, s)
    sc.set_hamiltonian(X)
    e, eed = np.real(sc.gs_energy()), np.real(sc.gs_energy(mode="ED"))
    assert e/s == pytest.approx(eed/s, abs=1e-6)
