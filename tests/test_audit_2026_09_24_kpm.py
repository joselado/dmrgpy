"""Regression tests for the kpm cluster of the 2026-09-24 audit.

Each test here locks in finding 3, 4 or 5 of
`docs/audit_2026_09_24_hole_hunt.md`, which records the original symptom,
the reproduction that was executed and the reviewer's analysis. Finding 6
of the same cluster is documentation only (the band-centre caveat in
`algebra/kpm.py::polynomials_for_broadening`), so it has no test here.

- #3: `get_distribution(mode="ED")` integrated to exactly 1/scale (0.1 at
  the default `scale=10`), because `edtk/distribution.py` undid only the
  1/pi of `dm_vivj_energy`'s pi/scale normalization, and with `xs=` its
  imaginary part was a copy of its real part. A 2-site chain on
  `itensor_version=3` got it without asking for ED.
- #4: every DMRG KPM route reconstructed the spectrum from n+2 Chebyshev
  moments (n+1 on the accelerated path at odd n) where the calibration
  and the ED route use n, so the DMRG line was narrower and taller than
  ED's by about 2/n: 7.4 per cent on the band-centre pole below at
  delta=0.2. The moments are now cut to n before the reconstruction.
- #5: a non-integer `kpm_n_scale` was rounded down silently on
  `"python"`, ED and `julia_live` (1.5 gave exactly 1x, anything below 1
  the 16-moment floor) where v2/v3 raised, and at `kpm_n_scale<=0` the two
  copies of the calibration disagreed. One validator now runs before
  dispatch on every route.

The chains are 2 to 4 sites. `itensor_version=3` cases skip themselves
when mpscpp3 is not compiled, and the `julia_live` ones when there is no
working Julia toolchain.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.algebra import kpm
from dmrgpy.kpmdmrg import dynamical_correlator_from_moments

from _helpers import julia_live_param


needs_v3 = pytest.mark.skipif(not cppext.available(3),
                              reason="mpscpp3 not compiled")
DMRG_BACKENDS = ["python", pytest.param(3, marks=needs_v3)]


# ---------------------------------------------------------------- helpers

def dimers(itensor_version="python"):
    """Two decoupled S=1/2 Heisenberg dimers, J=1. <Sz_0;Sz_0> has exactly
    one pole, at omega=1 with weight 1/4, and the many-body band is
    [E0,Emax] = [-1.5,0.5], so on the default bandwidth-centred window the
    pole sits at x=0, the band centre, where the calibration FWHM=2*delta
    is exact. Every quantity below is exact on it (emin/emax to 1e-12),
    so ED and DMRG have no reason to differ at all."""
    sc = spinchain.Spin_Chain(["S=1/2"] * 4, itensor_version=itensor_version)
    h = 0
    for (i, j) in [(0, 1), (2, 3)]:
        h = h + sc.Sx[i] * sc.Sx[j] + sc.Sy[i] * sc.Sy[j] + sc.Sz[i] * sc.Sz[j]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 20
    return sc


def heisenberg(n=4, itensor_version="python"):
    """Open S=1/2 Heisenberg chain, J=1."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 20
    return sc


def reconstruction_nodes(emin, emax, scale, n, xmax=0.9):
    """The energies dynamical_correlator_from_moments reconstructs on
    before interpolating onto `es` (its own 10n-point grid), kept to
    |x|<xmax, inside the ED route's valid window |x|<=0.95. At these
    energies the interpolation is exact, so a comparison against ED there
    sees the moments and the kernel only; anywhere else the linear
    interpolation off that grid costs about 5e-4 of the peak on this
    pole at delta=0.1, which is not what these tests are about."""
    xs = 0.99 * np.linspace(-1.0, 1.0, int(n * 10), endpoint=False)
    es = xs / scale + (emin + emax) / 2. - emin
    return es[np.abs(xs) < xmax]


# ------------------------------------------ #4: exactly n moments, DMRG

# delta=0.1 calibrates to n=52 on the dimers and delta=0.3 to n=17: one
# even and one odd count, since the accelerated recursion emits moments
# in pairs and used to return n+1 at odd n where the full one returned n+2
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
@pytest.mark.parametrize("delta,parity", [(0.1, 0), (0.3, 1)])
@pytest.mark.parametrize("accelerate", [True, False])
def test_dmrg_kpm_returns_exactly_the_calibrated_moment_count(
        itensor_version, delta, parity, accelerate):
    sc = dimers(itensor_version)
    sc.kpm_accelerate = accelerate
    mus, emin, emax, scale, n, _d = sc.get_dynamical_correlator_moments(
        name=[sc.Sz[0], sc.Sz[0]], delta=delta)
    assert n % 2 == parity   # the premise: this case exercises that parity
    assert n == kpm.polynomials_for_broadening(1.0 / scale, delta)
    assert len(mus) == n, \
        "len(mus)=%d, n=%d: the kernel would see the wrong N" % (len(mus), n)


@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_energy_truncated_kpm_returns_exactly_the_calibrated_moment_count(
        itensor_version):
    """The ground-state-anchored route is a separate session method on
    itensor_version=3 (kpm_dynamical_correlator_truncated), with its own
    moment loop, so it is pinned separately."""
    sc = dimers(itensor_version)
    sc.kpm_energy_truncate = True
    mus, emin, emax, scale, n, _d = sc.get_dynamical_correlator_moments(
        name=[sc.Sz[0], sc.Sz[0]], delta=0.1)
    assert len(mus) == n


@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_band_centre_pole_is_the_same_curve_on_dmrg_and_ed(itensor_version):
    """Before the cut the DMRG line was 0.952 of 2*delta wide against
    ED's 0.988 here, and max|ED-DMRG| was 4.6e-02 on a peak of 1.2202 at
    these very energies; with the same n moments on both routes the two
    are one curve (measured 1.4e-14 to 3.1e-14)."""
    delta = 0.1
    sc = dimers(itensor_version)
    name = [sc.Sz[0], sc.Sz[0]]
    mus, emin, emax, scale, n, _d = sc.get_dynamical_correlator_moments(
        name=name, delta=delta)
    es = reconstruction_nodes(emin, emax, scale, n)
    _x, yfm = dynamical_correlator_from_moments(mus, emin, emax, scale, n, es)
    _x, ydm = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                          name=name, delta=delta, es=es)
    _x, yed = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                                          name=name, delta=delta, es=es)
    yed = np.asarray(yed)
    assert np.max(np.real(yed)) == pytest.approx(1.2202, abs=1e-3)
    assert np.max(np.abs(np.asarray(yfm) - yed)) < 1e-6
    assert np.max(np.abs(np.asarray(ydm) - yed)) < 1e-6


@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_kpm_accelerate_does_not_change_the_curve_at_odd_n(itensor_version):
    """kpm_accelerate is a speed flag, but at odd n it used to change the
    moment count and so the curve (n+1 moments against n+2: a peak of
    0.5731 against 0.5849 on a 4-site staggered chain in the record)."""
    delta = 0.3   # n=17 on the dimers
    es = np.linspace(0.0, 2.0, 81)
    ys = []
    for accelerate in (True, False):
        sc = dimers(itensor_version)
        sc.kpm_accelerate = accelerate
        _x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                            name=[sc.Sz[0], sc.Sz[0]],
                                            delta=delta, es=es)
        ys.append(np.asarray(y))
    assert np.max(np.abs(ys[0] - ys[1])) < 1e-8


# ------------------------------------------ #4: exactly n moments, julia

@pytest.fixture
def spy_julia_reconstruction(monkeypatch):
    """Record what mpsjulialive's KPM hands the shared reconstruction,
    the only place its moments are visible (julia_live has no
    get_dynamical_correlator_moments)."""
    from dmrgpy.mpsjulialive import dynamics as jldyn
    seen = []
    orig = jldyn.dynamical_correlator_from_moments

    def spy(mus, emin, emax, scale, n, es, **kwargs):
        seen.append((len(mus), emin, emax, scale, n))
        return orig(mus, emin, emax, scale, n, es, **kwargs)
    monkeypatch.setattr(jldyn, "dynamical_correlator_from_moments", spy)
    return seen


@pytest.mark.parametrize("itensor_version", [julia_live_param()])
def test_julia_live_kpm_uses_exactly_the_calibrated_moment_count(
        itensor_version, spy_julia_reconstruction):
    delta = 0.1
    sc = dimers("python")
    sc.setup_julia()
    name = [sc.Sz[0], sc.Sz[0]]
    sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                delta=delta, es=np.linspace(0.0, 2.0, 11))
    nmus, emin, emax, scale, n = spy_julia_reconstruction[-1]
    assert nmus == n
    # and so it is ED's curve, at the reconstruction's own nodes
    es = reconstruction_nodes(emin, emax, scale, n)
    _x, yjl = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                          name=name, delta=delta, es=es)
    _x, yed = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                                          name=name, delta=delta, es=es)
    assert np.max(np.abs(np.asarray(yjl) - np.asarray(yed))) < 1e-6


# --------------------------------------- #3: ED get_distribution weight

def distribution_chain(n=4, itensor_version="python"):
    sc = heisenberg(n, itensor_version)
    return sc


def exact_distribution_weights(sc, X):
    """(eigenvalues, weights) of X's distribution in the ground state,
    P(x) = sum over the eigenvectors v of X at x of |<v|GS>|^2, from dense
    diagonalizations done here with numpy."""
    ed = sc.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    _e, vs = np.linalg.eigh(h)
    gs = vs[:, 0]
    xm = np.array(ed.MO2matrix(X).todense())
    lam, u = np.linalg.eigh(xm)
    w = np.abs(np.conjugate(u.T) @ gs) ** 2
    vals = np.unique(np.round(lam, 8))
    return vals, np.array([np.sum(w[np.abs(lam - v) < 1e-6]) for v in vals])


# Sz_0+Sz_1 has eigenvalues -1, 0, 1 and max|eig| = 1, so scale=None
# auto-scales to 2: the old weight 1/scale was 0.5 there, 0.1 at 10, and
# only a max|eig X| = 1/2 operator such as Sz_0 alone would have hidden it
@pytest.mark.parametrize("scale", [2.0, 10.0, None])
def test_ed_distribution_integrates_to_one(scale):
    sc = distribution_chain()
    X = sc.Sz[0] + sc.Sz[1]
    x, y = sc.get_distribution(mode="ED", X=X, delta=0.05, scale=scale)
    assert np.trapezoid(np.real(y), x) == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("scale", [2.0, 10.0])
def test_ed_distribution_peak_weights_match_the_exact_eigendecomposition(
        scale):
    """The weight of each line against P(x) itself; the reviewer measured
    0.0112/0.4774/0.0112 at scale=2, i.e. the exact
    0.022329/0.955342/0.022329 times 1/scale."""
    sc = distribution_chain()
    X = sc.Sz[0] + sc.Sz[1]
    vals, wexact = exact_distribution_weights(sc, X)
    assert np.allclose(vals, [-1.0, 0.0, 1.0])
    x, y = sc.get_distribution(mode="ED", X=X, delta=0.05, scale=scale)
    x, y = np.asarray(x), np.real(np.asarray(y))
    for v, w in zip(vals, wexact):
        s = np.abs(x - v) < 0.5
        assert np.trapezoid(y[s], x[s]) == pytest.approx(w, abs=5e-4)


def test_ed_distribution_on_requested_points_is_real_and_normalized():
    """With xs= the imaginary part used to be the real part read a second
    time (max|Im y| = 5.0 at scale=2 here); for a Hermitian X in its own
    ground state the distribution is real."""
    sc = distribution_chain()
    X = sc.Sz[0] + sc.Sz[1]
    xs = np.linspace(-1.5, 1.5, 3001)
    x, y = sc.get_distribution(mode="ED", X=X, delta=0.05, scale=2.0, xs=xs)
    y = np.asarray(y)
    assert np.max(np.abs(np.imag(y))) < 1e-12
    assert np.trapezoid(np.real(y), xs) == pytest.approx(1.0, abs=1e-3)


def test_two_site_v3_chain_gets_a_normalized_distribution_too():
    """mode.py routes a 2-site itensor_version=3 chain to ED (v3's two-site
    DMRG aborts below 3 sites, and without mpscpp3 it is ED anyway), so
    this chain reached the 1/scale weight with no mode= at all."""
    sc = distribution_chain(n=2, itensor_version=3)
    assert sc.get_mode() == "ED"
    x, y = sc.get_distribution(X=sc.Sz[0], delta=0.05)
    assert np.trapezoid(np.real(y), x) == pytest.approx(1.0, abs=1e-3)


# ---------------------------------------- #5: kpm_n_scale is validated

BAD_N_SCALES = [1.5, 0.5, 2.0, 0, -1, True]


@pytest.mark.parametrize("value", BAD_N_SCALES, ids=repr)
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_non_positive_integer_kpm_n_scale_raises_on_dmrg(itensor_version,
                                                         value):
    """1.5 used to give exactly 1x and 0.5 the 16-moment floor on
    "python" while v2/v3 raised pybind's own TypeError, and 0/-1 gave 16
    moments on "python" but 1x (61) on v3. match= keeps pybind's
    "incompatible function arguments" from passing for the check."""
    sc = heisenberg(4, itensor_version)
    sc.kpm_n_scale = value
    name = [sc.Sz[0], sc.Sz[0]]
    with pytest.raises((TypeError, ValueError), match="kpm_n_scale"):
        sc.get_dynamical_correlator_moments(name=name, delta=0.1)
    with pytest.raises((TypeError, ValueError), match="kpm_n_scale"):
        sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                    delta=0.1, es=np.linspace(0, 2, 5))


@pytest.mark.parametrize("value", BAD_N_SCALES, ids=repr)
def test_non_positive_integer_kpm_n_scale_raises_on_ed(value):
    sc = heisenberg(4, "python")
    sc.kpm_n_scale = value
    for kwargs in ({"submode": "KPM"}, {}):   # KPM is also ED's default
        with pytest.raises((TypeError, ValueError), match="kpm_n_scale"):
            sc.get_dynamical_correlator(mode="ED", name=[sc.Sz[0], sc.Sz[0]],
                                        delta=0.1, es=np.linspace(0, 2, 5),
                                        **kwargs)


@pytest.mark.parametrize("itensor_version", [julia_live_param()])
def test_non_positive_integer_kpm_n_scale_raises_on_julia_live(
        itensor_version):
    sc = heisenberg(4, "python")
    sc.setup_julia()
    for value in (1.5, 0):
        sc.kpm_n_scale = value
        with pytest.raises((TypeError, ValueError), match="kpm_n_scale"):
            sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                        name=[sc.Sz[0], sc.Sz[0]],
                                        delta=0.1, es=np.linspace(0, 2, 5))


def test_kpm_n_scale_is_only_checked_where_it_is_read():
    """The ED push validates for submode="KPM" alone, mirroring the DMRG
    side, where only the KPM route reaches the check: a submode that
    never reads kpm_n_scale does not fail on it."""
    sc = heisenberg(4, "python")
    sc.kpm_n_scale = 1.5
    for submode in ("ED", "INV"):
        _x, y = sc.get_dynamical_correlator(mode="ED", submode=submode,
                                            name=[sc.Sz[0], sc.Sz[0]],
                                            delta=0.1,
                                            es=np.linspace(0.1, 2, 5))
        assert np.all(np.isfinite(np.asarray(y)))


# 61 moments at delta=0.1 on the 4-site Heisenberg chain on every route
# before the validator (the reviewer's table), and 122 at kpm_n_scale=2
@pytest.mark.parametrize("value,expected", [(1, 61), (2, 122),
                                            (np.int64(2), 122)], ids=repr)
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_integer_kpm_n_scale_keeps_its_moment_count_on_dmrg(
        itensor_version, value, expected):
    sc = heisenberg(4, itensor_version)
    sc.kpm_n_scale = value
    mus, _emin, _emax, _scale, n, _d = sc.get_dynamical_correlator_moments(
        name=[sc.Sz[0], sc.Sz[0]], delta=0.1)
    assert n == expected
    assert len(mus) == expected


@pytest.mark.parametrize("value,expected", [(1, 61), (2, 122),
                                            (np.int64(2), 122)], ids=repr)
def test_integer_kpm_n_scale_keeps_its_moment_count_on_ed(monkeypatch,
                                                          value, expected):
    seen = []
    orig = kpm.polynomials_for_broadening

    def spy(*args, **kwargs):
        n = orig(*args, **kwargs)
        seen.append(n)
        return n
    monkeypatch.setattr(kpm, "polynomials_for_broadening", spy)
    sc = heisenberg(4, "python")
    sc.kpm_n_scale = value
    sc.get_dynamical_correlator(mode="ED", submode="KPM",
                                name=[sc.Sz[0], sc.Sz[0]], delta=0.1,
                                es=np.linspace(0, 2, 5))
    assert seen == [expected]


def test_validator_accepts_integers_only():
    assert kpm.validate_kpm_n_scale(3) == 3
    assert type(kpm.validate_kpm_n_scale(np.int64(3))) is int
    for bad in (1.5, 2.0, True, False, "2", None):
        with pytest.raises(TypeError, match="kpm_n_scale"):
            kpm.validate_kpm_n_scale(bad)
    for bad in (0, -1, np.int64(0)):
        with pytest.raises(ValueError, match="kpm_n_scale"):
            kpm.validate_kpm_n_scale(bad)
    # polynomials_for_broadening runs the same check, so a caller of the
    # calibration itself cannot get the old int() rounding either
    with pytest.raises(TypeError, match="kpm_n_scale"):
        kpm.polynomials_for_broadening(3.0, 0.05, n_scale=1.5)
    assert kpm.polynomials_for_broadening(3.0, 0.05, n_scale=2) == \
        2 * kpm.polynomials_for_broadening(3.0, 0.05)
