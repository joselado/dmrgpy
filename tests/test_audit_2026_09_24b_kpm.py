"""Regression tests for the kpm cluster of the 2026-09-24b audit.

Each test here locks in finding 2, 3, 4, 5 or 16 of
`docs/audit_2026_09_24b_hole_hunt.md`, which records the original symptom,
the reproduction that was executed and the reviewer's analysis.

- #2: since `765b537` `mode="ED"` KPM rescales like the DMRG routes, so
  below `kpm_scale=1/2` the ground state sits outside [-1,1], and the ED
  moment recursion had no divergence check: on the field chain below it
  returned a 0.80 spurious peak at 0.49 and an integral of 102.6 against a
  sum rule of 0.25 at 0.45. It now raises on the exact bound.
- #3: the DMRG guard `|mu_k| > 1e3*(||vi|| ||vj|| + 1)` let spectra up to
  109 times the true peak through just below 1/2, and with operators of
  norm 1e-2 it was an absolute threshold that never fired. Every backend
  now uses `1.5*||vi|| ||vj||`, and the accelerated loops check both
  moments of a step.
- #4: `mode="ED"` KPM ignored `kpm_energy_truncate`, returning the
  untruncated curve bit for bit; it now refuses it on the Hermitian KPM
  branch.
- #5: the ED `kpm_n_scale` check sat ahead of the Hermiticity branch, so a
  non-Hermitian KPM, which never reads the value, rejected 1.5 on ED
  while DMRG accepted it.
- #16: `i=`/`j=` next to an operator pair were dropped silently on every
  solver; they now raise.

The chains are 2 to 4 sites. `itensor_version=2`/`3` cases skip
themselves when the extension is not compiled, and the `julia_live` one
when there is no working Julia toolchain.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain

from _helpers import julia_live_param


needs_v2 = pytest.mark.skipif(not cppext.available(2),
                              reason="mpscpp2 not compiled")
needs_v3 = pytest.mark.skipif(not cppext.available(3),
                              reason="mpscpp3 not compiled")
DMRG_BACKENDS = ["python", pytest.param(3, marks=needs_v3),
                 pytest.param(2, marks=needs_v2)]
DIVERGING = "KPM moments diverging"
ES = np.linspace(-0.5, 4.0, 451)
DELTA = 0.1


def field_chain(itensor_version="python", field=0.2, nonherm=False, n=4):
    """Open S=1/2 Heisenberg chain plus field*Sz_0, the audit's chain.
    The field gives Sz_0|gs> an elastic weight <Sz_0>^2 = 0.0145 at E0,
    which is what grows once E0 leaves [-1,1]; without the field that
    weight is zero and the answer below 1/2 is exact."""
    np.random.seed(5)
    v = "python" if itensor_version == "julia_live" else itensor_version
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=v)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    h = h + field * sc.Sz[0]
    if nonherm:
        h = h + 0.1 * sc.Sx[n - 1] + 0.3j * sc.Sz[1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 20
    if itensor_version == "julia_live":
        sc.setup_julia()
    return sc


def kpm(sc, mode, name=None, **kwargs):
    if name is None:
        name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode=mode, submode="KPM", name=name,
                                       delta=DELTA, es=ES, **kwargs)
    return x, np.asarray(y)


class GroundStateSpy:
    """Records calls to the ED object's get_gs_array, so a test can pin
    that a precondition raised before any ground-state work."""

    def __init__(self, monkeypatch, sc):
        ed = sc.get_ED_obj()
        self.calls = []
        orig = ed.get_gs_array

        def spy(*args, **kwargs):
            self.calls.append(1)
            return orig(*args, **kwargs)
        monkeypatch.setattr(ed, "get_gs_array", spy)


# ------------------------------------------ #2: the ED moment recursion

@pytest.mark.parametrize("kpm_scale", [0.49, 0.45])
def test_ed_kpm_below_half_raises_where_it_returned_garbage(kpm_scale):
    """0.49 returned a 0.80 peak where the exact density is zero (sum rule
    0.2467), 0.45 an integral of 102.6 against 0.25."""
    sc = field_chain()
    sc.kpm_scale = kpm_scale
    with pytest.raises(RuntimeError, match=DIVERGING):
        kpm(sc, "ED")


def test_ed_kpm_below_half_still_exact_without_elastic_weight():
    """Without the field no pole outside the window carries weight, the
    moments stay at the exact bound and the answer is right, so the check
    must not reject it (an up-front raise at kpm_scale<1/2 would)."""
    sc = field_chain(field=0.0)
    sc.kpm_scale = 0.45
    x, y = kpm(sc, "ED")
    assert np.trapezoid(np.real(y), x) == pytest.approx(0.25, abs=1e-3)


def test_ed_kpm_at_the_default_scale_is_unchanged():
    sc = field_chain()
    x, y = kpm(sc, "ED")
    assert np.trapezoid(np.real(y), x) == pytest.approx(0.25, abs=1e-4)


# ------------------------------------------------ #3: the DMRG guards

@pytest.mark.parametrize("kpm_scale", [0.495, 0.49, 0.485])
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_dmrg_kpm_raises_in_the_band_the_old_guard_let_through(
        itensor_version, kpm_scale):
    """Silent on every backend before: max|y| 1.93, 5.67 and 29.5 against
    a true peak of about 0.82."""
    sc = field_chain(itensor_version)
    sc.kpm_scale = kpm_scale
    with pytest.raises(RuntimeError, match=DIVERGING):
        kpm(sc, "DMRG")


@pytest.mark.parametrize("itensor_version",
                         DMRG_BACKENDS + [julia_live_param()])
def test_dmrg_kpm_full_recursion_raises_too(itensor_version):
    """The two-vector loop at 0.49 (the accelerated one is the default for
    an auto pair). julia_live's guard is kpm.jl's `error()`, which reaches
    Python as juliacall.JuliaError rather than RuntimeError."""
    sc = field_chain(itensor_version)
    sc.kpm_scale = 0.49
    sc.kpm_accelerate = False
    expected = RuntimeError
    if itensor_version == "julia_live":
        from juliacall import JuliaError as expected
    with pytest.raises(expected, match=DIVERGING):
        kpm(sc, "DMRG")


@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_dmrg_guard_is_scale_invariant(itensor_version):
    """Operators scaled by 1e-2: ||vi|| ||vj|| = 2.5e-5, so the old +1 made
    the threshold an absolute 1e3 and the call returned 1.72e-2 against a
    sum rule of 2.5e-5."""
    sc = field_chain(itensor_version)
    sc.kpm_scale = 0.45
    A = 1e-2 * sc.Sz[0]
    with pytest.raises(RuntimeError, match=DIVERGING):
        kpm(sc, "DMRG", name=(A, A))


@pytest.mark.parametrize("kpm_scale", [0.7, 0.55])
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_dmrg_kpm_correct_runs_do_not_raise(itensor_version, kpm_scale):
    sc = field_chain(itensor_version)
    sc.kpm_scale = kpm_scale
    x, y = kpm(sc, "DMRG")
    assert np.trapezoid(np.real(y), x) == pytest.approx(0.25, abs=1e-3)


@pytest.mark.parametrize("accelerate", [True, False])
@pytest.mark.parametrize("itensor_version", DMRG_BACKENDS)
def test_dmrg_kpm_truncated_correct_runs_do_not_raise(itensor_version,
                                                      accelerate):
    """kpmmaxm=8 on 8 sites (full bond dimension 16) truncates every
    Chebyshev vector and still lands within a few per cent of ED, at a
    moment ratio of 1.0000 on v2/v3 unfixed; the cross pair (Sz_0,Sz_1)
    runs the two-vector loop, where the bound is ||vi|| ||vj|| and not
    |mu_0|. (A far harsher truncation, kpmmaxm=2 or 3 on 4 sites, drives
    the accelerated loop to a ratio of 2.6 to 10 with the spectrum 55 to
    680 per cent of its peak off, and that one raises now.)"""
    sc = field_chain(itensor_version, n=8)
    sc.kpm_scale = 0.55
    sc.kpm_accelerate = accelerate
    sc.kpmmaxm = 8
    for pair in ((0, 0), (0, 1)):
        name = (sc.Sz[pair[0]], sc.Sz[pair[1]])
        _x, y = kpm(sc, "DMRG", name=name)
        _x, ye = kpm(sc, "ED", name=name)
        err = np.max(np.abs(y - ye)) / np.max(np.abs(ye))
        assert err < 5e-2, (pair, err)


# ------------------------------- #4 and #5: the ED Hermitian KPM branch

def test_ed_kpm_refuses_energy_truncation_before_the_ground_state(
        monkeypatch):
    sc = field_chain()
    sc.kpm_energy_truncate = True
    spy = GroundStateSpy(monkeypatch, sc)
    with pytest.raises(NotImplementedError, match="kpm_energy_truncate"):
        kpm(sc, "ED")
    assert spy.calls == []


def test_two_site_v3_fallback_refuses_energy_truncation_by_name():
    """mode.py sends a 2-site v3 chain to ED even under mode="DMRG", so the
    message has to say how the call got there."""
    sc = field_chain(3, n=2)
    sc.kpm_energy_truncate = True
    with pytest.raises(NotImplementedError, match="fewer than 3"):
        kpm(sc, "DMRG")


def test_non_hermitian_ed_kpm_ignores_both_flags():
    """The non-Hermitian KPM reads neither kpm_energy_truncate nor
    kpm_n_scale (its count is n=), on either mode, so ED accepts both and
    returns the spectrum it returns without them."""
    es = np.linspace(0.0, 3.0, 13)
    kw = dict(mode="ED", submode="KPM", delta=0.2, es=es, E_max=10, n=60)
    sc = field_chain(nonherm=True)
    assert not sc.is_hermitian(sc.hamiltonian)
    _x, y1 = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **kw)
    sc = field_chain(nonherm=True)
    sc.kpm_n_scale = 1.5
    sc.kpm_energy_truncate = True
    _x, y2 = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **kw)
    np.testing.assert_allclose(np.asarray(y2), np.asarray(y1), rtol=0,
                               atol=1e-12)


@pytest.mark.parametrize("value", [1.5, 0])
def test_hermitian_ed_kpm_n_scale_raises_before_the_ground_state(
        monkeypatch, value):
    sc = field_chain()
    sc.kpm_n_scale = value
    spy = GroundStateSpy(monkeypatch, sc)
    with pytest.raises((TypeError, ValueError), match="kpm_n_scale"):
        kpm(sc, "ED")
    assert spy.calls == []


# ---------------------------------------------- #16: i=/j= with a pair

@pytest.mark.parametrize("sites", [dict(i=1), dict(j=1), dict(i=0, j=0)],
                         ids=["i=1", "j=1", "i=0,j=0"])
@pytest.mark.parametrize("mode", ["DMRG", "ED"])
def test_sites_next_to_an_operator_pair_raise(mode, sites):
    sc = field_chain()
    with pytest.raises(TypeError, match="i=/j="):
        kpm(sc, mode, **sites)


@pytest.mark.parametrize("mode", ["DMRG", "ED"])
def test_string_name_still_honours_the_sites(mode):
    sc = field_chain()
    _x, y11 = kpm(sc, mode, name=(sc.Sz[1], sc.Sz[1]))
    _x, y00 = kpm(sc, mode)
    _x, ys = kpm(sc, mode, name="ZZ", i=1, j=1)
    _x, ydef = kpm(sc, mode, name="ZZ")
    # the two sites differ on this chain (the field is on site 0), so the
    # comparison can tell them apart
    assert np.max(np.abs(y11 - y00)) > 0.1
    np.testing.assert_allclose(ys, y11, rtol=0, atol=1e-6)
    np.testing.assert_allclose(ydef, y00, rtol=0, atol=1e-6)
