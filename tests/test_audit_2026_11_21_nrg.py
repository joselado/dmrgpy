"""Regressions of the 2026-11-21 hole hunt, cluster nrg (finding 22 of
docs/audit_2026_11_21_hole_hunt.md): the NRG tail and the hand-over check
of nrgmps are relative to the band D.

wilsonchain.wilson_chain is a Lanczos on levels proportional to D, so the
chain of a band [-D, D] is D times the chain of [-1, 1] and breaks at the
same site whatever D is. The tail used to stop at an absolute scale of
1e-30 (116 sites at D = 1e5 and 120 at 1e8 at Lambda = 4, into the broken
end, against 99 at D = 1), and the hand-over check compared the head and
the tail to an absolute 1e-12, so a correct hand-over raised from D of
about 1e4 on, after the DMRG solve of the head. No test here solves DMRG.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import nrgmps, ysr

LAMBDA = 4.0


def _metal(D, ns=1, n=12, **kw):
    """A spin chain on a metal whose every energy is in units of D."""
    return ysr.spin_chain_on_metal(ns, n, JK=0.3*D, Js=0.5*D, D=D,
                                   Lambda=LAMBDA, itensor_version="python",
                                   **kw)


@pytest.mark.parametrize("D", [1e-3, 1.0, 1e5, 1e8])
def test_max_sites_is_the_same_at_every_band(D):
    """The tail length does not depend on D: 99 sites at Lambda = 4, where
    the absolute floor gave 89 at D = 1e-3, 116 at 1e5 and 120 at 1e8."""
    assert nrgmps.max_sites(LAMBDA, D) == 99
    for L in (2.0, 4.0, 6.0):
        assert nrgmps.max_sites(L, D) == nrgmps.max_sites(L)


@pytest.mark.parametrize("D", [1e-3, 1e5, 1e8])
def test_the_default_tail_is_the_unit_chain_times_D(D):
    """The default tail of a band D is the tail of the band 1 times D, site
    by site, and it ends before the broken sites: |eps|/t stays at
    roundoff, where the 116-site tail at D = 1e5 reached 2.3e8."""
    e1, t1 = nrgmps.tail_chain(_metal(1.0))
    eps, t = nrgmps.tail_chain(_metal(D))
    assert len(eps) == len(e1) == nrgmps.max_sites(LAMBDA)
    assert np.max(np.abs(t/D - t1)/np.abs(t1)) < 1e-12
    assert np.max(np.abs(eps[1:])/np.abs(t)) < 1e-10


def _absolute_check(eps, t, eps_h, t_h):
    """The hand-over check before the fix, an absolute 1e-12, as a
    verdict (True: it passes)."""
    return not (len(eps) < len(eps_h) or
                np.max(np.abs(eps[:len(eps_h)] - eps_h), initial=0.0) > 1e-12
                or np.max(np.abs(t[:len(t_h)] - t_h), initial=0.0) > 1e-12)


def _passes(eps, t, eps_h, t_h):
    try:
        nrgmps._check_prefix(eps, t, eps_h, t_h)
    except ValueError:
        return False
    return True


@pytest.mark.parametrize("D", [1e-4, 1.0, 1e5, 1e8])
def test_the_handover_check_is_relative_to_the_band(D):
    """A head and a tail of the same mesh pass at every D. A tail whose
    first sites are off by 1e-14 of themselves passes and one off by 1e-9
    is refused, at every D. At D = 1 the relative check gives the verdicts
    of the absolute one it replaces; that one refused the 1e-14 tail at
    D = 1e8 and passed the 1e-9 tail at D = 1e-4."""
    fc = _metal(D)
    eh, th = nrgmps.wilson_parameters(fc)
    et, tt = nrgmps.tail_chain(fc)
    assert _passes(et, tt, eh, th)
    scale = np.max(np.abs(th))
    for rel, ok in ((1e-14, True), (1e-9, False)):
        t2 = tt.copy()
        t2[:len(th)] = th*(1.0 + rel)
        e2 = et.copy()
        e2[:len(eh)] = eh + rel*scale
        assert _passes(et, t2, eh, th) == ok
        assert _passes(e2, tt, eh, th) == ok
        if D == 1.0:
            assert _absolute_check(et, t2, eh, th) == ok
            assert _absolute_check(e2, tt, eh, th) == ok
        if D == 1e8 and ok:
            assert not _absolute_check(et, t2, eh, th)
        if D == 1e-4 and not ok:
            assert _absolute_check(et, t2, eh, th)
    with pytest.raises(ValueError, match="does not start"):
        nrgmps._check_prefix(et[:5], tt[:4], eh, th)


def test_a_refused_handover_does_not_solve_the_head(monkeypatch):
    """head_tail_poles checks k0 and the tail before the DMRG solve of the
    head, so a refused hand-over raises without solving."""
    fc = _metal(1.0, n=6)
    eps, t = nrgmps.tail_chain(fc)
    bad = t.copy()
    bad[2] *= 1.0 + 1e-6

    def solve(*args, **kwargs):
        raise AssertionError("the head was solved before the hand-over "
                             "was checked")

    monkeypatch.setattr(nrgmps, "head_state", solve)
    with pytest.raises(ValueError, match="does not start"):
        nrgmps.head_tail_poles(fc, m=2, task=False, wilson_tail=(eps, bad))
    with pytest.raises(ValueError, match="k0"):
        nrgmps.head_tail_poles(fc, k0=5, task=False)


class _Reached(Exception):
    pass


def test_the_default_tail_hands_over_at_a_large_band(monkeypatch):
    """At D = 1e8 the default tail passes the hand-over check and the call
    goes on to the hand-over shell (with the head solve replaced by the
    unsolved chain, which the check does not read). Before, the absolute
    check refused it."""
    fc = _metal(1e8, n=6)
    monkeypatch.setattr(nrgmps, "head_state",
                        lambda c, Nf=None, Sz=None: (c, 0.0))

    def shell(*args, **kwargs):
        raise _Reached()

    monkeypatch.setattr(nrgmps, "handover_shell", shell)
    with pytest.raises(_Reached):
        nrgmps.head_tail_poles(fc, m=2, task=False, symmetric=False)


@pytest.mark.parametrize("L,n_tail,warns", [(4.0, 113, True),
                                            (4.0, 99, False),
                                            (2.0, 150, False),
                                            (2.0, 200, True)])
def test_an_explicit_tail_past_the_relative_floor_warns(L, n_tail, warns):
    """An explicit n_tail warns when the scale of the chain goes below
    OMEGA_MIN D (past 99 sites at Lambda = 4 and 199 at Lambda = 2), and
    not when it only goes past N_MAX; the default never warns."""
    fc = ysr.spin_chain_on_metal(1, 6, JK=0.3, Lambda=L,
                                 itensor_version="python")
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        eps, t = nrgmps.tail_chain(fc, n_tail)
        nrgmps.tail_chain(fc)
    assert len(eps) == n_tail
    ours = [w for w in record if "n_tail=" in str(w.message)]
    assert len(ours) == (1 if warns else 0)
    if warns:
        assert issubclass(ours[0].category, RuntimeWarning)
        assert "n_tail=%d" % n_tail in str(ours[0].message)


_X = np.logspace(-10, -1, 37)
_ES = np.concatenate([-_X[::-1], _X]) # in units of D


def _exact_run(D):
    """The reviewer's chain of finding 22: 3 spins, the metal under the
    middle one, head="exact", with the stop rule off so that the NRG runs
    to the end of the tail whatever the flow does."""
    fc = ysr.spin_chain_on_metal(3, 4, site=1, Js=0.5*D, JK=0.5*D,
                                 Lambda=LAMBDA, D=D, itensor_version="python")
    _, A, res = ysr.get_stm_spectrum(fc, D*_ES, submode="NRG", head="exact",
                                     keep=128, stop_tol=None,
                                     return_poles=True)
    return np.asarray(A), res[0]


@pytest.fixture(scope="module")
def unit_band():
    return _exact_run(1.0)


@pytest.mark.parametrize("D", [1e5, 1e8])
def test_the_exact_head_is_the_same_calculation_at_every_band(unit_band, D):
    """The NRG with the exact start shell on a band D is the one on the
    band 1 with every energy times D: the same number of shells, the
    shell scales times D, the same weights, and D A(D w) = A(w), since A
    is a density. Before, the tail at D = 1e5 had 116 sites and the
    curve moved by 7.3e-3 of its peak against a clean tail."""
    A1, r1 = unit_band
    A, r = _exact_run(D)
    assert r1.nshells == r.nshells == nrgmps.max_sites(LAMBDA)
    assert r1.stop == r.stop == "chain"
    np.testing.assert_allclose(r.omega_n/D, r1.omega_n, rtol=1e-12, atol=0.0)
    assert r.total == pytest.approx(r1.total, rel=1e-10)
    assert r.sum_rule == pytest.approx(r1.sum_rule, rel=1e-10)
    assert np.max(np.abs(D*A - A1)) < 1e-5*np.max(np.abs(A1))
