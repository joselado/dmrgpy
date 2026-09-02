"""Regression coverage for dmrgpy.atom's IETS entry points and for the
`dex` cutoff of the ED dynamical correlator.

Both were reported as silent-wrong-answer bugs:

* `atomtk/iets.py`'s `get_orbital_cotunneling`/`get_spinflip` accepted an
  `iorb` argument and then overwrote it with a hardcoded 0, so every call
  computed orbital 0 no matter what was requested -- for a d-shell atom
  with a crystal field that is a completely different spectrum, not a
  small perturbation.
* `edtk/dynamics.py`'s `dynamical_correlator_ED` treats every state below
  `dex` as an equally-weighted initial manifold, a hard step that makes
  the answer jump when a swept parameter moves a level across the cutoff.
  The equal-weight default is kept (changing it would silently change
  every existing caller's numbers); what is tested here is that the
  warning fires exactly in the cutoff-sensitive regime, and that the
  finite-temperature Lehmann route -- the recommended alternative --
  reproduces the degenerate-manifold average as T->0.

Everything runs in ED on a 5-orbital spinful atom (Hilbert space 4^5 =
1024), so the whole file is fast and deterministic.
"""
import warnings

import numpy as np
import pytest

from dmrgpy import atom
from dmrgpy.edtk.dynamics import check_dex_sensitivity

ORBS = ["dxy", "dyz", "dz2", "dxz", "dx2y2"]
# A *trigonal* crystal field: the a1g orbital (dz2) alone, the degenerate
# m=+-1 pair (dyz,dxz), and the degenerate m=+-2 pair (dxy,dx2y2). Some
# crystal field is essential -- with none, every orbital is equivalent and
# an ignored `iorb` is undetectable, which is exactly why the bug survived
# so long. Trigonal specifically, because the resulting pair degeneracy is
# a symmetry constraint that holds for *any* trigonal d-shell tij, giving
# a regression test that does not depend on these particular numbers
# (suggested by the caller who reported the bug, checked against their
# real VCl3 wannier tij as well as the synthetic one here).
TIJ = np.diag([0.6, 0.05, 0.0, 0.05, 0.6])
ES = np.linspace(0.0, 0.1, 151)


@pytest.fixture(scope="module")
def fc():
    return atom.generate_atom(orbs=ORBS, tij=TIJ, Ne=2, soc=0.03, U=4.0, J=0.5)


def _spinflip_reference(fc, iorb, es, delta=1e-3, dex=1e-3):
    """Hand-rolled spin-flip spectrum that respects `iorb` by construction,
    the independent reference get_spinflip is checked against."""
    total = 0.0
    for Si in [fc.Sx[iorb], fc.Sy[iorb], fc.Sz[iorb]]:
        v = fc.vev(Si, mode="ED")
        x, y = fc.get_dynamical_correlator(name=(Si - v, Si - v), es=es[es > 0.0],
                                           mode="ED", submode="ED",
                                           delta=delta, dex=dex)
        total = total + y
    return x, total


@pytest.mark.parametrize("iorb", [0, 2, 4])
def test_get_spinflip_uses_the_requested_orbital(fc, iorb):
    """get_spinflip(iorb=n) must equal a spin-flip built directly on
    orbital n -- it used to return orbital 0's spectrum for every n."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        x, y = atom.get_spinflip(fc, es=ES, iorb=iorb, delta=1e-3, dex=1e-3)
        xref, yref = _spinflip_reference(fc, iorb, ES)
    # get_spinflip returns the negative energies concatenated in front of
    # the positive ones; compare the positive half against the reference.
    assert y[len(yref):] == pytest.approx(yref, abs=1e-8)
    assert x[len(xref):] == pytest.approx(xref, abs=1e-12)


def test_spinflip_respects_the_trigonal_orbital_symmetry(fc):
    """The strongest form of the regression guard, and the one independent
    of any particular tij: in a trigonal field the two orbitals of each
    (m=+-1, m=+-2) pair are symmetry-equivalent and must give *identical*
    spectra, while inequivalent orbitals must give different ones. Code
    that ignores `iorb` returns the same spectrum for all five and so
    satisfies the equalities only trivially -- the inequalities are what
    catch it, and the equalities catch an `iorb` that is honoured but
    misrouted."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        y = {i: atom.get_spinflip(fc, es=ES, iorb=i, delta=1e-3, dex=1e-3)[1]
             for i in range(len(ORBS))}
    # degenerate pairs: identical up to ED numerics. The tolerance is
    # relative -- the spectra peak around 10, so an absolute 1e-10 would
    # be tighter than the diagonalization itself (the pairs agree to
    # ~2e-10 absolute here, i.e. ~2e-11 relative).
    assert y[0] == pytest.approx(y[4], rel=1e-8)  # dxy  == dx2y2 (m=+-2)
    assert y[1] == pytest.approx(y[3], rel=1e-8)  # dyz  == dxz   (m=+-1)
    # inequivalent orbitals: genuinely different
    assert y[0] != pytest.approx(y[1], rel=1e-3)  # m=+-2 vs m=+-1
    assert y[2] != pytest.approx(y[0], rel=1e-3)  # a1g   vs m=+-2
    assert y[2] != pytest.approx(y[1], rel=1e-3)  # a1g   vs m=+-1


@pytest.mark.parametrize("entry", [atom.get_spinflip, atom.get_orbital_cotunneling])
@pytest.mark.parametrize("iorb", [-1, 5, 9])
def test_out_of_range_orbital_raises(fc, entry, iorb):
    """An out-of-range orbital must raise rather than quietly becoming 0 --
    the failure mode this whole file exists to prevent."""
    with pytest.raises(ValueError):
        entry(fc, es=ES, iorb=iorb)


def test_orbital_cotunneling_loops_over_the_real_orbital_count(fc):
    """The output-orbital loop used to hardcode range(5); it must follow
    the chain's actual orbital count instead."""
    import inspect

    from dmrgpy.atomtk import iets
    src = inspect.getsource(iets.get_orbital_cotunneling)
    assert "range(5)" not in src
    assert "range(norb)" in src


@pytest.fixture(scope="module")
def cotunneling(fc):
    """get_orbital_cotunneling for every input orbital, computed once.

    This fixture is the point of the two tests below. Until it was added,
    nothing in the suite ever *ran* get_orbital_cotunneling: the two tests
    above call it only with an out-of-range iorb (which raises inside
    check_iorb, before a single operator is built) and read its source
    text. That is exactly how a regression reached a release with a green
    suite -- get_orbital_cotunneling is the one IETS entry point that
    builds its operators with toMPO() and passes them to
    get_dynamical_correlator's name=, and when name= resolution stopped
    accepting already-built operators the function raised ValueError on
    every call. get_spinflip passes MultiOperators and kept working, so
    the numeric tests above stayed green throughout. See
    tests/test_tompo_correlator_operators.py for the underlying contract;
    this covers the caller that depends on it.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return {i: atom.get_orbital_cotunneling(fc, es=ES, iorb=i,
                                                delta=1e-3, dex=1e-3)
                for i in range(len(ORBS))}


def test_orbital_cotunneling_runs_and_returns_a_real_spectrum(cotunneling):
    """The regression guard proper: it runs, and produces a spectrum
    rather than zeros or NaNs."""
    x, y = cotunneling[2]
    assert len(x) == len(y)
    assert np.all(np.isfinite(y))
    assert np.max(np.abs(y)) > 1e-6          # not silently empty
    # the returned grid is the negative energies concatenated in front of
    # the positive ones, so it is antisymmetric about 0 and sorted within
    # each half
    half = len(x)//2
    assert np.all(x[:half] < 0) and np.all(x[half:] > 0)
    assert np.allclose(np.sort(-x[:half]), np.sort(x[half:]))


def test_orbital_cotunneling_respects_the_trigonal_orbital_symmetry(cotunneling):
    """Same argument as the spin-flip version above, for the other entry
    point: in a trigonal field the two orbitals of each (m=+-1, m=+-2)
    pair must give identical cotunneling spectra, and inequivalent ones
    must differ. Independent of this particular tij."""
    y = {i: cotunneling[i][1] for i in cotunneling}
    assert y[0] == pytest.approx(y[4], rel=1e-8)   # dxy == dx2y2 (m=+-2)
    assert y[1] == pytest.approx(y[3], rel=1e-8)   # dyz == dxz   (m=+-1)
    assert y[2] != pytest.approx(y[0], rel=1e-3)   # a1g  vs m=+-2
    assert y[2] != pytest.approx(y[1], rel=1e-3)   # a1g  vs m=+-1


# ---------------------------------------------------------------------
# dex cutoff sensitivity
# ---------------------------------------------------------------------


def test_dex_warns_when_a_level_straddles_the_cutoff():
    """A level within a factor of a few of dex means the answer depends on
    where the cutoff was put -- warn there."""
    with pytest.warns(RuntimeWarning, match="dex"):
        check_dex_sensitivity(np.array([0.0, 0.0, 2e-3]), 1e-3)


def test_dex_is_silent_when_the_cutoff_sits_in_a_clean_gap():
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning fails the test
        check_dex_sensitivity(np.array([0.0, 0.0, 5e-2]), 1e-3)


def test_finite_T_reproduces_the_degenerate_dex_average(fc):
    """The recommended replacement for a too-coarse `dex`: at genuine
    degeneracy the Boltzmann-weighted sum must reproduce the equal-weight
    average as T->0, so switching to T is a smooth change of method and
    not a change of answer."""
    Si = fc.Sz[0]
    v = fc.vev(Si, mode="ED")
    name = (Si - v, Si - v)
    es = np.linspace(1e-4, 0.1, 101)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        _, y_dex = fc.get_dynamical_correlator(name=name, es=es, mode="ED",
                                               submode="ED", delta=1e-3, dex=1e-3)
        _, y_T = fc.get_dynamical_correlator(name=name, es=es, mode="ED",
                                             submode="ED", delta=1e-3, T=1e-6)
    assert y_T == pytest.approx(y_dex, abs=1e-8)
