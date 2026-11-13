"""The speedups of the TASK Kondo route (docs/kondo_speedup_plan.md, items 6
to 10, and the GPU rule on the residual of the block Krylov space). Every
test pins the numbers against the form of the code before the change, kept
here or reached through a switch, and the saving, as a count of what is no
longer built or solved, so that undoing the saving fails a test."""
import numpy as np
import pytest

from dmrgpy import cppext, spinchain, taskdmrg
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
from dmrgpy.pyitensor import task

from test_kondo_spectrum_task import EV, KW

VERSIONS = ["python", pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3), reason="v3 extension not compiled"))]
KB = 8.617333262e-5


def _chain(version="python", n=10, maxm=12):
    """A disordered XXZ chain in a staggered field, whose Krylov space is
    far from exhausted at the nkry used here."""
    rng = np.random.default_rng(5)
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    J = 1e-3 * rng.uniform(0.5, 1.5, n - 1)
    h = 0
    for i in range(n - 1):
        h = h + J[i] * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + 1.5 * sc.Sz[i] * sc.Sz[i + 1])
    for i in range(n):
        h = h + 6e-4 * (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 10, 1e-10
    sc.gs_energy()
    return sc


def _cold_T(sc, site, nkry):
    """A temperature at which the lowest pole the spin reaches has a
    Boltzmann factor of exp(-30)."""
    ks = TaskKondoSpectrum(sc, site, nkry=nkry)
    return np.min(ks.e[ks.e > 0.]) / (30. * KB)


def _count_builds(monkeypatch):
    calls = []
    build = task.spectrum_from_arrays

    def counted(*args, **kwargs):
        calls.append(kwargs.get("column", False))
        return build(*args, **kwargs)

    monkeypatch.setattr(task, "spectrum_from_arrays", counted)
    return calls


# ---------------------------------------------------------------- item 6


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("cold", [False, True])
def test_order_two_is_the_column_of_the_full_build_bit_for_bit(version, cold):
    sc = _chain(version)
    nkry = 6
    T = _cold_T(sc, 3, nkry) if cold else 0.
    full = TaskKondoSpectrum(sc, 3, nkry=nkry, T=T)
    col = TaskKondoSpectrum(sc, 3, nkry=nkry, T=T, order=2)
    assert col.column_only and not full.column_only
    assert full.tangent_dimension > full.dim  # the space is not exhausted
    assert np.array_equal(col.e, full.e)
    assert np.array_equal(col.residual, full.residual)
    for a, b in zip((col.Sx, col.Sy, col.Sz), (full.Sx, full.Sy, full.Sz)):
        assert np.array_equal(a[:, 0], b[:, 0])
        assert np.array_equal(a[0, :], b[0, :])
        assert not np.any(a[1:, 1:]) and np.any(b[1:, 1:])
    eV = np.linspace(-4e-3, 4e-3, 81)
    assert np.array_equal(cd.second_order_dIdV(col, eV, U=0.3),
                          cd.second_order_dIdV(full, eV, U=0.3))


def test_order_two_from_get_kondo_spectrum_is_the_old_route():
    """The old route built the full matrices for order=2 as well."""
    sc = _chain()
    eV = np.linspace(-4e-3, 4e-3, 81)
    kw = dict(site=3, Jrho_s=-0.1, U=0.3, T=0., omega0=20e-3, Gamma0=5e-5)
    _, new = sc.get_kondo_spectrum(eV, mode="DMRG", submode="TASK", nkry=6,
                                   order=2, **kw)
    old = cd.second_order_dIdV(TaskKondoSpectrum(sc, 3, nkry=6), eV, U=0.3)
    assert np.array_equal(new, old)


def test_order_two_applies_no_spin_operator_to_the_basis(monkeypatch):
    sc = _chain()
    applied = []
    apply_many = task.TangentSpace.apply_many

    def spy(self, op, T):
        applied.append(op[0] is self.Ws)
        return apply_many(self, op, T)

    monkeypatch.setattr(task.TangentSpace, "apply_many", spy)
    TaskKondoSpectrum(sc, 3, nkry=6, order=2)
    assert applied and all(applied)  # H only
    applied.clear()
    TaskKondoSpectrum(sc, 3, nkry=6)
    assert not all(applied)  # and the spin operators for the block


def test_the_column_and_the_full_build_are_kept_apart(monkeypatch):
    sc = _chain()
    calls = _count_builds(monkeypatch)
    col = TaskKondoSpectrum(sc, 3, nkry=6, order=2)
    full = TaskKondoSpectrum(sc, 3, nkry=6)
    col2 = TaskKondoSpectrum(sc, 3, nkry=6, order=2)
    full2 = TaskKondoSpectrum(sc, 3, nkry=6)
    assert calls == [True, False]  # one build of each, then both served
    assert np.any(full2.Sx[1:, 1:]) and not np.any(col2.Sx[1:, 1:])
    for a, b in ((col, col2), (full, full2)):
        for x, y in zip((a.Sx, a.Sy, a.Sz), (b.Sx, b.Sy, b.Sz)):
            assert np.array_equal(x, y)


def test_the_third_order_sums_refuse_a_column():
    sc = _chain()
    col = TaskKondoSpectrum(sc, 3, nkry=6, order=2)
    with pytest.raises(ValueError, match="column_only"):
        cd.third_order_kondo_dIdV(col, EV, -0.1)
    with pytest.raises(ValueError, match="column_only"):
        cd.third_order_potential_dIdV(col, EV, -0.1, 0.3)
    with pytest.raises(ValueError, match="order must be 2 or 3"):
        TaskKondoSpectrum(sc, 3, nkry=6, order=1)
