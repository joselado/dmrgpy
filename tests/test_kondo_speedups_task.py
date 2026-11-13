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


# ---------------------------------------------------------------- item 8


def _old_projected_spectrum(ts, ops, nkry):
    """projected_spectrum as it was before items 6 and 8: the full matrices,
    every operator applied as compress_mpo returns it, G complex and the
    Ritz transform Uh @ G @ U in complex arithmetic."""
    b = [ts.project_operator(W) for W, _ in ops]
    bd = [ts.project_operator(Wd) for _, Wd in ops]
    store, n, Hq = task.block_krylov(ts, b + bd, nkry)
    ws, U = np.linalg.eigh(Hq)
    Uh = U.conj().T

    def overlaps(v):
        return task.backend.to_host(task.backend.row_overlaps(store, v, n))[:n]

    mats = []
    for (W, Wd), bk, bdk in zip(ops, b, bd):
        op = ts.operator(task.compress_mpo(W))
        M = np.zeros((n + 1, n + 1), dtype=complex)
        M[0, 0] = ts.expectation(op)
        M[1:, 0] = Uh @ overlaps(ts.flat(bk))
        M[0, 1:] = (Uh @ overlaps(ts.flat(bdk))).conj()
        G = np.zeros((n, n), dtype=complex)
        nb = ts.batch_size(op)
        for j0 in range(0, n, nb):
            j1 = min(n, j0 + nb)
            ys = ts.flat_many(ts.apply_many(op, ts.unflat_many(store[j0:j1])))
            G[:, j0:j1] = task._overlap_columns(store, ys, n)
        M[1:, 1:] = Uh @ G @ U
        mats.append(M)
    return np.concatenate([[0.0], ws]), mats


def _spin_like_ops(seed, d=3, site=2):
    """Three Hermitian one-site operators as AutoMPO writes them: two real
    ones and one that is i times a real antisymmetric matrix, like S_y."""
    from test_task_batched_spectrum import _one_site_mpo
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(d, d))
    B = rng.normal(size=(d, d))
    C = rng.normal(size=(d, d))
    mats = (A + A.T, 1j * (B - B.T), C + C.T)
    return [(_one_site_mpo(O, site, d=d),) * 2 for O in mats]


def _rel(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(a))


@pytest.mark.parametrize("real", [True, False])
def test_the_matrices_are_those_of_the_complex_products(real):
    from test_task_batched_spectrum import _tensors
    _, Ms, Ws = _tensors(11, real=real)
    ts = task.TangentSpace(Ms, Ws)
    ops = _spin_like_ops(12)
    e0, old = _old_projected_spectrum(ts, ops, 4)
    e1, new = task.projected_spectrum(ts, ops, 4)
    assert np.array_equal(e0, e1)
    if not real:  # a complex basis keeps the complex code, bit for bit
        for a, b in zip(old, new):
            assert np.array_equal(a, b)
        return
    for a, b in zip(old, new):
        assert _rel(a, b) < 1e-14
    X, Y, Z = new
    assert not np.any(X.imag) and not np.any(Z.imag) and not np.any(Y.real)
    assert np.any(old[1].real)  # the old Sy had a real part of rounding

    class Spectrum:
        pass

    ks = Spectrum()
    ks.Sx, ks.Sy, ks.Sz = new
    assert cd._spin_triplet(ks)[1]


def test_on_a_real_basis_every_operator_is_applied_in_real_arithmetic(monkeypatch):
    from test_task_batched_spectrum import _tensors
    _, Ms, Ws = _tensors(11, real=True)
    ts = task.TangentSpace(Ms, Ws)
    seen = []
    apply_many = task.TangentSpace.apply_many

    def spy(self, op, T):
        if op[0] is not self.Ws:
            seen.append(any(np.iscomplexobj(W) for W in op[0]))
        return apply_many(self, op, T)

    monkeypatch.setattr(task.TangentSpace, "apply_many", spy)
    task.projected_spectrum(ts, _spin_like_ops(12), 4)
    assert seen and not any(seen)


def test_a_phase_is_split_off_only_from_a_phase_times_a_real_operator():
    from test_task_batched_spectrum import _dense
    ops = _spin_like_ops(13)
    Wy = task.compress_mpo(ops[1][0])
    assert np.iscomplexobj(Wy[0])
    phase, R = task.real_up_to_phase(Wy)
    assert phase in (1j, -1j) and not any(np.iscomplexobj(W) for W in R)
    assert _rel(_dense(Wy), phase * _dense(R)) < 1e-15
    # Sx + Sy is not a phase times a real operator
    mixed = task.compress_mpo([a + b for a, b in zip(ops[0][0], ops[1][0])])
    phase, R = task.real_up_to_phase(mixed)
    assert phase == 1.0 and all(r is m for r, m in zip(R, mixed))


@pytest.mark.skipif(not cppext.available(3), reason="v3 extension not compiled")
def test_a_real_export_gives_the_real_form_and_the_same_curve():
    """v3 exports a real state: Sy is then imaginary with exact zeros, the
    third-order sums run in real arithmetic, and the curve is the one of
    the complex products."""
    sc = _chain(3)
    ks = TaskKondoSpectrum(sc, 3, nkry=6)
    assert cd._spin_triplet(ks)[1]
    s = sc._session
    Ms = s.mps_site_arrays(s.gs_wavefunction())
    Ws = s.mpo_site_arrays(s.hamiltonian_mpo())
    ops = [(taskdmrg._v3_arrays(s, getattr(sc, o)[3].to_terms()),) * 2
           for o in ("Sx", "Sy", "Sz")]
    _, old = _old_projected_spectrum(task.TangentSpace(Ms, Ws), ops, 6)

    class Spectrum:
        pass

    ko = Spectrum()
    ko.__dict__.update(ks.__dict__)
    ko.Sx, ko.Sy, ko.Sz = old
    assert not cd._spin_triplet(ko)[1]
    eV = np.linspace(-4e-3, 4e-3, 81)

    def curve(k):
        kw = dict(omega0=20e-3, Gamma0=5e-5)
        return (cd.second_order_dIdV(k, eV, U=0.3)
                + cd.third_order_kondo_dIdV(k, eV, -0.1, **kw)
                + cd.third_order_potential_dIdV(k, eV, -0.1, 0.3, **kw))

    assert _rel(curve(ko), curve(ks)) < 1e-14
