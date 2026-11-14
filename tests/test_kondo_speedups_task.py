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
    np.random.seed(5)  # the random start of the solve, for one state every run
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
    # with the residual, as projected_spectrum calls it, so that Hq is the
    # same to the last bit and only the matrices are compared
    store, n, Hq, _ = task.block_krylov(ts, b + bd, nkry, residual=True)
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


# ---------------------------------------------------------------- item 9


def _old_block_krylov(ts, starts, nkry, breakdown=1e-10, residual=False):
    """task.block_krylov as it was before item 9: the H overlaps of every
    block from one product with the basis, and every Krylov candidate
    gauged once more before its two projection passes."""
    xp, backend = np, task.backend
    N = ts.flat_size()
    starts = task.real_starts(ts, starts)
    dtype = np.result_type(*[x.dtype for w in starts for x in w],
                           ts.As[0].dtype, *[W.dtype for W in ts.Ws])
    nmax = min(N, len(starts))
    store = xp.zeros((nmax, N), dtype=dtype)
    count = 0

    def add(block):
        nonlocal count, store
        for w in block:
            if count == nmax:
                break
            w = ts.flat(ts.gauge(w))
            scale = float(xp.linalg.norm(w))
            for rep in range(2):
                w = backend.project_rows(store, w, count)[1]
                w = ts.flat(ts.gauge(ts.unflat(w)))
            nw = float(xp.linalg.norm(w))
            if scale > 0.0 and nw > breakdown * scale:
                store = backend.setrow(store, count, w / nw)
                count += 1

    add(starts)
    first = store
    nmax = min(N, count * (nkry + 1))
    store = xp.zeros((nmax, N), dtype=dtype)
    for k in range(count):
        store = backend.setrow(store, k, first[k])
    done = 0
    Hu = np.zeros((nmax, nmax), dtype=complex)
    left = None
    for step in range(nkry + 1):
        n = count
        if n == done:
            break
        Hb = []
        nb = ts.batch_size()
        for j0 in range(done, n, nb):
            j1 = min(n, j0 + nb)
            hs = ts.flat_many(ts.apply_h_many(ts.unflat_many(store[j0:j1])))
            Hu[:n, j0:j1] = task._overlap_columns(store, hs, n)
            Hb.extend(ts.unflat(hs[j - j0]) for j in range(j0, j1))
        if step == nkry:
            if residual:
                out = []
                for h in Hb:
                    h = ts.flat(h)
                    for rep in range(2):
                        h = backend.project_rows(store, h, n)[1]
                    out.append(backend.to_host(h))
                out = np.array(out)
                left = (done, out.conj() @ out.T)
            done = n
            break
        done = n
        add(Hb)
    upper = np.triu(Hu[:count, :count])
    Hq = upper + np.triu(upper, 1).conj().T
    Hq = 0.5 * (Hq + Hq.conj().T)
    if not np.iscomplexobj(np.zeros(0, dtype=dtype)):
        Hq = Hq.real
    if residual:
        if left is None:
            left = (count, np.zeros((0, 0)))
        return store, count, Hq, left
    return store, count, Hq


def _starts(seed, real, k=3):
    from test_task_batched_spectrum import _tensors
    rand, Ms, Ws = _tensors(seed, real=real)
    ts = task.TangentSpace(Ms, Ws)
    return ts, [ts.gauge(ts.unflat(rand(ts.flat_size()))) for _ in range(k)]


@pytest.mark.parametrize("real", [True, False])
@pytest.mark.parametrize("nkry", [1, 4, 40])
@pytest.mark.parametrize("residual", [True, False])
def test_the_krylov_space_is_the_old_one(real, nkry, residual):
    """On a space of dimension 566, from three starts; the exhausted
    space is checked against ED in test_kondo_task_convergence.py."""
    ts, starts = _starts(21, real)
    old = _old_block_krylov(ts, starts, nkry, residual=residual)
    new = task.block_krylov(ts, starts, nkry, residual=residual)
    assert new[1] == old[1]  # the same number of vectors
    n, scale = new[1], np.max(np.abs(old[2]))
    assert np.max(np.abs(new[2] - old[2])) < 1e-13 * scale
    assert np.allclose(np.linalg.eigvalsh(new[2]), np.linalg.eigvalsh(old[2]),
                       rtol=0, atol=1e-13 * scale)
    if residual:
        assert new[3][0] == old[3][0]
        R0, R1 = old[3][1], new[3][1]
        assert R0.shape == R1.shape
        if R0.size:
            assert np.max(np.abs(R1 - R0)) < 1e-12 * np.max(np.abs(R0))


def test_the_converged_ritz_values_are_the_old_ones():
    """Once Ritz pairs converge, the next Krylov vectors are set by what
    rounding leaves of the residual, so the entries of Hq differ from the
    old code by much more than rounding (8.7e-3 of max|Hq| on the 20-site
    chain of the plan at nkry=100, and 2.7e-2 between the old code and
    itself with its starts scaled by 1+1e-15). What must not move is the
    converged spectrum."""
    sc = _chain()
    s = sc._session
    from dmrgpy.pyitensor.chain import _mps_arrays_lpr
    ts = task.TangentSpace(_mps_arrays_lpr(s.wf0), task.mpo_arrays(s.H))
    starts = [ts.project_operator(task.mpo_arrays(s._mpo(getattr(sc, o)[3].to_terms())))
              for o in ("Sx", "Sy", "Sz")]
    _, n0, H0, (j0, R0) = _old_block_krylov(ts, starts, 30, residual=True)
    _, n1, H1, _ = task.block_krylov(ts, starts, 30, residual=True)
    assert n1 == n0
    w0, U0 = np.linalg.eigh(H0)
    w1 = np.linalg.eigvalsh(H1)
    r0 = np.sqrt(np.abs(np.einsum("jf,jk,kf->f", U0[j0:].conj(), R0, U0[j0:])))
    conv = r0 < 1e-8 * np.max(np.abs(w0))
    assert conv.sum() > 5
    near = np.array([np.min(np.abs(w1 - x)) for x in w0[conv]])
    assert np.max(near) < 1e-13 * np.max(np.abs(H0))


def test_a_krylov_candidate_is_gauged_twice_and_not_overlapped_again(monkeypatch):
    ts, starts = _starts(22, True)
    counts = {"gauge": 0, "overlap": 0}
    gauge, overlap = task.TangentSpace.gauge, task._overlap_columns

    def g(self, X):
        counts["gauge"] += 1
        return gauge(self, X)

    def o(*a):
        counts["overlap"] += 1
        return overlap(*a)

    monkeypatch.setattr(task.TangentSpace, "gauge", g)
    monkeypatch.setattr(task, "_overlap_columns", o)
    store, n, Hq, _ = task.block_krylov(ts, starts, 4, residual=True)
    assert n == 5 * len(starts)  # every candidate of these tensors is kept
    assert counts["overlap"] == 0
    assert counts["gauge"] == 3 * len(starts) + 2 * (n - len(starts))


# ---------------------------------------------------------------- item 7


def _basis_dtypes(monkeypatch):
    seen = []
    build = task.block_krylov

    def spy(ts, starts, nkry, **kwargs):
        out = build(ts, starts, nkry, **kwargs)
        seen.append(np.iscomplexobj(out[0]))
        return out

    monkeypatch.setattr(task, "block_krylov", spy)
    return seen


def _heisenberg(n, field_y=0.):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + 1e-3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(h + 2e-4 * sc.Sz[0] + field_y * sc.Sy[1])
    return sc


def test_the_real_mpo_is_the_same_operator_and_a_complex_one_is_refused():
    from test_task_batched_spectrum import _dense
    sc = _heisenberg(5)
    sc.gs_energy()
    W = task.mpo_arrays(sc._session.H)
    assert any(np.any(w.imag) for w in W)  # S_y S_y on "python"
    Wr = task.real_mpo(W)
    assert Wr is not None and not any(np.iscomplexobj(w) for w in Wr)
    assert all(a.shape[3] <= b.shape[3] for a, b in zip(Wr, W))
    assert _rel(_dense(W), _dense(Wr)) < 1e-13
    sy = _heisenberg(5, field_y=3e-4)
    sy.gs_energy()
    assert task.real_mpo(task.mpo_arrays(sy._session.H)) is None


def test_a_python_state_of_a_real_h_runs_in_real_arithmetic(monkeypatch):
    sc = _chain()
    T = _cold_T(sc, 3, 200)
    seen = _basis_dtypes(monkeypatch)
    sc.__dict__.pop("_task_spectrum_cache", None)
    real = TaskKondoSpectrum(sc, 3, nkry=200, T=T)
    monkeypatch.setattr(taskdmrg, "REAL_GAUGE", False)
    cplx = TaskKondoSpectrum(sc, 3, nkry=200, T=T)
    assert seen == [False, True]
    assert 0. <= real.imaginary_dropped < taskdmrg.REAL_GAUGE_DROP
    assert cplx.imaginary_dropped is None
    assert cd._spin_triplet(real)[1] and not cd._spin_triplet(cplx)[1]
    width = np.max(np.abs(cplx.e))
    conv = cplx.residual < 1e-8 * width
    conv[0] = False  # the ground state itself
    assert conv.sum() > 3
    near = np.array([np.min(np.abs(real.e - x)) for x in cplx.e[conv]])
    assert np.max(near) < 1e-12 * width
    eV = np.linspace(-4e-3, 4e-3, 81)

    def curve(k):
        kw = dict(omega0=20e-3, Gamma0=5e-5)
        return (cd.second_order_dIdV(k, eV, U=0.3)
                + cd.third_order_kondo_dIdV(k, eV, -0.1, **kw)
                + cd.third_order_potential_dIdV(k, eV, -0.1, 0.3, **kw))

    assert _rel(curve(cplx), curve(real)) < 1e-6


def test_a_complex_hamiltonian_keeps_the_complex_route_bit_for_bit(monkeypatch):
    sc = _heisenberg(6, field_y=3e-4)
    seen = _basis_dtypes(monkeypatch)
    a = TaskKondoSpectrum(sc, 0, nkry=5)
    sc.__dict__.pop("_task_spectrum_cache", None)
    monkeypatch.setattr(taskdmrg, "REAL_GAUGE", False)
    b = TaskKondoSpectrum(sc, 0, nkry=5)
    assert seen == [True, True] and a.imaginary_dropped is None
    for x, y in zip((a.e, a.Sx, a.Sy, a.Sz), (b.e, b.Sx, b.Sy, b.Sz)):
        assert np.array_equal(x, y)


def test_a_complex_mixture_of_a_doublet_keeps_the_complex_route_bit_for_bit(
        monkeypatch):
    """A member of a degenerate manifold is in general a complex mixture
    of two real states, which no real gauge reaches."""
    sc = _heisenberg(5)
    sc.set_hamiltonian(sc.hamiltonian - 2e-4 * sc.Sz[0])  # a Kramers doublet
    sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 0.
    np.random.seed(7)  # one pair of members every run
    _, (m1, m2) = sc.get_excited_states(n=2, purify=True)
    w = m1 + 1j * m2
    sc.set_gs(w)
    seen = _basis_dtypes(monkeypatch)
    a = TaskKondoSpectrum(sc, 0, nkry=5)
    sc.__dict__.pop("_task_spectrum_cache", None)
    monkeypatch.setattr(taskdmrg, "REAL_GAUGE", False)
    b = TaskKondoSpectrum(sc, 0, nkry=5)
    assert seen == [True, True] and a.imaginary_dropped is None
    for x, y in zip((a.e, a.Sx, a.Sy, a.Sz), (b.e, b.Sx, b.Sy, b.Sz)):
        assert np.array_equal(x, y)


def test_the_real_gauge_measures_the_part_it_drops():
    """A real state in a complex gauge with a global phase is rotated back
    with nothing dropped; one with an imaginary part of 1e-3 is kept in
    its complex gauge."""
    from test_task_batched_spectrum import _tensors
    rng = np.random.default_rng(31)
    _, Ms, Ws = _tensors(30, real=True)
    cplx = [M.astype(complex) for M in Ms]
    for l in range(len(cplx) - 1):  # a complex gauge on every bond
        b = cplx[l].shape[2]
        X = rng.normal(size=(b, b)) + 1j * rng.normal(size=(b, b)) + 3 * np.eye(b)
        cplx[l] = np.tensordot(cplx[l], X, axes=(2, 0))
        cplx[l + 1] = np.tensordot(np.linalg.inv(X), cplx[l + 1], axes=(1, 0))
    cplx[0] = cplx[0] * np.exp(0.4j)
    Mr, Wr, dropped = taskdmrg.real_gauge_state(cplx, Ws)
    assert not any(np.iscomplexobj(M) for M in Mr) and dropped < 1e-7
    o = taskdmrg._mps_overlap
    assert abs(abs(o(Mr, Ms)) ** 2 / (o(Mr, Mr).real * o(Ms, Ms).real) - 1) < 1e-13
    bent = [M.copy() for M in cplx]
    bent[3] = bent[3] + 1e-3j * np.linalg.norm(bent[3]) * rng.normal(size=bent[3].shape) / np.sqrt(bent[3].size)
    Mb, Wb, gone = taskdmrg.real_gauge_state(bent, Ws)
    assert gone is None and all(a is b for a, b in zip(Mb, bent))


# ---------------------------------------------------------------- item 10


def _doublet(n=5, version="python"):
    """An odd Heisenberg chain, a Kramers doublet, at a bond dimension that
    holds it."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + 1e-3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 0.
    return sc


KW10 = dict(site=0, Jrho_s=-0.1, U=0.2, T=0., omega0=20e-3, Gamma0=5e-5,
            order=3, mode="DMRG", submode="TASK", nkry=10, n_gs=2)
EV10 = np.linspace(-4e-3, 4e-3, 41)


def _count_solves(sc):
    calls = []
    solve = sc.get_excited_states

    def counted(*a, **k):
        calls.append(1)
        return solve(*a, **k)

    sc.get_excited_states = counted
    return calls


def _call(sc, **kw):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sc.get_kondo_spectrum(EV10, **dict(KW10, **kw))[1]


@pytest.mark.parametrize("version", VERSIONS)
def test_a_repeat_call_reuses_the_manifold_and_its_spectra(monkeypatch, version):
    sc = _doublet(version=version)
    solves, builds = _count_solves(sc), _count_builds(monkeypatch)
    first = _call(sc)
    again = _call(sc, U=0.3)  # a scan over U at the same state
    assert len(solves) == 1 and len(builds) == 2
    assert np.array_equal(_call(sc), first)
    assert len(solves) == 1 and len(builds) == 2
    # the old behaviour, a solve and two builds per call, is what dropping
    # the entry gives back
    del sc._n_gs_manifold
    _call(sc)
    assert len(solves) == 2 and len(builds) == 4
    assert np.all(np.isfinite(again))


def test_a_new_hamiltonian_solves_the_manifold_again():
    sc = _doublet()
    solves = _count_solves(sc)
    _call(sc)
    sc.set_hamiltonian(sc.hamiltonian + 1e-5 * sc.Sx[1] * sc.Sx[2])
    _call(sc)
    assert len(solves) == 2


def test_other_solver_settings_solve_the_manifold_again():
    sc = _doublet()
    solves = _count_solves(sc)
    _call(sc)
    sc.maxm = 12
    _call(sc)
    assert len(solves) == 2


def test_the_key_holds_the_hamiltonian_and_the_settings():
    """set_hamiltonian() drops the entry, and the key does not rely on it:
    the same state under another H, n_gs or solver key is another key."""
    from dmrgpy.kondospectrumtk import kondomethod
    sc = _doublet()
    gs0 = sc.get_gs()
    key = kondomethod._manifold_key(sc, gs0, 2)
    assert kondomethod._manifold_key(sc, gs0, 2) == key
    assert kondomethod._manifold_key(sc, gs0, 3) != key
    sc.excited_gram_schmidt = True
    assert kondomethod._manifold_key(sc, gs0, 2) != key
    sc.excited_gram_schmidt = False
    sc.set_hamiltonian(sc.hamiltonian + 1e-5 * sc.Sz[1])
    sc.set_gs(gs0)
    sc.get_gs()
    assert kondomethod._manifold_key(sc, gs0, 2) != key


def test_a_copy_of_the_chain_solves_its_own_manifold():
    import copy
    sc = _doublet()
    first = _call(sc)
    cp = copy.deepcopy(sc)
    assert cp._n_gs_manifold is None
    solves = _count_solves(sc)
    out = _call(cp)
    assert np.max(np.abs(out - first)) < 1e-6 * np.max(np.abs(first))
    _call(sc)
    assert not solves  # the original still holds its own


# ---------------------------------------------------------------- the cache


def _raw_builds(monkeypatch):
    raw = []
    build = task.spectrum_from_arrays

    def spy(*args, **kwargs):
        out = build(*args, **kwargs)
        raw.append(out)
        return out

    monkeypatch.setattr(task, "spectrum_from_arrays", spy)
    return raw


def test_the_cache_keeps_a_real_basis_in_real_storage(monkeypatch):
    sc = _chain()
    raw = _raw_builds(monkeypatch)
    for k in range(2):  # built, then served
        ks = TaskKondoSpectrum(sc, 3, nkry=6)
        (entry,) = sc._task_spectrum_cache.values()
        assert [m[0] for m in entry[1]] == ["real", "imag", "real"]
        assert not any(np.iscomplexobj(m[1]) for m in entry[1])
        for M, R in zip((ks.Sx, ks.Sy, ks.Sz), raw[0][1]):
            assert np.iscomplexobj(M) and np.array_equal(M, R)
    assert len(raw) == 1
    assert 2 * sum(m[1].nbytes for m in entry[1]) == sum(R.nbytes for R in raw[0][1])


def test_the_cache_keeps_a_complex_basis_as_it_is(monkeypatch):
    sc = _heisenberg(6, field_y=3e-4)
    raw = _raw_builds(monkeypatch)
    ks = TaskKondoSpectrum(sc, 0, nkry=5)
    (entry,) = sc._task_spectrum_cache.values()
    assert [m[0] for m in entry[1]] == ["as is"] * 3
    for M, R in zip((ks.Sx, ks.Sy, ks.Sz), raw[0][1]):
        assert np.array_equal(M, R) and M is not R
