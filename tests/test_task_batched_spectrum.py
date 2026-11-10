"""The speed pass of the TASK Kondo spectrum (pyitensor/task.py,
taskdmrg.session_spectrum): the tangent-space applications on a stack of
vectors, the exact compression of an operator MPO, the real Krylov basis
from starts that are real up to a phase, and the spectrum kept on the
chain under a digest of the arrays it was computed from.

Each of them must leave the numbers where they were: the first three are
checked against the vector-by-vector code on random tensors, the last one
against a fresh computation.
"""

import numpy as np
import pytest

from dmrgpy import cppext, taskdmrg
from dmrgpy.pyitensor import task

from test_kondo_spectrum_task import EV, KW, impurity_chain

VERSIONS = ["python", pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3), reason="v3 extension not compiled"))]

L, D, DW = 7, 9, 4


def _tensors(seed, d=3, real=False):
    rng = np.random.default_rng(seed)

    def rand(*shape):
        x = rng.normal(size=shape)
        return x if real else x + 1j * rng.normal(size=shape)

    dims = [1] + [min(D, d ** min(l, L - l)) for l in range(1, L)] + [1]
    Ms = [rand(dims[l], d, dims[l + 1]) for l in range(L)]
    wd = [1] + [DW] * (L - 1) + [1]
    Ws = [rand(wd[l], d, d, wd[l + 1]) for l in range(L)]
    Ws = [W + W.conj().transpose(0, 2, 1, 3) for W in Ws]  # any operator will do
    return rand, Ms, Ws


def _one_site_mpo(O, site, d=3):
    """O on one site as AutoMPO writes it, with bond dimension 2: the
    channel 0 is "not placed yet", the channel 1 is "placed"."""
    eye = np.eye(d)
    Ws = []
    for l in range(L):
        W = np.zeros((2, d, d, 2), dtype=complex)
        W[0, :, :, 0] = eye
        W[1, :, :, 1] = eye
        if l == site:
            W[0, :, :, 1] = O
        Ws.append(W)
    Ws[0] = Ws[0][:1]
    Ws[-1] = Ws[-1][:, :, :, 1:]
    return Ws


def _dense(Ws):
    """The operator of an MPO as a matrix."""
    X = Ws[0][0]  # t s w
    for W in Ws[1:]:
        X = np.tensordot(X, W, axes=(X.ndim - 1, 0))
    X = X[..., 0]
    n = X.ndim // 2
    X = X.transpose(list(range(0, 2 * n, 2)) + list(range(1, 2 * n, 2)))
    dim = int(np.prod(X.shape[:n]))
    return X.reshape(dim, dim)


@pytest.mark.parametrize("real", [False, True])
def test_a_stack_of_vectors_is_applied_as_each_one_is(real):
    rand, Ms, Ws = _tensors(0, real=real)
    ts = task.TangentSpace(Ms, Ws)
    op = ts.operator(_one_site_mpo(rand(3, 3), 2))
    V = np.array([ts.flat(ts.gauge(ts.unflat(rand(ts.flat_size()))))
                  for _ in range(5)])
    for o in (op, (ts.Ws, ts.LH, ts.RH)):
        one = np.array([ts.flat(ts.apply(o, ts.unflat(v))) for v in V])
        many = ts.flat_many(ts.apply_many(o, ts.unflat_many(V)))
        assert np.max(np.abs(one - many)) < 1e-13 * np.max(np.abs(one))
    one = np.array([ts.flat(ts.apply_h(ts.unflat(v))) for v in V])
    many = ts.flat_many(ts.apply_h_many(ts.unflat_many(V)))
    assert np.max(np.abs(one - many)) < 1e-13 * np.max(np.abs(one))


def test_the_stack_size_follows_the_memory_budget(monkeypatch):
    _, Ms, Ws = _tensors(1)
    ts = task.TangentSpace(Ms, Ws)
    large = ts.batch_size()
    monkeypatch.setattr(task, "BATCH_BYTES", 1)
    assert ts.batch_size() == 1 < large


@pytest.mark.parametrize("site", [0, 3, L - 1])
def test_a_one_site_operator_compresses_to_a_product_operator(site):
    rand, _, _ = _tensors(2, d=2)
    Wo = _one_site_mpo(rand(2, 2), site, d=2)
    Wc = task.compress_mpo(Wo)
    assert all(W.shape[0] == 1 and W.shape[3] == 1 for W in Wc)
    full = _dense(Wo)
    assert np.max(np.abs(_dense(Wc) - full)) < 1e-13 * np.max(np.abs(full))


def test_compression_keeps_an_operator_that_needs_its_bonds():
    rand, _, Ws = _tensors(3, d=2)
    Ws = [W[:, :2, :2] for W in Ws]
    Wc = task.compress_mpo(Ws)
    full = _dense(Ws)
    assert np.max(np.abs(_dense(Wc) - full)) < 1e-12 * np.max(np.abs(full))
    assert max(W.shape[3] for W in Wc) > 1


def test_starts_that_are_real_up_to_a_phase_give_a_real_basis():
    rand, Ms, Ws = _tensors(4, real=True)
    ts = task.TangentSpace(Ms, Ws)
    v = [ts.gauge(ts.unflat(rand(ts.flat_size()))) for _ in range(3)]
    phases = (1.0, 1j, np.exp(0.7j))
    starts = [[p * x for x in w] for p, w in zip(phases, v)]
    store, n, Hq = task.block_krylov(ts, starts, 6)
    assert not np.iscomplexobj(store) and not np.iscomplexobj(Hq)
    # the same space as the one grown from the real vectors themselves
    _, m, Href = task.block_krylov(ts, v, 6)
    assert n == m
    assert np.allclose(np.linalg.eigvalsh(Hq), np.linalg.eigvalsh(Href),
                       rtol=0, atol=1e-10 * np.max(np.abs(Href)))


def test_a_start_that_is_not_a_phase_times_a_real_vector_stays_complex():
    rand, Ms, Ws = _tensors(5, real=True)
    ts = task.TangentSpace(Ms, Ws)
    a, b = (ts.gauge(ts.unflat(rand(ts.flat_size()))) for _ in range(2))
    w = [x + 1j * y for x, y in zip(a, b)]
    assert task.real_starts(ts, [w])[0] is w
    store, _, _ = task.block_krylov(ts, [w], 3)
    assert np.iscomplexobj(store)


def _count_builds(monkeypatch):
    calls = []
    build = task.spectrum_from_arrays

    def counted(*args, **kwargs):
        calls.append(1)
        return build(*args, **kwargs)

    monkeypatch.setattr(task, "spectrum_from_arrays", counted)
    return calls


@pytest.mark.parametrize("version", VERSIONS)
def test_a_scan_over_the_couplings_builds_one_spectrum(version, monkeypatch):
    sc = impurity_chain(version=version)
    calls = _count_builds(monkeypatch)
    kw = dict(KW, mode="DMRG", submode="TASK", nkry=50)
    _, first = sc.get_kondo_spectrum(EV, **kw)
    for J in (-0.05, -0.2):
        _, d = sc.get_kondo_spectrum(EV, **dict(kw, Jrho_s=J))
        _, ed = sc.get_kondo_spectrum(EV, **dict(KW, Jrho_s=J, mode="ED"))
        assert np.max(np.abs(d - ed)) < 1e-9 * np.max(np.abs(ed))
    _, again = sc.get_kondo_spectrum(EV, **kw)
    assert len(calls) == 1
    assert np.array_equal(first, again)


@pytest.mark.parametrize("version", VERSIONS)
def test_another_site_setting_or_hamiltonian_builds_again(version, monkeypatch):
    sc = impurity_chain(version=version)
    calls = _count_builds(monkeypatch)
    kw = dict(KW, mode="DMRG", submode="TASK", nkry=50)
    sc.get_kondo_spectrum(EV, **kw)
    sc.get_kondo_spectrum(EV, **dict(kw, site=1))
    assert len(calls) == 2
    sc.get_kondo_spectrum(EV, **dict(kw, nkry=20))
    assert len(calls) == 3
    sc.set_hamiltonian(sc.hamiltonian + 2e-4 * sc.Sz[0])
    _, d = sc.get_kondo_spectrum(EV, **kw)
    assert len(calls) == 4
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    assert np.max(np.abs(d - ed)) < 1e-9 * np.max(np.abs(ed))


@pytest.mark.parametrize("version", VERSIONS)
def test_one_vector_per_stack_gives_the_same_spectrum(version, monkeypatch):
    """The regime of a large bond dimension, where the memory budget cuts
    a Krylov block and the basis into stacks smaller than themselves."""
    monkeypatch.setattr(task, "BATCH_BYTES", 1)
    sc = impurity_chain(version=version)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **KW)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


@pytest.mark.parametrize("version", VERSIONS)
def test_a_complex_ground_state_keeps_its_complex_basis(version, monkeypatch):
    """A field along y makes H and the ground state complex, so no phase
    can be removed from the starts and the whole route stays complex."""
    seen = []
    build = task.block_krylov

    def spy(ts, starts, nkry, **kwargs):
        out = build(ts, starts, nkry, **kwargs)
        seen.append(np.iscomplexobj(out[0]))
        return out

    monkeypatch.setattr(task, "block_krylov", spy)
    sc = impurity_chain(version=version)
    sc.set_hamiltonian(sc.hamiltonian + 3e-4 * sc.Sy[0] + 1e-4 * sc.Sx[1])
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **KW)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **KW)
    assert seen == [True]
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


def test_the_stored_spectrum_is_not_handed_out(monkeypatch):
    sc = impurity_chain()
    taskdmrg.ground_state(sc)
    ops = [(sc.Sz[0].to_terms(), sc.Sz[0].get_dagger().to_terms())]
    e, mats, _ = taskdmrg.session_spectrum(sc, ops, 10)
    e[:] = 7.0
    mats[0][:] = 7.0
    e2, mats2, _ = taskdmrg.session_spectrum(sc, ops, 10)
    assert e2[0] == 0.0 and not np.any(mats2[0] == 7.0)


def test_the_chain_keeps_a_bounded_number_of_spectra(monkeypatch):
    sc = impurity_chain()
    monkeypatch.setattr(taskdmrg, "SPECTRUM_CACHE_SIZE", 2)
    taskdmrg.ground_state(sc)
    for site in range(4):
        ops = [(sc.Sz[site].to_terms(), sc.Sz[site].get_dagger().to_terms())]
        taskdmrg.session_spectrum(sc, ops, 5)
    assert len(sc._task_spectrum_cache) == 2
