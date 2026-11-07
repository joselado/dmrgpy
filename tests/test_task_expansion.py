"""TASK with expanded bonds (pyitensor/task.py, expand_bonds): every bond of
the state grows by at most k directions taken from the two-site discarded
space of B|GS> and of the first Lanczos step H P B|GS>, which enlarges the
one-site tangent space without changing the state. On short chains it
reaches the error of pad_bonds at a smaller tangent space, and
expand_bonds=None changes nothing."""
import warnings

import numpy as np
import pytest

from dmrgpy import taskdmrg, ysr
from dmrgpy.pyitensor import task
from dmrgpy.ysrtk import paritystate

_ES = np.linspace(-1.2, 1.2, 25)

# Measured on d22 (run short_all, 241 bias points): the S = 1 pair is off
# ED CVM by 1.3e-11 ("python") and 1.7e-11 (v3) at every k from 1 to 16.
PAIR_K = 4
PAIR_TOL = {"python": 1.3e-10, 3: 1.7e-10}
J0_TOL = {"python": 3.2e-13, 3: 2.8e-11}


def _random_mps(dims, chis, seed=3):
    rng = np.random.default_rng(seed)
    bonds = [1] + list(chis) + [1]
    return [rng.normal(size=(bonds[l], d, bonds[l + 1]))
            for l, d in enumerate(dims)]


def _random_mpo(dims, w, seed, cplx=False):
    rng = np.random.default_rng(seed)
    ws = [1] + [w]*(len(dims) - 1) + [1]
    out = []
    for l, d in enumerate(dims):
        X = rng.normal(size=(ws[l], d, d, ws[l + 1]))
        if cplx:
            X = X + 1j*rng.normal(size=X.shape)
        out.append(X)
    return out


def _state(Ms):
    psi = np.ones((1, 1))
    for M in Ms:
        psi = np.tensordot(psi, M, axes=(1, 0)).reshape(-1, M.shape[2])
    return psi.ravel()


@pytest.mark.parametrize("k", [1, 2, 5])
def test_expansion_keeps_the_state_and_the_gauge(k):
    """A random MPS, a random Hermitian MPO and a complex random operator:
    both gauge relations hold, the state, E0 and the real dtype do not
    change, and the dimension is the one of the shapes."""
    dims, chis = [2, 4, 3, 2, 3], [1, 3, 2, 2]
    Ms = _random_mps(dims, chis)
    Ws = [0.5*(W + W.transpose(0, 2, 1, 3)) for W in _random_mpo(dims, 3, 5)]
    Wb = _random_mpo(dims, 2, 7, cplx=True)
    t0 = task.TangentSpace(Ms, Ws)
    ts = task.TangentSpace(Ms, Ws, expand_bonds=k, expand_ops=[Wb])
    assert sum(ts.added) > 0
    assert max(ts.added) <= k
    bounds = task.expansion_bounds([a.shape for a in t0.As], k)
    assert all(0 <= n <= m for n, m in zip(ts.added, bounds))
    assert ts.dimension() == task.expanded_dimension(
        [a.shape for a in t0.As], ts.added) > t0.dimension()
    assert ts.As[0].dtype == np.float64
    assert ts.e0 == pytest.approx(t0.e0, abs=1e-12)
    for A, B in zip(ts.As, ts.Bs):
        assert A.shape == B.shape
        a, d, b = A.shape
        X = A.reshape(a*d, b)
        assert np.allclose(X.conj().T @ X, np.eye(b), atol=1e-12)
        Y = B.reshape(B.shape[0], -1)
        assert np.allclose(Y @ Y.conj().T, np.eye(B.shape[0]), atol=1e-12)
    ref = _state(t0.Cs[:1] + t0.Bs[1:])
    for l in range(len(Ms)):
        got = _state(ts.As[:l] + [ts.Cs[l]] + ts.Bs[l + 1:])
        assert np.max(np.abs(got - ref)) < 1e-12
    # the expanded space holds the unexpanded one: P b is unchanged
    Wo = _random_mpo(dims, 2, 9)
    b0 = _state_of_tangent(t0, t0.project_operator(Wo))
    b1 = _state_of_tangent(ts, ts.project_operator(Wo))
    assert np.linalg.norm(b1) >= np.linalg.norm(b0) - 1e-12
    assert abs(np.vdot(b0, b1) - np.vdot(b0, b0)) < 1e-10


def _state_of_tangent(ts, T):
    """The full vector sum_l A..A T_l B..B of a tangent vector."""
    return sum(_state(ts.As[:l] + [T[l]] + ts.Bs[l + 1:])
               for l in range(ts.L))


def test_the_hamiltonian_feed_keeps_the_state_and_the_gauge():
    """The feed P^{2perp} H|GS> (off by default) on a random MPS, which is
    not an eigenstate, so the feed is not rounding: it adds directions,
    and the gauge relations, the state and E0 hold as for the default
    feeds."""
    dims, chis = [2, 4, 3, 2, 3], [1, 3, 2, 2]
    Ms = _random_mps(dims, chis)
    Ws = [0.5*(W + W.transpose(0, 2, 1, 3)) for W in _random_mpo(dims, 3, 5)]
    t0 = task.TangentSpace(Ms, Ws)
    ts = task.TangentSpace(Ms, Ws)
    added = ts.expand(2, [], feeds=("hamiltonian",))
    assert 0 < sum(added) and max(added) <= 2
    assert ts.dimension() == task.expanded_dimension(
        [a.shape for a in t0.As], added)
    assert ts.e0 == pytest.approx(t0.e0, abs=1e-12)
    for A, B in zip(ts.As, ts.Bs):
        a, d, b = A.shape
        X = A.reshape(a*d, b)
        assert np.allclose(X.conj().T @ X, np.eye(b), atol=1e-12)
        Y = B.reshape(B.shape[0], -1)
        assert np.allclose(Y @ Y.conj().T, np.eye(B.shape[0]), atol=1e-12)
    ref = _state(t0.Cs[:1] + t0.Bs[1:])
    for l in range(len(Ms)):
        got = _state(ts.As[:l] + [ts.Cs[l]] + ts.Bs[l + 1:])
        assert np.max(np.abs(got - ref)) < 1e-12


def test_no_expansion_is_bit_for_bit():
    Ms = _random_mps([2, 4, 3, 2, 3], [1, 3, 2, 2])
    dims = [2, 4, 3, 2, 3]
    Ws = [0.5*(W + W.transpose(0, 2, 1, 3)) for W in _random_mpo(dims, 3, 5)]
    Wb = _random_mpo(dims, 2, 7)
    a = task.TangentSpace(Ms, Ws)
    b = task.TangentSpace(Ms, Ws, expand_bonds=None, expand_ops=[Wb])
    for x, y in zip((a.As, a.Bs, a.Cs, a.LH, a.RH),
                    (b.As, b.Bs, b.Cs, b.LH, b.RH)):
        for u, v in zip(x, y):
            assert np.array_equal(u, v)
    assert a.e0 == b.e0
    p = task.correlator_poles(Ms, Ws, Wb, nkry=20)
    q = task.correlator_poles(Ms, Ws, Wb, nkry=20, expand_bonds=None)
    for u, v in zip(p, q):
        if u is None:
            assert v is None
        elif isinstance(u, int):
            assert u == v
        else:
            assert np.array_equal(u, v)


@pytest.mark.parametrize("bad", [0, -2, True, 2.5, "full", "4"])
def test_expand_bonds_is_validated(bad):
    """expand_bonds is None or a positive int, expand_rounds a positive int
    that needs expand_bonds."""
    with pytest.raises(ValueError, match="expand_bonds"):
        task.check_expand_bonds(bad)
    fc = ysr.kondo_chain(2, J=0.6, Delta=0.1, itensor_version="python")
    with pytest.raises(ValueError, match="expand_bonds"):
        taskdmrg.get_poles(fc, name=("Sz", "Sz"), expand_bonds=bad)
    with pytest.raises(ValueError, match="expand_rounds"):
        task.check_expansion_keywords(None, 4, bad if bad != "full" else 0)
    with pytest.raises(ValueError, match="needs expand_bonds"):
        task.check_expansion_keywords(None, None, 2)


def test_padding_and_expansion_together_are_refused():
    fc = ysr.kondo_chain(2, J=0.6, Delta=0.1, itensor_version="python")
    with pytest.raises(ValueError, match="cannot be used together"):
        taskdmrg.get_poles(fc, name=("Sz", "Sz"), pad_bonds=4,
                           expand_bonds=4)
    with pytest.raises(ValueError, match="cannot be used together"):
        task.TangentSpace(_random_mps([2, 2], [2]),
                          [np.eye(2).reshape(1, 2, 2, 1)]*2,
                          pad_bonds="full", expand_bonds=1)


def test_an_expanded_tangent_space_too_large_is_refused(monkeypatch):
    """The bound is read from the shapes before any SVD: on the random
    chain k = 1 can reach a tangent space of 60 and k = 2 one of 82, so a
    cap of 65 lets the first through and refuses the second."""
    dims = [2, 4, 3, 2, 3]
    Ms = _random_mps(dims, [1, 3, 2, 2])
    Ws = [0.5*(W + W.transpose(0, 2, 1, 3)) for W in _random_mpo(dims, 3, 5)]
    Wb = _random_mpo(dims, 2, 7)
    monkeypatch.setattr(task, "MAX_PADDED_DIMENSION", 65)
    ts = task.TangentSpace(Ms, Ws, expand_bonds=1, expand_ops=[Wb])
    assert ts.dimension() == 60
    with pytest.raises(ValueError, match="tangent space of dimension up to 82"):
        task.TangentSpace(Ms, Ws, expand_bonds=2, expand_ops=[Wb])


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_promoted_sector_state_with_expanded_bonds_against_ed(itensor_version):
    """The degenerate pair S = 1, Dz = -0.05 on 4 bath sites of
    test_task_bond_padding.py, on the sector state. On d22 it is off ED CVM
    by 1.6e-3 unexpanded, by 2.4e-8 at pad_bonds=16 (a tangent space of
    1791), and by 1.3e-11 ("python") and 1.7e-11 (v3) expanded, where two
    of the four correlators grow from 562 to 658 ("python") and 621 (v3)
    and the other two have nothing to add."""
    tol = PAIR_TOL[itensor_version]
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, S="1", Dz=-0.05,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    _, ed = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    c, _ = paritystate.parity_ground_state(fc)
    _, plain = ysr.get_stm_spectrum(c, _ES, delta=0.08, U=0.3)
    assert np.max(np.abs(plain - ed)) > 5e-4
    _, dm = ysr.get_stm_spectrum(c, _ES, delta=0.08, U=0.3,
                                 expand_bonds=PAIR_K)
    assert np.max(np.abs(dm - ed)) < tol


def test_product_state_at_zero_coupling_with_expanded_bonds():
    """J = 0: a product of impurity and bath. Unexpanded TASK keeps 1/3 of
    the weight at U = 0 and 0.405 at U = 0.3 (the S^z channel). Expanded,
    the bond between the impurity and the bath gets the flipped impurity,
    and the whole weight comes back (missing below 5e-13 of it on d22, both
    backends, every k). The poles are not those of ED: the flipped
    impurity is carried on the first bond only, while the bath electron
    moves on through the chain, and the spectrum is off ED by 0.39 against
    0.27 unexpanded (run short_all). pad_bonds="full" is the cure there."""
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 15
    kept, missing, dims, plain = 0.0, 0.0, [], []
    for X in ysr.tunneling_operators(fc, U=0.3):
        for pair in ((X.get_dagger(), X), (X, X.get_dagger())):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                q = taskdmrg.get_poles(fc, name=pair, nkry=100)
            plain.append(q["dimension"])
            p = taskdmrg.get_poles(fc, name=pair, nkry=100, expand_bonds=2)
            kept += np.sum(np.real(p["weight"]))
            missing += abs(p["missing_weight"])
            dims.append(p["dimension"])
    assert missing < 1e-10*kept
    assert all(d > d0 for d, d0 in zip(dims, plain))


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_two_rounds_give_the_product_state_spectrum(itensor_version):
    """J = 0 with expand_rounds=2: the second round starts from the
    expanded state, so b and H b carry the flipped impurity onto the next
    bond, and the spectrum is that of ED (3.2e-14 "python" and 2.8e-12 v3 on
    d22 at k = 4 to 16, run rounds, at a tangent space of 515 against 1023
    for "full"; k = 2 needs three rounds)."""
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    _, ed = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _, dm = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3,
                                     expand_bonds=8, expand_rounds=2)
    assert np.max(np.abs(dm - ed)) < J0_TOL[itensor_version]
    # the route is exact, so the J = 0 "incomplete" warning stays silent
    assert not any("J = 0 on the chain" in str(w.message) for w in caught)


def test_stm_poles_pass_the_expansion_on():
    """zaverage.stm_poles takes expand_bonds and expand_rounds as
    get_stm_spectrum does: at J = 0 two rounds recover the whole weight, so
    neither the missing-weight warning of stm_poles nor the sum-rule warning
    of get_poles fires, while the unexpanded poles miss most of it."""
    from dmrgpy.ysrtk import zaverage
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 15
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        p0 = zaverage.stm_poles(fc, U=0.3, nkry=100)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p = zaverage.stm_poles(fc, U=0.3, nkry=100, expand_bonds=8,
                               expand_rounds=2)
    total = p["missing_weight"] + np.sum(np.abs(p["weight"]))
    assert p0["missing_weight"] > 0.5*total
    assert p["missing_weight"] < 1e-8*total
    assert not any(issubclass(w.category, RuntimeWarning) for w in caught)


def test_expansion_on_the_jax_backend_matches_numpy():
    """The SVDs and the assembly are xp calls: the poles of the S = 1 pair
    with expand_bonds on JAX are the NumPy ones, and the expanded isometries
    live on the backend. At nkry = 200 the Lanczos has converged and the two
    agree to 3e-15 (JAX on the CPU, run jaxdiag_cpu); at nkry = 60 its
    unconverged part differs by rounding, 2.5e-8 (CPU) and 2.6e-6 (GPU)."""
    jax = pytest.importorskip("jax", reason="the JAX backend needs jax")
    try:
        jax.devices()
    except RuntimeError as err:  # a GPU held by another process
        pytest.skip("JAX could not initialize a device: %s" % err)
    from dmrgpy.pyitensor import backend as bk
    from dmrgpy.pyitensor.chain import _mps_arrays_lpr
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, S="1", Dz=-0.05,
                         itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 15
    c, _ = paritystate.parity_ground_state(fc)
    from dmrgpy.groundstate import ground_state_on_session
    ground_state_on_session(c)
    s = c._session
    X = ysr.tunneling_operators(c, U=0.3)[0]
    args = (_mps_arrays_lpr(s.wf0), task.mpo_arrays(s.H),
            task.mpo_arrays(s._mpo(X.get_dagger().to_terms())))
    w0, wt0, _, dim0 = task.correlator_poles(*args, nkry=200, expand_bonds=4)
    assert dim0 > task.correlator_poles(*args, nkry=2)[3]  # it expanded
    try:
        try:
            bk.set_backend("jax")
            w1, wt1, _, dim1 = task.correlator_poles(*args, nkry=200,
                                                     expand_bonds=4)
            ts = task.TangentSpace(args[0], args[1], expand_bonds=4,
                                   expand_ops=[args[2]])
        except Exception as err:  # a device that runs out of memory
            if "RESOURCE_EXHAUSTED" in str(err) or "out of memory" in str(err):
                pytest.skip("the JAX device ran out of memory: %s" % err)
            raise
        assert isinstance(ts.As[2], jax.Array)
    finally:
        bk.set_backend("numpy")
    assert dim1 == dim0
    def curve(w, wt):
        return np.real(0.08/(np.pi*((_ES[:, None] - w[None, :])**2 + 0.08**2))
                       @ wt)
    assert np.max(np.abs(curve(w1, wt1) - curve(w0, wt0))) < 1e-10


@pytest.mark.parametrize("itensor_version", ["python", 3])
def test_expansion_does_not_move_a_converged_spectrum(itensor_version):
    """The chain of test_ysr.py::test_stm_spectrum_dmrg_against_ed, whose
    tangent space already holds the whole spectrum: the two-site weights of
    both feeds are rounding (below 1e-10 of their vectors on d22, run
    j0diag2), so nothing is added and the correlator is bit for bit the
    unexpanded one."""
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    X = ysr.tunneling_operators(fc, U=0.3)[0]
    pair = (X.get_dagger(), X)
    _, A, p = taskdmrg.dynamical_correlator(fc, es=_ES, delta=0.08, name=pair,
                                            return_poles=True, variance=False)
    _, B, q = taskdmrg.dynamical_correlator(fc, es=_ES, delta=0.08, name=pair,
                                            return_poles=True, variance=False,
                                            expand_bonds=4)
    assert p["dimension"] == q["dimension"]
    assert np.array_equal(A, B)
