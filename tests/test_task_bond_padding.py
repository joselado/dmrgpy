"""TASK with padded bonds (pyitensor/task.py, pad_bonds): every bond of the
state is padded with zero-weight directions before the recursion, which
enlarges the one-site tangent space without changing the state. The two
known failures of a small Schmidt rank, a promoted sector state of a
degenerate pair and the J = 0 product state, become exact on short chains
with pad_bonds="full", and pad_bonds=None changes nothing."""
import warnings

import numpy as np
import pytest

from dmrgpy import taskdmrg, ysr
from dmrgpy.pyitensor import task
from dmrgpy.ysrtk import paritystate, zaverage

_ES = np.linspace(-1.2, 1.2, 25)


def _random_mps(dims, chis, seed=3):
    rng = np.random.default_rng(seed)
    bonds = [1] + list(chis) + [1]
    return [rng.normal(size=(bonds[l], d, bonds[l + 1]))
            for l, d in enumerate(dims)]


def _state(Ms):
    psi = np.ones((1, 1))
    for M in Ms:
        psi = np.tensordot(psi, M, axes=(1, 0)).reshape(-1, M.shape[2])
    return psi.ravel()


@pytest.mark.parametrize("pad,expected", [("full", [2, 8, 6, 3]),
                                          (5, [2, 5, 5, 3]),
                                          (2, [2, 3, 2, 2])])
def test_padding_keeps_the_state_and_the_gauge(pad, expected):
    dims, chis = [2, 4, 3, 2, 3], [1, 3, 2, 2]
    Ms = _random_mps(dims, chis)
    assert task.structural_bounds(dims) == [2, 8, 6, 3]
    padded = task.pad_mps_bonds(Ms, pad)
    assert [M.shape[2] for M in padded[:-1]] == expected
    ref = _state(Ms)
    assert np.allclose(_state(padded), ref, atol=1e-14)
    As, Bs, Cs = task.mixed_gauge(Ms, pad_bonds=pad)
    for A, B in zip(As, Bs):
        a, d, b = A.shape
        X = A.reshape(a*d, b)
        assert np.allclose(X.conj().T @ X, np.eye(b), atol=1e-12)
        Y = B.reshape(B.shape[0], -1)
        assert np.allclose(Y @ Y.conj().T, np.eye(B.shape[0]), atol=1e-12)
    ref = ref/np.linalg.norm(ref)
    for l in range(len(Ms)):
        got = _state(As[:l] + [Cs[l]] + Bs[l + 1:])
        assert abs(abs(np.vdot(got, ref)) - 1.0) < 1e-12
    ts = task.TangentSpace(Ms, [np.eye(d).reshape(1, d, d, 1) for d in dims],
                           pad_bonds=pad)
    if pad == "full":
        assert ts.dimension() == int(np.prod(dims)) - 1


def test_no_padding_is_bit_for_bit():
    Ms = _random_mps([2, 4, 3, 2, 3], [1, 3, 2, 2])
    assert task.pad_mps_bonds(Ms, None) is Ms
    for x, y in zip(task.mixed_gauge(Ms), task.mixed_gauge(Ms, pad_bonds=None)):
        for a, b in zip(x, y):
            assert np.array_equal(a, b)


def test_a_padded_tangent_space_too_large_is_refused():
    """"full" on 21 sites of d = 4 asks for a tangent space of 4^21 - 1;
    the check reads only the shapes and raises before any array is made.
    An int pad that stays under the cap passes."""
    Ms = [np.ones((1, 4, 1)) for _ in range(21)]
    with pytest.raises(ValueError, match="tangent space of dimension %d"
                       % (4**21 - 1)):
        task.pad_mps_bonds(Ms, "full")
    with pytest.raises(ValueError, match="cap"):
        task.TangentSpace(Ms, [np.eye(4).reshape(1, 4, 4, 1)] * 21,
                          pad_bonds=1024)
    padded = task.pad_mps_bonds(Ms, 16)
    assert max(M.shape[2] for M in padded) == 16


@pytest.mark.parametrize("bad", [0, -2, True, 2.5, "Full", "max"])
def test_pad_bonds_is_validated(bad):
    with pytest.raises(ValueError, match="pad_bonds"):
        task.check_pad_bonds(bad)
    fc = ysr.kondo_chain(2, J=0.6, Delta=0.1, itensor_version="python")
    with pytest.raises(ValueError, match="pad_bonds"):
        taskdmrg.get_poles(fc, name=("Sz", "Sz"), pad_bonds=bad)


@pytest.mark.parametrize("S,Dz,n,itensor_version,before,tol",
                         [("1", -0.05, 4, "python", 8e-4, 3e-7),
                          ("3/2", -0.04, 3, 3, 4e-2, 1e-10)])
def test_promoted_sector_state_with_padded_bonds_against_ed(
        S, Dz, n, itensor_version, before, tol):
    """The degenerate pairs of test_ysr_parity.py at the chain's cutoff
    1e-12. Run on d22: without padding the tangent space is 562 (S=1) and
    205 (S=3/2) and the spectrum is off ED by 8.7e-4 and 4.8e-2; with
    "full" it is 4095 and the error 4.9e-8 ("python", S=1) and 8.9e-12
    (v3, S=3/2)."""
    fc = ysr.kondo_chain(n, J=0.6, Delta=0.1, S=S, Dz=Dz,
                         itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 64, 15
    _, ed = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    c, _ = paritystate.parity_ground_state(fc)
    _, plain = ysr.get_stm_spectrum(c, _ES, delta=0.08, U=0.3)
    assert np.max(np.abs(plain - ed)) > before
    _, dm = ysr.get_stm_spectrum(c, _ES, delta=0.08, U=0.3, pad_bonds="full")
    assert np.max(np.abs(dm - ed)) < tol
    X = ysr.tunneling_operators(c, U=0.3)[0]
    p = taskdmrg.get_poles(c, name=(X.get_dagger(), X), nkry=100,
                           pad_bonds="full")
    assert p["dimension"] == 4**6 - 1


def test_product_state_at_zero_coupling_with_full_padding():
    """J = 0: a product of impurity and bath. Without padding TASK keeps a
    part of the weight only; with "full" it is ED."""
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 15
    _, ed = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3, mode="ED",
                                 submode="CVM")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        plain = zaverage.stm_poles(fc, U=0.3, nkry=100)
    assert plain["missing_weight"] > 0.1
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        p = zaverage.stm_poles(fc, U=0.3, nkry=100, pad_bonds="full")
    assert [x for x in w if "stm_poles" in str(x.message)] == []
    assert p["missing_weight"] < 1e-10
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, dm = ysr.get_stm_spectrum(fc, _ES, delta=0.08, U=0.3,
                                     pad_bonds="full")
    assert [x for x in w if "J = 0 on the chain" in str(x.message)] == []
    assert np.max(np.abs(dm - ed)) < 1e-8


@pytest.mark.parametrize("itensor_version", ["python", 3])
@pytest.mark.parametrize("cutoff", [None, 0.0])
def test_padding_does_not_move_a_converged_spectrum(itensor_version, cutoff):
    """The chain of test_ysr.py::test_stm_spectrum_dmrg_against_ed. At the
    chain's cutoff its tangent space (515 in every run on d22) already
    holds the whole spectrum (ED to 2e-12); padded to 1023 it moves by
    6e-14 ("python") and 4e-12 (v3) on d22. At cutoff 0 the bonds reach
    their bound and padding changes nothing, bit for bit."""
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version=itensor_version)
    fc.maxm, fc.nsweeps = 60, 15
    if cutoff is not None:
        fc.cutoff = cutoff
    X = ysr.tunneling_operators(fc, U=0.3)[0]
    pair = (X.get_dagger(), X)
    a = taskdmrg.get_poles(fc, name=pair, nkry=100)
    b = taskdmrg.get_poles(fc, name=pair, nkry=100, pad_bonds="full")
    if cutoff is None:
        assert a["dimension"] < b["dimension"] == 1023
    else:
        assert (a["dimension"], b["dimension"]) == (1023, 1023)
    _, A = taskdmrg.dynamical_correlator(fc, es=_ES, delta=0.08, name=pair)
    _, B = taskdmrg.dynamical_correlator(fc, es=_ES, delta=0.08, name=pair,
                                         pad_bonds="full")
    if cutoff is None:
        assert np.max(np.abs(A - B)) < 1e-10
    else:
        assert np.array_equal(A, B)


def test_an_int_pad_inflates_the_pole_variances():
    """The variance is not an error bar under an int pad below the bound.
    On the converged 4-site chain, pad_bonds=4 leaves the spectrum within
    3e-13 of "full", but the largest variance of a pole of weight above
    1e-3 of the largest goes from 8e-21 to 1e-4 (run fixvar on d22); with
    "full" the two-site discarded space is empty and it is zero. Without
    a pad the variance is rounding, and its size depends on the random
    MPS start: 7e-21 to 3.4e-15 over seeds 0 to 7 on d22 (run
    fix_varseeds), and 1.1e-14 once in a full pytest run. So it is
    bounded at 1e-12, still six orders below the int pad."""
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 60, 15
    X = ysr.tunneling_operators(fc, U=0.3)[0]
    out = {}
    for pad in [None, 4, "full"]:
        _, A, p = taskdmrg.dynamical_correlator(fc, es=_ES, delta=0.08,
                                                name=(X.get_dagger(), X),
                                                return_poles=True,
                                                pad_bonds=pad)
        w = np.abs(p["weight"])
        out[pad] = (A, np.max(np.asarray(p["variance"])[w > 1e-3*w.max()]))
    assert np.max(np.abs(out[4][0] - out["full"][0])) < 1e-10
    assert out[None][1] < 1e-12
    assert out[4][1] > 1e-6
    assert out["full"][1] < 1e-25


def test_padding_on_the_jax_backend_matches_numpy():
    """The padding is xp.pad on the arrays of the backend: the J = 0 poles
    with "full" are the NumPy ones, and the padded isometries live on the
    backend."""
    jax = pytest.importorskip("jax", reason="the JAX backend needs jax")
    try:
        jax.devices()
    except RuntimeError as err:  # a GPU held by another process
        pytest.skip("JAX could not initialize a device: %s" % err)
    from dmrgpy.pyitensor import backend as bk
    from dmrgpy.pyitensor.chain import _mps_arrays_lpr
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 15
    fc.gs_energy()
    s = fc._session
    X = ysr.tunneling_operators(fc, U=0.3)[0]
    args = (_mps_arrays_lpr(s.wf0), task.mpo_arrays(s.H),
            task.mpo_arrays(s._mpo(X.to_terms())))
    w0, wt0, _, dim0 = task.correlator_poles(*args, nkry=60, pad_bonds="full")
    try:
        bk.set_backend("jax")
        w1, wt1, _, dim1 = task.correlator_poles(*args, nkry=60,
                                                 pad_bonds="full")
        As = task.mixed_gauge(args[0], pad_bonds="full")[0]
        assert isinstance(As[2], jax.Array)
    finally:
        bk.set_backend("numpy")
    assert dim1 == dim0 == 4**5 - 1
    def curve(w, wt):
        return np.real(0.08/(np.pi*((_ES[:, None] - w[None, :])**2 + 0.08**2))
                       @ wt)
    assert np.max(np.abs(curve(w1, wt1) - curve(w0, wt0))) < 1e-10
