"""itensor_version=3 against "python" for the public iMPS algebra of
Infinite_Many_Body_Chain (gs_method="vumps"): apply_operator (a unitary
flip and a non-unitary operator, checking .eta and vev/correlator on the
result), overlap and summps on the TFIM at D=1,2,3, and the finite window
of kpm_finite following the backend of the infinite chain. gs_method="idmrg"
on v3 must raise NotImplementedError (its cell cannot be exported)."""
import numpy as np
import pytest

from dmrgpy import cppext, infinitechain
from dmrgpy.infinitechain import _V3Uniform

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")


def _tfim(g, D, version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.etol = 1e-11
    ic.maxiter = 200
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_energy()
    return ic


def _observables(ic):
    # <Sx> is left out: at D=1 the TFIM state is a symmetry-broken product
    # state whose sign of <Sx> is a coin flip of the solver
    return np.array([ic.vev("Sz", 0)]
                    + [ic.correlator("Sz", 0, "Sz", r) for r in (1, 2, 3)]
                    + [ic.correlator("Sx", 0, "Sx", r) for r in (1, 2)])


@pytest.mark.parametrize("D", [1, 2, 3])
def test_apply_operator_unitary_flip(D):
    res = {}
    for v in (3, "python"):
        ic = _tfim(1.5, D, v)
        before = _observables(ic)
        flipped = ic.apply_operator(2 * ic.SxC[0])
        assert flipped.eta == pytest.approx(1.0, abs=1e-6)
        res[v] = (before, _observables(flipped))
        assert flipped.vev("Sz", 0) == pytest.approx(-before[0], abs=1e-6)
        with pytest.raises(RuntimeError):
            flipped.gs_energy()
    assert res[3][1] == pytest.approx(res["python"][1], abs=1e-6)


@pytest.mark.parametrize("D", [1, 2, 3])
def test_apply_operator_non_unitary(D):
    res = {}
    for v in (3, "python"):
        ic = _tfim(1.5, D, v)
        # the Sx term is sign-dependent at D=1 (see _observables)
        mix = 0.3 * (ic.SxC[0] if D > 1 else ic.get_operator("Id", 0))
        w = ic.apply_operator(ic.SzC[0] + mix)
        res[v] = (w.eta, _observables(w))
    assert res[3][0] == pytest.approx(res["python"][0], abs=1e-6)
    assert 0 < res[3][0] < 1
    assert res[3][1] == pytest.approx(res["python"][1], abs=1e-6)


def test_apply_operator_chained_on_v3():
    ic = _tfim(1.5, 2, 3)
    twice = ic.apply_operator(2 * ic.SxC[0]).apply_operator(2 * ic.SxC[0])
    assert twice.eta == pytest.approx(1.0, abs=1e-6)
    assert _observables(twice) == pytest.approx(_observables(ic), abs=1e-6)
    # the parent chain is untouched by deriving from it
    assert ic.vev("Sz", 0) == pytest.approx(twice.vev("Sz", 0), abs=1e-6)


def test_apply_operator_annihilating_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(2.0 * ic.SzC[0])  # all spins down
    ic.gs_energy()
    with pytest.raises(ValueError):
        ic.apply_operator(ic.get_operator("Sm", 0))
    with pytest.raises(ValueError):
        ic.apply_operator(ic.SxR[0])


@pytest.mark.parametrize("D", [1, 2, 3])
def test_overlap_matches_python(D):
    f = {}
    for v in (3, "python"):
        a, b = _tfim(1.5, D, v), _tfim(2.0, D, v)
        flipped = a.apply_operator(2 * a.SxC[0])
        # a and b overlap with a sign of <Sx> that is the solver's coin
        # flip at D=1, so that pair is compared from D=2 on
        f[v] = (a.overlap(a), a.overlap(flipped), a.overlap(b) if D > 1 else 0.0)
    assert f[3][0] == pytest.approx(1.0, abs=1e-8)
    assert f[3] == pytest.approx(f["python"], abs=1e-6)


def test_overlap_across_backends_raises():
    with pytest.raises(ValueError):
        _tfim(1.5, 1, 3).overlap(_tfim(1.5, 1, "python"))


def _scaled(ic, factor):
    s = ic._python_state("test")
    return ic._derived_chain(
        _V3Uniform(s.sites_uc, s.n_uc, s.D, s.d_g, s.AL * factor,
                   s.AR * factor, s.C * factor, factor**2), "test")


@pytest.mark.parametrize("D", [1, 2, 3])
def test_summps_dominant_branch(D):
    dom = _tfim(1.5, D, 3)
    other = _tfim(2.0, D, 3)
    summed = dom.summps(_scaled(other, 0.9))
    assert summed.eta == pytest.approx(1.0, abs=1e-6)
    assert _observables(summed) == pytest.approx(_observables(dom), abs=1e-6)
    with pytest.raises(Exception):
        dom.summps(other)  # two normalized states tie


def test_idmrg_on_v3_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.apply_operator(2 * ic.SxC[0])


def test_kpm_finite_window_follows_backend():
    spectra = {}
    for v in (3, "python"):
        ic = _tfim(1.5, 2, v)
        np.random.seed(1)
        es = np.linspace(-0.5, 8.0, 60)
        x, y = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=4,
                             window_chain_kwargs=dict(maxm=20, nsweeps=6),
                             delta=0.3, es=es)
        spectra[v] = np.real(np.asarray(y))
        with pytest.raises(TypeError):
            ic.kpm_finite("Sz", 0, "Sz", 0, n_window=4,
                          window_chain_kwargs=dict(itensor_version=3), es=es)
    a, b = spectra[3], spectra["python"]
    assert np.max(np.abs(a - b)) < 1e-3 * max(1.0, np.max(np.abs(b)))


def test_multisite_vumps_on_v3_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * 3, itensor_version=3)
    ic.gs_method = "vumps"
    ic.maxm = 2
    ic.set_hamiltonian(sum(ic.SzC[i] * ic.SzC[i + 1] for i in range(2))
                       + ic.SzC[2] * ic.SzR[0] + 0.5 * ic.SxC[0])
    with pytest.raises(NotImplementedError):
        ic.apply_operator(2 * ic.SxC[0])
