"""submode="TASK", the tangent-space Krylov correlator (arXiv:2510.07279).

When the ground-state bond dimension saturates the Hilbert space the
one-site tangent space is the whole orthogonal complement of the ground
state, so the projection TaSK makes is exact and the poles must reproduce
the Lehmann sum to Krylov precision, with every two-site variance zero.
Below saturation the variance must become nonzero, which is what makes it
an error estimate. The saturated chains set cutoff=0: the default cutoff
can drop a Schmidt state (measured on a 10-site Heisenberg chain, 31 of 32
on the middle bond), and the tangent space then misses a direction.

Every test runs on itensor_version="python" and 3; v3 exports its site
tensors to the same Python implementation (taskdmrg.session_poles),
and so does julia_live (mpsjulialive/arrays.jl), tested at the end of
this file through taskdmrg directly.
"""

import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, spinchain

from _helpers import julia_live_param

VERSIONS = ["python", pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3), reason="v3 extension not compiled"))]


def heisenberg(n=6, hz=0.3, maxm=64, version="python"):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
    h = h + hz * sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, 0.
    return sc


def complex_hopping_chain(n=5, seed=3, version="python"):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff = 64, 20, 0.
    return fc


def lehmann_curve(chain, A, B, es, delta):
    ed = chain.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    e, U = np.linalg.eigh(h)
    Ae = U.conj().T @ np.array(ed.MO2matrix(A).todense()) @ U
    Be = U.conj().T @ np.array(ed.MO2matrix(B).todense()) @ U
    M, D = Ae[0, :] * Be[:, 0], e - e[0]
    lor = delta / (np.pi * ((es[:, None] - D[None, :]) ** 2 + delta ** 2))
    return lor @ M


ES = np.linspace(-0.5, 5., 300)


@pytest.mark.parametrize("version", VERSIONS)
def test_saturated_hermitian_pair_matches_lehmann_and_has_zero_variance(version):
    sc = heisenberg(version=version)
    A = sc.Sz[2]
    _, y = sc.get_dynamical_correlator(submode="TASK", name=(A, A), es=ES,
                                       delta=0.05)
    ref = lehmann_curve(sc, A, A, ES, 0.05)
    assert np.max(np.abs(y - ref)) < 1e-10
    _, _, p = sc.get_dynamical_correlator(submode="TASK", name=(A, A),
                                          es=ES, delta=0.05, return_poles=True)
    assert p["dimension"] == 2 ** 6 - 1
    assert np.max(p["variance"]) < 1e-20


@pytest.mark.parametrize("version", VERSIONS)
def test_elastic_term_carries_the_ground_state_weight(version):
    """<Sz_0> != 0 in the field, so the n=GS term at w=0 is nonzero and
    the sum rule closes only with it."""
    sc = heisenberg(hz=0.6, version=version)
    from dmrgpy import taskdmrg
    p = taskdmrg.get_poles(sc, name=(sc.Sz[0], sc.Sz[0]))
    assert abs(p["elastic"]) > 1e-3
    assert abs(np.sum(p["weight"]) + p["elastic"] - 0.25) < 1e-10


@pytest.mark.parametrize("version", VERSIONS)
def test_saturated_complex_off_diagonal_pair_matches_lehmann(version):
    fc = complex_hopping_chain(version=version)
    A, B = fc.Cdag[0], fc.C[2]
    _, y = fc.get_dynamical_correlator(submode="TASK", name=(A, B), es=ES,
                                       delta=0.1)
    ref = lehmann_curve(fc, A, B, ES, 0.1)
    assert np.max(np.abs(ref.imag)) > 1e-2  # a genuinely complex M_n
    assert np.max(np.abs(y - ref)) < 1e-8


@pytest.mark.parametrize("version", VERSIONS)
def test_truncated_bond_dimension_shows_up_in_the_variance(version):
    sc = heisenberg(n=8, maxm=4, version=version)
    from dmrgpy import taskdmrg
    p = taskdmrg.get_poles(sc, name=(sc.Sz[3], sc.Sz[3]), variance=True)
    assert p["dimension"] < 2 ** 8 - 1
    big = np.abs(p["weight"]) > 1e-3
    assert np.all(p["variance"][big] > 1e-8)


def test_v3_conserved_sector_matches_the_dense_chain():
    """A conserved-sector state exports with its QN blocks scattered into
    dense site tensors, and the Sz correlator stays in the sector."""
    if not cppext.available(3):
        pytest.skip("v3 extension not compiled")
    dense = heisenberg(n=8, hz=0., version=3)
    sector = heisenberg(n=8, hz=0., version=3)
    sector.set_conserved_sector(Sz=0)
    A = dense.Sz[3]
    _, yd = dense.get_dynamical_correlator(submode="TASK", name=(A, A),
                                           es=ES, delta=0.05)
    _, ys = sector.get_dynamical_correlator(submode="TASK",
            name=(sector.Sz[3], sector.Sz[3]), es=ES, delta=0.05)
    assert np.max(np.abs(ys - yd)) < 1e-8 * np.max(np.abs(yd))


@pytest.mark.parametrize("version", [2])
def test_other_backends_raise_before_solving(version):
    if version == 2 and not cppext.available(2):
        pytest.skip("extension not compiled")
    sc = spinchain.Spin_Chain(["S=1/2"] * 4, itensor_version=version)
    sc.set_hamiltonian(sc.Sz[0] * sc.Sz[1])
    with pytest.raises(NotImplementedError, match="python"):
        sc.get_dynamical_correlator(submode="TASK", name="ZZ")
    assert not sc.computed_gs


# ------------------------------------------------------------ julia_live

JULIA = [julia_live_param()]


def _pole_curve(p, es, delta):
    ws = np.concatenate([p["omega"], [0.]])
    wt = np.concatenate([p["weight"], [p["elastic"]]])
    return (delta / (np.pi * ((es[:, None] - ws[None, :]) ** 2 + delta ** 2))) @ wt


@pytest.mark.parametrize("version", JULIA)
def test_julia_saturated_poles_match_python_and_lehmann(version):
    """A saturated 6-site chain: the tangent space is the whole orthogonal
    complement on both backends, so the julia_live poles are the python
    ones and the Lehmann ones, for a Hermitian pair, for (S+, S-) and for
    an off-diagonal pair, with the sum rule closed."""
    from dmrgpy import taskdmrg
    jl = heisenberg(version=version)
    py = heisenberg(version="python")
    for pair in (lambda c: (c.Sz[0], c.Sz[0]), lambda c: (c.Sx[0] + 1j * c.Sy[0], c.Sx[0] - 1j * c.Sy[0]),
                 lambda c: (c.Sz[1], c.Sz[4])):
        pj = taskdmrg.get_poles(jl, name=pair(jl), variance=True)
        pp = taskdmrg.get_poles(py, name=pair(py))
        assert pj["dimension"] == 2 ** 6 - 1
        assert abs(pj["missing_weight"]) < 1e-6
        assert np.max(pj["variance"]) < 1e-16
        yj, yp = _pole_curve(pj, ES, 0.05), _pole_curve(pp, ES, 0.05)
        ref = lehmann_curve(py, *pair(py), ES, 0.05)
        assert np.max(np.abs(yj - yp)) < 1e-8
        assert np.max(np.abs(yj - ref)) < 1e-8
        wj = np.sort(pj["omega"][np.abs(pj["weight"]) > 1e-6])
        wp = np.sort(pp["omega"][np.abs(pp["weight"]) > 1e-6])
        assert len(wj) == len(wp) and np.allclose(wj, wp, atol=1e-8)
    _, y = taskdmrg.dynamical_correlator(jl, name=(jl.Sz[2], jl.Sz[2]),
                                         es=ES, delta=0.05)
    ref = lehmann_curve(py, py.Sz[2], py.Sz[2], ES, 0.05)
    assert np.max(np.abs(y - ref)) < 1e-8


@pytest.mark.parametrize("version", JULIA)
def test_julia_fermion_pair_matches_lehmann(version):
    """A pair of hoppings with complex amplitudes on a saturated chain:
    the Jordan-Wigner strings of the Julia MPOs of B and A^dagger, and a
    complex M_n. Parity-even, since julia_live builds no MPO of a single
    C yet."""
    from dmrgpy import taskdmrg
    fc = complex_hopping_chain(version=version)
    A, B = fc.Cdag[0] * fc.C[3], fc.Cdag[1] * fc.C[4]
    _, y = taskdmrg.dynamical_correlator(fc, name=(A, B), es=ES, delta=0.1)
    ref = lehmann_curve(fc, A, B, ES, 0.1)
    assert np.max(np.abs(ref.imag)) > 1e-2
    assert np.max(np.abs(y - ref)) < 1e-8


@pytest.mark.parametrize("version", JULIA)
def test_julia_padding_and_expansion(version):
    """Below saturation (maxm=2 on 6 sites) pad_bonds="full" closes the
    sum rule of the state and expand_bonds grows the tangent space, as on
    python."""
    from dmrgpy import taskdmrg
    sc = heisenberg(maxm=2, version=version)
    A = sc.Sz[0]
    p0 = taskdmrg.get_poles(sc, name=(A, A))
    pf = taskdmrg.get_poles(sc, name=(A, A), pad_bonds="full")
    pe = taskdmrg.get_poles(sc, name=(A, A), expand_bonds=2)
    assert pf["dimension"] == 2 ** 6 - 1
    assert pe["dimension"] > p0["dimension"]
    assert abs(pf["missing_weight"]) < 1e-8
    assert abs(np.sum(pf["weight"]) + pf["elastic"] - 0.25) < 1e-8


@pytest.mark.parametrize("version", JULIA)
def test_julia_kondo_task_spectrum_matches_python(version):
    """The TASK Kondo spectrum (taskdmrg.session_spectrum, fed with terms)
    on julia_live against python on a saturated chain."""
    from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
    sj = TaskKondoSpectrum(heisenberg(version=version), 0)
    sp = TaskKondoSpectrum(heisenberg(version="python"), 0)
    assert len(sj.e) == len(sp.e)
    assert np.allclose(np.sort(sj.e), np.sort(sp.e), atol=1e-8)
