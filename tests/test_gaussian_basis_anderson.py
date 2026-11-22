"""anderson_chain(basis="gaussian"), set_gaussian_basis and
GaussianBasis.natural_basis: the rotation is exact for any reference, the
chain keeps its basis so that get_kondo_screening measures the original
bath, and the natural orbitals of a measured state keep the spins apart."""

import numpy as np
import pytest

from dmrgpy import fermionchain, wilsonchain


def _cpp3():
    from dmrgpy import cppext
    return pytest.mark.skipif(not cppext.available(3),
                              reason="the v3 extension is not compiled")


@pytest.mark.parametrize("kw", [
    dict(ed=-0.5, U=1.0, tD=0.3),                       # symmetric
    dict(ed=-0.2, U=1.0, tD=0.3, B=0.1),                # off symmetry, field
    dict(ed=-0.2, U=1.0, tD=0.3, reference_level=-0.2),  # the bare level
    dict(ed=[-0.4, -0.1], U=0.6, tD=[0.3, 0.2]),         # two levels
])
def test_spectrum_is_the_one_of_the_site_basis(kw):
    n = 3
    site = wilsonchain.anderson_chain(n, Lambda=3.0, **{k: v for k, v in
                                      kw.items() if k != "reference_level"})
    gauss = wilsonchain.anderson_chain(n, Lambda=3.0, basis="gaussian", **kw)
    assert gauss.gaussian_basis is not None and site.gaussian_basis is None
    e0 = site.get_excited(n=4, mode="ED")
    e1 = gauss.get_excited(n=4, mode="ED")
    assert np.allclose(e0, e1, atol=1e-9)


def test_kondo_screening_is_the_one_of_the_site_basis():
    kw = dict(ed=-0.3, U=1.0, tD=0.3, Lambda=3.0, B=0.05)
    site = wilsonchain.anderson_chain(5, **kw)
    gauss = wilsonchain.anderson_chain(5, basis="gaussian", **kw)
    a = site.get_kondo_screening(mode="ED")
    b = gauss.get_kondo_screening(mode="ED")
    for key in ("M_imp", "M_env", "M_n", "schmidt"):
        assert np.allclose(a[key], b[key], atol=1e-9), key


def test_kondo_screening_of_two_levels_per_channel():
    kw = dict(ed=[-0.5, -0.3], U=0.8, tD=[0.3, 0.2], Lambda=3.0)
    site = wilsonchain.anderson_chain(3, **kw)
    gauss = wilsonchain.anderson_chain(3, basis="gaussian", **kw)
    a = site.get_kondo_screening(mode="ED")
    b = gauss.get_kondo_screening(mode="ED")
    for key in ("M_imp", "M_env", "schmidt"):
        assert np.allclose(np.array(a[key], dtype=float),
                           np.array(b[key], dtype=float), atol=1e-9), key


@_cpp3()
def test_deep_wilson_chain_screening_on_dmrg():
    """A deep chain, where a lost low-energy hopping would not show in the
    energy: the screened moment (M_env = -M_imp), the Schmidt values and
    the cloud M_n of the site basis. On this chain the Gaussian basis needs
    the same maxm as the site basis (1e-6 at maxm=64 for both, 3e-5 for it
    at maxm=32), so both run at maxm=64."""
    kw = dict(ed=-0.25, U=0.5, tD=0.2, Lambda=4.0, itensor_version=3)
    site = wilsonchain.anderson_chain(21, **kw)
    gauss = wilsonchain.anderson_chain(21, basis="gaussian", **kw)
    for fc in (site, gauss):
        fc.maxm, fc.nsweeps = 64, 16
    assert gauss.gs_energy() == pytest.approx(site.gs_energy(), abs=1e-6)
    a, b = site.get_kondo_screening(), gauss.get_kondo_screening()
    assert b["M_env"] == pytest.approx(-b["M_imp"], abs=1e-4)
    assert b["M_imp"] == pytest.approx(a["M_imp"], abs=1e-4)
    assert np.allclose(b["schmidt"], a["schmidt"], atol=1e-3)
    assert np.allclose(b["M_n"], a["M_n"], atol=1e-3)


def test_set_hamiltonian_drops_the_basis():
    fc = wilsonchain.anderson_chain(3, basis="gaussian")
    H = fc.site_hamiltonian
    fc.set_hamiltonian(H)
    assert fc.gaussian_basis is None and fc.site_hamiltonian is None
    with pytest.raises(ValueError, match="keeps no Hamiltonian"):
        fc.set_gaussian_basis(fc.get_gaussian_basis(
            h=np.diag([0.0, 1.0, 2.0, 3.0]), nf=4, fixed=[0]))


def test_bad_keywords_raise():
    with pytest.raises(ValueError, match="basis must be one of"):
        wilsonchain.anderson_chain(3, basis="natural")
    with pytest.raises(ValueError, match="needs basis"):
        wilsonchain.anderson_chain(3, reference_level=0.0)
    fc = wilsonchain.anderson_chain(3, basis="gaussian")
    with pytest.raises(ValueError, match="fixed orbitals"):
        fc.get_kondo_screening(imp=1, mode="ED")


def _complex_chain(n):
    h = np.zeros((n, n), dtype=complex)
    for i in range(n-1):
        h[i, i+1] = -np.exp(0.3j*i)
        h[i+1, i] = np.conj(h[i, i+1])
    h[2, 5], h[5, 2] = 0.4j, -0.4j
    fc = fermionchain.Fermionic_Chain(n)
    H = sum(h[i, j]*fc.Cdag[i]*fc.C[j] for i in range(n) for j in range(n)
            if h[i, j] != 0)
    return fc, h, H


def test_measured_correlation_matrix_is_the_reference_at_U0():
    """The convention of original_correlation_matrix: at U=0 the state of
    the rotated chain, taken back to the original modes, is the
    reference, and its natural basis is the same basis."""
    fc, h, H = _complex_chain(8)
    b = fc.get_gaussian_basis(h=h, nf=4, fixed=[0])
    fc.set_gaussian_basis(b, H)
    cm = b.original_correlation_matrix(mode="ED")
    reference = b.V.conj().T @ b.correlation_matrix() @ b.V
    assert np.max(np.abs(cm-reference)) < 1e-12
    b2 = b.natural_basis(mode="ED")
    assert np.max(np.abs(b2.V-b.V)) < 1e-10


def test_natural_basis_is_exact_and_keeps_the_spins_apart():
    fc = wilsonchain.anderson_chain(3, ed=-0.3, U=1.0, tD=0.3, Lambda=3.0,
                                    B=0.05, basis="gaussian")
    e0 = fc.gs_energy(mode="ED")
    b = fc.gaussian_basis.natural_basis(mode="ED")
    assert not b.mixes_spins()
    from dmrgpy import gaussianbasis
    assert b.eigval_cutoff == gaussianbasis.NATURAL_EIGVAL_CUTOFF
    fc.set_gaussian_basis(b)
    assert fc.gaussian_basis is b
    assert fc.gs_energy(mode="ED") == pytest.approx(e0, abs=1e-10)


@_cpp3()
def test_natural_basis_of_a_dmrg_state_drops_the_spin_mixing():
    fc = wilsonchain.anderson_chain(9, ed=-0.25, U=0.5, tD=0.2, Lambda=3.0,
                                    basis="gaussian", itensor_version=3)
    fc.maxm, fc.nsweeps = 16, 8
    e0 = fc.gs_energy()  # 1.3e-4 above the energy at maxm=128, at 8 and 30 sweeps
    b = fc.gaussian_basis.natural_basis()
    assert b.spin_mixing > 1e-12  # 1.4e-4 measured: the spins would mix
    assert not b.mixes_spins()
    fc.set_gaussian_basis(b)
    assert fc.gs_energy() < e0-2e-5  # 5.2e-5 above it, at the same maxm
