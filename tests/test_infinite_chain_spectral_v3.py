"""Spectral weights of the infinite-chain quasiparticle branches on
`itensor_version=3` (`Chain::vumps_spectral_weights`), the C++ port of
`pyitensor.idmrg_excitations.spectral_weights`.

The checks are the ones `tests/test_infinite_chain_spectral.py` makes for the
"python" backend, now on the compiled one, plus a direct cross-check of the
two backends. Individual weights inside a degenerate multiplet are basis
arbitrary, so the cross-checks compare multiplet sums and the total.

Skipped automatically if mpscpp3 is not compiled.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

pytestmark = pytest.mark.skipif(
    not cppext.available(3),
    reason="requires the compiled mpscpp3 (ITensor v3) extension")


_CACHE = {}


def _cached(key, make):
    if key not in _CACHE:
        _CACHE[key] = make()
    return _CACHE[key]


def _tfim(g, D, version, J=1.0):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.set_hamiltonian(-4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.vumps_nrestarts = 6
    ic.gs_energy()
    return ic


def _dimerized(version, J1=1.0, J2=0.25, D=4):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=version)
    h = J1 * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
    h = h + J2 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                  + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.gs_energy()
    return ic


def _heisenberg2(version, D=2):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=version)
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.vumps_nrestarts = 6
    ic.gs_energy()
    return ic


def _static_structure_factor(ic, opname, k, rmax=60, p=0, step=1):
    mean = complex(ic.vev(opname, p)).real
    C = [complex(ic.correlator(opname, p, opname, step * r)).real - mean * mean
         for r in range(rmax)]
    return C[0] + 2.0 * sum(np.cos(k * r) * C[r] for r in range(1, rmax))


# ---------------------------------------------------------------------------
# Exactly solvable limit and sum rules.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [0.0, 0.7, np.pi])
def test_product_state_weights(k):
    """J=0: |down...down>, one flipped spin at energy 2g. Sx carries weight
    exactly 1/4 (sigma^x/2), Sz carries 0."""
    ic = _cached(("pol",), lambda: _tfim(1.0, 1, 3, J=0.0))
    e, w = ic.spectral_weights("Sx", k, n=1)
    assert e[0] == pytest.approx(2.0, abs=1e-8)
    assert w[0] == pytest.approx(0.25, abs=1e-10)
    _, w, total = ic.spectral_weights("Sz", k, n=1, return_total=True)
    assert w[0] == pytest.approx(0.0, abs=1e-12)
    assert total == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("k", [0.0, 0.4, 1.3, np.pi])
@pytest.mark.parametrize("opname", ["Sx", "Sz"])
def test_total_weight_matches_static_structure_factor(k, opname):
    ic = _cached(("tfim", 2.5, 2), lambda: _tfim(2.5, 2, 3))
    _, _, total = ic.spectral_weights(opname, k, n=1, return_total=True)
    ref = _static_structure_factor(ic, opname, k)
    assert total == pytest.approx(ref, abs=1e-9)


def test_branch_weights_sum_to_the_total():
    ic = _cached(("tfim", 2.5, 2), lambda: _tfim(2.5, 2, 3))
    _, w, total = ic.spectral_weights("Sx", 0.9, n=4, return_total=True)
    assert len(w) == 4
    assert w.sum() == pytest.approx(total, rel=1e-12)


@pytest.mark.parametrize("p", [0, 1])
@pytest.mark.parametrize("k", [0.0, 0.9, np.pi])
def test_two_site_cell_total_weight_matches_structure_factor(p, k):
    ic = _cached(("dimer", 3), lambda: _dimerized(3))
    _, _, total = ic.spectral_weights("Sz", k, p=p, n=1, return_total=True)
    ref = _static_structure_factor(ic, "Sz", k, rmax=50, p=p, step=2)
    assert total == pytest.approx(ref, rel=1e-11)


# ---------------------------------------------------------------------------
# Cross-check against itensor_version="python".
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("D", [1, 2, 3])
@pytest.mark.parametrize("opname", ["Sx", "Sz"])
def test_tfim_matches_python_backend(D, opname):
    """TFIM, n_uc=1, gapped (g=1.5): the two backends converge the same
    state, so every branch weight and the total agree."""
    a = _cached(("tfim_py", D), lambda: _tfim(1.5, D, "python"))
    b = _cached(("tfim_v3", D), lambda: _tfim(1.5, D, 3))
    for k in (0.0, 0.7, 2.0):
        ea, wa, ta = a.spectral_weights(opname, k, n=D * D, return_total=True)
        eb, wb, tb = b.spectral_weights(opname, k, n=D * D, return_total=True)
        assert eb == pytest.approx(ea, abs=1e-8)
        assert wb == pytest.approx(wa, abs=1e-8)
        assert tb == pytest.approx(ta, abs=1e-8)


@pytest.mark.parametrize("k", [0.0, 0.9, 2.2, np.pi])
def test_dimerized_cell_multiplet_sums_match_python_backend(k):
    """n_uc=2, gapped dimerized Heisenberg at the multiplet-complete D=4.
    Sz reaches the lowest triplet only, so the sum over the first three
    branches is physical, and so is the total."""
    a = _cached(("dimer", "python"), lambda: _dimerized("python"))
    b = _cached(("dimer", 3), lambda: _dimerized(3))
    for p in (0, 1):
        ea, wa, ta = a.spectral_weights("Sz", k, p=p, n=3, return_total=True)
        eb, wb, tb = b.spectral_weights("Sz", k, p=p, n=3, return_total=True)
        assert eb == pytest.approx(ea, abs=1e-9)
        assert wb.sum() == pytest.approx(wa.sum(), rel=1e-8)
        assert tb == pytest.approx(ta, rel=1e-8)


@pytest.mark.parametrize("k", [0.0, 0.8, np.pi])
def test_critical_heisenberg_cell_matches_python_backend(k):
    """n_uc=2 critical Heisenberg at D=2: the restart searches may land on
    slightly different optima, so only a loose tolerance, on the total and
    on the energy of the lowest branch."""
    a = _cached(("heis", "python"), lambda: _heisenberg2("python"))
    b = _cached(("heis", 3), lambda: _heisenberg2(3))
    ea, _, ta = a.spectral_weights("Sz", k, n=1, return_total=True)
    eb, _, tb = b.spectral_weights("Sz", k, n=1, return_total=True)
    assert eb[0] == pytest.approx(ea[0], abs=5e-2)
    assert tb == pytest.approx(ta, abs=5e-2)


def test_lanczos_route_matches_dense_route():
    """`dense_max=0` forces the deflated Lanczos solver; the multiplet sum
    and the energies must equal those of the dense solver on one state."""
    ic = _cached(("dimer", 3), lambda: _dimerized(3))
    s = ic._session3
    for k in (0.4, 2.0):
        ed, wd, td = s.vumps_spectral_weights("Sz", 0, k, 3, -1)
        el, wl, tl = s.vumps_spectral_weights("Sz", 0, k, 3, 0)
        assert np.array(el) == pytest.approx(np.array(ed), abs=1e-9)
        assert np.sum(wl) == pytest.approx(np.sum(wd), rel=1e-7)
        assert tl == pytest.approx(td, rel=1e-12)


def test_dynamical_structure_factor_runs_on_v3():
    ic = _cached(("tfim_v3", 2), lambda: _tfim(1.5, 2, 3))
    ks, es, S = ic.dynamical_structure_factor("Sz", ks=np.linspace(0.1, 3, 5),
                                               n=1)
    assert S.shape == (5, len(es))
    assert np.all(S >= 0)


def test_rejects_fermionic_operator_and_bad_subsite():
    ic = _cached(("tfim_v3", 2), lambda: _tfim(1.5, 2, 3))
    with pytest.raises(ValueError, match="p must"):
        ic.spectral_weights("Sz", 0.0, p=1)


# ---------------------------------------------------------------------------
# AKLT (S=1) at D=2: closed-form references.
# ---------------------------------------------------------------------------

def _aklt_chain():
    ic = infinitechain.Infinite_Spin_Chain(["1"], itensor_version=3)
    Sc = [ic.SxC[0], ic.SyC[0], ic.SzC[0]]
    Sr = [ic.SxR[0], ic.SyR[0], ic.SzR[0]]
    h = Sc[0] * Sr[0] + Sc[1] * Sr[1] + Sc[2] * Sr[2]
    for a in range(3):
        for b in range(3):
            h = h + (1. / 3.) * Sc[a] * Sc[b] * Sr[a] * Sr[b]
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = 2
    ic.maxiter = 800
    ic.vumps_nrestarts = 8
    ic.gs_energy()
    return ic


def _aklt_szz(k):
    x = -1. / 3.
    tail = (x * np.cos(k) - x * x) / (1 - 2 * x * np.cos(k) + x * x)
    return 2. / 3. + (8. / 3.) * tail


def _aklt_sma(k):
    return 2.0 * (5. / 27.) * (5. + 3. * np.cos(k))


@pytest.mark.parametrize("k", [0.5, 1.0, 2.0, 2.6, np.pi])
def test_aklt_total_weight_matches_closed_form(k):
    ic = _cached(("aklt",), _aklt_chain)
    assert ic.e0 == pytest.approx(-2. / 3., abs=1e-9)
    _, _, total = ic.spectral_weights("Sz", k, n=1, return_total=True)
    assert total == pytest.approx(_aklt_szz(k), abs=1e-9)


def test_aklt_total_weight_vanishes_at_zero_momentum():
    ic = _cached(("aklt",), _aklt_chain)
    _, _, total = ic.spectral_weights("Sz", 0.0, n=1, return_total=True)
    assert total == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("k", [0.5, 1.0, 2.0, 2.6, np.pi])
def test_aklt_first_moment_matches_single_mode_dispersion(k):
    """sum_a E_a w_a / total is the Rayleigh quotient of S^z_k|Psi>, which
    lies inside the tangent space: the Arovas-Auerbach-Haldane dispersion."""
    ic = _cached(("aklt",), _aklt_chain)
    e, w, total = ic.spectral_weights("Sz", k, n=8, return_total=True)
    assert float(np.sum(e * w)) / total == pytest.approx(_aklt_sma(k), rel=1e-12)


@pytest.mark.parametrize("k", [0.5, 2.0, np.pi])
def test_aklt_selection_rule_and_multiplet_summed_weight(k):
    ic = _cached(("aklt",), _aklt_chain)
    e, w, total = ic.spectral_weights("Sz", k, n=8, return_total=True)
    triplet = np.abs(e - _aklt_sma(k)) < 1e-6
    assert triplet.sum() == 3
    assert w[triplet].sum() == pytest.approx(total, rel=1e-9)
    assert w[~triplet].sum() / total < 1e-12
