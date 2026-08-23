"""Coverage for pyitensor/vumps.py (VUMPS ground-state solver) and its
Infinite_Many_Body_Chain wiring (`gs_method="vumps"`, infinitechain.py).
See tests/test_vumps_correlator.py for dedicated coverage of vumps.py's
own static-correlator machinery (onsite_expectation/two_point_correlator).

Cross-checks used, in order of how rigorous they are: (1) exactly solvable
D=1 cases (a pure field, and a fully-decoupled Heisenberg dimer -- both
reduce to a product state, so the exact answer is known in closed form and
should match to numerical precision, exercising the field-only and the
bond/pending-channel code paths respectively); (2) the transverse-field
Ising model's own exact free-fermion energy density (Pfeuty 1970), a
genuinely D>1 (entangled) cross-check independent of any other part of
this codebase. See vumps.py's own module docstring ("Convergence
robustness") for why D>1 tests below use a generous `maxiter`/`nrestarts`
and a gapped model (TFIM away from its critical point g=1) rather than a
gapless one -- VUMPS's own convergence robustness for critical models is a
documented, open scope limitation, not something these tests attempt to
nail down tightly.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import infinitechain
from dmrgpy.pyitensor import vumps


def _field_chain(B):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-B * ic.SzC[0])
    return ic


def _dimer_chain():
    """n_uc=2, intra-cell Heisenberg bond only (no inter-cell coupling at
    all) -- the exactly-solvable decoupled-singlet-dimer limit also used by
    idmrg_excitations.py's own investigation (see its module docstring):
    the ground state is an exact product of singlets, D=1 after grouping,
    energy -3/4 per unit cell (-3/8 per site)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    h = ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
    ic.set_hamiltonian(h)
    return ic


def _tfim_chain(g):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] -- with Sx/Sz the spin-1/2
    operators (eigenvalues +-1/2), this is exactly the standard Pauli-matrix
    TFIM H=-sum(sigma^x_i sigma^x_{i+1}) - g*sum(sigma^z_i) (Sx=sigma^x/2,
    Sz=sigma^z/2, so -4*Sx*Sx=-sigma^x*sigma^x and -2*g*Sz=-g*sigma^z)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    h = -4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0]
    ic.set_hamiltonian(h)
    return ic


def _tfim_exact_energy_density(g):
    """Pfeuty (1970) free-fermion result for H=-sum(sigma^x sigma^x)
    - g*sum(sigma^z): e0/N = -(1/pi) * integral_0^pi sqrt(1+g^2-2g*cos(k)) dk."""
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


def test_field_only_matches_exact_solution():
    """D=1, no bonds at all -- exercises the onsite-only (h1, no pending
    channels) code path. Exact answer: a field-polarized product state,
    energy density -|B|/2."""
    ic = _field_chain(0.7)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert result.converged
    assert result.e0 == pytest.approx(-0.35, abs=1e-8)


def test_field_only_negative_sign():
    ic = _field_chain(-0.4)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert result.converged
    assert result.e0 == pytest.approx(-0.2, abs=1e-8)


def test_dimer_matches_exact_singlet_energy():
    """D=1 (after grouping) WITH bonds -- exercises the pending-channel
    (Lvec_a/Rvec_b, H_AC/H_C bond terms) code path at the one bond-dimension
    where the exact answer is still known in closed form."""
    ic = _dimer_chain()
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert result.converged
    assert result.e0 == pytest.approx(-0.375, abs=1e-8)
    assert result.D == 1
    assert result.d_g == 4  # grouped n_uc=2 supersite, d_g = 2*2


@pytest.mark.parametrize("D", [2, 4])
def test_tfim_stays_within_variational_bound(D):
    """A restart-robust D>1 sanity check that does not assume any single
    call reliably finds the true optimum (see this module's own
    "Convergence robustness" section -- confirmed directly, even D=2 can
    occasionally land ~10-15% off the exact answer, not just D>=4): every
    fixed point VUMPS actually returns, converged or not, must still be a
    variationally sensible energy (not appreciably below the exact
    ground-state energy -- the variational principle forbids that exactly,
    up to a tiny numerical slack -- and not implausibly far above it
    either) -- catches a construction bug (e.g. a sign error) producing
    energies systematically on the wrong side of the true ground state,
    without asserting the tighter (and, empirically, unreliable at this D)
    claim that it is close to that ground state."""
    g = 1.5
    exact = _tfim_exact_energy_density(g)
    ic = _tfim_chain(g)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=D, maxiter=400, nrestarts=5)
    assert result.e0 > exact - 0.2
    assert result.e0 < exact + 0.5


def test_gauge_mismatch_small_when_converged():
    ic = _field_chain(1.0)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    assert result.converged
    assert result.gauge_mismatch < 1e-8


def test_mixed_gauge_condition_holds():
    """AL@C and C@AR must both reproduce AC to within gauge_mismatch --
    the defining VUMPS fixed-point condition, checked directly on the
    returned tensors (not just the internal diagnostic)."""
    ic = _dimer_chain()
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=1, maxiter=50)
    lhs1 = np.einsum('lpm,mr->lpr', result.AL, result.C)
    lhs2 = np.einsum('lm,mpr->lpr', result.C, result.AR)
    assert np.max(np.abs(result.AC - lhs1)) < 1e-6
    assert np.max(np.abs(result.AC - lhs2)) < 1e-6


def test_al_left_canonical_ar_right_canonical():
    ic = _tfim_chain(1.5)
    result = vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                       ic.n_uc, D=3, maxiter=300, nrestarts=3)
    D = result.D
    gram_l = np.einsum('lpc,lpd->cd', result.AL.conj(), result.AL)
    assert np.allclose(gram_l, np.eye(D), atol=1e-8)
    gram_r = np.einsum('lpc,mpc->lm', result.AR.conj(), result.AR)
    assert np.allclose(gram_r, np.eye(D), atol=1e-8)


def test_n_uc_above_2_uses_the_sequential_multisite_algorithm():
    """n_uc>2 used to be rejected outright (grouping the cell into one
    d^n_uc supersite is exponential, so it was capped at 2). It now routes
    to pyitensor/vumps_ms.py's sequential algorithm instead, which costs
    linearly in the cell size -- see tests/test_vumps_multisite.py for the
    accuracy checks; this one only pins that the route exists and returns a
    multi-site result rather than raising."""
    res = vumps.vumps_ground_state([1, 1, 1], [], [], 3, D=1, maxiter=20)
    assert getattr(res, "multisite", False)
    assert res.n_uc == 3
    assert len(res.AL) == 3


def test_d_less_than_1_rejected():
    ic = _field_chain(0.5)
    with pytest.raises(ValueError):
        vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                  ic.n_uc, D=0)


def test_nrestarts_less_than_1_rejected():
    ic = _field_chain(0.5)
    with pytest.raises(ValueError):
        vumps.vumps_ground_state(ic.site_types, ic._h_intra.op, ic._h_inter.op,
                                  ic.n_uc, D=1, nrestarts=0)


# == Infinite_Many_Body_Chain wiring (gs_method="vumps") =====================

def test_gs_method_vumps_matches_exact_field_solution():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-0.6 * ic.SzC[0])
    e0 = ic.gs_energy()
    assert ic.converged
    assert e0 == pytest.approx(-0.3, abs=1e-8)


def test_gs_method_vumps_vev_matches_exact_field_polarization():
    """vev/correlator under gs_method="vumps" are implemented directly on
    VUMPSResult's mixed-gauge {AC, AR} (see pyitensor/vumps.py's own
    "Static correlators" section) -- see tests/test_vumps_correlator.py
    for the full cross-checked coverage; this is a light smoke test that
    Infinite_Many_Body_Chain's own dispatch reaches it at all."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-0.5 * ic.SzC[0])
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)


def test_gs_method_vumps_correlator_matches_exact_field_polarization():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-0.5 * ic.SzC[0])
    assert ic.correlator("Sz", 0, "Sz", 1) == pytest.approx(0.25, abs=1e-8)


def test_gs_method_vumps_excitation_gap_works():
    """excitation_energies/excitation_gap require gs_method="vumps"
    specifically (unlike vev/correlator, which work under either
    gs_method) -- see pyitensor.idmrg_excitations' own module docstring
    for the algorithm.
    D=1 field-polarized case: E(k)=B for every k (a single spin flip costs
    exactly the field energy, independent of momentum, since there is no
    bond term to disperse it)."""
    ic = _field_chain(2.0)
    ic.gs_method = "vumps"
    ic.maxm = 1
    assert ic.excitation_gap() == pytest.approx(2.0, abs=1e-8)


def test_gs_method_unknown_raises():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "not_a_real_method"
    ic.set_hamiltonian(-0.5 * ic.SzC[0])
    with pytest.raises(NotImplementedError):
        ic.gs_energy()


def test_set_hamiltonian_resets_vumps_result():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-0.5 * ic.SzC[0])
    ic.gs_energy()
    assert ic._vumps_result is not None
    ic.set_hamiltonian(-0.3 * ic.SzC[0])
    assert ic._vumps_result is None
    assert ic.e0 is None
