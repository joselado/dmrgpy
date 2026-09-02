"""Spectral weights of the infinite-chain quasiparticle branches --
`Infinite_Many_Body_Chain.spectral_weights` /
`dynamical_structure_factor`, i.e. the exact delta-peak residues each
tangent-space excitation branch contributes to S(k,w) in the
thermodynamic limit (pyitensor.idmrg_excitations._spectral_source_vector).

Three independent references are used here, deliberately, because a
spectral weight has no single golden value to pin:

- an **exactly solvable product state** (the J=0 transverse-field Ising
  chain, a D=1 fully polarized state whose only excitation is one flipped
  spin): sigma^x/sigma^y must carry weight exactly 1 at every momentum,
  sigma^z exactly 0 -- the latter also being the check that the *connected*
  subtraction at k=0 happens on its own, since <sigma^z> = -1 there;
- the **static sum rule**: a one-site operator applied to a uniform MPS
  lands exactly inside the variational tangent space, so the total weight
  over every branch is the per-site connected static structure factor
  `sum_r e^{ikr}(<O_0 O_r> - <O>^2)` -- compared here against a real-space
  sum of `Infinite_Many_Body_Chain.correlator`, machinery that shares
  nothing with the excitation ansatz;
- the **f-sum rule**: the first moment `sum_a E_a w_a(k)` must equal
  `(1/2)<[O_k^dagger,[H,O_k]]>`, which for O=sigma^z on the TFIM is
  `4J(<sigma^x_0 sigma^x_1> - cos(k) <sigma^y_0 sigma^y_1>)` in closed
  form. Unlike the static sum rule this one involves the excitation
  *energies* as well as the weights, so it tests H_eff(k) and the source
  vector together, and it converges with the ground state's own bond
  dimension (2e-4 at D=2, 1e-7 at D=4 -- measured).
"""

import numpy as np
import pytest

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations


def _tfim(g, J=1.0, D=2):
    """H = -J*sum(sigma^x sigma^x) - g*sum(sigma^z), converged with VUMPS
    at bond dimension D (Sx=sigma^x/2, so -4*J*Sx*Sx = -J*sigma^x sigma^x
    and -2*g*Sz = -g*sigma^z -- the same convention tests/test_vumps.py's
    own `_tfim_chain` uses)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-4.0 * J * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.gs_energy()
    return ic


_CACHE = {}


def _tfim_cached(g, J=1.0, D=2):
    """VUMPS convergence is by far the dominant cost in this file and the
    chains here are read-only afterwards, so converge each (g,J,D) once."""
    key = (g, J, D)
    if key not in _CACHE:
        _CACHE[key] = _tfim(g, J=J, D=D)
    return _CACHE[key]


def _static_structure_factor(ic, opname, k, rmax=60, p=0, step=1):
    """sum_r e^{ikr} (<O_0 O_r> - <O>^2) summed over r in (-inf, inf) from
    the real-space correlator, `step` physical sites per unit cell -- the
    independent reference for spectral_weights' own `return_total`.

    Uses C(-r) = C(r), true for every (reflection-symmetric, real) model
    used in this file. Needs a *gapped* chain: on a critical one C(r)
    decays algebraically and this sum does not converge."""
    mean = complex(ic.vev(opname, p)).real
    C = [complex(ic.correlator(opname, p, opname, step * r)).real - mean * mean
         for r in range(rmax)]
    return C[0] + 2.0 * sum(np.cos(k * r) * C[r] for r in range(1, rmax))


# ---------------------------------------------------------------------------
# Exactly solvable limit: J=0, a polarized product state (D=1).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [0.0, 0.7, np.pi])
@pytest.mark.parametrize("opname", ["Sx", "Sy"])
def test_product_state_spin_flip_weight_is_one(k, opname):
    """With J=0 the ground state is the exact product state |down...down>
    (sigma^z=-1 everywhere) and the only excitation is a single flipped
    spin, at energy 2g for every momentum. sigma^x and sigma^y each flip
    that spin with unit amplitude, so the weight is exactly 1 -- Sx/Sy are
    sigma^x/2, hence 1/4 here."""
    ic = _tfim_cached(1.0, J=0.0, D=1)
    e, w = ic.spectral_weights(opname, k, n=1)
    assert e[0] == pytest.approx(2.0, abs=1e-8)
    assert w[0] == pytest.approx(0.25, abs=1e-10)


@pytest.mark.parametrize("k", [0.0, 0.7, np.pi])
def test_product_state_sz_weight_is_zero(k):
    """sigma^z is diagonal in the product state's own basis, so it cannot
    reach the one-magnon state at all: weight exactly 0. At k=0 this is
    also the check that the *connected* subtraction is automatic -- the
    disconnected part is <sigma^z>^2 = 1, not small, and would completely
    swamp the answer if it were not projected out (see
    idmrg_excitations._spectral_resolvent's own docstring)."""
    ic = _tfim_cached(1.0, J=0.0, D=1)
    _, w, total = ic.spectral_weights("Sz", k, n=1, return_total=True)
    assert w[0] == pytest.approx(0.0, abs=1e-12)
    assert total == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Static sum rule against an independent real-space correlator sum.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [0.0, 0.4, 1.3, np.pi])
@pytest.mark.parametrize("opname", ["Sx", "Sz"])
def test_total_weight_matches_static_structure_factor(k, opname):
    ic = _tfim_cached(2.5, D=2)
    _, _, total = ic.spectral_weights(opname, k, n=1, return_total=True)
    ref = _static_structure_factor(ic, opname, k)
    assert total == pytest.approx(ref, abs=1e-9)


def test_branch_weights_sum_to_the_total():
    """`return_total` is defined as the sum over EVERY branch, so asking
    for all of them must reproduce it exactly -- this is what makes
    `weights.sum()/total` a meaningful "how much of the response is in
    these branches" fraction."""
    ic = _tfim_cached(2.5, D=2)
    env = ic._get_excitation_environment()
    dim = env.D * env.D * (env.d_g - 1)
    _, w, total = ic.spectral_weights("Sx", 0.9, n=dim, return_total=True)
    assert len(w) == dim
    assert w.sum() == pytest.approx(total, rel=1e-12)


def test_single_mode_fraction_separates_one_particle_from_continuum():
    """A physics check the sum rules alone do not make: in the
    paramagnetic TFIM sigma^x is parity-odd and creates exactly ONE
    quasiparticle, so the lowest branch carries essentially all of its
    weight. sigma^z is parity-even and cannot reach that branch at all --
    an exact selection rule, and one no sum rule would catch, since the
    weight it forbids on the lowest branch simply moves to the higher ones
    (the ansatz's own approximation to the two-particle continuum) without
    changing any total. Measured at ~1e-21, i.e. exactly zero."""
    ic = _tfim_cached(2.5, D=2)
    _, wx, tx = ic.spectral_weights("Sx", 0.9, n=1, return_total=True)
    _, wz, tz = ic.spectral_weights("Sz", 0.9, n=1, return_total=True)
    assert wx[0] / tx > 0.99
    assert wz[0] / tz < 1e-12
    assert tz > 1e-3      # ... while the total sigma^z response is not zero


# ---------------------------------------------------------------------------
# f-sum rule: the first moment, in closed form.
# ---------------------------------------------------------------------------

def _f_sum_rule(ic, k, J=1.0):
    """(1/2)<[O_{-k},[H,O_k]]> for O=Sz on H = -J*sum(sigma^x sigma^x)
    - g*sum(sigma^z). Only the bond term contributes (the field commutes
    with sigma^z); evaluating the double commutator gives
    4J(<sigma^x_0 sigma^x_1> - cos(k)<sigma^y_0 sigma^y_1>) for O=sigma^z,
    i.e. a quarter of that for O=Sz=sigma^z/2 -- and `correlator` returns
    <Sx Sx> = <sigma^x sigma^x>/4, so the two factors of 4 cancel and the
    expression below is written directly in Sx/Sy correlators."""
    XX = complex(ic.correlator("Sx", 0, "Sx", 1)).real   # = <sigma^x sigma^x>/4
    YY = complex(ic.correlator("Sy", 0, "Sy", 1)).real
    return 4.0 * J * (XX - np.cos(k) * YY)


@pytest.mark.parametrize("k", [0.0, 0.4, 1.3, 2.2])
def test_first_moment_matches_f_sum_rule(k):
    """Exact in the converged-ground-state limit, so the tolerance is set
    by the bond dimension rather than by machine precision: measured 2e-4
    relative at D=2, 1e-7 at D=4."""
    ic = _tfim_cached(2.5, D=4)
    env = ic._get_excitation_environment()
    dim = env.D * env.D * (env.d_g - 1)
    e, w = ic.spectral_weights("Sz", k, n=dim)
    assert float(np.sum(e * w)) == pytest.approx(_f_sum_rule(ic, k), rel=1e-5)


# ---------------------------------------------------------------------------
# Two-site unit cell (grouped supersite).
# ---------------------------------------------------------------------------

def _dimerized_chain(J1=1.0, J2=0.25, D=4):
    """Dimerized Heisenberg chain, two sites per cell, converged with
    VUMPS at bond dimension `D`.

    D=4 is a physical choice, not a comfortable one. This Hamiltonian is
    SU(2)-symmetric, so the Schmidt spectrum across a cell bond comes in
    degenerate multiplets: measured at D=8 it is

        1.000e+00 | 3.355e-02  3.355e-02  3.355e-02 | 4.664e-05 (x3) | ...

    i.e. a singlet then triplets, so the only truncations that keep an
    invariant subspace are D = 1, 4, 7, ... Cutting a triplet in half
    leaves the retained space non-invariant and VUMPS cannot settle in
    it: at D=3 the converged *energy* itself varied by 3e-6 between runs
    (2e-5 at D=2), and the spectral-weight check below came out at 2e-5
    typical and 6.2e-5 worst over 30 runs against a 1e-4 tolerance --
    close enough to fail intermittently, which it did. At D=4 the same 30
    checks land at 7e-15 and the energy is reproducible to 1.7e-16.

    So do not lower this back to 3 to save time (it saves none -- see the
    cache below) and do not raise it to 6, which cuts the *second*
    triplet and is 1e-7 rather than machine precision. 1, 4 and 7 are the
    safe values; D=1 drops the inter-cell entanglement this file's
    two-site-cell tests exist to exercise. Same lesson as
    tests/test_vumps_excitations_v3.py's
    test_excitation_accuracy_degrades_with_redundant_bond_dimension: ask
    for the bond dimension the state actually wants."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    h = J1 * (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1])
    h = h + J2 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                  + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = D
    ic.gs_energy()
    return ic


_DIMER_CACHE = {}


def _dimerized_cached(J1=1.0, J2=0.25, D=4):
    """Converge the dimerized chain once, for the same reason
    _tfim_cached exists: VUMPS dominates this file's runtime and the
    chains are read-only afterwards. This also removes five of the six
    independent random VUMPS starts the parametrization used to take."""
    key = (J1, J2, D)
    if key not in _DIMER_CACHE:
        _DIMER_CACHE[key] = _dimerized_chain(J1=J1, J2=J2, D=D)
    return _DIMER_CACHE[key]


@pytest.mark.parametrize("p", [0, 1])
@pytest.mark.parametrize("k", [0.0, 0.9, np.pi])
def test_two_site_cell_total_weight_matches_structure_factor(p, k):
    """n_uc=2: momentum is per unit cell, so the real-space reference sums
    over CELL separations (2 physical sites each) with both operators on
    the same sub-site p. The dimerized chain is gapped, so the sum
    converges.

    The tolerance used to be 1e-4 relative, on the reasoning that it was
    set by how well VUMPS had converged the grouped d_g=4 supersite
    rather than by the spectral-weight arithmetic. That was true but
    avoidable, and it made this test flaky: the residual sat at 2e-5
    typical / 6.2e-5 worst, close enough to 1e-4 to trip. The cause was
    the bond dimension splitting an SU(2) triplet, not anything about
    this check -- see _dimerized_chain. At the multiplet-complete D=4 the
    agreement is 7e-15 over 30 runs, so this is now pinned three orders
    inside that, tight enough to actually catch a wrong weight."""
    ic = _dimerized_cached()
    _, _, total = ic.spectral_weights("Sz", k, p=p, n=1, return_total=True)
    ref = _static_structure_factor(ic, "Sz", k, rmax=50, p=p, step=2)
    assert total == pytest.approx(ref, rel=1e-11)


def test_dimerized_cell_bond_spectrum_is_multiplet_degenerate():
    """Pins the reason _dimerized_chain runs at D=4 rather than 3.

    The Hamiltonian is SU(2)-symmetric, so the Schmidt values across a
    cell bond come in degenerate multiplets and only D = 1, 4, 7, ...
    keep an invariant subspace. If this ever stops holding -- a changed
    J1/J2, a symmetry-breaking term, a different grouping -- the D=4
    choice above stops being principled and the tight tolerance there
    will start failing for a reason nobody remembers. Check it here
    instead, where the failure names itself."""
    ic = _dimerized_cached(D=8)
    env = ic._get_excitation_environment()
    s = np.linalg.svd(np.array(env.C), compute_uv=False)
    s = s / s[0]
    assert s[0] == pytest.approx(1.0)                    # the cell singlet
    assert s[1] == pytest.approx(s[2], rel=1e-8)         # a triplet, not
    assert s[1] == pytest.approx(s[3], rel=1e-8)         # three singlets
    assert s[3] > 100 * s[4]                             # and a clear gap
    # so truncating at 2 or 3 splits that triplet, while 1 and 4 do not
    assert s[1] / s[0] < 0.1


# ---------------------------------------------------------------------------
# Symmetry, and the dense-vs-iterative linear algebra.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("k", [0.4, 1.7])
def test_weights_are_symmetric_under_k_to_minus_k(k):
    """A reflection-symmetric Hamiltonian has S(k,w)=S(-k,w), so the sign
    convention documented on spectral_weights can only relabel k, never
    change a number."""
    ic = _tfim_cached(2.5, D=2)
    e_p, w_p = ic.spectral_weights("Sx", k, n=1)
    e_m, w_m = ic.spectral_weights("Sx", -k, n=1)
    assert e_p[0] == pytest.approx(e_m[0], abs=1e-9)
    assert w_p[0] == pytest.approx(w_m[0], abs=1e-9)


def test_iterative_resolvent_matches_dense(monkeypatch):
    """`_spectral_resolvent` switches from an LU-factored dense solve to
    GMRES above `_RESOLVENT_DENSE_MAX`; forcing the threshold to 0 must not
    move any number. Both caches are cleared so the second pass really
    rebuilds its resolvents rather than reusing the dense ones."""
    ic = _tfim_cached(2.5, D=2)
    ks = [0.0, 0.8, np.pi]
    dense = [ic.spectral_weights("Sx", k, n=1, return_total=True) for k in ks]

    env = ic._get_excitation_environment()
    env.resolvent_cache.clear()
    env._spectral_resolvent_cache = {}
    monkeypatch.setattr(idmrg_excitations, "_RESOLVENT_DENSE_MAX", 0)
    iterative = [ic.spectral_weights("Sx", k, n=1, return_total=True) for k in ks]
    env.resolvent_cache.clear()
    env._spectral_resolvent_cache = {}

    for (e_d, w_d, t_d), (e_i, w_i, t_i) in zip(dense, iterative):
        assert e_i[0] == pytest.approx(e_d[0], abs=1e-9)
        assert w_i[0] == pytest.approx(w_d[0], abs=1e-9)
        assert t_i == pytest.approx(t_d, abs=1e-9)


# ---------------------------------------------------------------------------
# S(k,w) grid, and the argument checks.
# ---------------------------------------------------------------------------

def test_structure_factor_grid_integrates_to_the_branch_weight():
    """Each broadened peak is a normalized Lorentzian, so integrating
    S(k,w) over w at fixed k returns that momentum's summed branch weight.
    The tolerance is the Lorentzian tail the finite energy window cuts off,
    ~(delta/pi)*(1/(E-e_lo) + 1/(e_hi-E)) -- 1e-3 relative for the window
    used here, which is why the window is deliberately far wider than the
    band rather than snug around it."""
    ic = _tfim_cached(2.5, D=2)
    ks = np.linspace(-np.pi, np.pi, 7)
    es = np.linspace(-20.0, 35.0, 11001)
    ks_out, es_out, S = ic.dynamical_structure_factor(
        "Sx", ks=ks, energies=es, delta=0.05, n=1)
    assert S.shape == (len(ks), len(es))
    assert np.allclose(ks_out, ks) and np.allclose(es_out, es)
    for i, k in enumerate(ks):
        _, w = ic.spectral_weights("Sx", k, n=1)
        assert np.trapezoid(S[i], es) == pytest.approx(w[0], rel=3e-3)


def test_structure_factor_default_grid_covers_the_band():
    ic = _tfim_cached(2.5, D=2)
    ks, es, S = ic.dynamical_structure_factor(
        "Sx", ks=np.linspace(-np.pi, np.pi, 5), delta=0.1)
    assert S.shape == (5, 200)
    band = [ic.spectral_weights("Sx", k, n=1)[0][0] for k in ks]
    assert es[0] < min(band) and es[-1] > max(band)


def test_rejects_out_of_range_subsite():
    ic = _tfim_cached(2.5, D=2)
    with pytest.raises(ValueError, match="p must be in"):
        ic.spectral_weights("Sz", 0.0, p=1)


def test_rejects_fermionic_operator():
    """A parity-odd operator's Jordan-Wigner string would have to be closed
    at infinity -- the same reason correlator() rejects an odd-parity
    pair."""
    ic = infinitechain.Infinite_Many_Body_Chain([0], itensor_version="python")
    Cd, C = ic.get_operator("Cdag", 0, "C"), ic.get_operator("C", 0, "C")
    CR = ic.get_operator("C", 0, "R")
    CdR = ic.get_operator("Cdag", 0, "R")
    ic.set_hamiltonian(-(Cd * CR + CdR * C) + 0.5 * ic.get_operator("N", 0, "C"))
    ic.gs_method = "vumps"
    ic.maxm = 2
    with pytest.raises(ValueError, match="fermionic"):
        ic.spectral_weights("C", 0.0)


def test_rejects_itensor_version_3():
    """The mpscpp3 port of the excitation ansatz returns energies only --
    deliberately, the same way it has not picked up the iterative
    eigensolver (docs/idmrg_improvement_plan.md)."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 5.0 * ic.SzC[0])
    ic.gs_method = "vumps"
    with pytest.raises(NotImplementedError, match="itensor_version"):
        ic.spectral_weights("Sz", 0.0)


# ---------------------------------------------------------------------------
# AKLT (S=1): a D=2-exact reference with closed-form answers.
#
# The AKLT ground state IS a bond-dimension-2 MPS exactly, so at D=2 VUMPS
# converges to the exact state and every quantity below has a closed form --
# which makes this the only reference in this file that tests the spectral
# weights on a d>2 (spin-1) chain against something other than a numerical
# cross-check. It also exercises a genuinely SU(2)-symmetric model, where the
# branches come in degenerate multiplets: at D=2 the D^2*(d_g-1) = 8 branches
# split into a magnon TRIPLET and a quintuplet at every momentum, and the
# two cross near k~0.9. Both facts shape what a test here may assert (see
# `spectral_weights`' own docstring).
# ---------------------------------------------------------------------------

def _aklt_chain():
    """H = sum_j [ S_j.S_{j+1} + (1/3)(S_j.S_{j+1})^2 ], the AKLT point.

    Note the biquadratic term is the FULL (S.S)^2 = sum_{ab} S^a_i S^b_i
    S^a_j S^b_j, not the diagonal-only `sum_a (S^a_i S^a_j)^2` that
    `examples/spin_models/bilinear_biquadratic` uses -- only the former has
    the exact valence-bond-solid ground state this test is built on."""
    ic = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
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


_AKLT = []


def _aklt_cached():
    if not _AKLT:
        _AKLT.append(_aklt_chain())
    return _AKLT[0]


def _aklt_szz(k):
    """The exact AKLT static structure factor sum_r e^{ikr}<S^z_0 S^z_r>,
    summed in closed form from <S_0.S_r> = 4(-1/3)^r (r>=1), <S^z S^z> =
    <S.S>/3 by isotropy, and <S^z_0 S^z_0> = S(S+1)/3 = 2/3."""
    x = -1. / 3.
    tail = (x * np.cos(k) - x * x) / (1 - 2 * x * np.cos(k) + x * x)
    return 2. / 3. + (8. / 3.) * tail


def _aklt_sma(k):
    """Arovas, Auerbach & Haldane's single-mode dispersion for the AKLT
    chain, PRL 60, 531 (1988): w(k) = (5/27)(5 + 3 cos k) for the PROJECTOR
    normalization H = sum_j P^(2)_{j,j+1}. `_aklt_chain` uses H = S.S +
    (1/3)(S.S)^2 = 2*sum_j P^(2) + const, hence the factor 2."""
    return 2.0 * (5. / 27.) * (5. + 3. * np.cos(k))


def test_aklt_ground_state_is_exact_at_bond_dimension_two():
    """Guards the premise every other AKLT test rests on."""
    ic = _aklt_cached()
    assert ic.e0 == pytest.approx(-2. / 3., abs=1e-9)
    for r in range(1, 5):
        assert complex(ic.correlator("Sz", 0, "Sz", r)).real == pytest.approx(
            (4. / 3.) * (-1. / 3.) ** r, abs=1e-9)


@pytest.mark.parametrize("k", [0.5, 1.0, 2.0, 2.6, np.pi])
def test_aklt_total_weight_matches_closed_form_structure_factor(k):
    """The sum rule against an analytic answer rather than against another
    dmrgpy code path -- the strongest form of the check in this file."""
    ic = _aklt_cached()
    _, _, total = ic.spectral_weights("Sz", k, n=1, return_total=True)
    assert total == pytest.approx(_aklt_szz(k), abs=1e-9)


def test_aklt_total_weight_vanishes_at_zero_momentum():
    """S^zz(0) = 0 exactly: sum_j S^z_j is conserved, so the k=0 operator
    cannot connect the ground state to anything. The closed form gives 0
    there too, and this is also the sharpest possible test of the automatic
    connected subtraction -- nothing survives it at all."""
    ic = _aklt_cached()
    _, _, total = ic.spectral_weights("Sz", 0.0, n=1, return_total=True)
    assert total == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("k", [0.5, 1.0, 2.0, 2.6, np.pi])
def test_aklt_first_moment_matches_arovas_auerbach_haldane(k):
    """`sum_a E_a w_a / total` is the Rayleigh quotient of the single-mode
    state S^z_k|Psi>, and that state lies exactly inside the tangent space,
    so this must reproduce AAH's published SMA dispersion exactly -- an
    analytic check that pins the excitation ENERGIES and the weights
    together, on a model where the ground state carries no bond-dimension
    error at all."""
    ic = _aklt_cached()
    e, w, total = ic.spectral_weights("Sz", k, n=8, return_total=True)
    assert float(np.sum(e * w)) / total == pytest.approx(_aklt_sma(k), rel=1e-7)


@pytest.mark.parametrize("k", [0.5, 2.0, np.pi])
def test_aklt_spin_selection_rule_and_multiplet_summed_weight(k):
    """S^z is a spin-1 operator, so it reaches the magnon triplet and
    nothing else: the quintuplet's weight is ~1e-23, and the triplet's three
    (individually basis-arbitrary, see `spectral_weights`' docstring)
    weights sum to the entire total. Includes k=0.5, which is below the
    multiplet crossing -- there the lowest branch is the *quintuplet*, so
    `n=1` alone would report a weight of zero for a momentum whose total
    response is perfectly finite."""
    ic = _aklt_cached()
    e, w, total = ic.spectral_weights("Sz", k, n=8, return_total=True)
    triplet = np.abs(e - _aklt_sma(k)) < 1e-6
    assert triplet.sum() == 3
    assert w[triplet].sum() == pytest.approx(total, rel=1e-9)
    assert w[~triplet].sum() / total < 1e-12
