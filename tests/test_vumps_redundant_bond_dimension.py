"""VUMPS asked for more bond dimension than the state needs.

A gapped model's exact ground state often needs far fewer than the
requested `maxm` directions -- a field-polarized chain needs exactly one.
The extra directions then carry no Schmidt weight, and the state's
transfer matrix picks up a decoupled unimodular block, i.e. a DEGENERATE
dominant eigenvalue.

That is benign, but it looks identical to the one thing
`_check_dominant_eigenvalue_nondegenerate` (and its C++ counterpart
`vx_check_perron_nondegenerate`) exists to reject: a "cat state", two
branches with matched *nonzero* weight, where no single dominant fixed
point is meaningful. `itensor_version=3` rejected the benign case too, so
`gs_energy()` raised "every attempt at D=... failed" for any polarized
chain at `maxm>1` -- on the sequential solver always, and on the grouped
one whenever an iteration happened to land exactly on the degeneracy
(measured at a second eigenvalue of 0.99996, just outside the guard's own
1e-9, i.e. it survived by luck rather than by design).

Both environment builders now fall back to the fixed points the state
itself names, `C C^dag` and `C^dag C`. See `Chain::vx_bond_fixed_points`
for why that is the right element of the degenerate subspace in both the
benign and the pathological case.

Separately but reachable through the same models: `Chain::vms_ground_state`
had no D-ramp warm start at all (its `reuse` test compared the previous
rung's tensor size against the *new* D, so it could never hold), which
made every rung start from pure noise and land in exactly the redundant
configuration above. `Chain::vms_grow_init` is the fix.

Everything here has an exact answer, so a wrong fixed point cannot pass.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

BACKENDS = ["python"] + ([3] if cppext.available(3) else [])

FIELD = 4.0
J = 0.7
# Fully polarized: <Sz> = 1/2 on every site, so <Sz Sz> = 1/4 at every
# separation and the energy density is -FIELD/2 + J/4 per site.
EXACT_E = -FIELD / 2.0 + J / 4.0


def _polarized(n_uc, reach, D, backend):
    """-FIELD sum Sz + J sum Sz_i Sz_{i+reach} on an n_uc-site cell.

    `reach > n_uc` is what routes the chain to the sequential solver at
    small `n_uc`; `reach == 1` on a 1- or 2-site cell keeps it on the
    grouped one. Both are exercised below, because the guard is in both
    environment builders.
    """
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                            itensor_version=backend)
    ic.maxm, ic.maxiter, ic.etol = D, 300, 1e-12
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        h = h - FIELD * ic.SzC[i]
        k = i + reach
        other = (ic.SzC[k] if k < n_uc
                 else ic.get_operator("Sz", k % n_uc, group=k // n_uc))
        h = h + J * ic.SzC[i] * other
    ic.set_hamiltonian(h)
    return ic


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("D", [2, 4])
def test_grouped_solver_tolerates_redundant_bond_dimension(backend, D):
    """A reach-1 chain on a 1-site cell: the GROUPED path, at two bond
    dimensions the exact (D=1) state does not need."""
    ic = _polarized(1, 1, D, backend)
    assert ic.gs_energy() == pytest.approx(EXACT_E, abs=1e-9)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("n_uc,reach", [(1, 2), (2, 3), (3, 1)])
def test_sequential_solver_tolerates_redundant_bond_dimension(backend, n_uc, reach):
    """The SEQUENTIAL path, reached three ways -- a coupling past the cell
    on a 1-site and on a 2-site cell, and a cell longer than 2 sites.

    On `itensor_version=3` every one of these raised before: at D=1 the
    ramp never grows and the answer was already right, and at D=2 the
    first rung produced the redundant state that the guard rejected."""
    ic = _polarized(n_uc, reach, 2, backend)
    assert ic.gs_energy() == pytest.approx(EXACT_E, abs=1e-9)


@pytest.mark.parametrize("backend", BACKENDS)
def test_observables_are_exact_on_a_redundant_state(backend):
    """The energy alone would not catch a wrongly-chosen fixed point that
    happens to preserve it -- and picking an arbitrary branch out of a
    degenerate subspace is exactly the failure mode that would. `<Sz>` and
    the correlator at several separations are read off the same converged
    state, and are exact for the polarized chain."""
    ic = _polarized(1, 2, 2, backend)
    ic.gs_energy()
    assert ic.vev("Sz", 0).real == pytest.approx(0.5, abs=1e-9)
    for r in (1, 2, 3, 5):
        assert ic.correlator("Sz", 0, "Sz", r).real == pytest.approx(0.25, abs=1e-9)


def _aklt(n_uc, D, backend):
    """H = sum_j [ S_j.S_{j+1} + (1/3)(S_j.S_{j+1})^2 ] on an n_uc-site cell
    of S=1 sites -- the AKLT point, whose exact ground state is a
    bond-dimension-TWO MPS with energy density -2/3.

    The biquadratic term is the FULL (S.S)^2 = sum_{ab} S^a_i S^b_i S^a_j
    S^b_j, not the diagonal-only form; only that has the valence-bond-solid
    ground state (same fixture as tests/test_infinite_chain_spectral.py).
    """
    ic = infinitechain.Infinite_Spin_Chain(["1"] * n_uc,
                                            itensor_version=backend)
    h = 0
    for i in range(n_uc):
        Sc = [ic.SxC[i], ic.SyC[i], ic.SzC[i]]
        j = i + 1
        Sr = ([ic.SxC[j], ic.SyC[j], ic.SzC[j]] if j < n_uc
              else [ic.get_operator(o, j % n_uc, group=j // n_uc)
                    for o in ("Sx", "Sy", "Sz")])
        for a in range(3):
            h = h + Sc[a] * Sr[a]
        for a in range(3):
            for b in range(3):
                h = h + (1. / 3.) * Sc[a] * Sc[b] * Sr[a] * Sr[b]
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter, ic.etol = D, 800, 1e-10
    ic.vumps_nrestarts = 3
    return ic


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("n_uc", [1, 3])
@pytest.mark.parametrize("D", [3, 4])
def test_aklt_is_exact_above_its_own_bond_dimension(backend, n_uc, D):
    """The case that actually DISCRIMINATES between the two ways to resolve
    a degenerate dominant fixed point.

    On a polarized chain every element of the degenerate subspace is the
    ground state (the transfer matrix factorizes as `a (x) U` with `U` a
    bond unitary), so any choice gives the exact energy and the test above
    cannot tell `C C^dag` apart from an eigensolver's arbitrary pick. AKLT
    can: its exact state is genuinely entangled at D=2, so at D=3 or D=4
    the redundant directions sit alongside a real two-dimensional bond, and
    a wrongly-selected branch would carry its own energy. Both `n_uc` values
    are here because they take different code paths (1 -> grouped, 3 ->
    sequential), and both backends because `vumps_ms._cell_fixed_points`
    makes the arbitrary pick while the C++ takes `C C^dag`: -2/3 on both is
    what says the two agree where they could have differed."""
    assert _aklt(n_uc, D, backend).gs_energy() == pytest.approx(-2. / 3., abs=1e-7)


@pytest.mark.skipif(3 not in BACKENDS, reason="mpscpp3 not compiled")
def test_v3_sequential_ramp_actually_warm_starts():
    """The D-ramp on the sequential solver has to CARRY the previous rung
    forward, not merely visit it.

    Pinned through a consequence rather than by reaching into the C++: a
    cold-started rung at redundant D lands on a state whose dominant fixed
    point is degenerate, so before `vms_grow_init` this raised for every
    D>=2 while D=1 answered exactly. Asserting the sequence rather than one
    point is what distinguishes "the ramp works" from "D=2 happened to get
    lucky"."""
    for D in (1, 2, 4, 8):
        ic = _polarized(3, 1, D, 3)
        assert ic.gs_energy() == pytest.approx(EXACT_E, abs=1e-9), D
