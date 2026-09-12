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
point is meaningful. Rejecting the benign case too makes `gs_energy()`
raise "every attempt at D=... failed" for a polarized chain at `maxm>1`,
and every solver here has been through that:

* `itensor_version=3` raised on the sequential solver always, and on the
  grouped one whenever an iteration happened to land exactly on the
  degeneracy (measured at a second eigenvalue of 0.99996, just outside the
  guard's own 1e-9, i.e. it survived by luck rather than by design).
* `itensor_version="python"` raised intermittently on the GROUPED path
  even after its sequential one was fixed -- 7 of 20 runs of the D=4 case
  below, every failure traced to `vumps._environments`' two fixed-point
  calls and to nowhere else in that module. It is 0 of 20 now.

All four environment builders (`vumps.py` and `vumps_ms.py` on the Python
side, `Chain::vx_*`'s grouped and sequential halves on the C++ one) now
prefer the fixed points the state itself names -- `C C^dag` for the AL
transfer's right one and `conj(C^dag C)` for the AR transfer's left one,
in this codebase's `X[ket, bra]` ordering -- whenever those reproduce
themselves under the transfer map: an exact algebraic identity in mixed
canonical gauge, and so a yes/no test rather than a tuned threshold, with
the guarded eigensolver used only when it does not hold. See
`Chain::vx_bond_fixed_points` and `vumps._transfer_fixed_points` for why
that is the right element of the degenerate subspace in both the benign
and the pathological case, and for why a threshold on `C`'s own weight
spectrum (the shape tried first) is not.

The C++ half of that sentence was not true until 2026-09-12: both
`Chain::` builders reached `vx_bond_fixed_points` only from a `catch`,
i.e. as a FAILURE fallback, so an arbitrary element of a merely
NEAR-degenerate subspace -- which trips no guard -- was accepted and its
energy returned. `Chain::vx_choose_fixed_point` is the shared
bond-candidate-first selection that closed it, and
`docs/known_issue_v3_vumps_variational_floor.md` is the record of what it
looked like while open. Fixing it also surfaced a second, latent defect
in `vx_bond_fixed_points` itself: its LEFT candidate was `C^dag C` where
the index ordering needs the conjugate, invisible for as long as nothing
measured the candidate and for as long as the only models it was
exercised on (a polarized chain, AKLT) had a real symmetric `C^dag C`.

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
@pytest.mark.parametrize("D", [2, 4, 6, 8])
def test_grouped_solver_tolerates_redundant_bond_dimension(backend, D):
    """A reach-1 chain on a 1-site cell: the GROUPED path, at four bond
    dimensions the exact (D=1) state does not need.

    What this pins is a fallback, not luck. `vumps._transfer_fixed_points`
    hands `_environments` the fixed points `C` names whenever they
    reproduce themselves, so the degenerate dominant eigenvalue this state
    always has never has to be resolved by an eigensolver at all. Before
    that, D=4 here raised in 7 of 20 runs on `itensor_version="python"`
    (and this test was correspondingly flaky in the suite); it is 0 of 20
    now.

    `itensor_version=3` reached the same ordering later (see the module
    docstring). In between, this test was itself intermittently red on
    that backend -- not raising, but returning an energy past the `abs=
    1e-9` below asserted here, in 2 of 80 runs at D=6 and 3 of 80 at D=8.
    Both are 0 of 80 with `Chain::vx_choose_fixed_point`, as is D=4, which
    had been 0 of 80 throughout.

    Both asserts matter: the value, because a wrongly-chosen element of
    the degenerate subspace carries its own energy, and the one-sided
    bound, because the observed failure mode of an unguarded eigensolver
    was an energy BELOW the exact variational minimum (-2.30 and -4382
    against -1.825), which a two-sided tolerance alone would also catch but
    which is worth naming as the thing being excluded."""
    e = _polarized(1, 1, D, backend).gs_energy()
    assert e == pytest.approx(EXACT_E, abs=1e-9), D
    assert e >= EXACT_E - 1e-9, D


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("D", [2, 4, 6, 8])
@pytest.mark.parametrize("n_uc,reach", [(1, 2), (2, 3), (3, 1)])
def test_sequential_solver_tolerates_redundant_bond_dimension(backend, n_uc,
                                                              reach, D):
    """The SEQUENTIAL path, reached three ways -- a coupling past the cell
    on a 1-site and on a 2-site cell, and a cell longer than 2 sites -- at
    four bond dimensions the exact (D=1) state does not need.

    On `itensor_version=3` every one of these raised before: at D=1 the
    ramp never grows and the answer was already right, and from D=2 up the
    first rung produced the redundant state that the guard rejected.

    Swept over D rather than checked at D=2 alone because the redundancy
    grows with the gap between the requested bond dimension and the
    state's own (which is 1 here): D=2 leaves one surplus direction, D=8
    leaves seven, and a fallback that only happened to work for the first
    of those would not be the fix. Both asserts are the same pair as in the
    grouped test above -- the exact value, and the one-sided variational
    bound that an arbitrary element of the degenerate subspace violates.

    `itensor_version="python"` satisfies both at every (n_uc, reach, D)
    here, worst case 3.6e-15 over 6 runs each. `itensor_version=3` does
    too, worst case 8.5e-14.

    The (1, 2) cell at D>2 carried a narrow, non-strict `xfail` on
    `backend == 3` until 2026-09-12, for the C++ half of the same
    wrong-fixed-point defect -- `Chain::vms_environments` reached
    `vx_bond_fixed_points` only from a `catch`, so a near-degenerate
    eigenvector that did not trip the guard was accepted silently and the
    energy came back BELOW the exact variational minimum (measured 6 of 30
    runs at D=4, worst 6.1e-04). `vx_choose_fixed_point` put both C++
    builders on the same bond-candidate-first ordering the two Python ones
    already had; that cell is 0 of 30 at every one of D=4, 6 and 8
    afterwards, and the marker is gone rather than loosened."""
    e = _polarized(n_uc, reach, D, backend).gs_energy()
    assert e == pytest.approx(EXACT_E, abs=1e-9), (n_uc, reach, D)
    assert e >= EXACT_E - 1e-9, (n_uc, reach, D)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("n_uc,reach", [(1, 1), (1, 2)])
def test_variational_bound_holds_over_repeated_runs(backend, n_uc, reach):
    """The same one-sided bound as the two tests above, but over REPEATED
    runs of the same cell -- which is what it takes to catch this defect
    rather than to catch it 20% of the time.

    Nothing here is reproducible run to run: the solver starts from a
    random MPS (`vumps_random_init`, and `pyitensor`'s own equivalent), so
    whether a given run lands on the degenerate fixed point is a property
    of that run's start. The C++ excursions this file's own history
    records were 6 of 30 runs on the (1, 2) cell and 13 of 80 on the
    (1, 1) one -- so a single-run test, which is what the two tests above
    are, had a ~20% and ~16% chance of seeing it. Eight runs raise that to
    ~83% and ~74%, which is why this exists as its own test instead of
    tightening those.

    Both cells are here because the two solvers reach their environments
    through different code: (1, 1) is reach-1 on a one-site cell, so the
    GROUPED builder, and (1, 2) has a coupling past the cell, so the
    SEQUENTIAL one. `D=4` is the cheapest bond dimension above the exact
    state's own (which is 1) at which the C++ half was measured failing.

    The bound is the assert that matters here, not the value: a converged
    run can legitimately sit slightly ABOVE the exact minimum, while
    nothing legitimate sits below it."""
    for run in range(8):
        e = _polarized(n_uc, reach, 4, backend).gs_energy()
        assert e >= EXACT_E - 1e-9, (n_uc, reach, run, e)


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
    sequential), and both backends because the four environment builders
    reach `C C^dag` by four separate implementations of the same rule
    (`vumps._transfer_fixed_points`, `vumps_ms._cell_fixed_points`, and the
    grouped and sequential halves of the C++ `vx_*`): -2/3 on all of them is
    what says they agree where an arbitrary pick would have made them
    differ."""
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
