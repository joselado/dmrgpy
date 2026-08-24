"""Ground-state bond-dimension ramp (Many_Body_Chain.bond_ramp).

The ramp changes only the DMRG *sweep schedule* -- it grows the bond
dimension geometrically over the first (at most half of the) sweeps
instead of running every sweep at the full self.maxm, and decays the
noise term once the ramp is done (see manybodychain.py's own comment,
mpscppN/chain_session.h's Chain::make_sweeps_ramped() and
pyitensor/chain.py's _make_sweeps_ramped()). It must therefore not
change any result, only the time taken to get it.

These tests target the two ways that could silently go wrong -- both
produce a plausible-looking wrong number rather than a crash:

* the ramp truncating a *warm* starting wavefunction (an MPS handed in
  by set_initial_wf, or simply the previous gs_energy() call's own
  converged solution) back down to the ramp's starting bond dimension;
* gs_energy_single()'s Hamiltonian send-cache not keying on the ramp
  settings, so flipping self.bond_ramp between two otherwise identical
  gs_energy() calls returns the *cached* energy from the previous
  schedule instead of re-solving.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain

TOL = 1e-6

BACKENDS = [2, 3, "python"]


def heisenberg(n=8, itensor_version=3, maxm=30, nsweeps=12):
    """Uniform S=1/2 Heisenberg chain, small enough to cross-check ED."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1]
        h = h + sc.Sy[i] * sc.Sy[i + 1]
        h = h + sc.Sz[i] * sc.Sz[i + 1]
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("version", BACKENDS)
def test_ramped_ground_state_matches_ed(version):
    """The ramp is a scheduling change: the ground-state energy must
    still be the exact one, on every backend."""
    sc = heisenberg(itensor_version=version)
    sc.bond_ramp = True
    e_dmrg = sc.gs_energy()
    e_ed = heisenberg(itensor_version=version).gs_energy(mode="ED")
    assert e_dmrg == pytest.approx(e_ed, abs=TOL)


@pytest.mark.parametrize("version", BACKENDS)
def test_ramp_and_flat_schedules_agree(version):
    """Ramped and un-ramped schedules must converge to the same energy
    *and* the same state (checked through a local observable, which a
    merely-close energy does not guarantee)."""
    sc_ramp = heisenberg(itensor_version=version)
    sc_ramp.bond_ramp = True
    e_ramp = sc_ramp.gs_energy()

    sc_flat = heisenberg(itensor_version=version)
    sc_flat.bond_ramp = False
    e_flat = sc_flat.gs_energy()

    assert e_ramp == pytest.approx(e_flat, abs=TOL)
    for i in range(sc_ramp.ns):
        assert sc_ramp.vev(sc_ramp.Sz[i]).real == pytest.approx(
            sc_flat.vev(sc_flat.Sz[i]).real, abs=1e-5)


@pytest.mark.parametrize("version", BACKENDS)
def test_toggling_the_ramp_invalidates_the_energy_cache(version):
    """gs_energy_single() only re-sends the Hamiltonian to the session
    when its cache key changed, and the session keeps its own energy
    cache across a skipped re-send. The ramp settings are part of that
    key (groundstate.py's ramp_key()) -- without them, flipping
    bond_ramp between two bare gs_energy() calls would hand back the
    previous schedule's cached energy, silently, forever.
    """
    sc = heisenberg(itensor_version=version)
    sc.bond_ramp = False
    e_flat = sc.gs_energy()

    from dmrgpy.groundstate import ramp_key
    key_before = ramp_key(sc)
    sc.bond_ramp = True
    assert ramp_key(sc) != key_before  # the cache key really did change

    sc.computed_gs = False
    sc.skip_dmrg_gs = False
    e_ramp = sc.gs_energy()
    assert e_ramp == pytest.approx(e_flat, abs=TOL)


@pytest.mark.parametrize("version", BACKENDS)
def test_warm_start_is_not_truncated_by_the_ramp(version):
    """Re-running gs_energy() on an already-converged chain must not
    make the energy worse. The ramp's early sweeps run far below the
    target bond dimension, so ramping down from a warm, converged state
    would truncate it -- the ramp floors its starting bond dimension at
    whatever the incoming wavefunction already carries (see
    Chain::gs_energy()'s floor_dim / _make_sweeps_ramped()).
    """
    sc = heisenberg(itensor_version=version, maxm=40)
    sc.bond_ramp = True
    e_first = sc.gs_energy()

    # force a real second solve, warm-started from the converged state
    sc.computed_gs = False
    sc.skip_dmrg_gs = False
    sc._session_ham_cache = None
    e_second = sc.gs_energy()

    # variational: a warm restart may only improve the energy, never
    # regress it beyond numerical noise
    assert e_second <= e_first + TOL
    assert e_second == pytest.approx(e_first, abs=TOL)


@pytest.mark.parametrize("version", BACKENDS)
def test_ramp_disabled_still_works(version):
    """bond_ramp=False must reproduce the original flat schedule."""
    sc = heisenberg(itensor_version=version)
    sc.bond_ramp = False
    e_ed = heisenberg(itensor_version=version).gs_energy(mode="ED")
    assert sc.gs_energy() == pytest.approx(e_ed, abs=TOL)


@pytest.mark.parametrize("version", BACKENDS)
def test_ramp_start_above_maxm_is_a_flat_schedule(version):
    """A ramp whose starting bond dimension already exceeds the target
    has nothing to ramp; it must degrade to the flat schedule rather
    than producing a schedule that overshoots self.maxm."""
    sc = heisenberg(itensor_version=version, maxm=20)
    sc.bond_ramp = True
    sc.bond_ramp_start = 200  # >> maxm
    e_ed = heisenberg(itensor_version=version).gs_energy(mode="ED")
    assert sc.gs_energy() == pytest.approx(e_ed, abs=TOL)


def test_python_backend_ramp_schedule_shape():
    """The schedule itself: the bond dimension grows at (nearly) every
    ramping sweep, starts at bond_ramp_start, reaches the full target
    exactly when the ramp ends and never overshoots it, and the noise
    term is off from there on (so the final, converged sweeps are
    noise-free).

    Checked on the pure-Python backend, whose Sweeps object is directly
    inspectable; mpscppN/chain_session.h's make_sweeps_ramped() is a
    line-by-line counterpart of the same function.
    """
    from dmrgpy.pyitensor.chain import Chain

    ch = Chain([2] * 8)
    ch.set_sweep_params(120, 16, 1e-12, 1e-7)
    ch.set_bond_ramp(True, 10, 0.5, 0.1)
    sweeps = ch._make_sweeps_ramped(16, 120, 0)

    dims = [sweeps.at(i)[0] for i in range(16)]
    noise = [sweeps.at(i)[2] for i in range(16)]

    assert dims[0] == 10
    assert dims == sorted(dims)  # non-decreasing
    assert max(dims) == 120  # never overshoots the target...
    assert dims[16 // 2] == 120  # ...and reaches it at the half-way point
    assert all(d == 120 for d in dims[8:])
    # the point of a fraction-based ramp: the whole budget is used, so
    # most ramping sweeps really do run below the target
    assert sum(1 for d in dims if d < 120) == 8
    assert noise[0] == pytest.approx(1e-7)
    assert all(nz == 0.0 for nz in noise[8:])
    # noise decays monotonically across the ramp
    assert noise[:8] == sorted(noise[:8], reverse=True)


@pytest.mark.parametrize("fraction,expected_cheap", [(0.25, 4), (0.5, 8),
                                                     (0.75, 12)])
def test_ramp_fraction_controls_how_much_of_the_schedule_ramps(
        fraction, expected_cheap):
    """bond_ramp_fraction is the fraction of the sweeps spent below the
    target bond dimension."""
    from dmrgpy.pyitensor.chain import Chain

    ch = Chain([2] * 8)
    ch.set_sweep_params(120, 16, 1e-12, 1e-7)
    ch.set_bond_ramp(True, 10, fraction, 0.1)
    dims = [ch._make_sweeps_ramped(16, 120, 0).at(i)[0] for i in range(16)]
    assert sum(1 for d in dims if d < 120) == expected_cheap
    assert dims[-1] == 120  # the schedule always ends at the full target


def test_warm_start_floor_disables_the_early_truncation():
    """The ramp floor, directly: a schedule built with floor_dim=64 must
    never dip below 64, while the same schedule with floor_dim=0 does."""
    from dmrgpy.pyitensor.chain import Chain

    ch = Chain([2] * 8)
    ch.set_sweep_params(120, 16, 1e-12, 1e-7)
    ch.set_bond_ramp(True, 10, 0.5, 0.1)

    cold = [ch._make_sweeps_ramped(16, 120, 0).at(i)[0] for i in range(16)]
    warm = [ch._make_sweeps_ramped(16, 120, 64).at(i)[0] for i in range(16)]

    assert min(cold) < 64
    assert min(warm) >= 64


def test_inhomogeneous_hubbard_ramp_matches_flat():
    """A spatially inhomogeneous Heisenberg-Hubbard chain -- the model
    examples/groundstate/bond_dimension_ramp uses at n=30, here at a
    size the test suite can afford. Site-dependent hopping, Hubbard U
    and Heisenberg exchange, so no two bonds are equivalent and the
    ground state has real spatial structure for the ramp to reproduce.
    """
    n = 6

    def build(ramp):
        fc = fermionchain.Spinful_Fermionic_Chain(n)
        h = 0
        for i in range(n - 1):
            t = 1.0 + 0.3 * np.cos(2 * np.pi * i / 7.0)
            J = 0.4 + 0.3 * np.sin(2 * np.pi * i / 5.0)
            h = h + t * fc.Cdagup[i] * fc.Cup[i + 1]
            h = h + t * fc.Cdagdn[i] * fc.Cdn[i + 1]
            h = h + J * (fc.Sx[i] * fc.Sx[i + 1] + fc.Sy[i] * fc.Sy[i + 1]
                         + fc.Sz[i] * fc.Sz[i + 1])
        for i in range(n):
            U = 2.0 + 1.5 * np.cos(2 * np.pi * i / 11.0)
            h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
        h = h + h.get_dagger()
        h = 0.5 * h
        fc.maxm = 40
        fc.nsweeps = 14
        fc.bond_ramp = ramp
        fc.set_hamiltonian(h)
        return fc

    fc_ramp, fc_flat = build(True), build(False)
    e_ramp, e_flat = fc_ramp.gs_energy(), fc_flat.gs_energy()
    assert e_ramp == pytest.approx(e_flat, abs=1e-5)

    d_ramp = np.array(fc_ramp.get_density()).real
    d_flat = np.array(fc_flat.get_density()).real
    assert np.max(np.abs(d_ramp - d_flat)) < 1e-4
