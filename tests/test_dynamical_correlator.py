"""Functional coverage for Many_Body_Chain.get_dynamical_correlator,
distilled from examples/hubbard_chain and examples/dynamical_correlator/
dynamical_correlator_minimal down to a small, fast, deterministic
system. Uses the ED "INV" submode (the same one examples/hubbard_chain
compares its DMRG result against), which is exact and fast for small
systems -- no KPM polynomial-order tuning needed."""
import numpy as np
import pytest

from dmrgpy import spinchain

DELTA = 0.05

# 4-site Heisenberg chain: ground state is a non-degenerate singlet, and
# the first excited state is a 3-fold degenerate triplet at this exact
# gap above it (E1-E0, both golden values from test_excited_states.py's
# test_excited_states_energies_and_ordering). Any local dynamical
# correlator built from a single-site operator (e.g. <Sz_0(t)Sz_0(0)>)
# must show its dominant resonance here.
HEISENBERG_4_GAP = 0.658919


def _heisenberg_chain(n=4):
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc


from _helpers import setup_backend as _setup_backend


def test_dynamical_correlator_peaks_at_excitation_gap():
    """4-site Heisenberg chain: the ground state (singlet) has a 3-fold
    degenerate triplet excited state at gap 0.658919 above it (see
    test_excited_states.py). The <Sz_0(t) Sz_0(0)> dynamical correlator
    must show a resonance there and be small away from it -- checked
    both qualitatively (peak >> off-resonance values) and against golden
    regression values at each frequency."""
    n = 4
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)

    name = (sc.Sz[0], sc.Sz[0])
    es = np.array([0.0, 0.658919, 1.5])
    x, y = sc.get_dynamical_correlator(mode="ED", submode="INV", name=name,
                                        es=es, delta=DELTA)

    expected = np.array([0.0066913407688090855, 1.045617999495007, 0.06866722208702186])
    assert y.real == pytest.approx(expected, abs=1e-6)

    # the resonance is the dominant feature by a wide margin
    assert y.real[1] > 10 * y.real[0]
    assert y.real[1] > 10 * y.real[2]


# Frequency window bracketing HEISENBERG_4_GAP with 0.01 spacing -- fine
# enough to locate the peak to within half a grid step, coarse enough to
# stay fast (submode="CVM" solves one linear system per frequency point).
_PEAK_ES = np.linspace(0.3, 1.0, 71)


@pytest.mark.parametrize("itensor_version", [2, 3, "python"])
def test_cvm_dynamical_correlator_peak_matches_exact_gap(itensor_version):
    """CVM (cvm.py, a Correction Vector Method solved via conjugate
    gradient) solves a resolvent linear system at each frequency, which
    for an exact ground state is numerically exact rather than a
    controlled approximation -- confirmed in
    examples/dynamical_correlator_VS_ED to match ED pointwise to
    ~1e-13..1e-14. So its peak should land on the exact gap to within
    the frequency grid spacing, on all three DMRG backends."""
    sc = _heisenberg_chain()
    _setup_backend(sc, itensor_version)
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="CVM", name=name,
                                        es=_PEAK_ES, delta=DELTA)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.02)


@pytest.mark.parametrize("itensor_version", [2, 3, "python"])
def test_kpm_dynamical_correlator_peak_matches_exact_gap(itensor_version):
    """KPM (kpmdmrg.py, Chebyshev moment expansion) is a genuine
    controlled approximation and does not reproduce ED's exact lineshape
    point-by-point (confirmed in examples/dynamical_correlator_VS_ED:
    raising the moment count doesn't shrink that pointwise gap, since
    KPM's own effective resolution changes with it). What it must still
    get right is the location of the physical excitation energy: its
    dominant peak should land near the exact gap, on all three DMRG
    backends."""
    sc = _heisenberg_chain()
    _setup_backend(sc, itensor_version)
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                        es=_PEAK_ES, delta=DELTA)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.05)


@pytest.mark.parametrize("itensor_version", [2, 3, "python"])
def test_ex_dynamical_correlator_peak_matches_exact_gap(itensor_version):
    """EX (dcex.py) builds the correlator from a small number of
    explicitly-computed DMRG excited states (Lehmann sum over a
    Lagrange-multiplier-penalty excited-state search, see excited.py/
    chain_session.h's Chain::excited_states) rather than KPM's Chebyshev
    expansion or CVM's resolvent linear solve. Like CVM, for a small
    system where the excited-state search essentially recovers the exact
    spectrum, its peak should land on the exact gap to within the
    frequency grid spacing, consistently across all three DMRG backends
    -- this is the cross-backend consistency check for submode="EX"
    (previously untested, see dcex.py/excited.py)."""
    sc = _heisenberg_chain()
    _setup_backend(sc, itensor_version)
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="EX", name=name,
                                        es=_PEAK_ES, delta=DELTA, nex=6)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.02)


def test_ex_dynamical_correlator_peak_matches_kpm_and_cvm():
    """Benchmark submode="EX" directly against the other two DMRG
    dynamical-correlator submodes (KPM and CVM) on the same chain/
    operator/frequency grid, rather than only against the analytic gap:
    all three must locate their dominant peak at the same frequency,
    since they are three different numerical routes to the same
    physical spectral function."""
    sc = _heisenberg_chain()
    name = (sc.Sz[0], sc.Sz[0])
    peaks = {}
    for submode, kwargs in [("EX", {"nex": 6}), ("KPM", {}), ("CVM", {})]:
        x, y = sc.get_dynamical_correlator(mode="DMRG", submode=submode,
                                            name=name, es=_PEAK_ES,
                                            delta=DELTA, **kwargs)
        x, y = np.array(x), np.array(y).real
        peaks[submode] = x[np.argmax(y)]
    assert peaks["EX"] == pytest.approx(peaks["KPM"], abs=1e-9)
    assert peaks["EX"] == pytest.approx(peaks["CVM"], abs=1e-9)


@pytest.mark.parametrize("itensor_version", [2, 3, "python"])
def test_tdz_dynamical_correlator_peak_matches_exact_gap(itensor_version):
    """TDZ (tdz.py, complex-time evolution + perturbative real-axis
    reconstruction, arXiv:2311.10909) should locate its dominant peak at
    the exact gap to within the dt=0.1 (default) time-discretization
    error, consistently across all three DMRG backends: TDVP
    (itensor_version 3 or "python") and the MPO-Taylor fallback
    (itensor_version=2, which has no TDVP -- see tdz.py's
    _advance_complex_time_step). Confirmed directly (see this module's
    own dev notes) that all three backends land on the identical peak
    here, i.e. the TDVP-vs-Taylor-MPO choice does not itself introduce a
    cross-backend discrepancy at this alpha0/n_max/dt."""
    sc = _heisenberg_chain()
    _setup_backend(sc, itensor_version)
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=name,
                                        es=_PEAK_ES, delta=DELTA)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.03)
