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


def _setup_backend(sc, itensor_version):
    if itensor_version == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(itensor_version)


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
