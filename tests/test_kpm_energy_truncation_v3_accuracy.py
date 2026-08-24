"""End-to-end coverage of ITensor v3's native KPM energy truncation
(mpscpp3/chain_session.h's kpm_dynamical_correlator_truncated(), Holzner
et al., PRB 83, 195115 (2011), Sec. III-B) through the public
Many_Body_Chain.get_dynamical_correlator(submode="KPM") API -- the C++
v3 counterpart of test_kpm_energy_truncation_accuracy.py's pyitensor
coverage. See test_kpm_energy_truncation_v3.py for coverage of the
lower-level scaled_hamiltonian_gs_anchored()/kpm_energy_truncate()
primitives themselves.

test_narrow_kpm_scale_without_truncation_raises below exercises
mpscpp3's check_kpm_moment() divergence guard directly -- this used to
hard-abort the whole process (SIGABRT) instead of raising a catchable
Python exception (ITensor's Error() macro prints and calls abort()
unconditionally, never actually throwing, despite ITError being a real
std::runtime_error-derived type), a real, pre-existing bug found while
building this test file, fixed by changing check_kpm_moment to
`throw ITError(...)` directly (see chain_session.h's own comment on that
one call site, and test_kpm_divergence_guard_catchable.py for a focused
regression test of the fix across all three backends). Confirmed the
same bug existed in mpscpp2 too and was fixed there identically -- both
backends' check_kpm_moment() are otherwise unrelated to energy
truncation itself.

Same 4-site Heisenberg chain and exact gap as
test_dynamical_correlator.py's KPM/CVM/EX peak tests (golden value
0.658919 from test_excited_states.py) and the pyitensor accuracy test.
"""
import numpy as np
import pytest

from dmrgpy import cppext, spinchain

pytestmark = pytest.mark.skipif(
    not cppext.available(3), reason="requires the compiled mpscpp3 (ITensor v3) extension")

DELTA = 0.05
HEISENBERG_4_GAP = 0.658919
_PEAK_ES = np.linspace(0.3, 1.0, 71)

# See test_kpm_energy_truncation_accuracy.py's own module docstring for
# why 0.65 (not something much smaller): narrow enough to require energy
# truncation under the ground-state-anchored rescaling, but not so
# narrow that it starts cutting into this small system's own correlator
# weight (confirmed directly, mirroring the pyitensor findings -- this
# system's own physics doesn't depend on which backend computes it).
_NARROW_KPM_SCALE = 0.65


def _heisenberg_chain(n=4):
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    sc.setup_cpp(3)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc


def test_narrow_kpm_scale_without_truncation_raises():
    """Baseline: today's guard rail. Enabling kpm_energy_truncate alone
    switches to the ground-state-anchored rescaling convention (see
    module docstring) but this test forces the truncation sweeps
    themselves off (kpm_truncate_nsweeps=0), isolating that the
    convention switch by itself does not avoid the leakage --
    chain_session.h's check_kpm_moment must still catch it with a
    catchable RuntimeError (see module docstring for the fix that made
    this safe to test at all) -- this is the failure real energy
    truncation (the next test) exists to fix, not a new behavior being
    introduced here."""
    sc = _heisenberg_chain()
    sc.kpm_scale = _NARROW_KPM_SCALE
    sc.kpm_energy_truncate = True
    sc.kpm_truncate_nsweeps = 0
    name = (sc.Sz[0], sc.Sz[0])
    with pytest.raises(RuntimeError, match="KPM moments diverging"):
        sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                     es=_PEAK_ES, delta=DELTA)


def test_narrow_kpm_scale_with_truncation_matches_exact_gap():
    """A window narrow enough to need energy truncation must still
    locate the exact excitation gap once truncation is actually
    enabled."""
    sc = _heisenberg_chain()
    sc.kpm_scale = _NARROW_KPM_SCALE
    sc.kpm_energy_truncate = True
    sc.kpm_truncate_dK = 30
    sc.kpm_truncate_nsweeps = 10
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                                        es=_PEAK_ES, delta=DELTA)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.05)


def test_truncated_v3_matches_untruncated_python_backend():
    """Cross-backend check: the same narrow kpm_scale, truncation
    enabled, on both the native v3 port and the pyitensor port of the
    same algorithm, must locate the same peak (within the frequency
    grid's own spacing -- the two backends' independent numerics can
    land on an adjacent grid point even when the underlying spectral
    functions agree closely) -- confirming the C++ port and the Python
    port of Sec. III-B agree, not just that each independently looks
    reasonable."""
    sc_v3 = _heisenberg_chain()
    sc_v3.kpm_scale = _NARROW_KPM_SCALE
    sc_v3.kpm_energy_truncate = True
    sc_v3.kpm_truncate_dK = 30
    sc_v3.kpm_truncate_nsweeps = 10
    _, y_v3 = sc_v3.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                              name=(sc_v3.Sz[0], sc_v3.Sz[0]),
                                              es=_PEAK_ES, delta=DELTA)
    peak_v3 = _PEAK_ES[np.argmax(np.array(y_v3).real)]

    sc_py = spinchain.Spin_Chain(["S=1/2" for _ in range(4)])
    sc_py.setup_python()
    h = 0
    for i in range(3):
        h = h + sc_py.Sx[i] * sc_py.Sx[i + 1] + sc_py.Sy[i] * sc_py.Sy[i + 1] + sc_py.Sz[i] * sc_py.Sz[i + 1]
    sc_py.set_hamiltonian(h)
    sc_py.kpm_scale = _NARROW_KPM_SCALE
    sc_py.kpm_energy_truncate = True
    sc_py.kpm_truncate_dK = 30
    sc_py.kpm_truncate_nsweeps = 10
    _, y_py = sc_py.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                              name=(sc_py.Sz[0], sc_py.Sz[0]),
                                              es=_PEAK_ES, delta=DELTA)
    peak_py = _PEAK_ES[np.argmax(np.array(y_py).real)]

    assert peak_v3 == pytest.approx(peak_py, abs=0.02)
    assert peak_v3 == pytest.approx(HEISENBERG_4_GAP, abs=0.05)


def test_energy_truncation_matches_ed_qualitative_shape():
    """Cross-check against the exact ED "INV" spectral function (same
    ground truth test_dynamical_correlator_peaks_at_excitation_gap
    uses): the resonance must be the largest of the three sampled
    points, matching ED's own qualitative shape (see the pyitensor
    accuracy test's docstring for why this is qualitative, not a strict
    dominance-factor bound, on this small a system)."""
    sc_ed = _heisenberg_chain()
    name_ed = (sc_ed.Sz[0], sc_ed.Sz[0])
    es = np.array([0.0, HEISENBERG_4_GAP, 1.5])
    _, y_ed = sc_ed.get_dynamical_correlator(mode="ED", submode="INV",
                                              name=name_ed, es=es, delta=DELTA)

    sc = _heisenberg_chain()
    sc.kpm_scale = _NARROW_KPM_SCALE
    sc.kpm_energy_truncate = True
    sc.kpm_truncate_dK = 30
    sc.kpm_truncate_nsweeps = 10
    name = (sc.Sz[0], sc.Sz[0])
    _, y_kpm = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                            name=name, es=es, delta=DELTA)

    y_ed, y_kpm = y_ed.real, y_kpm.real
    assert y_ed[1] > 10 * y_ed[0]
    assert y_ed[1] > 10 * y_ed[2]
    assert y_kpm[1] > y_kpm[0]
    assert y_kpm[1] > y_kpm[2]


def test_six_site_truncated_spectra_agree_across_backends():
    """The strongest form of the cross-backend check above, on a chain
    big enough for the paper's recommended dK=30 to be smaller than the
    local effective Hamiltonians it is applied to (a 4-site chain's local
    spaces are only 16-dimensional, so dK=30 spans them completely and
    both ports trivially agree there): compare the *whole* reconstructed
    spectrum, not just its peak.

    The two implementations run the identical algorithm on the identical
    rescaled Hamiltonian, so they agree far beyond what this asserts
    (measured 6e-12 relative L2 difference); the tolerance is loose only
    to leave room for the two backends' independent DMRG ground states.

    Regression test for a real bug: pyitensor's _local_krylov_projection
    orthogonalized each new Krylov vector against the accepted ones only
    once, so its basis lost orthogonality at dK=30 while mpscpp3's port
    (two passes) stayed correct -- the two backends' spectra then
    disagreed at the O(1) level, with the pyitensor peak landing at
    ~1.06 J instead of ~0.48 J."""
    es = np.linspace(0.2, 1.6, 80)
    curves = {}
    for version in ("python", 3):
        np.random.seed(0)
        sc = spinchain.Spin_Chain(["S=1/2" for _ in range(6)])
        if version == "python":
            sc.setup_python()
        else:
            sc.setup_cpp(version=3)
        h = 0
        for i in range(5):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        sc.set_hamiltonian(h)
        sc.kpm_scale = _NARROW_KPM_SCALE
        sc.kpm_energy_truncate = True
        sc.kpm_truncate_dK = 30
        sc.kpm_truncate_nsweeps = 2
        _, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                            name=(sc.Sz[0], sc.Sz[0]),
                                            es=es, delta=DELTA)
        curves[version] = np.array(y).real

    ref = curves[3]
    diff = np.linalg.norm(curves["python"] - ref) / np.linalg.norm(ref)
    assert diff < 1e-3, "backends disagree by %.3e in relative L2" % diff
