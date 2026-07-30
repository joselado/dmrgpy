"""End-to-end coverage of KPM energy truncation (Holzner et al., PRB 83,
195115 (2011), Sec. III-B; src/dmrgpy/pyitensor/kpm_energy_truncation.py)
through the public Many_Body_Chain.get_dynamical_correlator(submode="KPM")
API -- see test_kpm_energy_truncation.py for coverage of energy_truncate()
itself, independent of the Chebyshev recursion it plugs into here.

Same 4-site Heisenberg chain and exact gap as
test_dynamical_correlator.py's KPM/CVM/EX peak tests (golden value
0.658919 from test_excited_states.py) so this is directly comparable to
that file's existing (kpm_scale=0.7, no truncation) KPM coverage.

Enabling kpm_energy_truncate also switches chain.py's
kpm_dynamical_correlator() from _scaled_hamiltonian() (which centers the
KPM window on the full many-body bandwidth's *midpoint*) to
_scaled_hamiltonian_gs_anchored() (which anchors it at the ground state
instead, per the paper's own Eq. (21b)) -- confirmed directly while
building this test: with the bandwidth-midpoint convention, the ground
state itself (which sits at the spectrum's *edge*, not its middle) falls
outside the window before the window is ever narrow enough to be useful,
so the two rescaling conventions need separate kpm_scale intuitions; see
_scaled_hamiltonian_gs_anchored()'s own docstring.

Also confirmed directly: on a system this small, a 4-site chain's local
operator (Sz on one site) has real spectral weight spread across a good
fraction of the *entire* many-body bandwidth -- there is no large
separation between the correlator's own spectral width and the full
bandwidth the way there would be for a longer chain (see the paper's own
Sec. II.B on why Ws<<W pays off for extensive systems). So this test
picks a *modestly* narrowed kpm_scale (already enough to require energy
truncation) rather than an aggressively small one: pushed further, the
window starts cutting into real physical weight of this particular
correlator and the reconstructed lineshape degrades (confirmed directly:
kpm_scale=0.5 already visibly distorts the lineshape here), which is an
expected small-system limitation of the technique, not a bug in the
truncation mechanism itself (that mechanism's correctness is covered
independently, and more rigorously, in test_kpm_energy_truncation.py).
"""
import numpy as np
import pytest

from dmrgpy import spinchain

DELTA = 0.05
HEISENBERG_4_GAP = 0.658919
_PEAK_ES = np.linspace(0.3, 1.0, 71)

# See module docstring: kpm_scale=0.65 already needs energy truncation
# once the ground-state-anchored convention is in effect (confirmed
# directly: 0.55-0.7 all raise without real truncation), but is not so
# narrow that it starts cutting into this correlator's own real spectral
# weight (confirmed directly: 0.6-0.7 reproduce the exact gap; 0.5
# already visibly distorts the lineshape, see module docstring).
_NARROW_KPM_SCALE = 0.65


def _heisenberg_chain(n=4):
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    sc.setup_python()
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
    convention switch by itself does not avoid the leakage -- chain.py's
    _check_kpm_moment must still catch it, exactly as it does today for
    any accidentally-too-tight kpm_scale. This is the failure real energy
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
    """With energy truncation actually enabled, the same narrow
    kpm_scale that raises above must instead compute cleanly and still
    locate the exact excitation gap."""
    np.random.seed(0)
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


def test_energy_truncation_matches_ed_qualitative_shape():
    """Cross-check against the exact ED "INV" spectral function (same
    ground truth test_dynamical_correlator_peaks_at_excitation_gap uses)
    at the resonance and at two off-resonance points. Only a qualitative
    check (resonance is the largest of the three, matching ED's own
    shape) rather than ED's own ">10x dominant" bar: see module docstring
    for why a narrowed, truncation-requiring window measurably broadens
    this particular small-system correlator's lineshape relative to the
    safe default."""
    sc_ed = _heisenberg_chain()
    name_ed = (sc_ed.Sz[0], sc_ed.Sz[0])
    es = np.array([0.0, HEISENBERG_4_GAP, 1.5])
    _, y_ed = sc_ed.get_dynamical_correlator(mode="ED", submode="INV",
                                              name=name_ed, es=es, delta=DELTA)

    np.random.seed(0)
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
