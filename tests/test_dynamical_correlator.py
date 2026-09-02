"""Functional coverage for Many_Body_Chain.get_dynamical_correlator,
distilled from examples/hubbard_chain and examples/dynamical_correlator/
dynamical_correlator_minimal down to a small, fast, deterministic
system. Uses the ED "INV" submode (the same one examples/hubbard_chain
compares its DMRG result against), which is exact and fast for small
systems -- no KPM polynomial-order tuning needed."""
import numpy as np
import pytest

from dmrgpy import cppext, spinchain

DELTA = 0.05

# itensor_version=2 exercises TDZ's real MPO-Taylor-fallback path (see
# test_tdz_dynamical_correlator_peak_matches_exact_gap's own docstring) --
# when mpscpp2 isn't compiled, mode.py's own DMRG->ED fallback (see
# CLAUDE.md) silently redirects here instead, and ED has no TDZ
# implementation at all (a DMRG/TDVP-only complex-time-evolution method),
# so this needs the same skip other multi-backend test files already use
# for an optional/uncompiled extension (e.g. test_nh_dmrg.py,
# test_four_point_correlator.py), not just for itensor_version=3.
_TDZ_VERSIONS = [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2),
        reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    pytest.param("python", id="python"),
]

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


@pytest.mark.parametrize("itensor_version", [2, 3])
def test_ex_dynamical_correlator_peak_matches_exact_gap(itensor_version):
    """EX (dcex.py) builds the correlator from a small number of
    explicitly-computed DMRG excited states (Lehmann sum over a
    Lagrange-multiplier-penalty excited-state search, see excited.py/
    chain_session.h's Chain::excited_states) rather than KPM's Chebyshev
    expansion or CVM's resolvent linear solve. Like CVM, for a small
    system where the excited-state search essentially recovers the exact
    spectrum, its peak should land on the exact gap to within the
    frequency grid spacing -- this is the cross-backend consistency check
    for submode="EX" (previously untested, see dcex.py/excited.py).
    "python" is intentionally excluded here: pyitensor's excited-state
    search has no seeded start, and this test was intrinsically flaky
    on that backend (~10-30% failure rate, e.g. peak at ~0.77 vs the
    asserted 0.658919+-0.02, reproduced on an unmodified tree, i.e. not
    a regression) -- see nex=6 in excited.py's overlap-penalty search.
    The [2]/[3] variants have never been observed to flake."""
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


@pytest.mark.parametrize("itensor_version", _TDZ_VERSIONS)
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


@pytest.mark.parametrize("damping", ["exp", "gaussian", "parzen"])
def test_td_window_options_preserve_peak_position(damping):
    """submode="TD" (timedependent.py::dynamical_correlator) windows the
    raw time-domain correlator before the FFT via a selectable `damping`
    taper (see timedependent.py::_damping_window and
    docs/td_dynamical_correlator_sharpening_plan.md): "exp" (the
    pre-existing default, an exact Lorentzian broadening), "gaussian"
    (faster-decaying tail), and "parzen" (compact-support taper,
    suppresses truncation ringing). All three only change the *shape* of
    the peak, not its location -- each should still land on the exact
    gap to within the same tolerance the TDZ test above uses (it shares
    this same FFT tail via _fourier_transform_correlator)."""
    sc = _heisenberg_chain()
    sc.setup_python()  # cheap, always-available backend for this check
    name = (sc.Sz[0], sc.Sz[0])
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name,
                                        es=_PEAK_ES, delta=DELTA,
                                        damping=damping)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.03)


def test_td_linear_prediction_sharpens_peak():
    """Linear-prediction extrapolation (dynamicstk/linearprediction.py,
    enabled via dynamical_correlator's predict=True) extends the raw C(t)
    before windowing, so the same real TDVP simulation should yield a
    narrower resonance than the plain windowed FFT -- the entire point of
    the technique (see docs/td_dynamical_correlator_sharpening_plan.md).
    Checked on the same 4-site Heisenberg gap used throughout this file:
    both the plain and linear-prediction spectra must still peak at the
    exact gap, but the linear-prediction one must have a strictly
    narrower full width at half max."""
    sc = _heisenberg_chain()
    sc.setup_python()
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(0.3, 1.0, 281)  # finer grid to resolve FWHM reliably
    delta = DELTA

    x0, y0 = sc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name,
                                          es=es, delta=delta, predict=False)
    x1, y1 = sc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name,
                                          es=es, delta=delta, predict=True,
                                          lp_order=15, lp_extend_factor=8)
    x0, y0 = np.array(x0), np.abs(y0)
    x1, y1 = np.array(x1), np.abs(y1)

    peak0, peak1 = x0[np.argmax(y0)], x1[np.argmax(y1)]
    assert peak0 == pytest.approx(HEISENBERG_4_GAP, abs=0.03)
    assert peak1 == pytest.approx(HEISENBERG_4_GAP, abs=0.03)

    def fwhm(x, y):
        i0 = np.argmax(y)
        half = y[i0]/2
        above = x[y >= half]
        return above.max() - above.min()

    assert fwhm(x1, y1) < fwhm(x0, y0)


def test_finite_T_dynamical_correlator_independent_of_es_window():
    """get_dynamical_correlator(mode="ED", submode="ED", T=...) (the
    finite-temperature Lehmann sum, edtk/dynamics.py's
    dynamical_correlator_finite_T) must return the same values at a given
    frequency regardless of how far the requested `es` window extends,
    since the *initial*-state Boltzmann sum is physically independent of
    which frequencies the caller happens to be asking for.

    Regression test for a real bug (caught by code review): an earlier
    version cropped the initial-state index space to `es`'s own window
    too (not just the final-state space), which silently dropped any
    thermally-populated state above the window's max energy from the
    numerator while the partition function Z (from the full spectrum)
    still counted it in the denominator -- a systematic under-count
    growing with how much Boltzmann weight got left out. At T=2.0 on this
    4-site chain, an es window capped at 1.0 only captures ~41% of the
    total Boltzmann weight (the bulk of the spectrum sits above that),
    so the old bug would have inflated every value here by roughly
    1/0.41 =~ 2.4x relative to a window wide enough to capture ~100%."""
    n = 4
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)

    name = (sc.Sz[0], sc.Sz[0])
    T = 2.0
    es_narrow = np.linspace(-1.0, 1.0, 41)
    es_wide = np.linspace(-1.0, 6.0, 141)

    x1, y1 = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name,
                                          T=T, es=es_narrow, delta=0.05)
    x2, y2 = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name,
                                          T=T, es=es_wide, delta=0.05)
    y2_on_x1 = np.interp(x1, x2, y2.real)

    assert y1.real == pytest.approx(y2_on_x1, abs=1e-6)


# ---------------------------------------------------------------------
# The ED Lehmann sum is exact: no eigenstate truncation, and therefore
# no dependence on where the spectrum sits or how wide `es` is.
#
# edtk/dynamics.py used to crop the final-state sum to `emu < max(es)`,
# comparing absolute eigenvalues against a frequency measured from the
# ground state. Symptoms, all reproduced before the fix: a constant added
# to H changed the answer (1.1e-3 at H+3*Id on this chain) though a
# constant can move no pole; a spectrum sitting entirely above max(es)
# cropped itself away and raised numba's "zero-size array to reduction"
# from inside dynamical_sum; and the answer depended on the width of the
# requested window. Both dynamical_correlator_ED and
# dynamical_correlator_finite_T carried it.
# ---------------------------------------------------------------------


def _shifted_chain(shift, n=6):
    from dmrgpy import multioperator
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h + shift * multioperator.identity())
    return sc


# +10 and +500 both cropped the entire spectrum away before the fix, so
# these are crash regressions as well as numerical ones.
@pytest.mark.parametrize("shift", [1.0, 3.0, -3.0, 10.0, 500.0, -500.0])
@pytest.mark.parametrize("T", [None, 0.5])
def test_ed_correlator_invariant_under_a_constant_shift_of_H(shift, T):
    """H -> H + c*Id shifts every eigenvalue by c and every *excitation*
    energy by nothing, so a correlator -- whose poles live at excitation
    energies -- must be completely unchanged. Covers the T=0 sum and the
    finite-temperature one, which had the same crop."""
    es = np.linspace(0.05, 3.0, 120)

    def run(c):
        kwargs = dict(name=None, mode="ED", submode="ED", es=es, delta=0.05)
        sc = _shifted_chain(c)
        kwargs["name"] = (sc.Sz[0], sc.Sz[0])
        if T is not None:
            kwargs["T"] = T
        return sc.get_dynamical_correlator(**kwargs)[1].real

    base = run(0.0)
    scale = np.max(np.abs(base))
    assert scale > 1e-6  # the comparison would be vacuous otherwise
    # the residual is float cancellation from adding +-500 to eigenvalues
    # of order 1, not truncation; it grows with |shift| and stays ~1e-12
    assert np.max(np.abs(run(shift) - base)) / scale < 1e-10


@pytest.mark.parametrize("T", [None, 0.5])
def test_ed_correlator_independent_of_the_es_window_at_fixed_frequency(T):
    """The T=0 counterpart of the finite-T window test above: widening
    `es` must not change the value returned at a frequency both windows
    contain. The crop made max(es) a truncation threshold, so it did."""
    sc = _shifted_chain(0.0)
    name = (sc.Sz[0], sc.Sz[0])
    kw = dict(name=name, mode="ED", submode="ED", delta=0.05)
    if T is not None:
        kw["T"] = T
    # a shared step, so the narrow grid is an exact subset of the wide one
    # and the two are compared point for point rather than through an
    # interpolation whose own error would swamp what is being tested
    step = 0.025
    x1, y1 = sc.get_dynamical_correlator(es=0.05 + step*np.arange(40), **kw)
    x2, y2 = sc.get_dynamical_correlator(es=0.05 + step*np.arange(320), **kw)
    assert x1 == pytest.approx(x2[:len(x1)], abs=1e-12)  # subset, as intended
    assert y1.real == pytest.approx(y2.real[:len(y1)], abs=1e-9)


def test_ed_correlator_matches_a_textbook_lehmann_sum():
    """Against an independent reference written straight from the
    spectral representation -- full eigendecomposition, every final state,
    no dmrgpy code in the sum itself. This is what "submode='ED' is the
    exact Lehmann sum" is supposed to mean, and what the truncation
    quietly made untrue."""
    from dmrgpy.algebra import algebra
    from dmrgpy.edtk import dynamics as edyn

    sc = _shifted_chain(0.0)
    es, delta = np.linspace(0.05, 3.0, 120), 0.05
    y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), mode="ED",
                                    submode="ED", es=es, delta=delta)[1]

    ed = sc.get_ED_obj()
    a0 = edyn.EDOperator(sc.Sz[0], ed).SO
    emu, vs = algebra.eigh(ed.get_hamiltonian())
    ex = emu - emu.min()
    Am = np.conjugate(vs.T) @ a0 @ vs
    ref = np.zeros(len(es), dtype=complex)
    for j in range(len(ex)):  # non-degenerate ground state: i = 0 only
        ref += Am[0, j] * Am[j, 0] * (1.0 / (es + 1j * delta - ex[j])
                                      - 1.0 / (es - 1j * delta - ex[j]))
    ref = -ref.imag / (2 * np.pi)
    assert np.max(np.abs(y.real - ref)) / np.max(np.abs(ref)) < 1e-12


def test_realify_spin_terms_matches_verbatim_ampo():
    """The Sy -> (S+ - S-)/(2i) rewrite every mpscpp3 AutoMPO build goes
    through (mo_terms.h's build_ampo(), see docs/documentation.md 4.8a) is a
    pure representation change: it exists so a real Hamiltonian written the
    textbook dmrgpy way -- Sx.Sx + Sy.Sy + Sz.Sz -- yields a *real*-valued
    MPO instead of a complex one (ITensor picks the MPO's element type from
    the AutoMPO coefficients alone, and never notices that the `Sy` site
    matrix itself is purely imaginary), which makes every downstream
    applyMPO/sum/inner run in real rather than complex arithmetic -- ~3.4x
    on the L=30 KPM dynamical correlator.

    So it must not move any number. This pins that on the two dynamical-
    correlator submodes it was written for, plus the ground state and a
    genuinely-complex observable, by running each one twice with the
    rewrite enabled and disabled (the process-global
    set_realify_spin_terms() escape hatch) and comparing:

    * <Sy_i> alone has an *odd* number of Sy factors -- a genuinely
      imaginary operator -- so should_realify() must decline to rewrite it
      at all, and it must stay imaginary rather than silently becoming
      real;
    * <Sy_i Sy_j> is even, so it *is* rewritten, and must agree with ED.
    """
    cppext = pytest.importorskip("dmrgpy.cppext")
    if not cppext.available(3):
        pytest.skip("mpscpp3 extension not compiled")
    from dmrgpy.mpscpp3 import _dmrgcpp

    n = 6
    es = np.linspace(-1.0, 6.0, 80)

    def run():
        sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
        h = 0
        for i in range(n - 1):
            h = h + (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                     + 0.7 * sc.Sz[i] * sc.Sz[i + 1])
        for i in range(n):
            h = h + 0.3 * sc.Sx[i] + 0.11 * sc.Sz[i]
        sc.set_hamiltonian(h)
        name = (sc.Sz[2], sc.Sz[2])
        return dict(
            e0=sc.gs_energy(),
            e0_ed=sc.gs_energy(mode="ED"),
            syy=sc.vev(sc.Sy[2] * sc.Sy[3]),
            syy_ed=sc.vev(sc.Sy[2] * sc.Sy[3], mode="ED"),
            sy=sc.vev(sc.Sy[2]),
            kpm=sc.get_dynamical_correlator(name=name, submode="KPM",
                                             es=es, delta=0.3)[1],
            td=sc.get_dynamical_correlator(name=name, submode="TD",
                                            es=es, delta=0.3)[1],
        )

    try:
        _dmrgcpp.set_realify_spin_terms(True)
        assert _dmrgcpp.get_realify_spin_terms()
        on = run()
        _dmrgcpp.set_realify_spin_terms(False)
        off = run()
    finally:
        _dmrgcpp.set_realify_spin_terms(True)  # restore the default

    # the rewritten path is still exactly the same physics as ED
    assert on["e0"] == pytest.approx(on["e0_ed"], abs=1e-8)
    assert on["syy"].real == pytest.approx(on["syy_ed"].real, abs=1e-6)

    # ...and identical to what the verbatim (complex-MPO) build produces
    assert on["e0"] == pytest.approx(off["e0"], abs=1e-8)
    assert on["syy"] == pytest.approx(off["syy"], abs=1e-6)
    assert on["kpm"] == pytest.approx(off["kpm"], abs=1e-6)
    assert on["td"] == pytest.approx(off["td"], abs=1e-6)

    # an odd Sy count is left alone: <Sy> stays a (vanishing here, but
    # genuinely imaginary-typed) quantity on both paths
    assert abs(on["sy"].real) < 1e-8 and abs(off["sy"].real) < 1e-8
    assert on["sy"] == pytest.approx(off["sy"], abs=1e-6)
