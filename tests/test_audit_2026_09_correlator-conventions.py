"""Regression tests for the 2026-09 audit's dynamical-correlator findings.

Three findings from `docs/audit_2026_09_hole_hunt.md`:

* **#5** -- the same submode name computed different quantities under
  `mode="ED"` and `mode="DMRG"`, and three of the routes sat off the
  convention the rest of them share.
* **#11** -- `submode="CVM_explicit"` returned exactly twice the spectral
  function it is documented to share with `submode="CVM"`, through an
  `np.abs()` that also discarded a phase.
* **#29** -- `kpmdmrg.general_kpm_moments`'s bare `raise` for a missing
  `X=`, surfacing as "RuntimeError: No active exception to reraise".

THE CONVENTION THESE TESTS PIN is stated once, in
`src/dmrgpy/dynamics.py`'s module docstring: every submode, on every
backend, returns the *complex* Lehmann density

    C_AB(w) = sum_n M_n * delta/(pi*((w-D_n)^2+delta^2)) = i(G^R-G^A)/(2pi)

with M_n = <GS|A|n><n|B|GS> and D_n = E_n-E_0. The durable statement of
it is the kernel-independent sum rule

    integral dw C_AB(w) = sum_n M_n = <GS|A B|GS>,

which no golden number can go stale against, and which the other common
convention -(1/pi) Im G^R cannot satisfy: its dispersive term
-Im(M_n)(w-D_n)/(pi*D) has principal-value tails reaching arbitrarily far
outside any finite window, and it integrates to Re<GS|A B|GS> plus a
window-dependent leftover (measured 0.13145 against the exact
0.11052176-0.27135291j on the chain used here).

Everything below is tested on an operator pair with a genuinely COMPLEX
M_n -- A = Cdag_0, B = C_2 on a chain with complex hoppings, i.e. an
ordinary off-diagonal Green's function. That is the point: for the
Hermitian pair A = B^dagger which every other test in tests/ uses, M_n is
real and the two conventions coincide exactly, so the whole family of
bugs is invisible. `test_hermitian_pair_*` below pins that coincidence
too, since a change *there* would mean something was broken rather than
fixed.
"""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain


# ---------------------------------------------------------------- helpers

def complex_hopping_chain(n=4, itensor_version=3, seed=3):
    """Spinless fermions with a random complex Hermitian hopping matrix
    plus a nearest-neighbour interaction. The complex hoppings are what
    make <GS|A|n><n|B|GS> complex for an off-diagonal pair; without them
    H is real, its eigenvectors can be chosen real, and every convention
    agrees."""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    return fc


def lehmann(chain, A, B):
    """(D_n, M_n): excitation energies and complex Lehmann weights
    <GS|A|n><n|B|GS>, from a full dense diagonalization done here with
    numpy -- deliberately no dmrgpy code in the reference itself."""
    ed = chain.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    emu, vs = np.linalg.eigh(h)
    U = np.array(vs)
    Uh = np.conjugate(U.T)
    Ae = Uh @ np.array(ed.MO2matrix(A).todense()) @ U
    Be = Uh @ np.array(ed.MO2matrix(B).todense()) @ U
    return emu - emu[0], Ae[0, :] * Be[:, 0]


def density(D, M, es, delta):
    """The house convention evaluated exactly from (D_n, M_n)."""
    return np.array([np.sum(M * (delta / np.pi) / ((w - D) ** 2 + delta ** 2))
                     for w in es])


ES = np.linspace(-1.0, 6.0, 40)   # pointwise comparison window
DELTA = 0.15

# Wide window / sharp broadening, for the sum rule: the integral only
# converges to <A B> once the Lorentzian tails are inside the window.
ES_WIDE = np.linspace(-12.0, 18.0, 1500)
DELTA_WIDE = 0.05


# ----------------------------------------------- #5: one convention, all
#                                                  submodes, both modes

@pytest.mark.parametrize("mode,submode", [
    ("ED", "ED"),      # was Re[C_AB]: the Im(M_n) term was dropped outright
    ("ED", "INV"),     # was already on the convention
    ("ED", "CVM"),     # was already on the convention
    ("ED", "ROOTN"),   # was -(1/pi) Im G^R
    ("DMRG", "CVM"),   # was -(1/pi) Im G^R
])
def test_every_submode_returns_the_complex_lehmann_density(mode, submode):
    """Pointwise against an exact numpy Lehmann sum on an off-diagonal
    pair. Before the fix the five non-KPM rows split into two families
    3e-1 and 6e-1 apart, larger than the correlator's own peak (0.426)."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    D, M = lehmann(fc, A, B)
    ref = density(D, M, ES, DELTA)
    _x, y = fc.get_dynamical_correlator(mode=mode, submode=submode,
                                        name=[A, B], es=ES, delta=DELTA)
    y = np.asarray(y, dtype=np.complex128)
    # The ED routes are exact here; mode="DMRG" submode="CVM" is limited
    # by its own conjugate-gradient tolerance (self.cvm_tol, 1e-5 by
    # default), which is what the looser bound tracks -- both are five
    # orders of magnitude below the 3e-1/6e-1 the two conventions differ
    # by on this pair. submode="KPM" is deliberately absent: it is on the
    # same convention but with a Chebyshev, not Lorentzian, kernel, so
    # only the sum rule below can pin it quantitatively.
    tol = 2e-5 if mode == "DMRG" else 1e-6
    assert np.max(np.abs(y - ref)) < tol, \
        "%s/%s is off the complex Lehmann density by %.3e" \
        % (mode, submode, np.max(np.abs(y - ref)))


@pytest.mark.parametrize("mode,submode", [
    ("ED", "ED"), ("ED", "INV"), ("ED", "CVM"), ("ED", "ROOTN"),
    ("ED", "KPM"), ("DMRG", "KPM"),
])
def test_sum_rule(mode, submode):
    """integral dw C_AB(w) = <GS|A B|GS>, complex. Kernel-independent, so
    this covers KPM (whose kernel is a Chebyshev expansion, not a
    Lorentzian) on exactly the same footing as the resolvent submodes,
    and it is the reason the complex density -- not -(1/pi) Im G^R -- is
    the house convention: only the density satisfies it."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    _D, M = lehmann(fc, A, B)
    exact = np.sum(M)   # = <GS|A B|GS>, completeness
    _x, y = fc.get_dynamical_correlator(mode=mode, submode=submode,
                                        name=[A, B], es=ES_WIDE,
                                        delta=DELTA_WIDE)
    got = np.trapezoid(np.asarray(y, dtype=np.complex128), ES_WIDE)
    # 2e-3 absolute: the Lorentzian submodes lose ~6e-4 to the finite
    # window, which is the only error left once the convention is right.
    assert abs(got - exact) < 2e-3, \
        "%s/%s integrates to %s, not <A B>=%s" % (mode, submode, got, exact)


def test_submode_ED_keeps_the_imaginary_part():
    """The sharpest single statement of #5: submode="ED", the exact
    Lehmann sum everything else is validated against, returned
    -out.imag/(2*pi) of a complex accumulator, i.e. Re(M_n) only. Its
    imaginary part was identically zero where the exact answer's peaks
    at 0.57."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    _x, y = fc.get_dynamical_correlator(mode="ED", submode="ED",
                                        name=[A, B], es=ES, delta=DELTA)
    assert np.max(np.abs(np.asarray(y).imag)) > 0.5


def test_ed_and_dmrg_agree_on_the_same_submode_name():
    """mode="ED" and mode="DMRG" must not mean two different observables
    under one submode name: max|CVM(ED)-CVM(DMRG)| was 0.5746 against a
    correlator peak of 0.4257."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    _x, y_ed = fc.get_dynamical_correlator(mode="ED", submode="CVM",
                                           name=[A, B], es=ES, delta=DELTA)
    _x, y_dmrg = fc.get_dynamical_correlator(mode="DMRG", submode="CVM",
                                             name=[A, B], es=ES, delta=DELTA)
    assert np.max(np.abs(np.asarray(y_ed) - np.asarray(y_dmrg))) < 2e-5


def test_finite_temperature_sum_keeps_the_imaginary_part():
    """dynamical_correlator_finite_T shares dynamical_sum's kernel and
    carried the identical `-out.imag` defect. At a temperature far below
    the gap (0.122 here) it must reproduce the T=0 answer exactly,
    imaginary part included."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    _x, y0 = fc.get_dynamical_correlator(mode="ED", submode="ED",
                                         name=[A, B], es=ES, delta=DELTA)
    _x, yT = fc.get_dynamical_correlator(mode="ED", submode="ED", T=1e-3,
                                         name=[A, B], es=ES, delta=DELTA)
    assert np.max(np.abs(np.asarray(yT).imag)) > 0.5
    assert np.max(np.abs(np.asarray(yT) - np.asarray(y0))) == pytest.approx(
        0.0, abs=1e-6)


def test_dmrg_rootn_is_on_the_convention_too():
    """The last route that was still returning -(1/pi) Im G^R:
    `rootndmrg.rootn_correction_vector` ran the fractional-resolvent
    recursion once, at +i*eta, where `algebra/rootn.py`'s ED
    implementation of the same submode name runs it at both signs and
    returns i(G^R-G^A)/(2*pi). Measured on this pair before the fix,
    mode="DMRG" submode="ROOTN" sat 5.720e-01 from the density (peak
    0.6135) and 4.8e-08 from the old convention, i.e. one submode name
    meant two different observables depending on `mode=`; it is now
    8.2e-07 from the density, which is the MPS/Lanczos error the
    mode="ED" route (7.9e-16) does not pay.

    It gets its own test rather than a row in
    test_every_submode_returns_the_complex_lehmann_density because it is
    the expensive one: N sequential Lanczos subspaces of dimension nkry,
    each step a truncated MPO application, now twice over -- hence the
    reduced N=6/nkry=16 here."""
    fc = complex_hopping_chain()
    A, B = fc.Cdag[0], fc.C[2]
    D, M = lehmann(fc, A, B)
    ref = density(D, M, ES, DELTA)
    _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="ROOTN",
                                        name=[A, B], es=ES, delta=DELTA,
                                        N=6, nkry=16)
    assert np.max(np.abs(np.asarray(y) - ref)) < 1e-5


# ----------------------------------- A = B^dagger: nothing may change
#                                     (and #11's factor of 2)

def staggered_heisenberg(n=6, itensor_version=3):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n): h = h + 0.3 * (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 20
    return sc


# Golden values recorded BEFORE the convention change, from the audit's
# own reproduction of finding #11 (6-site staggered Heisenberg, A=B=Sz_0,
# delta=0.3, es=linspace(0.2,4,8)). Every submode except CVM_explicit
# already produced these; CVM_explicit produced exactly twice them. They
# are pinned here precisely because A = B^dagger is the case that must
# NOT have moved: M_n = |<n|B|GS>|^2 is real and non-negative there, so
# the house convention and -(1/pi) Im G^R are the same number.
_HERMITIAN_GOLDEN = np.array([0.1084749, 0.0460046, 0.086561, 0.0525879])


@pytest.mark.parametrize("mode,submode", [
    ("ED", "ED"), ("ED", "CVM"), ("ED", "ROOTN"), ("ED", "EX"),
    ("DMRG", "CVM"), ("DMRG", "ROOTN"), ("DMRG", "CVM_explicit"),
])
def test_hermitian_pair_is_unchanged(mode, submode):
    """For A = B^dagger the fix must be a no-op -- a change here would
    mean something got broken, not fixed. CVM_explicit is the exception:
    it was returning 2x these values (#11), through an np.abs() that
    would also have flipped the sign of any negative weight."""
    sc = staggered_heisenberg()
    es = np.linspace(0.2, 4, 8)
    kw = {}
    if submode == "EX": kw["nex"] = 12
    _x, y = sc.get_dynamical_correlator(mode=mode, submode=submode,
                                        name=(sc.Sz[0], sc.Sz[0]),
                                        es=es, delta=0.3, **kw)
    y = np.asarray(y, dtype=np.complex128)
    # EX is a truncated excited-state expansion, the others are exact or
    # CG-converged on this tiny chain
    tol = 1e-3 if submode == "EX" else 1e-5
    assert np.max(np.abs(y[:4].real - _HERMITIAN_GOLDEN)) < tol, \
        "%s/%s moved off the pre-fix Hermitian-pair values: %s" \
        % (mode, submode, y[:4].real)
    assert np.max(np.abs(y.imag)) < 1e-5, \
        "a Hermitian pair has a real spectral density"


def test_cvm_explicit_matches_cvm():
    """The cross-check submode has to agree with what it cross-checks;
    the two were never compared on amplitude before, only on shape."""
    sc = staggered_heisenberg()
    es = np.linspace(0.2, 4, 8)
    name = (sc.Sz[0], sc.Sz[0])
    _x, y_cvm = sc.get_dynamical_correlator(submode="CVM", name=name,
                                            es=es, delta=0.3)
    _x, y_exp = sc.get_dynamical_correlator(submode="CVM_explicit",
                                            name=name, es=es, delta=0.3)
    assert np.max(np.abs(np.asarray(y_exp) - np.asarray(y_cvm))) < 1e-5


def test_cvm_explicit_names_its_own_restriction():
    """Its A^dagger == B guard was a bare `raise`: "RuntimeError: No
    active exception to reraise", after a bare print()."""
    fc = complex_hopping_chain()
    with pytest.raises(NotImplementedError) as exc:
        fc.get_dynamical_correlator(submode="CVM_explicit",
                                    name=[fc.Cdag[0], fc.C[2]],
                                    es=ES, delta=DELTA)
    assert "CVM_explicit" in str(exc.value)


# ------------------------------------------------- #29: the missing X=

def test_get_distribution_without_X_names_the_argument():
    """`if X is None: raise` in kpmdmrg.general_kpm_moments gave
    "RuntimeError: No active exception to reraise", naming neither the
    function nor the argument -- and get_distribution()'s own docstring
    documents no required argument at all."""
    sc = staggered_heisenberg(n=4)
    with pytest.raises(TypeError) as exc:
        sc.get_distribution(es=np.linspace(-1, 1, 5))
    assert "X=" in str(exc.value)


def test_kpm_moments_wfa_wfb_without_X_names_the_argument():
    """The same bare `raise` sat in the sibling entry point of the same
    file, reached with the two wavefunctions supplied directly."""
    from dmrgpy import kpmdmrg
    sc = staggered_heisenberg(n=4)
    wf = sc.get_gs()
    with pytest.raises(TypeError) as exc:
        kpmdmrg.kpm_moments_wfa_wfb(sc, wfa=wf, wfb=wf)
    assert "X=" in str(exc.value)
