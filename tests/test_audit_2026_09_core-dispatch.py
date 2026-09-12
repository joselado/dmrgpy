"""Regression tests for the core-dispatch findings of the 2026-09 audit.

Each test here locks in one finding from `docs/audit_2026_09_hole_hunt.md`
(findings #6, #14, #15, #16, #21, #23, #26, #27, #28), which records the
original symptom, the reproduction that was executed and the reviewer's
analysis. They share a shape: a dispatch decision taken without the
information that should inform it -- a `**kwargs` nobody consumes, a
module reading `self.mode` instead of asking `mode.py`, a type test that
went stale, an `else` doing two jobs -- see documentation.md 4.10.

Chains are deliberately tiny: ED is exact at this size and is the
reference wherever one is needed.
"""

import numpy as np
import pytest

from dmrgpy import spinchain, cppext


# ---------------------------------------------------------------- helpers

def heisenberg(n=4, itensor_version=3, maxm=30, nsweeps=20):
    """Uniform S=1/2 Heisenberg chain."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.SS(i, i + 1)
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


needs_v3 = pytest.mark.skipif(not cppext.available(3),
                              reason="the ITensor v3 extension is not compiled")


# ------------------------------------------------- #6 gs_energy_fluctuation

def test_gs_energy_fluctuation_forwards_mode():
    """`mode=` used to be swallowed by an unconsumed `**kwargs`.

    The body was `e=self.vev(h); e2=self.vev(h,npow=2)` with nothing
    forwarded, so `gs_energy_fluctuation(mode="ED")` returned the DMRG
    number byte for byte, while setting `sc.mode="ED"` on the same chain
    returned a different one -- proof the kwarg never reached the
    dispatcher.
    """
    sc = heisenberg(n=4)
    by_kwarg = sc.gs_energy_fluctuation(mode="ED")
    sc2 = heisenberg(n=4)
    sc2.mode = "ED"
    by_attribute = sc2.gs_energy_fluctuation()
    assert by_kwarg == pytest.approx(by_attribute, abs=1e-10)
    # and on an exact ED eigenstate the fluctuation is zero
    assert abs(by_kwarg) < 1e-6


@needs_v3
def test_gs_energy_fluctuation_on_the_automatic_ed_fallback():
    """A 2-site itensor_version=3 chain is routed to ED by mode.py with
    nobody opting in (v3's two-site dmrg() aborts below 3 sites).  The
    exact ground state of the 2-site Heisenberg model is an eigenstate at
    E0=-0.75, so the fluctuation must be ~0; it used to come back as
    1.1456439237389602 = sqrt(|E0-E0^2|), i.e. <H^2> answered with <H>.
    """
    sc = heisenberg(n=2, itensor_version=3)
    assert sc.get_mode() == "ED"
    assert sc.gs_energy_fluctuation() == pytest.approx(0.0, abs=1e-6)


def test_gs_energy_fluctuation_rejects_npow():
    """It sets the power itself; accepting npow= would silently collide."""
    sc = heisenberg(n=4)
    with pytest.raises(TypeError):
        sc.gs_energy_fluctuation(npow=3)


# --------------------------------------------------------- #14 exponential

def test_exponential_takes_the_dmrg_path_for_a_two_site_hamiltonian():
    """`mpsalgebra.exponential` gated on the *symbolic*
    MultiOperator.is_hermitian(), which false-rejects every two-site term
    (simplify() does not know that get_dagger()'s factor-order reversal
    is a no-op across sites).  Both branches failed and control fell into
    an uncontrolled 2-term Taylor truncation: measured 1.5%, 16% and 247%
    relative error at z = 0.25, 0.5 and 1.0 on this chain.

    The reference is ED's own exp(z*H) on the same state, which is the
    quantity the user guide documents (`e^{h}|psi>`); note the DMRG path
    additionally computed e^{-z*H} for real z, a sign flip invisible for
    as long as this branch was unreachable.
    """
    sc = heisenberg(n=4)
    h = sc.get_hamiltonian()
    assert not h.is_hermitian()       # the symbolic test still says False...
    assert sc.is_hermitian(h)         # ...while the numerical one is right
    wf = sc.get_gs()
    wf_ed = sc.get_gs(mode="ED")
    for z in [0.25, 0.5, 1.0]:
        ref = np.vdot(wf_ed.v, sc.exponential(z * h, wf_ed).v).real
        got = sc.overlap(wf, sc.exponential(z * h, wf)).real
        assert got == pytest.approx(ref, rel=1e-5)


def test_exponential_imaginary_step_is_unchanged():
    """The tau sign fix must not move the purely-imaginary-dt callers
    (timeevolution.evolve_WF), whose tau was already correct: the old
    formula complex(-dt.real,dt.imag) negated only the real part.
    exp(1j*t*H)|gs> = exp(1j*t*E0)|gs> for a ground state, so the overlap
    with |gs> is a pure phase of unit modulus and known argument.

    exponential_dmrg is called directly here, the way evolve_WF calls it
    -- going through exponential() instead would take the anti-Hermitian
    branch, which is a different code path (and is covered above)."""
    from dmrgpy.mpsalgebra import exponential_dmrg
    sc = heisenberg(n=4)
    h = sc.get_hamiltonian()
    e0 = sc.gs_energy()
    wf = sc.get_gs()
    t = 0.3
    out = exponential_dmrg(sc, h, wf, dt=1j * t, nt0=200)
    got = sc.overlap(wf, out)
    assert got == pytest.approx(np.exp(1j * t * e0), abs=1e-4)


def test_exponential_of_an_antihermitian_operator():
    """The anti-Hermitian branch passed dt=-1j alongside the operator
    -1j*h, i.e. exp((-1j)*(-1j)*h) = exp(-h): the same sign flip as the
    Hermitian branch, one composition further along."""
    sc = heisenberg(n=4)
    wf = sc.get_gs()
    wf_ed = sc.get_gs(mode="ED")
    K = 1j * 0.4 * (sc.Sz[0] * sc.Sz[1] + sc.Sz[1] * sc.Sz[2])  # K^dag = -K
    ref = np.vdot(wf_ed.v, sc.exponential(K, wf_ed).v)
    got = sc.overlap(wf, sc.exponential(K, wf))
    assert got == pytest.approx(ref, abs=1e-6)


def test_exponential_rejects_a_non_hermitian_operator_on_dmrg():
    """The old `else` printed "Warning, using 3rd order taylor expansion
    mode" and returned an unconverged number; there is no convergent DMRG
    route for this case, so it must raise."""
    sc = heisenberg(n=4)
    wf = sc.get_gs()
    A = sc.Sx[0] + 1j * sc.Sy[0]   # S+, neither Hermitian nor anti-Hermitian
    with pytest.raises(NotImplementedError):
        sc.exponential(A, wf)


# ------------------------------- #15 the automatic fallback and self.mode

@needs_v3
def test_algebra_primitives_follow_the_automatic_ed_fallback():
    """On a 2-site v3 chain mode.py answers with ED, but overlap/aMb read
    `self.mode` (still None) and took their DMRG branch anyway, dying with
    "'State' object has no attribute 'cpp_handle'" -- and random_state()
    handed back an MPS on the very same chain whose get_gs() returned a
    State, so no overlap between the two was possible at all.
    """
    sc = heisenberg(n=2, itensor_version=3)
    gs = sc.get_gs()
    assert sc.overlap(gs, gs) == pytest.approx(1.0, abs=1e-8)
    assert sc.aMb(gs, sc.get_hamiltonian(), gs) == pytest.approx(-0.75, abs=1e-8)
    # get_gs and random_state must agree about what kind of object this
    # chain holds
    assert type(sc.random_state()) is type(gs)


@needs_v3
def test_get_rdm_names_the_ed_routing_instead_of_an_attributeerror():
    """get_rdm has no ED implementation and consulted neither self.mode
    nor get_mode(), so every ED route reached session-only code.  It now
    refuses by name, the way get_distribution_moments does, and accepts
    mode= rather than rejecting it as an unexpected keyword."""
    sc = heisenberg(n=2, itensor_version=3)
    with pytest.raises(NotImplementedError):
        sc.get_rdm(i=0)
    with pytest.raises(NotImplementedError):
        sc.get_rdm(i=0, mode="ED")


# -------------------------------------- #16 the stale np.ndarray type test

def test_applyoperator_and_summps_work_on_the_ed_backend():
    """Both tested `type(wf)==np.ndarray` for their ED branch, but no ED
    route has produced a bare ndarray since EDchain.get_gs started
    returning an edtk.edchain.State: the branches were dead and every ED
    call fell into a bare `raise`, i.e. "RuntimeError: No active
    exception to reraise".  applyinverse, ten lines below in the same
    file, already tested State.
    """
    from dmrgpy.edtk.edchain import State
    sc = heisenberg(n=4)
    sc.mode = "ED"
    wf = sc.get_gs()
    assert isinstance(wf, State)
    A = sc.Sz[0]
    assert isinstance(sc.applyoperator(A, wf), State)
    assert isinstance(sc.summps(wf, wf), State)
    assert isinstance(sc.scale_mps(2.0, wf), State)
    # <gs|Sz0|gs> either way round, as a numerical check that the ED
    # branch computes the right thing and not merely the right type
    assert sc.overlap(wf, sc.applyoperator(A, wf)) == \
        pytest.approx(sc.vev(A), abs=1e-8)


def test_algebra_primitives_reject_an_unknown_wavefunction_type():
    """The bare `raise` reported "No active exception to reraise"."""
    sc = heisenberg(n=4)
    with pytest.raises(TypeError):
        sc.applyoperator(sc.Sz[0], "not a wavefunction")


# ------------------------------------------ #21 the documented mode= surface

def test_every_mps_algebra_primitive_accepts_mode():
    """The user guide's §2 table prefaces these with "Each takes the same
    mode=/**kwargs as the rest of the API".  Four of them raised TypeError
    on mode= (applyinverse -- which forwarded it into applyinverse_dmrg,
    whose only kwargs are delta/maxn -- scale_mps, operator_norm and
    is_zero_operator on top of it).
    """
    sc = heisenberg(n=4)
    wf = sc.get_gs()
    A = sc.Sz[0] * sc.Sz[1]
    assert sc.overlap(wf, wf, mode="DMRG") == pytest.approx(1.0, abs=1e-6)
    assert sc.aMb(wf, A, wf, mode="DMRG") == pytest.approx(sc.vev(A), abs=1e-6)
    sc.applyoperator(A, wf, mode="DMRG")
    sc.summps(wf, wf, mode="DMRG")
    sc.applyinverse(A + 2.0, wf, mode="DMRG")
    sc.scale_mps(2.0, wf, mode="DMRG")
    sc.operator_norm(A, mode="DMRG")
    assert sc.is_zero_operator(A - A, mode="DMRG")
    sc.trace(A, mode="DMRG")
    sc.exponential(0.1 * A, wf, mode="DMRG")


def test_mode_disagreeing_with_the_wavefunction_is_refused():
    """These primitives take the backend from the wavefunction's type, so
    a mode= naming the other one used to be dropped on the floor and the
    call silently ran the wrong backend."""
    sc = heisenberg(n=4)
    wf = sc.get_gs()               # an MPS
    with pytest.raises(TypeError):
        sc.applyoperator(sc.Sz[0], wf, mode="ED")


# --------------------------------------------------- #23 the dead exchange

def test_get_hamiltonian_without_set_hamiltonian_raises_by_name():
    """Spin_Chain.get_hamiltonian's "conventional way" fallback iterated
    self.exchange/self.fields, which no surviving builder populates --
    they are permanently the integer 0 -- so it died with "TypeError:
    'int' object is not iterable", and so did gs_energy_fluctuation(),
    which calls get_hamiltonian() unconditionally.
    """
    sc = spinchain.Spin_Chain(["S=1/2"] * 4)
    with pytest.raises(ValueError):
        sc.get_hamiltonian()
    with pytest.raises(ValueError):
        sc.gs_energy_fluctuation()


def test_get_hamiltonian_takes_self():
    """Many_Body_Chain.get_hamiltonian was `def get_hamiltonian():` --
    no self -- so every chain class that did not override it (all of them
    but Spin_Chain) raised TypeError on the plain call."""
    from dmrgpy import fermionchain
    fc = fermionchain.Fermionic_Chain(4)
    h = fc.Cdag[0] * fc.C[1] + fc.Cdag[1] * fc.C[0]
    fc.set_hamiltonian(h)
    assert fc.get_hamiltonian() is not None


# ----------------------------------------------- #26 the list-vs-ndarray n=1

@pytest.mark.parametrize("mode", ["DMRG", "ED"])
def test_get_excited_states_n1_returns_an_array(mode):
    """excited.py had two n==1 short circuits: the non-Hermitian one
    returned np.array([e0]), the Hermitian one a bare Python list.  Every
    other (n,mode) returns an ndarray, and get_gs_manifold's
    es[np.abs(es-e0)<tol] raises TypeError on a list.
    """
    sc = heisenberg(n=4)
    es, ws = sc.get_excited_states(n=1, mode=mode)
    assert isinstance(es, np.ndarray)
    assert len(sc.get_gs_manifold(n=1, mode=mode)) >= 1


# ------------------------------------------------ #27 the off-by-one bond guard

@pytest.mark.parametrize("itensor_version", [3, "python"])
def test_get_bond_entropy_rejects_an_out_of_range_site(itensor_version):
    """The guard read `b>self.ns` where bonds run 1..ns-1, so
    get_bond_entropy(wf, ns-1, ns) reached ITensor, whose own check calls
    abort(): the whole process died with an uncatchable SIGABRT and a core
    dump.

    NOTE: if that guard ever regresses on a C++ backend this test does not
    fail, it SIGABRTs the pytest session -- which is inherent to what it
    pins.
    """
    if itensor_version == 3 and not cppext.available(3):
        pytest.skip("the ITensor v3 extension is not compiled")
    sc = heisenberg(n=4, itensor_version=itensor_version, maxm=20, nsweeps=10)
    wf = sc.get_gs()
    # the valid bonds still answer
    assert sc.get_bond_entropy(wf, 1, 2) > 0.0
    with pytest.raises(IndexError):
        sc.get_bond_entropy(wf, sc.ns - 1, sc.ns)
    with pytest.raises(IndexError):
        sc.get_bond_entropy(wf, -1, 0)
    with pytest.raises(ValueError):
        sc.get_bond_entropy(wf, 0, 2)   # not adjacent


# ------------------------------------------------------ #28 unvalidated mode

def test_a_mistyped_chain_mode_raises_instead_of_returning_none():
    """resolve_mode returned self.mode unchecked (the `if mode in
    ["ED","DMRG"]` test only ever saw the *call argument*), so
    `sc.mode = "ed"` made get_gs() fall off the end of its if/elif and
    return None silently, while gs_energy() reported "RuntimeError: No
    active exception to reraise".
    """
    sc = heisenberg(n=4)
    sc.mode = "ed"
    with pytest.raises(ValueError):
        sc.get_gs()
    with pytest.raises(ValueError):
        sc.gs_energy()


def test_a_mistyped_mode_kwarg_names_the_valid_options():
    sc = heisenberg(n=4)
    with pytest.raises(ValueError):
        sc.gs_energy(mode="dmrg")


@needs_v3
def test_a_mistyped_chain_mode_is_caught_even_behind_a_fallback():
    """The validation has to happen *before* mode.py's own fallbacks can
    return "ED": on a 2-site v3 chain every one of them returns early, so
    checking where self.mode is read would let the typo through exactly
    where a wrong solver is hardest to notice."""
    sc = heisenberg(n=2, itensor_version=3)
    sc.mode = "Ed"
    with pytest.raises(ValueError):
        sc.get_mode()
