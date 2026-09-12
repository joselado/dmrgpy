"""Regression tests for the leftovers of the 2026-09 audit fix round.

Five items that the first round of fix agents could not reach because the
files were outside their lane:

  1. audit #12 -- `Parafermionic_Chain.get_dynamical_correlator` shadowed
     the base method and branched on the raw `mode=` argument, so an
     enforced `self.mode="ED"` was ignored and the call went to DMRG
     anyway (SIGABRT on itensor_version=3). The override is now deleted.
  2. `multioperator.obj2MO`'s bare `raise`, reachable from public API.
  3. `entropytk/correlationentropy.py`'s bare `raise` on a non-fermionic
     chain (already fixed by the first round; pinned here).
  4. `pychainwrapper.old2ampo`, a dead rebuild of the permanently-zero
     `self.exchange`/`self.fields`, and its one dead call site.
  5. `Many_Body_Chain.setup_cpp`/`setup_python`/`setup_julia` left the
     chain half-switched when `initialize()` raised.

Chains are deliberately tiny -- ED is exact at this size, and item 1's
pre-fix symptom on itensor_version=3 was a process abort, so the dispatch
tests run on `itensor_version="python"`, where the same mis-dispatch was a
catchable RuntimeError and the dispatch code being pinned is identical.
"""

import numpy as np
import pytest

from dmrgpy import spinchain, bosonchain, parafermionchain, cppext
from dmrgpy import multioperator, pychainwrapper


needs_v2 = pytest.mark.skipif(not cppext.available(2),
                              reason="the ITensor v2 extension is not compiled")
needs_v3 = pytest.mark.skipif(not cppext.available(3),
                              reason="the ITensor v3 extension is not compiled")


# ---------------------------------------------------------------- helpers

def heisenberg(n=4, itensor_version="python", maxm=30, nsweeps=20):
    """Uniform S=1/2 Heisenberg chain."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.SS(i, i + 1)
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


def parafermion(n=4, Z=3, itensor_version="python"):
    """Z_N parafermion chain with a hopping and a transverse field."""
    sc = parafermionchain.Parafermionic_Chain(n, Z=Z,
                                              itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sig[i] * sc.Sigd[i + 1]
    for i in range(n):
        h = h + 0.4 * sc.Tau[i]
    h = h + h.get_dagger()
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 64, 25
    return sc


# ------------------------------- #12 Parafermionic_Chain dynamical correlator

def test_parafermion_dynamical_correlator_honours_enforced_ed_mode():
    """audit #12: `self.mode="ED"` used to be ignored by this class.

    `Parafermionic_Chain` overrode `get_dynamical_correlator` with a body
    branching on the raw `mode=` argument, never consulting
    `self.get_mode()`. So a chain put into ED mode still dispatched to
    DMRG -- where the Hamiltonian had never been pushed into the session,
    which on itensor_version=3 reached ITensor's `Error()`/`abort()` and
    killed the whole process. The implicit call must now return exactly
    what the explicit `mode="ED"` call returns.
    """
    sc = parafermion()
    sc.mode = "ED"
    assert sc.get_mode() == "ED"
    es = np.linspace(-1.0, 1.0, 6)
    kw = dict(name=(sc.Tau[1], sc.Tau[1]), es=es, delta=0.3)
    x0, y0 = sc.get_dynamical_correlator(mode="ED", **kw)
    x1, y1 = sc.get_dynamical_correlator(**kw)  # no mode= : self.mode wins
    assert x1 == pytest.approx(x0, abs=1e-12)
    assert y1 == pytest.approx(y0, abs=1e-12)


def test_parafermion_dynamical_correlator_resolves_string_names():
    """audit #12, second facet: the override never resolved `name=`.

    The documented string form is resolved by
    `Many_Body_Chain.get_dynamical_correlator` through
    `operatornames.str2MO`; the override bypassed that, so a string went
    several frames deep and died inside `EDOperator` with "takes a
    MultiOperator or another EDOperator, got str". A parafermion chain has
    no Sx/Sy/Sz, so the right answer for name="ZZ" is a ValueError naming
    the operator -- not the opaque type error, and not a crash.
    """
    sc = parafermion()
    sc.mode = "ED"
    with pytest.raises(ValueError, match="not available on this chain"):
        sc.get_dynamical_correlator(name="ZZ", i=1, j=1,
                                    es=np.linspace(-1.0, 1.0, 4), delta=0.3)


def test_parafermion_class_no_longer_shadows_the_base_dispatcher():
    """The override contained nothing the base method does not do."""
    assert ("get_dynamical_correlator"
            not in parafermionchain.Parafermionic_Chain.__dict__)


def test_parafermion_commutation_failure_reports_what_failed():
    """`test_commutation`'s three bare `raise`s said "No active exception
    to reraise"; the diagnostic existed only in a `print`."""
    sc = parafermion(n=3)
    sc.test()  # the real relations hold, so this must pass
    bad = parafermion(n=3)
    # strip the Tau string off Chi: bare Sig operators on different sites
    # commute instead of satisfying the Z_N relation
    bad.Chi = [bad.Sig[i] for i in range(bad.ns)]
    bad.Chid = [o.get_dagger() for o in bad.Chi]
    np.random.seed(0)  # the site pairs it checks are drawn at random
    with pytest.raises(AssertionError, match="commutation test failed"):
        bad.test(ntries=200)


# -------------------------------------------------------- obj2MO bare raise

def test_obj2MO_reports_what_it_was_given():
    """`else: raise` -> "RuntimeError: No active exception to reraise".

    Reachable from public API: `gs_energy_generalized(A=...)`, `vev()`,
    `mpsalgebra` and `infinitechain` all funnel a user-supplied operator
    through `obj2MO`.
    """
    for bad in ["Sz", np.zeros((2, 2)), {"Sz": 0}]:
        with pytest.raises(TypeError, match="obj2MO"):
            multioperator.obj2MO(bad)


def test_obj2MO_still_accepts_its_three_documented_forms():
    mo = multioperator.obj2MO([["Sz", 0]])
    assert type(mo) == multioperator.MultiOperator
    assert multioperator.obj2MO(mo) is mo
    assert multioperator.obj2MO(2.0).op[0][0] == pytest.approx(2.0)


def test_multioperator_scalar_algebra_errors_name_the_operand():
    """The same bare-`raise` pattern at the other sites in that file."""
    mo = multioperator.obj2MO([["Sz", 0]])
    with pytest.raises(TypeError, match="divide"):
        mo / "a"
    with pytest.raises(TypeError, match="numpy array"):
        mo * np.zeros(3)
    with pytest.raises(TypeError, match="multiply_scalar"):
        mo.multiply_scalar("a")
    # the valid forms are untouched
    assert (2.0 * mo).op[0][0] == pytest.approx(2.0)
    assert (mo / 2).op[0][0] == pytest.approx(0.5)


# ----------------------------------------- correlation matrix on a spin chain

def test_correlation_matrix_without_operators_on_a_spin_chain():
    """`print("Unrecognized type",...)` + bare `raise` -> RuntimeError.

    A non-fermionic chain has no C/Cdag to build a default single-particle
    operator set from; say that, and say `operators=` can be passed.
    """
    sc = heisenberg()
    with pytest.raises(ValueError, match="operators="):
        sc.get_correlation_matrix()


# ----------------------------------------------------- pychainwrapper.old2ampo

def test_old2ampo_is_gone():
    """It rebuilt a Hamiltonian from `self.exchange`/`self.fields`, which
    `set_exchange()`/`set_fields()` (both removed) were the only things
    that ever populated -- so `for c in self.exchange` could only die with
    "TypeError: 'int' object is not iterable", and the field half
    referenced a bare undefined name `fields` on top of that."""
    assert not hasattr(pychainwrapper, "old2ampo")


def test_get_full_hamiltonian_still_works_and_says_when_there_is_no_H():
    sc = heisenberg()
    h = sc.get_full_hamiltonian()
    assert h.shape == (2 ** sc.ns, 2 ** sc.ns)
    e0 = np.linalg.eigvalsh(np.array(h.todense()))[0]
    assert e0 == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-8)
    empty = spinchain.Spin_Chain(["S=1/2"] * 4)  # no set_hamiltonian()
    with pytest.raises(RuntimeError, match="no Hamiltonian"):
        empty.get_full_hamiltonian()


# ------------------------------------------------- setup_cpp half-switch

@needs_v3
def test_failed_backend_switch_leaves_the_chain_on_its_old_backend():
    """`setup_cpp` assigned `self.itensor_version` before `initialize()`.

    `Bosonic_Chain.initialize()` refuses itensor_version=2 for any local
    dimension other than 4 (ITensor v2 knows only BosonFourSite and would
    abort the process). The chain was then left with
    `itensor_version == 2` and the *v3* session still in `_session`:
    anything keying on the session kept working, anything keying on the
    version -- `tevol_method`'s "TDVP only on 3", `__deepcopy__`'s
    `cppext.get_backend(self.itensor_version).Chain(...)` -- diverged from
    it.
    """
    bc = bosonchain.Bosonic_Chain(3, maxnb=[6, 6, 6], itensor_version=3)
    h = 0
    for i in range(2):
        h = h + bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i]
    for i in range(3):
        h = h + 0.3 * bc.N[i]
    bc.set_hamiltonian(h)
    e0 = bc.gs_energy(mode="DMRG")
    session = bc._session
    with pytest.raises(ValueError, match="itensor_version=2"):
        bc.setup_cpp(version=2)
    assert bc.itensor_version == 3          # not left on the version that failed
    assert bc._session is session           # and still holding its own session
    assert bc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)


@needs_v3
def test_failed_backend_switch_under_a_conserved_sector():
    """The mirror image: `sites.py::initialize` assigns `_session` first
    and only then applies the conserved sector, which refuses a backend
    with no quantum numbers -- so here it is the *session* that has
    already been overwritten when the raise happens."""
    sc = heisenberg(itensor_version=3)
    sc.set_conserved_sector(Sz=0)
    e0 = sc.gs_energy(mode="DMRG")
    session = sc._session
    with pytest.raises(Exception):
        sc.setup_cpp(version=2)
    assert sc.itensor_version == 3
    assert sc._session is session
    assert sc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)


@needs_v3
def test_successful_backend_switch_is_unaffected():
    """A switch that works must still work, on every backend present."""
    sc = heisenberg(itensor_version=3)
    e0 = sc.gs_energy(mode="ED")
    assert sc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)
    sc.setup_python()
    assert sc.itensor_version == "python"
    assert sc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)
    if cppext.available(2):
        sc.setup_cpp(version=2)
        assert sc.itensor_version == 2
        assert sc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)
    sc.setup_cpp(version=3)
    assert sc.itensor_version == 3
    assert sc.gs_energy(mode="DMRG") == pytest.approx(e0, abs=1e-6)
