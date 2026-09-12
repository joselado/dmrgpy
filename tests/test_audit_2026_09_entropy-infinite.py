"""Regression tests for the 2026-09 audit's entropy/parity/infinite-chain
findings (#17, #22, #31, #35, #36).

`docs/audit_2026_09_hole_hunt.md` records the symptom, the reproduction and
the reviewer's analysis for each. They share a shape with the 2026-08 file
next door: a call either crashed where an equivalent one worked, or told
the user something that is not what the code does.

Chains are deliberately tiny -- ED is exact at this size and is the
reference wherever one is needed.
"""

import inspect

import numpy as np
import pytest

from dmrgpy import fermionchain, infinitechain
from dmrgpy.entropytk import correlationentropy as ce


# ---------------------------------------------------------------- helpers

def interacting_chain(n=4, itensor_version=3, mode=None):
    """Spinless hopping chain with a nearest-neighbour interaction, so the
    correlation matrix is not the trivial free-fermion one."""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    for i in range(n - 1):
        h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 20, 8
    if mode is not None: fc.mode = mode
    return fc


class _FakeWF:
    """The two attributes _default_dmmode actually reads off a state."""

    def __init__(self, MBO, cpp_handle=None):
        self.MBO = MBO
        self.cpp_handle = cpp_handle


# ----------------------------------------- #17: sector + ED, the dmmode default

def test_sector_correlation_matrix_works_on_the_ED_backend():
    """The 2026-08 audit's fix for its own finding #11 defaulted a
    sector-mode chain's dmmode to "full", which is
    wf.MBO._session.correlation_matrix(...) -- session-only. That was
    written when a sector-mode chain could not be answered by ED at all;
    d62a306 then gave ED its own sector implementation, and from there the
    default died with AttributeError: 'MBFermion' object has no attribute
    '_session' on every ED-answered sector chain. Enabling a sector broke
    a call that works without one."""
    fc = interacting_chain(mode="ED")
    ref = np.asarray(fc.get_correlation_matrix())  # no sector, works
    fc2 = interacting_chain(mode="ED")
    fc2.set_conserved_sector(Nf=2)
    got = np.asarray(fc2.get_correlation_matrix())  # default dmmode
    explicit = np.asarray(fc2.get_correlation_matrix(dmmode="explicit"))
    # the default must agree with the backend-agnostic route on the same
    # chain, which is the strongest statement that does not assume the
    # global ground state lies in Nf=2
    assert got == pytest.approx(explicit, abs=1e-8)
    # and everything layered on it must survive too (entanglement.py's
    # get_correlation_eigenvalues/_entropy went down with it)
    assert np.all(np.isfinite(fc2.get_correlation_eigenvalues()))
    assert np.isfinite(fc2.get_correlation_entropy())
    # this particular ground state does lie in Nf=2, so the sector and the
    # unconstrained answer coincide -- the oracle the audit used
    assert got == pytest.approx(ref, abs=1e-6)


@pytest.mark.parametrize("itensor_version", [3, "python"])
def test_sector_correlation_matrix_still_works_on_the_DMRG_backends(itensor_version):
    """The other half of the same default: a sector chain answered by DMRG
    keeps the session route ("full"), on both backends that have quantum
    numbers -- the fix must not have moved those off it."""
    from dmrgpy import cppext
    if not cppext.available(itensor_version):
        pytest.skip("itensor_version=%s not available" % (itensor_version,))
    fc = interacting_chain(itensor_version=itensor_version)
    fc.set_conserved_sector(Nf=2)
    wf = fc.get_gs()
    assert ce._default_dmmode(fc, wf, "electron") == "full"
    got = np.asarray(fc.get_correlation_matrix())
    ed = np.asarray(interacting_chain(itensor_version=itensor_version,
                                      mode="ED").get_correlation_matrix())
    assert got == pytest.approx(ed, abs=1e-6)


def test_sector_dmmode_default_is_resolved_from_the_state_not_the_backend():
    """The point of the fix: the default is chosen from what the
    wavefunction handed in can actually do, so it cannot go stale again
    when another backend becomes able to reach this code."""
    fc = interacting_chain()
    # no sector: the cheapest route, unchanged
    assert ce._default_dmmode(fc, _FakeWF(fc), "electron") == "fast"
    fc.set_conserved_sector(Nf=2)
    # a sector chain measured on an ED state (no handle, and an MBO with
    # no session at all) must not pick the session-only route
    assert ce._default_dmmode(fc, _FakeWF(object()), "electron") == "explicit"
    # a state with a live session that can answer directly -> "full"
    from dmrgpy import cppext
    if not cppext.available(3):
        pytest.skip("no compiled mpscpp3 extension, so no session to test")
    assert ce._default_dmmode(fc, _FakeWF(fc, cpp_handle=object()),
                              "electron") == "full"
    # the C++ correlation_matrix has no Nambu form, so a session is not
    # enough there either
    assert ce._default_dmmode(fc, _FakeWF(fc, cpp_handle=object()),
                              "Nambu") == "explicit"


# ------------------------------------- #22: a typo'd mode string names the options

def test_mistyped_mode_strings_name_the_valid_options():
    """Each of these ended its if/elif chain in a bare `raise` outside any
    except block, so a typo came back as "RuntimeError: No active
    exception to reraise" -- naming neither the argument nor its accepted
    values."""
    fc = interacting_chain(n=4)
    wf = fc.get_gs()
    with pytest.raises(ValueError, match="ctmode"):
        wf.get_four_correlation_tensor(ctmode="sweeep")
    with pytest.raises(ValueError, match="dmmode"):
        fc.get_correlation_matrix(dmmode="fasst")
    with pytest.raises(ValueError, match="fpmode"):
        wf.get_fermionic_parity(fpmode="fulll")
    # `basis` is the same documented enumeration one branch away, and was
    # worse than a bad message: a typo fell into the `else` and silently
    # returned the electron-basis matrix
    with pytest.raises(ValueError, match="basis"):
        fc.get_correlation_matrix(basis="nambu")
    # the valid values still work
    assert np.isfinite(complex(wf.get_fermionic_parity(fpmode="full")).real)


def test_a_mistyped_dmmode_is_rejected_before_the_ground_state_solve():
    """The validation sits ahead of get_gs(), so a misspelling does not
    cost a full DMRG solve before it is reported."""
    fc = interacting_chain(n=4)
    assert getattr(fc, "computed_gs", False) is False
    with pytest.raises(ValueError, match="dmmode"):
        fc.get_correlation_matrix(dmmode="fasst")
    assert getattr(fc, "computed_gs", False) is False


# --------------------------- #31: idmrg handles any reach, whatever the docstring said

@pytest.mark.parametrize("itensor_version", ["python"])
def test_idmrg_handles_a_coupling_past_one_unit_cell(itensor_version):
    """get_operator's docstring said gs_method="idmrg" is reach-1 only and
    raises past it. It is not and it does not -- the growth loop carries
    one pending channel per site of a term's reach, and the answer is
    exact. Only the tangent-space excitation ansatz is reach-1."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"],
                                           itensor_version=itensor_version)
    # H = -Sz + 0.5 Sz_i Sz_{i+1} + 0.25 Sz_i Sz_{i+2} on a 1-site cell:
    # the polarized state is exact at -1*0.5 + 0.5*0.25 + 0.25*0.25
    h = -ic.SzC[0] + 0.5 * ic.SzC[0] * ic.get_operator("Sz", 0, group=1) \
        + 0.25 * ic.SzC[0] * ic.get_operator("Sz", 0, group=2)
    ic.set_hamiltonian(h)
    ic.maxm = 4
    ic.gs_method = "idmrg"
    assert ic.gs_energy() == pytest.approx(-0.3125, abs=1e-6)
    assert complex(ic.vev("Sz", 0)).real == pytest.approx(0.5, abs=1e-6)
    # the one caller that genuinely is reach-1 still says so
    with pytest.raises(NotImplementedError):
        ic.excitation_energies(0.0)


def test_get_operator_docstring_matches_what_idmrg_does():
    doc = infinitechain.Infinite_Many_Body_Chain.get_operator.__doc__
    assert "BOTH ground-state methods" in doc
    # the sentence that was wrong: idmrg listed among what raises for reach>1
    assert "`gs_method=\"idmrg\"` and\n        `excitation_energies`" not in doc


# ------------------------------- #35: the four-point docstring names the real order

def test_four_correlation_tensor_docstring_lists_the_real_resolver_order():
    """The public docstring described a three-way choice ("sweep", then
    "full", then "explicit") that predates both "batched" (now the first
    thing the resolver tries) and "fold"."""
    doc = ce.get_four_correlation_tensor.__doc__
    for name in ("batched", "sweep", "fold", "full", "explicit"):
        assert '"%s"' % name in doc
    # and the order claimed is the order implemented
    order = [doc.index('"%s"' % n) for n in
             ("batched", "sweep", "fold", "full", "explicit")]
    assert order == sorted(order)
    fc = interacting_chain(itensor_version="python")
    assert ce._four_correlation_tensor_default_ctmode(_FakeWF(fc)) == "batched"


# ------------------------------------------- #36: the VUMPS restart default is 4

def test_vumps_nrestarts_default_is_four():
    """documentation.md reported its VUMPS timing table "at the default
    nrestarts=6"; the default is 4 in all three places and has never been
    6 (the measurement itself was run at 6, explicitly). Pinned here so
    the prose and the code cannot drift apart again unnoticed."""
    from dmrgpy.pyitensor import vumps, vumps_ms
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    assert ic.vumps_nrestarts == 4
    assert inspect.signature(
        vumps.vumps_ground_state).parameters["nrestarts"].default == 4
    assert inspect.signature(
        vumps_ms.ground_state).parameters["nrestarts"].default == 4
