"""Regression tests for the "ED backend and model classes" findings of the
2026-09 cross-backend audit (``docs/audit_2026_09_hole_hunt.md``).

Covered here: #4 (ED boson occupation projectors were ``|k><0|+|k><1|+...``
instead of ``|k><k|``), #6's ED half (``EDchain.vev`` swallowed ``npow=``),
#18 (a non-default ``maxnb`` on ``itensor_version=2`` aborted the process),
#19 (``Bosonic_Chain``/``Parafermionic_Chain`` never cached their ED
object), #25's code half (``SpinBoson_Chain`` ignored its own ``maxnb=``)
and #30 (these constructors rejected ``itensor_version=``).

The audit entry for each one records the original symptom, the
reproduction and the reviewer's analysis. Chains are deliberately tiny --
ED is exact at this size and is the reference wherever one is needed.
"""

import numpy as np
import pytest

from dmrgpy import bosonchain, parafermionchain, spinchain, spinfermionchain
from dmrgpy import cppext
from dmrgpy.pyboson import boson


# ---------------------------------------------------------------- helpers

def driven_boson_chain(n=3, maxnb=None, **kwargs):
    """Bosonic chain whose Hamiltonian does NOT conserve the boson number.

    The drive term ``0.5*(A_i+Adag_i)`` is what makes finding #4 visible:
    with a number-conserving Hamiltonian the ground state has a definite
    total N, every off-diagonal ``|n><m|`` piece of the broken projector
    has exactly zero expectation value, and ED agrees with DMRG to 1e-9
    even with the bug in place.
    """
    if maxnb is None: maxnb = [4] * n
    bc = bosonchain.Bosonic_Chain(n, maxnb=maxnb, **kwargs)
    h = 0
    for i in range(n - 1):
        h = h + bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i]
    for i in range(n):
        h = h + 0.3 * bc.N[i] * bc.N[i] + 0.5 * (bc.A[i] + bc.Adag[i])
    bc.set_hamiltonian(h)
    bc.maxm, bc.nsweeps = 30, 20
    return bc, h


def heisenberg(n=4, **kwargs):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, **kwargs)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 20
    return sc, h


def parafermion_chain(n=4, **kwargs):
    pf = parafermionchain.Parafermionic_Chain(n, **kwargs)
    h = 0
    for i in range(n - 1):
        h = h + pf.Chi[i] * pf.Chid[i + 1]
    h = h + h.get_dagger()
    pf.set_hamiltonian(h)
    pf.maxm, pf.nsweeps = 20, 10
    return pf, h


# ------------------------------------------------- #4: boson projectors

def test_ed_boson_occupation_operators_are_projectors():
    """`N<k>` on an ED boson site must be |k><k|: Hermitian, idempotent,
    and summing over k to the identity. It used to be sum_m |k><m| (the
    whole row k), because `op[n] = 1.0` on a (d,d) array assigns a row."""
    def dense(M):  # one2many hands back a sparse matrix
        return np.asarray(M.todense()) if hasattr(M, "todense") else np.asarray(M)
    b = boson.BosonChain([3, 3])
    dim = 9
    total = np.zeros((dim, dim), dtype=np.complex128)
    for k in range(3):
        M = dense(b.operators[("N" + str(k), 0)])
        assert np.allclose(M, M.conj().T), "N%d is not Hermitian" % k
        assert np.allclose(M @ M, M), "N%d is not a projector" % k
        total = total + M
    assert np.allclose(total, np.eye(dim))


def test_boson_occupation_probabilities_ed_matches_dmrg():
    """<D[i][k]> must be a probability on every backend: in [0,1],
    summing to 1 over k, with sum_k k*<D[i][k]> == <N_i>. Under ED it
    used to give negative entries summing to 0.058 on this chain."""
    bc, h = driven_boson_chain(n=3, maxnb=[4, 4, 4])
    e_ed = bc.gs_energy(mode="ED")
    e_dmrg = bc.gs_energy(mode="DMRG")
    assert e_ed == pytest.approx(e_dmrg, abs=1e-6)  # the Hamiltonian agrees
    for mode in ("ED", "DMRG"):
        ps = [bc.vev(bc.D[0][k], mode=mode).real for k in range(4)]
        n0 = bc.vev(bc.N[0], mode=mode).real
        for k, p in enumerate(ps):
            assert -1e-8 <= p <= 1 + 1e-8, "P(n=%d) = %g under %s" % (k, p, mode)
        assert sum(ps) == pytest.approx(1.0, abs=1e-6)
        assert sum(k * p for k, p in enumerate(ps)) == pytest.approx(n0, abs=1e-6)


# ------------------------------------------------------- #6: vev(npow=)

def test_ed_vev_honours_npow():
    """vev(h, npow=n, mode="ED") must return <H^n>. The ED ground state
    is an exact eigenstate, so <H^n> = E0^n identically -- a reference
    that needs no DMRG. npow used to fall into **kwargs and be dropped,
    so every power returned <H>."""
    sc, h = heisenberg(n=4)
    e0 = sc.gs_energy(mode="ED")
    assert sc.vev(h, npow=0, mode="ED") == pytest.approx(1.0, abs=1e-6)
    for n in (1, 2, 3):
        assert sc.vev(h, npow=n, mode="ED") == pytest.approx(e0 ** n, abs=1e-6)
        assert sc.vev(h, npow=n, mode="DMRG") == pytest.approx(e0 ** n, abs=1e-6)


def test_ed_vev_npow_on_an_operator_identity():
    """(Sz_i Sz_j)^2 = 1/16 identically for S=1/2, so <A^2> = 0.0625 in
    ANY state -- independent of convergence, and of which degenerate
    ground state the solver happens to land on."""
    sc, h = heisenberg(n=4)
    A = sc.Sz[0] * sc.Sz[2]
    for mode in ("ED", "DMRG"):
        assert sc.vev(A, npow=2, mode=mode) == pytest.approx(0.0625, abs=1e-6)


def test_gs_energy_fluctuation_is_zero_on_the_ed_route():
    """gs_energy_fluctuation() is sqrt(|<H^2>-<H>^2|); with npow dropped
    it became sqrt(|<H>-<H>^2|), which on this chain returned ~2.06 for
    a state that is an exact eigenvector of H."""
    sc, h = heisenberg(n=4)
    sc.mode = "ED"  # the attribute route -- the mode= kwarg is not consumed
    assert abs(sc.gs_energy_fluctuation()) < 1e-6


def test_gs_energy_fluctuation_on_the_automatic_ns_lt_3_fallback():
    """The route nobody opts into: mode.py falls back to ED for
    itensor_version=3 below 3 sites, so a 2-site chain on the *default*
    backend returned 1.1456 = sqrt(|E0-E0^2|) for an exact eigenstate."""
    sc, h = heisenberg(n=2, itensor_version=3)
    assert abs(sc.gs_energy_fluctuation()) < 1e-6


def test_ed_vev_rejects_a_negative_power():
    sc, h = heisenberg(n=4)
    with pytest.raises(ValueError):
        sc.vev(h, npow=-1, mode="ED")


# ----------------------------------------------- #18: maxnb on v2 aborts

def test_non_default_maxnb_on_itensor_version_2_raises():
    """mpscpp2/get_sites.h only knows the site-type code 104 (ITensor's
    BosonFourSite) and calls ITensor's Error() -- i.e. abort() -- for
    anything else, so this used to kill the interpreter with SIGABRT
    before any Python code could see it. Both routes onto the
    combination must raise instead."""
    with pytest.raises(ValueError):
        bosonchain.Bosonic_Chain(3, maxnb=[3, 3, 3], itensor_version=2)
    bc = bosonchain.Bosonic_Chain(3, maxnb=[3, 3, 3])  # default backend, fine
    with pytest.raises(ValueError):
        bc.setup_cpp(version=2)


@pytest.mark.skipif(not cppext.available(2), reason="mpscpp2 not compiled")
def test_maxnb_4_still_works_on_itensor_version_2():
    """The boundary is sharp: 4 levels is exactly what v2 does support,
    and must keep working (and agreeing with ED)."""
    bc, h = driven_boson_chain(n=3, maxnb=[4, 4, 4])
    e_ed = bc.gs_energy(mode="ED")
    bc.set_hamiltonian(h)
    bc.setup_cpp(version=2)
    assert bc.gs_energy(mode="DMRG") == pytest.approx(e_ed, abs=1e-6)


# ------------------------------------------------------- #19: ED caching

@pytest.mark.parametrize("build", [driven_boson_chain, parafermion_chain])
def test_ed_object_is_cached(build):
    """Every chain class must follow Many_Body_Chain's has_ED_obj/ED_obj
    protocol: without it each mode="ED" call rebuilt the sparse operator
    dictionary and redid the whole ground-state solve."""
    chain, h = build()
    assert chain.get_ED_obj() is chain.get_ED_obj()


@pytest.mark.parametrize("build", [driven_boson_chain, parafermion_chain])
def test_ed_cache_is_invalidated_by_a_new_hamiltonian(build):
    """The cache must not survive set_hamiltonian() -- that is the shape
    of the 2026-08 audit's own finding #2, and caching without checking
    it would simply have recreated that bug somewhere new. Checked
    empirically (a new object AND the right number), not by trusting
    restart()."""
    chain, h = build()
    obj = chain.get_ED_obj()
    e1 = chain.gs_energy(mode="ED")
    h2 = h + 0.7 * chain.N[0] * chain.N[0]
    chain.set_hamiltonian(h2)
    assert chain.get_ED_obj() is not obj
    e2 = chain.gs_energy(mode="ED")
    assert abs(e2 - e1) > 1e-6  # the new Hamiltonian really is different
    fresh, _ = build()
    fresh.set_hamiltonian(h2)
    assert fresh.gs_energy(mode="ED") == pytest.approx(e2, abs=1e-6)


# --------------------------------------------- #25: SpinBoson_Chain maxnb

def test_spinboson_chain_honours_maxnb():
    """maxnb= used to be accepted and dropped: the caller got 4-level
    boson sites (code 104) whatever they asked for, i.e. a silently
    wrong Hilbert space."""
    sb = bosonchain.SpinBoson_Chain(["B", "S=1/2"], maxnb=[6, None])
    assert sb.sites == [106, 2]
    assert sb.maxnb == [6, None]
    assert len(sb.D[0]) == 6 and len(sb.D[1]) == 0
    # "B<k>" says the same thing in the label
    assert bosonchain.SpinBoson_Chain(["B6", "S=1/2"]).sites == [106, 2]
    # and the default is unchanged
    assert bosonchain.SpinBoson_Chain(["B", "S=1/2"]).sites == [104, 2]


def test_spinboson_chain_rejects_a_bad_site_label():
    """The label the user guide used to give ("boson") died with
    `RuntimeError: No active exception to reraise` from a bare raise."""
    with pytest.raises(ValueError):
        bosonchain.SpinBoson_Chain(["boson", "S=1/2"])
    with pytest.raises(ValueError):  # n= is redundant, not silently dropped
        bosonchain.SpinBoson_Chain(["B", "S=1/2"], n=3)


def test_spinboson_chain_with_maxnb_agrees_across_dmrg_backends():
    """The wired maxnb has to actually build on both DMRG backends --
    SpinBoson_Chain's ED backend is an unfinished stub, so a
    v3-vs-"python" cross-check is the available reference."""
    def build(version):
        sb = bosonchain.SpinBoson_Chain(["B", "S=1/2", "S=1/2"],
                                        maxnb=[6, None, None],
                                        itensor_version=version)
        h = 0.5 * sb.Sz[1] + 0.5 * sb.Sz[2] + 0.3 * sb.Sx[1] * sb.Sx[2] \
            + 1.0 * sb.Adag[0] * sb.A[0] \
            + 0.4 * (sb.A[0] + sb.Adag[0]) * sb.Sx[1]
        sb.set_hamiltonian(h)
        sb.maxm, sb.nsweeps = 40, 20
        return sb
    e_python = build("python").gs_energy(mode="DMRG")
    assert build(3).gs_energy(mode="DMRG") == pytest.approx(e_python, abs=1e-6)


# ------------------------------------------- #30: itensor_version= kwarg

@pytest.mark.parametrize("build", [
    lambda v: bosonchain.Bosonic_Chain(3, maxnb=[6, 6, 6], itensor_version=v),
    lambda v: bosonchain.SpinBoson_Chain(["B", "S=1/2"], itensor_version=v),
    lambda v: parafermionchain.Parafermionic_Chain(4, itensor_version=v),
    lambda v: spinfermionchain.Spin_Fermion_Hamiltonian(["S=1/2"] * 2,
                                                        itensor_version=v),
])
def test_constructors_accept_itensor_version(build):
    """deb6bf8 gave the fermionic subclasses a **kwargs passthrough; these
    four were left with fixed signatures, so the kwarg the user guide
    tells boson users to pass was a TypeError."""
    assert build("python").itensor_version == "python"
