"""Regression coverage for itensor_version="julia_live" (mpsjulialive/),
targeting the specific behaviors fixed across two rounds of code review on
the mpo-algebra-tdvp-gse branch: CVM's maxm restoration, the non-Hermitian
dispatch guard, and TDVP's per-step renormalization target. Before this
file, julia_live had zero persisted test coverage anywhere in tests/ --
every prior validation of these fixes was done with ad hoc, non-committed
scripts, so a future change could silently regress any of them with
nothing to catch it.

Requires a working juliacall/Julia toolchain (see
mpsjulialive/juliasession.py); the whole module is skipped if importing it
fails, mirroring how tests/ already skips itensor_version 2/3 when the
corresponding compiled extension isn't available (see cppext.available()
in test_nh_dmrg.py etc.)."""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain, timedependent

try:
    from dmrgpy.mpsjulialive import juliasession as _juliasession  # noqa: F401
    _JULIA_AVAILABLE = True
    _JULIA_UNAVAILABLE_REASON = ""
except Exception as _e:  # pragma: no cover - environment dependent
    _JULIA_AVAILABLE = False
    _JULIA_UNAVAILABLE_REASON = str(_e)

pytestmark = pytest.mark.skipif(not _JULIA_AVAILABLE,
        reason="requires a working juliacall/Julia toolchain: %s"
                % _JULIA_UNAVAILABLE_REASON)

DELTA = 0.05
HEISENBERG_4_GAP = 0.658919  # see test_dynamical_correlator.py


def _heisenberg_chain(n=4):
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    sc.setup_julia()
    sc.maxm = 20
    sc.nsweeps = 4
    return sc


def test_cvm_restores_maxm():
    """cvm.py::_cvm_sweep_params temporarily raises self.maxm to
    self.cvm_maxm for the CG solve and must restore it afterward, even on
    julia_live where this is a plain attribute (not a C++ session call) --
    the first code-review round found this wasn't restored at all."""
    sc = _heisenberg_chain()
    maxm0 = sc.maxm
    name = (sc.Sz[0], sc.Sz[0])
    sc.get_dynamical_correlator(mode="DMRG", submode="CVM", name=name,
            es=np.array([0.5]), delta=DELTA)
    assert sc.maxm == maxm0


def test_non_hermitian_blocks_kpm_but_not_excited_states():
    """The julia_live non-Hermitian guard in dynamics.py must only block
    KPM/CVM/TDZ (which assume a Hermitian spectrum), not EX -- EX's own
    non-Hermitian path (excited_states_non_hermitian, via
    mpsalgebra.mpsarnoldi) is itensor_version-agnostic and already works.
    An earlier fix placed this guard before submode dispatch, incorrectly
    blocking EX too (second code-review round)."""
    n = 4
    sc = _heisenberg_chain(n=n)
    h_nh = sc.hamiltonian + 1j * sum(sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h_nh)
    assert not sc.is_hermitian(sc.hamiltonian)

    name = (sc.Sz[0], sc.Sz[0])
    with pytest.raises(NotImplementedError):
        sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name,
                es=np.array([0.5]), delta=DELTA)

    # The underlying non-Hermitian excited-states machinery itself must
    # work on julia_live (this is what the dispatch fix was meant to
    # unblock).
    es, wfs = sc.get_excited_states(n=3, purify=False)
    assert len(es) == 3
    assert len(wfs) == 3

    # Going through the top-level dispatch, submode="EX" must not be
    # rejected by the julia_live-specific guard (NotImplementedError).
    # dcex.py has a separate, pre-existing kwarg-forwarding bug (its own
    # scale= default leaks into excited_states_non_hermitian/mpsarnoldi,
    # which don't accept it) that still breaks the full round trip on
    # every backend, not just julia_live -- unrelated to this guard and
    # out of scope here, so any *other* exception is accepted.
    try:
        sc.get_dynamical_correlator(mode="DMRG", submode="EX", name=name,
                es=np.array([0.5]), delta=DELTA, nex=3)
    except NotImplementedError:
        pytest.fail("submode='EX' must not be blocked by the julia_live "
                     "non-Hermitian guard (KPM/CVM/TDZ-only)")
    except Exception:
        pass


def test_evolution_aba_preserves_input_norm():
    """evolve_and_measure_tdvp (mpsjulialive/tdvp.jl) must renormalize
    every step back to the trajectory's own starting norm, not force it
    to 1 -- evolution_ABA() feeds it wfA = A*wf0, generally not
    unit-norm. Sz[0]^2 = (1/4)*Identity exactly for a spin-1/2 site, so
    ||Sz[0]*wf0||^2 == 0.25 exactly regardless of wf0, giving an exact,
    backend-independent expected value: with B left as the identity
    (evolution_ABA's default), the returned correlator <psi(t)|psi(t)>
    must stay at 0.25 throughout, not collapse to 1 (the bug the second
    code-review round found and this regression-tests)."""
    sc = _heisenberg_chain(n=3)
    wf0 = sc.get_gs()
    A = sc.Sz[0]
    wfA = A * wf0
    norm0 = wfA.dot(wfA).real
    assert norm0 == pytest.approx(0.25, abs=1e-8)

    ts, cs = timedependent.evolution_ABA(sc, A=A, mode="DMRG", wf=wf0,
            nt=5, dt=0.05)
    assert cs.real == pytest.approx(norm0 * np.ones(len(ts)), abs=1e-6)
    assert cs.imag == pytest.approx(np.zeros(len(ts)), abs=1e-6)


def test_kpm_peak_matches_exact_gap():
    """End-to-end sanity check that the KPM optimization (apply_op()
    instead of applyoperator() in the Chebyshev recursion and the seed
    vectors, plus summps()'s cutoff fix) didn't change the physical
    result: the dominant peak must still land on the exact 4-site
    Heisenberg gap."""
    sc = _heisenberg_chain()
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(0.3, 1.0, 71)
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
            name=name, es=es, delta=DELTA)
    x, y = np.array(x), np.array(y).real
    peak = x[np.argmax(y)]
    assert peak == pytest.approx(HEISENBERG_4_GAP, abs=0.05)


def _quench_from_product_state(tevol_method, gse_sweeps, n=6):
    """Spin-1/2 chain whose starting state is the exact *product* ground
    state of a pure staggered field (bond dimension 1), then quenched to
    Heisenberg. The product-state start is the point: one-site TDVP
    conserves bond dimension exactly, so without the global subspace
    expansion it structurally cannot follow this quench -- see
    examples/time_evolution/tdvp_gse_julia_VS_ED_time_evolution."""
    sc = spinchain.Spin_Chain([2 for _ in range(n)])
    sc.setup_julia()
    sc.tevol_method = tevol_method
    sc.tdvp_gse_sweeps = gse_sweeps
    sc.tdvp_gse_krylov_order = 3
    sc.tdvp_gse_cutoff = 1e-10
    h0, h1 = 0, 0
    for i in range(n):
        h0 = h0 + (-1)**i * sc.Sz[i]
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    wf_ed = sc.get_gs(mode="ED")
    sc.set_hamiltonian(h1)
    return sc, wf, wf_ed


def _max_diff_vs_ed(sc, wf, wf_ed, nt=40, dt=0.05):
    op = sc.Sz[0]
    _ts, sz = timedependent.evolve_and_measure(sc, operator=op, nt=nt,
            dt=dt, wf=wf)
    _ts_ed, sz_ed = timedependent.evolve_and_measure(sc, operator=op, nt=nt,
            dt=dt, wf=wf_ed, mode="ED")
    return np.max(np.abs(sz - sz_ed))


def test_tdvp_gse_matches_ed():
    """tevol_method="TDVP_GSE" on julia_live (mpsjulialive/tdvp.jl's
    evolve_and_measure_tdvp_gse: ITensorMPS.jl's own
    expand(...; alg="global_krylov") + nsite=1 tdvp) must track exact
    diagonalization, same as the plain two-site TDVP path."""
    sc, wf, wf_ed = _quench_from_product_state("TDVP_GSE", gse_sweeps=40)
    assert _max_diff_vs_ed(sc, wf, wf_ed) < 1e-4


def test_tdvp_gse_expansion_is_what_makes_it_work():
    """The companion negative control for the test above, and the only
    thing that proves the expansion actually ran rather than the call
    silently falling through to the two-site path: the *same* one-site
    integrator with tdvp_gse_sweeps=0 must fail badly on this quench,
    since one-site TDVP alone cannot grow a bond-dimension-1 state."""
    sc, wf, wf_ed = _quench_from_product_state("TDVP_GSE", gse_sweeps=0)
    assert _max_diff_vs_ed(sc, wf, wf_ed) > 1e-2


def test_unsupported_tevol_method_raises():
    """julia_live implements tevol_method "TDVP"/"TDVP_GSE" only. A
    request for "TEBD" (or the legacy "MPO" path) must raise rather than
    silently running plain TDVP instead -- a silent fallback would make a
    backend-comparison script compare the wrong integrators without
    saying so."""
    sc = _heisenberg_chain()
    sc.tevol_method = "TEBD"
    with pytest.raises(NotImplementedError):
        timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=3,
                dt=0.05, wf=sc.get_gs())


def test_nhdmrg_works_with_mode_ed_set():
    """chain.nhdmrg() must not care that self.mode=="ED" is set.

    The julia_live attempt functions build their own random start, and
    routing that through the chain-level Many_Body_Chain.random_state()
    made it honor self.mode and hand back an edtk State (no .jlmps),
    which blew up several frames later as an opaque AttributeError inside
    the juliacall call. The session backends build their random start
    inside their own session and are unaffected, so julia_live was the
    only backend where setting mode="ED" turned nhdmrg() into a crash
    rather than just running the DMRG solver (NH-DMRG has no ED
    implementation for mode= to select)."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_julia()
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + 0.6 * fc.Cdag[i + 1] * fc.C[i]
    for i in range(n):
        h = h + 1j * 0.7 * (-1)**i * fc.Cdag[i] * fc.C[i]
    fc.set_hamiltonian(h)
    fc.mode = "ED"
    e, psil, psir = fc.nhdmrg()   # must not raise AttributeError
    assert abs(psil.dot(psir) - 1.0) < 1e-6


def test_tdz_honors_tevol_method():
    """tdz.py's julia_live branch used to return before any
    tevol_method inspection, so submode="TDZ" silently ran plain two-site
    TDVP while the real-time entry points raised for the very same
    setting -- two calls on one chain disagreeing about which integrator
    "TEBD" or "TDVP_GSE" means. It now goes through the same
    _check_tevol_method guard."""
    from dmrgpy.tdz import _advance_complex_time_step
    sc = _heisenberg_chain()
    wf = sc.get_gs()
    Hop = sc.toMPO(sc.hamiltonian)
    sc.tevol_method = "TEBD"
    with pytest.raises(NotImplementedError):
        _advance_complex_time_step(sc, Hop, wf, 0.05 - 0.01j)
    # ...and the two implemented methods both work through the same path
    for method in ("TDVP", "TDVP_GSE"):
        sc.tevol_method = method
        out = _advance_complex_time_step(sc, Hop, wf, 0.05 - 0.01j, do_gse=True)
        assert out.dot(out).real > 0.0
