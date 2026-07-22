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

from dmrgpy import spinchain, timedependent

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
