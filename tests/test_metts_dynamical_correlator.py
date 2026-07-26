"""Regression coverage for Many_Body_Chain.metts_dynamical_correlator()
(real-time finite-temperature dynamical correlators C_AB(t)=<A(t)B>_T via
the dynamical METTS algorithm, Wang, McClarty, Dankova, Honecker & Wietek,
"Spectroscopy and complex-time correlations using minimally entangled
typical thermal states", arXiv:2405.18484, Sec. II / Algorithm 1).

Validated directly against an exact ED reference: a full Boltzmann-weighted
Lehmann sum over eigenstates, computed inline here in the time domain
(the same physics as edtk/dynamics.py's dynamical_correlator_finite_T,
which returns the frequency-domain spectral function instead -- computing
the time-domain sum directly here avoids needing an FFT round-trip for
this cross-check).

METTS is a Monte Carlo method: exact agreement isn't expected, only
agreement within a generous multiple of the reported (Markov-correlated,
so likely optimistic) standard error -- same convention as
test_metts_vev.py.
"""
import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.edtk.edchain import EDOperator

from _helpers import setup_backend

VERSIONS = [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
]


def _heisenberg_field_chain(n, version, J=1.0, B=0.3):
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    setup_backend(sc, version)
    h = 0
    for i in range(n - 1):
        h = h + J * sc.Sx[i] * sc.Sx[i + 1] + J * sc.Sy[i] * sc.Sy[i + 1] + J * sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + B * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc, h


def _ed_time_domain_reference(sc, name, T, ts):
    """Exact C_AB(t)=<A(t)B>_T via a full Boltzmann-weighted Lehmann sum
    over every ED eigenstate -- A,B used exactly as given (no dagger),
    same convention edtk/dynamics.py's dynamical_correlator_finite_T
    uses (both built from the same Uh@op@U matrix-element machinery)."""
    edobj = sc.get_ED_obj()
    emu, vs = edobj.get_diagonalized_hamiltonian()
    a0 = EDOperator(name[0], edobj).SO
    b0 = EDOperator(name[1], edobj).SO
    e0 = np.min(emu)
    ex = emu - e0
    beta = 1.0 / T
    w = np.exp(-beta * ex)
    w = w / w.sum()
    U = np.array(vs)
    Uh = np.conjugate(U.T)
    A = Uh @ a0 @ U
    B = Uh @ b0 @ U
    out = np.zeros(len(ts), dtype=complex)
    for k, t in enumerate(ts):
        phase = np.exp(1j * (ex[:, None] - ex[None, :]) * t)
        out[k] = np.sum(w[:, None] * phase * A * B.T)
    return out


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_dynamical_correlator_matches_ed_lehmann_sum(version):
    """4-site Heisenberg+field chain: C_AB(t) for A=B=Sz[0] from dynamical
    METTS must agree with the exact ED Lehmann sum (within a generous
    multiple of METTS's own reported standard error) at every sampled
    time, not just at t=0."""
    n = 4
    sc, h = _heisenberg_field_chain(n, version, B=0.4)
    T = 1.0
    name = (sc.Sz[0], sc.Sz[0])
    nt, dt = 10, 0.2

    ts, means, stderrs = sc.metts_dynamical_correlator(
        name, T, nt=nt, dt=dt, nsamples=200, nwarmup=30,
        dbeta_half_step=0.08, seed=2024, njobs=4 if version == "python" else 1)

    ref = _ed_time_domain_reference(sc, name, T, ts)

    diff = np.abs(means - ref)
    tol = np.maximum(5 * stderrs, 0.03)  # 5-sigma headroom, floored for the near-zero-stderr t=0 point
    assert np.all(diff <= tol), list(zip(ts, diff, tol))


def test_metts_dynamical_correlator_t0_is_static_vev():
    """At t=0, v_i(0)=B|psi_i> and w_i(0)=|psi_i>, so C^i(0)=<psi_i|A B|psi_i>
    -- for A=B=Sz[0] on a spin-1/2 site this is exactly <(Sz[0])^2>=1/4
    regardless of temperature or sampling (Sz[0]^2 = Identity/4 for
    spin-1/2), a model-independent sanity check on the v/w bookkeeping
    at the very first measurement."""
    sc, h = _heisenberg_field_chain(3, "python", B=0.4)
    T = 0.8
    name = (sc.Sz[0], sc.Sz[0])
    ts, means, stderrs = sc.metts_dynamical_correlator(
        name, T, nt=1, dt=0.1, nsamples=20, nwarmup=10, seed=1)
    assert means[0].real == pytest.approx(0.25, abs=1e-8)
    assert means[0].imag == pytest.approx(0.0, abs=1e-8)


def test_metts_dynamical_correlator_requires_supported_backend():
    """Same backend restriction as metts_vev (see its own test) -- mpscpp2
    has no TDVP module at all."""
    sc, h = _heisenberg_field_chain(2, "python")
    sc.itensor_version = 2
    with pytest.raises(NotImplementedError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1)


def test_metts_dynamical_correlator_v3_rejects_short_chain():
    """Same ITensor v3 two-site-TDVP-aborts-below-3-sites guard as
    metts_vev's own test -- Python-side guard only, no compiled extension
    needed to exercise it."""
    sc, h = _heisenberg_field_chain(2, "python")
    sc.itensor_version = 3
    with pytest.raises(NotImplementedError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_dynamical_correlator_rejects_non_positive_temperature(version):
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), -1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1)
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 0.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_dynamical_correlator_rejects_empty_basis_ops(version):
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1, basis_ops=())


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_dynamical_correlator_rejects_non_positive_nsamples(version):
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=0, nwarmup=1)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_dynamical_correlator_rejects_non_positive_nt(version):
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=0, dt=0.1,
                                       nsamples=2, nwarmup=1)


def test_metts_dynamical_correlator_njobs_rejects_v3():
    """Same reason as metts_vev's own njobs>1 guard: mpscpp3's Chain is a
    single live in-process C++ session with no per-worker copy a
    multiprocessing pool could hand out."""
    sc, h = _heisenberg_field_chain(3, 3)
    with pytest.raises(NotImplementedError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=4, nwarmup=1, njobs=2)


def test_metts_dynamical_correlator_tdvp_cutoff_maxdim_reach_python_backend():
    """tdvp_cutoff/tdvp_maxdim must actually be reachable end-to-end from
    Many_Body_Chain.metts_dynamical_correlator down to
    pyitensor.metts.metts_dynamical_correlator (not silently dropped
    somewhere in the vevtk/mettsdynamicalcorrelator.py <-> pyitensor.chain
    <-> pyitensor.metts call chain, which used to leave them dead -- no
    caller could ever set them to anything but their None-default,
    confirmed by code review). Exercised via a value tight enough to
    force truncation (tdvp_maxdim=1, well below what a real evolution
    would need) so the call path is genuinely exercised, not just
    accepted and ignored; only checks the call succeeds and returns
    sane-shaped output, not the resulting numerical accuracy."""
    sc, h = _heisenberg_field_chain(3, "python", B=0.4)
    ts, means, stderrs = sc.metts_dynamical_correlator(
        (sc.Sz[0], sc.Sz[0]), 1.0, nt=3, dt=0.1, nsamples=3, nwarmup=1,
        seed=1, tdvp_cutoff=1e-6, tdvp_maxdim=1)
    assert len(ts) == len(means) == len(stderrs) == 3


def test_metts_dynamical_correlator_v3_rejects_tdvp_cutoff_maxdim():
    """v3's Chain::metts_dynamical_correlator has no separate real-time
    cutoff/maxdim knob at all -- passing either must raise rather than
    silently being ignored."""
    sc, h = _heisenberg_field_chain(3, 3)
    with pytest.raises(NotImplementedError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1, tdvp_cutoff=1e-6)
    with pytest.raises(NotImplementedError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=2, nwarmup=1, tdvp_maxdim=10)


def test_metts_dynamical_correlator_rejects_non_positive_njobs():
    sc, h = _heisenberg_field_chain(3, "python")
    with pytest.raises(ValueError):
        sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), 1.0, nt=2, dt=0.1,
                                       nsamples=4, nwarmup=1, njobs=0)
