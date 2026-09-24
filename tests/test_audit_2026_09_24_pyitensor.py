"""Regression tests for the pyitensor cluster of the 2026-09-24 audit.

Each test here locks in finding 15 or finding 16 of
`docs/audit_2026_09_24_hole_hunt.md`, which records the original symptom,
the reproduction that was executed and the reviewer's analysis.

- #15: above `_DENSE_EIG_MAX`, `excitation_energies(k, n>=2)` was one
  ARPACK call for all n values from one constant start. A single Krylov
  space holds one direction of a degenerate eigenspace, so the second copy
  of an exactly degenerate level came back as the next distinct level, a
  genuine eigenpair no residual check catches. The iterative path at n>=2
  is now one deflated Lanczos run per value, each from a fresh generic
  start (`_lowest_iterative_deflated`, the C++ port's vx_lanczos_lowest).
  The tests force it with `_DENSE_EIG_MAX=0` on the n_uc=2 Heisenberg
  cell at maxm=4, dim=48, where the old call misses in about 0.1 s.
- #16: under `set_pad_bonds(K)`, tevol_method="TDVP_GSE" ran one-site
  TDVP on the padded manifold: qr_split completed the padded zero
  directions into live ones, and at K=maxm the Krylov expansion had no
  room to add any. The one-site route is now exempt from padding the way
  the MPO is, so a padded run follows the unpadded one.

The pre-fix constructions are rebuilt in-test (the single ARPACK call; the
exemption bypassed by monkeypatching), so each test also shows that it
would have caught the defect. Everything runs on itensor_version="python".
"""

import contextlib

import numpy as np
import pytest

from dmrgpy import infinitechain, spinchain, timedependent
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import chain as chainmod
from dmrgpy.pyitensor import gse as gsemod
from dmrgpy.pyitensor import idmrg_excitations as ie
from dmrgpy.pyitensor.mpsalgebra import _link_at
from dmrgpy.pyitensor.svd import Spectrum
from dmrgpy.pyitensor.svd import svd as _real_svd


# ------------------------------------------------------ finding 15 helpers

@pytest.fixture(scope="module")
def heisenberg_cell():
    """The record's cell at maxm=4: critical S=1/2 Heisenberg on n_uc=2.
    Its H_eff(k) has an exactly degenerate lowest pair at k!=0 (the two
    transverse members of the magnon triplet of the finite-D state), with
    the longitudinal member split off just above it. VUMPS here is not
    reproducible run to run, which does not matter: every comparison below
    is between two solvers on the same environment."""
    np.random.seed(5)
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.gs_method = "vumps"
    ic.maxm = 4
    ic.vumps_nrestarts = 3
    ic.gs_energy()
    env = ic._get_excitation_environment()
    assert env.D * env.D * (env.d_g - 1) == 48
    return ic, env


def _dense(monkeypatch, env, k, n):
    monkeypatch.setattr(ie, "_DENSE_EIG_MAX", 10 ** 9)
    return np.asarray(ie.excitation_energies(env, k, n=n))


def _iterative(monkeypatch, env, k, n):
    monkeypatch.setattr(ie, "_DENSE_EIG_MAX", 0)
    return np.asarray(ie.excitation_energies(env, k, n=n))


def _single_arpack_call(env, k, n):
    """The pre-fix iterative path, verbatim: one eigsh call for all n
    values from the constant start vector."""
    from scipy.sparse.linalg import LinearOperator, eigsh
    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    dim = Dx * D
    op = LinearOperator((dim, dim), dtype=complex,
                        matvec=lambda x: ie._h_eff_action(k, x.reshape(Dx, D), env).reshape(-1))
    v0 = np.ones(dim, dtype=complex) / np.sqrt(dim)
    w, V = eigsh(op, k=n, which="SA", v0=v0, tol=ie._ITERATIVE_EIG_TOL)
    return np.sort(w) - env.lam_AC


# ----------------------------------------------------------- finding 15

@pytest.mark.parametrize("k", [0.37, 1.0])
def test_iterative_path_returns_the_degenerate_pair_with_multiplicity(heisenberg_cell, monkeypatch, k):
    _, env = heisenberg_cell
    dense = _dense(monkeypatch, env, k, 3)
    # the premise: an exactly degenerate lowest pair, and a third level
    # close above it but distinct, which is where one Krylov space misses
    assert abs(dense[1] - dense[0]) < 1e-8, dense
    assert dense[2] - dense[1] > 1e-6, dense
    for n in (2, 3):
        got = _iterative(monkeypatch, env, k, n)
        assert got == pytest.approx(dense[:n], abs=1e-10), (k, n, got, dense)


def test_single_arpack_call_misses_the_copy_on_this_cell(heisenberg_cell):
    """What the test above would have caught: the pre-fix call returns
    the next distinct level in place of the degenerate copy for at least
    one (k, n) here (at k=0.37, n=2 and k=1.0, n=3 in every run measured,
    off by 1.3e-4 and 0.10)."""
    _, env = heisenberg_cell
    misses = 0
    for k in (0.37, 1.0):
        dense = np.asarray(ie._lowest_dense(k, env, 3)[0]) - env.lam_AC
        for n in (2, 3):
            misses += np.max(np.abs(_single_arpack_call(env, k, n) - dense[:n])) > 1e-6
    assert misses >= 1


@pytest.mark.parametrize("k", [0.37, 1.0])
def test_n1_iterative_path_is_the_single_call_bit_for_bit(heisenberg_cell, monkeypatch, k):
    """n=1 keeps the one ARPACK call from the constant start, unchanged."""
    _, env = heisenberg_cell
    assert np.array_equal(_iterative(monkeypatch, env, k, 1), _single_arpack_call(env, k, 1))


def test_non_ascending_runs_are_refused_and_dense_answers(heisenberg_cell, monkeypatch):
    """A run landing below the previous one means an earlier run missed a
    level, so the whole iterative answer is refused rather than sorted into
    shape, and excitation_energies falls back to the dense path."""
    _, env = heisenberg_cell
    k = 0.37
    dense = _dense(monkeypatch, env, k, 3)
    real_run = ie._deflated_lanczos_run
    calls = []

    def run_that_undershoots(act, v0, found, niter, residual_tol):
        val, vec, top = real_run(act, v0, found, niter, residual_tol)
        calls.append(val)
        if found.shape[1] == 1:
            val = val - 0.05       # pretend the second run found a lower level
        return val, vec, top

    monkeypatch.setattr(ie, "_deflated_lanczos_run", run_that_undershoots)
    assert ie._lowest_iterative(k, env, 3) is None
    got = _iterative(monkeypatch, env, k, 3)
    assert got == pytest.approx(dense, abs=1e-12)
    assert len(calls) >= 2


def test_spectral_weights_multiplet_sums_match_dense(heisenberg_cell, monkeypatch):
    """The consumer: within the degenerate pair the split is basis-
    arbitrary, but its sum, the third branch and the returned fraction of
    the total are not, and they now agree with the dense path."""
    ic, _ = heisenberg_cell
    k = 0.37
    for op in ("Sx", "Sz"):
        monkeypatch.setattr(ie, "_DENSE_EIG_MAX", 10 ** 9)
        e_d, w_d, tot_d = ic.spectral_weights(op, k, p=0, n=3, return_total=True)
        monkeypatch.setattr(ie, "_DENSE_EIG_MAX", 0)
        e_i, w_i, tot_i = ic.spectral_weights(op, k, p=0, n=3, return_total=True)
        assert np.asarray(e_i) == pytest.approx(np.asarray(e_d), abs=1e-10)
        assert tot_i == pytest.approx(tot_d, abs=1e-12)
        assert w_i[0] + w_i[1] == pytest.approx(w_d[0] + w_d[1], abs=1e-8)
        assert w_i[2] == pytest.approx(w_d[2], abs=1e-8)


# ------------------------------------------------------ finding 16 helpers

N_SITES, K_PAD = 10, 4     # K = maxm, below the rank the evolved state needs
NT, DT = 20, 0.05
NT_LONG = 40               # where the sweeps=3 gap has grown past roundoff


@contextlib.contextmanager
def _no_suspension(suspend=True):
    yield


class _BackendWithoutExemption:
    """backend as chain.py sees it, with pad_bonds_suspended a no-op."""

    def __getattr__(self, name):
        if name == "pad_bonds_suspended":
            return _no_suspension
        return getattr(bk, name)


def _svd_reporting_the_padded_rank(*args, **kwargs):
    """gse.py's svd, with the spectrum stretched to the (padded) bond
    dimension, so `_gse_bond_step` counts the padded rows as rank the way
    it did before the fix."""
    U, S, V, spec = _real_svd(*args, **kwargs)
    width = S.inds[0].dim
    probs = np.concatenate([spec.eigs(), np.zeros(width - len(spec.eigs()))])
    return U, S, V, Spectrum(spec.singular_values, probs, spec.truncerr)


@pytest.fixture
def pre_fix(monkeypatch):
    """Rebuild the pre-fix one-site route: no suspension, no strip, and the
    padded bond dimension taken as the rank in the expansion."""
    def apply():
        monkeypatch.setattr(chainmod, "_bk", _BackendWithoutExemption())
        monkeypatch.setattr(chainmod, "_strip_bond_padding", lambda psi: psi)
        monkeypatch.setattr(gsemod, "svd", _svd_reporting_the_padded_rank)
    return apply


def _neel_chain(pad, n=N_SITES):
    """XXZ chain in a weak field, prepared in the Neel product state as the
    ground state of a staggered field; returns (chain, quench Hamiltonian).
    A product state is the case that separates the two routes: its true
    bond dimension is 1 everywhere, far below the pad width."""
    np.random.seed(11)
    bk.set_pad_bonds(pad)
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    sc.tevol_method = "TDVP_GSE"
    sc.maxm = K_PAD
    sc.nsweeps = 8
    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + 0.7 * sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h1 = h1 + 0.1 * sc.Sz[i]
    sc.set_hamiltonian(h0)
    sc.get_gs()
    return sc, h1


def _evolve(pad, sweeps, nt=NT):
    """<Sz_0>(t) after the Neel quench, through the public evolve_and_measure
    (the session's evolve_and_measure_tdvp_gse), plus the bond dimensions
    of the caller's starting wavefunction after the run."""
    try:
        sc, h1 = _neel_chain(pad)
        sc.tdvp_gse_sweeps = sweeps
        wf = sc.get_gs()
        sc.set_hamiltonian(h1)
        out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=DT,
                                               wf=wf, return_wf=True)
        wf_bonds = [_link_at(wf.cpp_handle, i, i + 1).dim for i in range(1, N_SITES)]
        return np.real(np.asarray(out[1])), wf_bonds
    finally:
        bk.set_pad_bonds(None)


def _quench(pad, sweeps, n=8):
    """C(t) = <A gs| exp(-iHt) A|gs> on the Neel state through the session's
    quench_tdvp_gse, called directly (the frequency stage of
    timedependent.py is not what is being pinned here)."""
    try:
        sc, h1 = _neel_chain(pad, n=n)
        A = sc.Sz[0] + sc.Sx[1]
        c, _ = sc._session.quench_tdvp_gse(h1.to_terms(), A.to_terms(), A.to_terms(),
                                           30, DT, sweeps, 3, 1e-8)
        return np.asarray(c)
    finally:
        bk.set_pad_bonds(None)


# ----------------------------------------------------------- finding 16

def test_padded_one_site_tdvp_follows_the_unpadded_trajectory(pre_fix):
    """tdvp_gse_sweeps=0, no expansion at all: one-site TDVP from a product
    state conserves bond dimension 1, so <Sz_0> stays at -1/2 unpadded.
    Padded, it used to leave it (by 0.49 over 40 steps), because qr_split
    had turned the padded zeros into live directions."""
    unpadded, _ = _evolve(None, 0)
    padded, wf_bonds = _evolve(K_PAD, 0)
    assert np.max(np.abs(padded - unpadded)) < 1e-12
    # the caller's own wavefunction keeps its padding: the strip ran on a copy
    assert wf_bonds == [K_PAD] * (N_SITES - 1)
    pre_fix()
    before, _ = _evolve(K_PAD, 0)
    assert np.max(np.abs(before - unpadded)) > 1e-2


def test_padded_tdvp_gse_follows_the_unpadded_trajectory(pre_fix):
    """The default tdvp_gse_sweeps=3, over 40 steps (over 20 both routes
    still agree to 6e-9). The two runs evolve the same state in different
    gauges, and the expansion's truncation (cutoff, and the cap at maxm)
    is gauge-sensitive at this level, so the agreement is to ~1e-7 rather
    than to roundoff: 1.7e-7, measured, against 4.1e-6 before the fix,
    both well inside TDVP's own 9e-5 error against ED on this quench."""
    unpadded, _ = _evolve(None, 3, nt=NT_LONG)
    padded, _ = _evolve(K_PAD, 3, nt=NT_LONG)
    assert np.max(np.abs(padded - unpadded)) < 1e-6
    pre_fix()
    before, _ = _evolve(K_PAD, 3, nt=NT_LONG)
    assert np.max(np.abs(before - unpadded)) > 1e-6


@pytest.mark.parametrize("sweeps", [0, 3])
def test_unpadded_runs_are_unchanged_by_the_exemption(pre_fix, sweeps):
    """With padding off the exemption and the true-rank count are no-ops,
    bit for bit."""
    fixed, _ = _evolve(None, sweeps)
    pre_fix()
    before, _ = _evolve(None, sweeps)
    assert np.array_equal(fixed, before)


def test_padded_quench_tdvp_gse_follows_the_unpadded_correlator(pre_fix):
    """The same exemption on quench_tdvp_gse, the route the TD dynamical
    correlator takes; before the fix this correlator moved by 0.33."""
    unpadded = _quench(None, 0)
    assert np.max(np.abs(_quench(K_PAD, 0) - unpadded)) < 1e-12
    pre_fix()
    assert np.max(np.abs(_quench(K_PAD, 0) - unpadded)) > 1e-2


def _expanded_bonds(pad, direct):
    """Bond dimensions after one Krylov expansion of Sz_0|Neel>, through
    Chain.global_subspace_expand, or with direct=True through gse.py
    itself with padding left on (no exemption), which is what the true-rank
    count in gse._gse_bond_step keeps right."""
    try:
        sc, h1 = _neel_chain(pad)
        s = sc._session
        H = sc.toMPO(h1).cpp_handle
        psi = s.apply_pure_operator(sc.toMPO(sc.Sz[0]).cpp_handle, sc.wf0.cpp_handle)
        if direct:
            out = gsemod.global_subspace_expand(H, psi, 3, 1e-8, bond_maxdim=K_PAD)
        else:
            out = s.global_subspace_expand(H, psi, 3, 1e-8, 0)
        return [_link_at(out, i, i + 1).dim for i in range(1, N_SITES)]
    finally:
        bk.set_pad_bonds(None)


def test_krylov_expansion_sees_the_true_bond_dimension_under_padding(pre_fix):
    """At K=maxm every padded bond already sat at the cap, so the expansion
    added nothing anywhere and returned [K]*(n-1), the padding itself. It
    now grows the bonds exactly as it does unpadded, both through the
    session (padding suspended) and when gse.py is called directly with
    padding on (the true-rank count alone)."""
    unpadded = _expanded_bonds(None, False)
    assert unpadded != [K_PAD] * (N_SITES - 1)
    assert max(unpadded) <= K_PAD
    assert _expanded_bonds(K_PAD, False) == unpadded
    assert _expanded_bonds(K_PAD, True) == unpadded
    pre_fix()
    assert _expanded_bonds(K_PAD, False) == [K_PAD] * (N_SITES - 1)
    assert _expanded_bonds(None, False) == unpadded
