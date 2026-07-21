"""Regression coverage for the 4-point fermionic correlator tensor
<Cdag_i C_j Cdag_k C_l> (mps.State.get_four_correlation_tensor(),
entropytk/correlationentropy.py), distilled from a manual ED-vs-DMRG
verification across itensor_version in (2, 3, "python") down to small,
fast (n=4) systems (examples/four_correlation_tensor/main.py is the
original, larger, non-assert exploratory version of this check).

Each model below was picked to have a comfortably large gap (>0.5)
between its ground and first excited state (checked once with ED at
authoring time). This matters: a manual sweep across random Hamiltonians
found that whenever the ground state is *nearly* degenerate (gap
~1e-4 or below), DMRG's tiny residual energy error (~1e-5, well within
its own convergence tolerance) can correspond to a several-percent
wavefunction admixture with the nearly-degenerate partner state --
energy is second-order insensitive to wavefunction error near a
variational minimum, but off-diagonal multi-point observables are not,
so the 4-point tensor can differ from the true (ED) one by O(0.1-0.5)
even though the ground-state energy looks converged. That is a real
numerical-accuracy property of near-degenerate spectra, not a bug in
the correlator code -- confirmed by all three DMRG backends still
agreeing with *each other* to ~1e-16 in that regime, i.e. they
consistently converge to the same (slightly wrong) state together. The
models here are deliberately kept well-gapped to give a clean,
non-flaky correctness signal instead.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, cppext

from _helpers import setup_backend

DMRG_TOL = 1e-4

VERSIONS = [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2),
        reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    pytest.param("python", id="python"),
]

N = 4  # chain length: small enough to be fast, large enough for a real tensor


def hopping_interaction_chain(version):
    """Complex NN hopping + NN density-density interaction. Well-gapped
    (gap ~1.22 at these parameters, checked with ED at authoring time)."""
    fc = fermionchain.Fermionic_Chain(N)
    setup_backend(fc, version)
    hops = [0.9 + 0.4j, 0.5 + 0.7j, 1.1 + 0.2j]
    h = 0
    for i in range(N - 1):
        h = h + hops[i] * fc.Cdag[i] * fc.C[i + 1]
        h = h + 0.7 * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    h = h + h.get_dagger()
    fc.maxm = 30
    fc.nsweeps = 10
    fc.set_hamiltonian(h)
    return fc, h


def free_fermion_chain(version):
    """Uniform NN hopping, no interaction. Well-gapped (gap ~0.62)."""
    fc = fermionchain.Fermionic_Chain(N)
    setup_backend(fc, version)
    h = 0
    for i in range(N - 1):
        h = h + 1.0 * fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    fc.maxm = 30
    fc.nsweeps = 10
    fc.set_hamiltonian(h)
    return fc, h


def staggered_field_chain(version):
    """Complex NN hopping plus a staggered on-site field breaking
    translational symmetry. Well-gapped (gap ~1.03)."""
    fc = fermionchain.Fermionic_Chain(N)
    setup_backend(fc, version)
    h = 0
    for i in range(N - 1):
        h = h + (1.0 + 0.3j) * fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    for i in range(N):
        h = h + ((-1) ** i) * 0.8 * fc.N[i]
    fc.maxm = 30
    fc.nsweeps = 10
    fc.set_hamiltonian(h)
    return fc, h


MODELS = {
    "hopping_interaction": hopping_interaction_chain,
    "free_fermions": free_fermion_chain,
    "staggered_field": staggered_field_chain,
}


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("model", MODELS.keys())
def test_four_correlation_tensor_matches_ed(model, version):
    build = MODELS[model]

    fc_ed, h = build(version)
    wf_ed = fc_ed.get_gs(mode="ED")
    ct_ed = wf_ed.get_four_correlation_tensor()

    fc_dmrg, _ = build(version)
    wf_dmrg = fc_dmrg.get_gs(mode="DMRG")
    ct_dmrg = wf_dmrg.get_four_correlation_tensor(ctmode="full", accelerate=True)

    assert np.max(np.abs(ct_dmrg - ct_ed)) == pytest.approx(0.0, abs=DMRG_TOL)


@pytest.mark.parametrize("version", VERSIONS)
def test_four_correlation_tensor_ctmode_and_accelerate_agree(version):
    """ctmode="full" (native AutoMPO, ITensor's own fermionic JW) and
    ctmode="explicit" (Python loop over vev() of MultiOperator products,
    multioperatortk's own JW) are two independent implementations of the
    same tensor -- they should agree with each other regardless of
    whether they also happen to agree with ED. Likewise accelerate=True
    (exploiting the <Cdag_i C_j Cdag_k C_l>^dagger = Cdag_l C_k Cdag_j C_i
    Hermitian symmetry to only compute half the tensor) must reproduce
    the non-accelerated result exactly."""
    fc, h = hopping_interaction_chain(version)
    wf = fc.get_gs(mode="DMRG")

    ct_full_accel = wf.get_four_correlation_tensor(ctmode="full", accelerate=True)
    ct_full_noaccel = wf.get_four_correlation_tensor(ctmode="full", accelerate=False)
    ct_explicit = wf.get_four_correlation_tensor(ctmode="explicit")

    assert np.max(np.abs(ct_full_accel - ct_full_noaccel)) == pytest.approx(0.0, abs=1e-8)
    assert np.max(np.abs(ct_full_accel - ct_explicit)) == pytest.approx(0.0, abs=1e-8)
