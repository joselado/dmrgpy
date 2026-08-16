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


# -- Regression coverage for the local-fold repeated-index path
# (pyitensor/chain.py's _four_pt_fill_repeated). Those entries used to go
# through a per-tuple AutoMPO+to_mpo+inner build, which measured at 96% of
# four_correlation_tensor_sweep's total runtime at n=12 despite its comment
# calling it "subdominant". They are now folded locally, which needs a
# fermionic sign HTerm.resolve() does not supply -- see
# _four_pt_site_sort_sign. These tests pin both the sign and the scale.

def _fermion_chain_for_fold(n, maxm=20):
    from dmrgpy import fermionchain
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    for i in range(n - 1):
        h = h + 0.4 * fc.N[i] * fc.N[i + 1]
    h = h + 0.23 * fc.N[0]  # break particle-hole symmetry
    fc.set_hamiltonian(h)
    fc.maxm = maxm
    fc.gs_energy()
    return fc


def test_sweep_repeated_index_entries_match_the_slow_paths():
    """Every entry of the fast tensor must equal the always-correct 'full'
    and 'explicit' per-tuple builds. Those two are far too slow to use in
    anger but are exact oracles, and they cover the repeated-index entries
    the local fold now produces."""
    import numpy as np
    wf = _fermion_chain_for_fold(6).get_gs()
    fast = wf.get_four_correlation_tensor(ctmode="sweep")
    full = wf.get_four_correlation_tensor(ctmode="full")
    explicit = wf.get_four_correlation_tensor(ctmode="explicit")
    assert np.abs(fast - full).max() < 1e-12
    assert np.abs(fast - explicit).max() < 1e-12
    # and specifically the repeated-index entries, so this cannot pass by
    # the pairwise-distinct sweep alone being right
    n = 6
    rep = [(i, j, k, l)
           for i in range(n) for j in range(n) for k in range(n) for l in range(n)
           if len({i, j, k, l}) < 4]
    assert len(rep) == n**4 - n * (n - 1) * (n - 2) * (n - 3)
    worst = max(abs(fast[t] - full[t]) for t in rep)
    assert worst < 1e-12


def test_sweep_repeated_index_sign_is_not_dropped():
    """A tuple whose sites need an ODD reordering permutation to reach site
    order, e.g. (i,j,k,l)=(0,0,2,1). HTerm.resolve() composes same-site
    factors but drops the cross-site fermionic sign, so a fold that forgets
    _four_pt_site_sort_sign returns exactly the negative here -- checked
    against the independent 'full' oracle rather than a pinned number."""
    import numpy as np
    from dmrgpy.pyitensor.chain import _four_pt_site_sort_sign
    assert _four_pt_site_sort_sign((0, 0, 2, 1)) == -1
    assert _four_pt_site_sort_sign((0, 0, 1, 2)) == 1
    wf = _fermion_chain_for_fold(5).get_gs()
    fast = wf.get_four_correlation_tensor(ctmode="sweep")
    full = wf.get_four_correlation_tensor(ctmode="full")
    t = (0, 0, 2, 1)
    assert abs(fast[t] - full[t]) < 1e-12
    assert abs(fast[t]) > 1e-6          # non-trivial, so the sign is testable
    assert abs(fast[t] + full[t]) > 1e-6  # would vanish if the sign flipped


def test_sweep_is_raw_not_normalized():
    """Both halves of the tensor must be the raw <wf|Op|wf>. A wf scaled by
    c must scale every entry by |c|^2*... -- concretely by c^2 for real c,
    uniformly across distinct- and repeated-index entries alike. An earlier
    version of the sweep normalized one half only."""
    import numpy as np
    wf = _fermion_chain_for_fold(6).get_gs()
    base = wf.get_four_correlation_tensor(ctmode="sweep")
    scaled = (wf * 2.5).get_four_correlation_tensor(ctmode="sweep")
    assert np.abs(scaled - 6.25 * base).max() < 1e-10


def test_repeated_tuple_enumeration_is_exact():
    """_four_pt_repeated_tuples must yield each non-pairwise-distinct tuple
    exactly once -- it replaced an O(n^4) scan-and-skip, so an enumeration
    bug would silently leave entries at zero."""
    from dmrgpy.pyitensor.chain import _four_pt_repeated_tuples
    for n in (4, 6, 9):
        got = list(_four_pt_repeated_tuples(n))
        assert len(got) == len(set(got))
        assert len(got) == n**4 - n * (n - 1) * (n - 2) * (n - 3)
        assert all(len(set(t)) < 4 for t in got)


def test_four_pt_sign_agrees_with_hterm():
    """`_four_pt_site_sort_sign` is a second implementation of the
    anticommutation sign `autompo.HTerm.add()` already applies while
    insertion-sorting a term's factors. Pin them against each other over
    every site pattern a four-point tuple can have, so the simplified form
    (valid because all four factors are fermionic) cannot drift from the
    general one."""
    import itertools
    from dmrgpy.pyitensor.autompo import HTerm
    from dmrgpy.pyitensor.chain import _four_pt_site_sort_sign

    n = 4
    for i, j, k, l in itertools.product(range(n), repeat=4):
        ht = HTerm(1.0)
        for name, site in (("Cdag", i + 1), ("C", j + 1),
                           ("Cdag", k + 1), ("C", l + 1)):
            ht.add(name, site)
        # HTerm folds the sign into its coefficient while sorting
        assert ht.coef.real == _four_pt_site_sort_sign((i, j, k, l)), (i, j, k, l)
        assert abs(ht.coef.imag) < 1e-15
        # and the sorted factor order is what resolve() then consumes
        assert [s for _nm, s in ht.ops] == sorted([i + 1, j + 1, k + 1, l + 1])
