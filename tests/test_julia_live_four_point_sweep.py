"""ctmode="sweep" four-point tensor <Cdag_i C_j Cdag_k C_l> on julia_live
(mpsjulialive/fourpoint.jl). The oracle is ctmode="explicit" on the same
Julia state, which builds one MPO per tuple with ITensors' own Jordan-Wigner
threading; maxm is above the exact bond dimension, so the oracle does not
truncate. ED checks that the state itself is converged."""
import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy.entropytk.correlationentropy import \
    _four_correlation_tensor_default_ctmode

from _helpers import julia_live_param, setup_backend

pytestmark = pytest.mark.parametrize("version", [julia_live_param()])

ED_TOL = 1e-8


def spinless_chain(version, n=6):
    """Complex hopping, density-density interaction, a field on one end:
    gapped, no particle-hole symmetry, complex amplitudes."""
    fc = fermionchain.Fermionic_Chain(n)
    setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + (0.9 + 0.4j) * fc.Cdag[i] * fc.C[i + 1]
        h = h + 0.7 * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    h = h + h.get_dagger() + 0.23 * fc.N[0]
    fc.maxm = 40
    fc.nsweeps = 20
    fc.set_hamiltonian(h)
    return fc


def hubbard_chain(version, n=3):
    """Interleaved spinful chain (2n spinless sites), Hubbard U and a
    Zeeman field, so up and down differ."""
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagdn[i] * fc.Cdn[i + 1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 2.0 * (fc.Nup[i] - 0.5) * (fc.Ndn[i] - 0.5)
        h = h + 0.3 * fc.Sz[i] + 0.1 * (i + 1) * fc.Ntot[i]
    fc.maxm = 40
    fc.nsweeps = 20
    fc.set_hamiltonian(h)
    return fc


@pytest.mark.parametrize("build", [spinless_chain, hubbard_chain],
                         ids=["spinless", "hubbard"])
def test_sweep_matches_explicit_and_ed(build, version):
    fc = build(version)
    wf = fc.get_gs()
    assert _four_correlation_tensor_default_ctmode(wf) == "sweep"
    sweep = wf.get_four_correlation_tensor(ctmode="sweep")
    explicit = wf.get_four_correlation_tensor(ctmode="explicit")
    assert np.abs(sweep - explicit).max() < 1e-8
    # repeated indices with an odd reordering sign, singled out
    for t in [(0, 0, 2, 1), (2, 1, 0, 0), (1, 3, 1, 0), (3, 3, 3, 3)]:
        assert abs(sweep[t] - explicit[t]) < 1e-8
    noacc = wf.get_four_correlation_tensor(ctmode="sweep", accelerate=False)
    assert np.abs(sweep - noacc).max() < 1e-12
    # the default ctmode is the sweep
    assert np.abs(wf.get_four_correlation_tensor() - sweep).max() < 1e-12
    ed = build(version).get_gs(mode="ED").get_four_correlation_tensor()
    assert np.abs(sweep - ed).max() < ED_TOL


def test_sweep_is_the_raw_expectation_value(version):
    """Not normalized: scaling the state by 2 scales the tensor by 4."""
    wf = spinless_chain(version, n=5).get_gs()
    ct = wf.get_four_correlation_tensor(ctmode="sweep")
    ct2 = (2.0 * wf).get_four_correlation_tensor(ctmode="sweep")
    assert np.abs(ct2 - 4.0 * ct).max() < 1e-10


def test_sweep_larger_chain_against_ed(version):
    fc = spinless_chain(version, n=8)
    ct = fc.get_gs().get_four_correlation_tensor(ctmode="sweep")
    ed = spinless_chain(version, n=8).get_gs(mode="ED")
    assert np.abs(ct - ed.get_four_correlation_tensor()).max() < ED_TOL


def test_sweep_matches_the_python_sweep(version):
    """Same model on the "python" backend: two different engines and two
    independent DMRG solves of a gapped state."""
    ct = spinless_chain(version).get_gs().get_four_correlation_tensor(
        ctmode="sweep")
    ref = spinless_chain("python").get_gs().get_four_correlation_tensor(
        ctmode="sweep")
    assert np.abs(ct - ref).max() < ED_TOL
