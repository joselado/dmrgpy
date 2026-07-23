"""Regression coverage for the 4-point fermionic correlator tensor
<Cdag_i C_j Cdag_k C_l> (mps.State.get_four_correlation_tensor(),
entropytk/correlationentropy.py) on
fermionchain.Spinful_Fermionic_Chain_Native, mirroring
test_four_point_correlator.py's coverage of the spinless
Fermionic_Chain.

Spinful_Fermionic_Chain_Native exposes flat self.C/self.Cdag (mode
2*i = up, 2*i+1 = down at physical site i, see the class docstring in
fermionchain.py), so entropytk/correlationentropy.py's ctmode="explicit"
path -- the only one that works for native (Electron/Hubbard-type)
sites, see below -- needs no changes at all to cover this class; these
tests exercise that wiring plus cross-check the resulting tensor,
index for index, against the interleaved Spinful_Fermionic_Chain's own
(identical-convention) tensor.

ctmode="full" is NOT available for this class: it calls
Chain::four_correlation_tensor (mpscpp3/chain_session.h), which builds
its AutoMPO terms with the literal "Cdag"/"C" operator names -- not
defined on ITensor's ElectronSite (only Cup/Cdn/Cdagup/Cdagdn are), so
it would raise from ITensor's own site.op() dispatch. This mirrors why
plain ctmode="full" also isn't meaningful for a Cup/Cdn-only chain in
general -- there is no single unified fermion flavor to ask ITensor's
AutoMPO to auto-thread a plain "C"/"Cdag" Jordan-Wigner string through.

Each model was checked with ED at authoring time to have either a
comfortable gap or (for the spin-symmetric ones) an explicit symmetry-
breaking field, following test_spinful_fermion_chain.py's precedent:
a plain spin-up/down-symmetric spinful Hamiltonian has an *exactly*
2-fold degenerate ground state (verified: gap ~1e-16 without the field
term below), and an unconstrained DMRG search converging to either
member of that subspace makes the 4-point tensor (unlike the energy)
differ by O(0.1-0.5) from ED's arbitrary choice of the same subspace --
see test_four_point_correlator.py's module docstring for the general
phenomenon this reproduces.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain

DMRG_TOL = 1e-5

N = 3  # orbitals -> 2*N = 6 flat modes: small enough to be fast, large
       # enough for a real (non-trivial) 4-index tensor


def hopping_hubbard_chain():
    """Complex NN spin-conserving hopping + on-site Hubbard U, plus a
    small field lifting the up/down ground-state degeneracy (gap ~0.34
    at these parameters, checked with ED at authoring time)."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(N)
    h = 0
    for i in range(N - 1):
        h = h + (0.9 + 0.2j) * fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + (0.9 + 0.2j) * fc.Cdagdn[i] * fc.Cdn[i + 1]
    for i in range(N):
        h = h + 1.3 * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()
    h = h + 0.3 * (fc.Nup[0] - fc.Ndn[0])
    fc.maxm = 30
    fc.nsweeps = 12
    fc.set_hamiltonian(h)
    return fc


def free_fermion_chain():
    """Uniform NN spin-conserving hopping, no interaction, plus a
    degeneracy-lifting field (gap ~0.45 at N=3, checked with ED --
    needs a bigger field than the other two models here to reach a
    comfortable gap, since the plain free-fermion spectrum at this size
    is closer to degenerate to begin with)."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(N)
    h = 0
    for i in range(N - 1):
        h = h + 1.0 * fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + 1.0 * fc.Cdagdn[i] * fc.Cdn[i + 1]
    h = h + h.get_dagger()
    h = h + 1.0 * (fc.Nup[0] - fc.Ndn[0])
    fc.maxm = 30
    fc.nsweeps = 12
    fc.set_hamiltonian(h)
    return fc


def staggered_field_chain():
    """Complex NN spin-conserving hopping plus a staggered on-site
    charge field breaking translational symmetry (which, unlike the two
    models above, is enough on its own to leave a non-degenerate ground
    state -- checked with ED, gap ~1.1)."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(N)
    h = 0
    for i in range(N - 1):
        h = h + (1.0 + 0.3j) * fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + (1.0 + 0.3j) * fc.Cdagdn[i] * fc.Cdn[i + 1]
    h = h + h.get_dagger()
    for i in range(N):
        h = h + ((-1) ** i) * 0.8 * fc.Ntot[i]
    fc.maxm = 30
    fc.nsweeps = 12
    fc.set_hamiltonian(h)
    return fc


MODELS = {
    "hopping_hubbard": hopping_hubbard_chain,
    "free_fermions": free_fermion_chain,
    "staggered_field": staggered_field_chain,
}


@pytest.mark.parametrize("model", MODELS.keys())
def test_four_correlation_tensor_matches_ed(model):
    build = MODELS[model]

    fc_ed = build()
    wf_ed = fc_ed.get_gs(mode="ED")
    ct_ed = wf_ed.get_four_correlation_tensor()

    fc_dmrg = build()
    fc_dmrg.setup_cpp(3)
    wf_dmrg = fc_dmrg.get_gs(mode="DMRG")
    ct_dmrg = wf_dmrg.get_four_correlation_tensor(ctmode="explicit")

    assert np.max(np.abs(ct_dmrg - ct_ed)) == pytest.approx(0.0, abs=DMRG_TOL)


def test_four_correlation_tensor_matches_interleaved_chain():
    """Spinful_Fermionic_Chain_Native's self.C/self.Cdag use exactly the
    same flat (2*i=up, 2*i+1=down) convention as
    Spinful_Fermionic_Chain's own, so the two classes' 4-point tensors
    must agree index for index for the same physical Hamiltonian --
    cross-checked here against Spinful_Fermionic_Chain's own
    ctmode="full" (C++-accelerated, ITensor's native AutoMPO Jordan-
    Wigner) result, an independent implementation from the ctmode=
    "explicit" path both classes otherwise share.
    """
    U = 1.3

    def build(chain_cls):
        fc = chain_cls(N)
        h = 0
        for i in range(N - 1):
            h = h + (0.9 + 0.2j) * fc.Cdagup[i] * fc.Cup[i + 1]
            h = h + (0.9 + 0.2j) * fc.Cdagdn[i] * fc.Cdn[i + 1]
        for i in range(N):
            h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
        h = h + h.get_dagger()
        h = h + 0.3 * (fc.Nup[0] - fc.Ndn[0])
        fc.maxm = 30
        fc.nsweeps = 12
        fc.set_hamiltonian(h)
        fc.setup_cpp(3)
        return fc

    fc_native = build(fermionchain.Spinful_Fermionic_Chain_Native)
    ct_native = fc_native.get_gs(mode="DMRG").get_four_correlation_tensor(
            ctmode="explicit")

    fc_double = build(fermionchain.Spinful_Fermionic_Chain)
    ct_double = fc_double.get_gs(mode="DMRG").get_four_correlation_tensor(
            ctmode="full", accelerate=True)

    assert np.max(np.abs(ct_native - ct_double)) == pytest.approx(0.0, abs=DMRG_TOL)
