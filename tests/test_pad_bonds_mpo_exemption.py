"""`backend.set_pad_bonds` must pad MPS bonds and leave MPO bonds alone.

Padding exists to stop the engine minting a fresh array shape every time a
bond dimension changes, so that XLA compiles one kernel per operation
instead of one per bond dimension (see backend.set_pad_bonds). An MPS bond
really does change from sweep to sweep. An *operator's* does not: the
Hamiltonian MPO is built once by mpobuilder.to_mpo and keeps its bond
dimension for the whole run, so padding it cannot stabilize a shape that
was never moving -- it only inflates the MPO bond `w` that every
environment tensor carries and that the two-site matvec's dominant
O(chi^3 d^2 w) term is linear in.

It was not free. Before mpscontainer._Chain._pad_bonds existed, to_mpo's
own compression sweep went through the same svd() and was padded with
everything else: measured here, `set_pad_bonds(60)` took this Hamiltonian's
MPO from bond dimension 8 to 60. On a 6 GB GTX 1060 that turned a padded
ground state at maxm=60 into a single 1.77 GiB allocation the card could
not serve, XLA fell back to a slower plan, and the solve did not finish in
40 minutes against 6.4 s on one CPU core; at maxm=30, where it did fit, it
cost 1.84x (27.5 s against 15.0 s). On a large-memory device the same
padding is only an invisible constant factor, which is why it survived the
original port -- so this test pins the exemption on the host, where it
needs no GPU to check.

Runs on NumPy: the rule is a property of the engine, not of the array
library, and svd() applies padding identically on both.
"""

import pytest

from dmrgpy import spinchain
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor.tensor import commonIndex


PAD = 24


@pytest.fixture
def pad_restored():
    """Process-wide state, so put it back whatever the test does."""
    yield
    bk.set_pad_bonds(None)


def _solved_chain(pad, n=10, maxm=PAD, nsweeps=3):
    bk.set_pad_bonds(pad)
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):     # nearest neighbour
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n - 2):     # and next-nearest, so the MPO bond is > 5
        h = h + sc.Sx[i] * sc.Sx[i + 2] + sc.Sy[i] * sc.Sy[i + 2] \
            + sc.Sz[i] * sc.Sz[i + 2]
    sc.set_hamiltonian(h)
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    energy = sc.gs_energy()
    return sc, energy


def _mpo_dims(sc):
    H = sc._session.H
    return [max(ind.dim for ind in H.A(k).inds)
            for k in range(1, H.length() + 1)]


def _mps_bonds(sc):
    psi = sc.get_gs().cpp_handle
    return [commonIndex(psi.A(k), psi.A(k + 1)).dim
            for k in range(1, psi.length())]


def test_padding_does_not_touch_the_hamiltonian_mpo(pad_restored):
    unpadded, _ = _solved_chain(None)
    padded, _ = _solved_chain(PAD)
    assert _mpo_dims(padded) == _mpo_dims(unpadded)
    # and the exemption is only interesting because the MPO is smaller than
    # the pad width -- otherwise padding would have been a no-op anyway.
    assert max(_mpo_dims(unpadded)) < PAD


def test_padding_still_freezes_every_mps_bond(pad_restored):
    """The other half of the rule: what padding is *for* still happens."""
    unpadded, _ = _solved_chain(None)
    padded, _ = _solved_chain(PAD)
    assert set(_mps_bonds(padded)) == {PAD}
    # the unpadded chain is the interesting comparison: its bonds ramp up
    # from the chain edges, which is exactly the shape churn padding removes
    assert len(set(_mps_bonds(unpadded))) > 1


def test_padding_changes_no_energy(pad_restored):
    """Exact in representation: the appended singular values are zeros."""
    _, e_unpadded = _solved_chain(None)
    _, e_padded = _solved_chain(PAD)
    assert e_padded == pytest.approx(e_unpadded, abs=1e-6)
