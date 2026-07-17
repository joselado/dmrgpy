"""Functional coverage for MPS.get_site_entropy/get_mutual_information,
distilled from examples/entanglement_entropy/entropy_dimer and
examples/mutual_information down to small, fast systems.

Unlike the other test files here, entanglement entropy/mutual
information have no ED implementation at all (MPS.get_site_entropy()
raises if self.MBO's wavefunction isn't a real MPS -- see mps.py): ED
ground states are plain State objects with no notion of a bond to take
an entropy across. These tests therefore need at least one compiled
C++ DMRG backend and are skipped otherwise; wherever both v2 and v3 are
compiled, they also double as a v2-vs-v3 cross-check.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy import cppext

pytestmark = pytest.mark.skipif(
    not (cppext.available(2) or cppext.available(3)),
    reason="entanglement entropy/mutual information need a compiled C++ DMRG backend",
)

DMRG_TOL = 1e-6


def _available_versions():
    return [v for v in (2, 3) if cppext.available(v)]


def test_bell_pair_site_entropy_is_log2():
    """Two-site Heisenberg dimer with no symmetry-breaking field: the
    exact ground state is the singlet (Sx0.Sx1+...), a maximally
    entangled Bell pair, so the entanglement entropy of either site with
    the rest of the system is exactly ln(2).

    v2 only: itensor_version=3 crashes hard ("LocalOp is default
    constructed", an ITensor v3 internal check in
    itensor/mps/localop.h) for any exactly-2-physical-site chain,
    independent of physics type -- a genuine mpscpp3 bug, confirmed
    also for spinless-fermion dimers (see test_fermion_chain.py). 3+
    sites is unaffected (see test_mutual_information_... below, and
    test_long_chain.py's 30-site chains)."""
    n = 2
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = sc.Sx[0] * sc.Sx[1] + sc.Sy[0] * sc.Sy[1] + sc.Sz[0] * sc.Sz[1]

    if not cppext.available(2):
        pytest.skip("only itensor_version=3 is compiled, which crashes on 2-site chains")

    sc.set_hamiltonian(h)
    sc.setup_cpp(2)
    wf = sc.get_gs(mode="DMRG")
    s0 = wf.get_site_entropy(0)
    assert s0 == pytest.approx(np.log(2), abs=DMRG_TOL)


def test_mutual_information_decreases_with_weaker_coupling():
    """4-site Heisenberg chain (examples/mutual_information): the
    mutual information between site 1 and site 2 must decrease as the
    coupling between them (c01) is turned down, since that coupling is
    exactly the bond straddling the (1,2) partition being measured."""
    n = 4
    spins = ["S=1/2" for _ in range(n)]

    def get_mi(version, c01):
        sc = spinchain.Spin_Chain(spins)
        h = 0
        for i in range(n - 1):
            c = 1.0 if i != 1 else c01
            h = h + c * sc.Sx[i] * sc.Sx[i + 1]
            h = h + c * sc.Sy[i] * sc.Sy[i + 1]
            h = h + c * sc.Sz[i] * sc.Sz[i + 1]
        sc.set_hamiltonian(h)
        sc.setup_cpp(version)
        wf = sc.get_gs(mode="DMRG")
        return wf.get_mutual_information(1, 2)

    for version in _available_versions():
        mi_weak = get_mi(version, 0.1)
        mi_strong = get_mi(version, 1.0)
        assert mi_weak < mi_strong
