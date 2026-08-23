"""Setting the Hamiltonian from an already-built MPO, not a term list.

`set_hamiltonian` normally takes a MultiOperator -- a symbolic sum of
products of named site operators, which the backend turns into an MPO via
ITensor's AutoMPO. It also accepts a `StaticOperator`, i.e. an MPO that has
already been built and possibly combined with the MPO algebra that class
exposes (`*` for products, `+` for compressed sums, scalar scaling).

Why that route exists: some operators are cheap as an MPO and expensive as
a term list. The motivating case is a total-spin penalty
`g*S^2_total = g*((sum Sx)^2 + (sum Sy)^2 + (sum Sz)^2)`, used to project
onto a spin sector in a ground-state-only solver. Symbolically that is
O(n^2) terms -- measured at 2030 terms for 12 spinful orbitals against 302
for the Hamiltonian itself -- while as MPOs it is three squares of
extensive one-body operators, each of small bond dimension.

The checks here are all "the two routes agree", because that is the whole
contract: an MPO-built Hamiltonian must give the same physics as the
symbolic one it is equal to. Nothing here asserts that the MPO route is
FASTER, and it should not -- AutoMPO compresses a term list well enough
that the term-count blow-up largely does not reach the solve (6.7x the
terms cost 1.4x the sweep time on the model above).
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import fermionchain
from dmrgpy import spinchain

pytestmark = pytest.mark.skipif(
    not cppext.available(3),
    reason="set_hamiltonian_mpo is implemented for itensor_version=3 only")


def _spin_chain(n=6):
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    sc.itensor_version = 3
    sc.maxm = 20
    sc.nsweeps = 12
    return sc


def _heisenberg(sc):
    h = 0
    n = sc.ns
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
            + sc.Sz[i] * sc.Sz[i + 1]
    return h


def test_mpo_hamiltonian_reproduces_the_symbolic_one():
    """The same Hamiltonian, set symbolically and as a prebuilt MPO, must
    give the same ground-state energy."""
    sc = _spin_chain()
    h = _heisenberg(sc)
    sc.set_hamiltonian(h)
    e_terms = sc.gs_energy()

    sc2 = _spin_chain()
    h2 = _heisenberg(sc2)
    sc2.set_hamiltonian(sc2.toMPO(h2))
    e_mpo = sc2.gs_energy()
    assert e_mpo == pytest.approx(e_terms, abs=1e-6)


def test_mpo_sum_matches_the_symbolic_sum():
    """H + g*Sz_total, built as an MPO SUM of two independently built MPOs,
    against the same thing written symbolically."""
    g = 0.3
    sc = _spin_chain()
    h = _heisenberg(sc)
    sz = 0
    for i in range(sc.ns):
        sz = sz + sc.Sz[i]
    sc.set_hamiltonian(h + g * sz)
    e_terms = sc.gs_energy()

    sc2 = _spin_chain()
    h2 = _heisenberg(sc2)
    sz2 = 0
    for i in range(sc2.ns):
        sz2 = sz2 + sc2.Sz[i]
    sc2.set_hamiltonian(sc2.toMPO(h2) + g * sc2.toMPO(sz2))
    e_mpo = sc2.gs_energy()
    assert e_mpo == pytest.approx(e_terms, abs=1e-6)


def test_total_spin_penalty_via_mpo_products():
    """The motivating case: S^2_total assembled as a sum of MPO SQUARES,
    which is what makes this route worth having. Checked against the
    symbolic S^2 rather than against a known eigenvalue, since the claim is
    equality of the two constructions."""
    g = 1.0
    sc = _spin_chain(n=4)
    h = _heisenberg(sc)
    sx = sy = sz = 0
    for i in range(sc.ns):
        sx = sx + sc.Sx[i]
        sy = sy + sc.Sy[i]
        sz = sz + sc.Sz[i]
    sc.set_hamiltonian(h + g * (sx * sx + sy * sy + sz * sz))
    e_terms = sc.gs_energy()

    sc2 = _spin_chain(n=4)
    h2 = _heisenberg(sc2)
    sx2 = sy2 = sz2 = 0
    for i in range(sc2.ns):
        sx2 = sx2 + sc2.Sx[i]
        sy2 = sy2 + sc2.Sy[i]
        sz2 = sz2 + sc2.Sz[i]
    S2 = (sc2.toMPO(sx2) * sc2.toMPO(sx2)
          + sc2.toMPO(sy2) * sc2.toMPO(sy2)
          + sc2.toMPO(sz2) * sc2.toMPO(sz2))
    sc2.set_hamiltonian(sc2.toMPO(h2) + g * S2)
    e_mpo = sc2.gs_energy()
    assert e_mpo == pytest.approx(e_terms, abs=1e-6)


def test_mpo_hamiltonian_matches_exact_diagonalization():
    """Against ED rather than against the other DMRG route, so a shared
    misconception between the two symbolic/MPO paths cannot pass."""
    sc = _spin_chain(n=6)
    h = _heisenberg(sc)
    sc.set_hamiltonian(h)
    e_ed = sc.gs_energy(mode="ED")

    sc2 = _spin_chain(n=6)
    h2 = _heisenberg(sc2)
    sc2.set_hamiltonian(sc2.toMPO(h2))
    assert sc2.gs_energy() == pytest.approx(e_ed, abs=1e-6)


def test_a_fermionic_chain_works_too():
    """Jordan-Wigner strings are threaded when the term list is built, so an
    MPO-set Hamiltonian must carry them identically."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.itensor_version = 3
    fc.maxm, fc.nsweeps = 20, 12
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    fc.set_hamiltonian(h)
    e_terms = fc.gs_energy()

    fc2 = fermionchain.Fermionic_Chain(n)
    fc2.itensor_version = 3
    fc2.maxm, fc2.nsweeps = 20, 12
    h2 = 0
    for i in range(n - 1):
        h2 = h2 + fc2.Cdag[i] * fc2.C[i + 1] + fc2.Cdag[i + 1] * fc2.C[i]
    fc2.set_hamiltonian(fc2.toMPO(h2))
    assert fc2.gs_energy() == pytest.approx(e_terms, abs=1e-6)
