"""pyitensor's MPO builder: exactness, minimal bond dimension, and the two
bugs the rewrite fixed.

`mpobuilder.to_mpo` used to build one exact bond-dimension-1 MPO per term
and compress the block-diagonal concatenation of all of them. That reached
the right answer at the right final bond dimension, but only by way of an
intermediate whose bond dimension was the *term count*, which made the
build O(L^4): at L=100 it was 95% of a `gs_energy()` call (239 s, against
~13 s for the DMRG itself). It is now assembled directly as a finite-state
machine over the terms' partial products -- see mpobuilder's own docstring
and docs/pip_install_and_pyitensor_performance_plan.md.

Since that is a rewrite of the thing every single pyitensor operator is
built from -- the Hamiltonian, but equally every vev/correlator/time-
evolution vertex -- the check that matters is exactness against an
independent reference, on a zoo of term shapes rather than one model.
`AutoMPO.dense_matrix()` is that reference: it Kronecker-multiplies the
same per-site matrices without ever forming an MPO, so it shares no code
with either builder.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import contract_many

# site type codes, as manybodychain.py's callers use them (see siteset.py)
SPIN_HALF, SPIN_ONE, FERMION, ELECTRON = 2, 3, 0, 1


def mpo_dense(mpo, sites):
    """The MPO contracted out to a plain matrix, in the same (out,in),
    site-1-most-significant convention AutoMPO.dense_matrix() uses."""
    n = mpo.length()
    T = contract_many([mpo.A(i) for i in range(1, n + 1)])
    si = [sites.si(i) for i in range(1, n + 1)]
    arr = np.asarray(T.transpose_to([i.prime(1) for i in si] + si))
    dim = int(np.prod([i.dim for i in si]))
    return arr.reshape(dim, dim)


def bond_dims(mpo):
    return [mpo.A(i).inds[-1].dim for i in range(1, mpo.length())]


def heisenberg_terms(n):
    return [(1.0, [(op, i), (op, i + 1)])
            for i in range(1, n) for op in ("Sx", "Sy", "Sz")]


# Every term shape the builder has to get right: nearest-neighbour and
# long-range, gaps in the support, several factors on one site, complex
# coefficients, mixed local dimensions, a bare coefficient with no
# operators at all, and -- the ones with a Jordan-Wigner string threaded
# through them -- spinless and spinful fermions, including a term of odd
# fermion parity whose string runs to the end of the chain.
TERM_ZOO = [
    ("nn heisenberg", [SPIN_HALF] * 5, heisenberg_terms(5)),
    ("long range with a gap", [SPIN_HALF] * 5, [(0.7, [("Sz", 1), ("Sz", 5)])]),
    ("single site field", [SPIN_HALF] * 5, [(0.3, [("Sz", 3)])]),
    ("two factors on one site", [SPIN_HALF] * 5, [(1.0, [("S+", 2), ("S-", 2)])]),
    ("complex coefficient", [SPIN_HALF] * 5, [(0.5 + 0.25j, [("Sx", 1), ("Sy", 4)])]),
    ("two identical terms", [SPIN_HALF] * 5,
     [(1.0, [("Sz", 2), ("Sz", 3)]), (2.0, [("Sz", 2), ("Sz", 3)])]),
    ("bare identity term", [SPIN_HALF] * 5, [(1.5, [("Id", 1)])]),
    ("mixed local dimensions", [SPIN_HALF, SPIN_ONE, SPIN_HALF, SPIN_ONE],
     [(1.0, [("Sz", 1), ("Sz", 2)]), (1.0, [("S+", 2), ("S-", 3)]),
      (0.4, [("Sz", 4)])]),
    ("spinless fermion hopping", [FERMION] * 5,
     [(1.0, [("Cdag", i), ("C", i + 1)]) for i in range(1, 5)]
     + [(1.0, [("Cdag", i + 1), ("C", i)]) for i in range(1, 5)]),
    ("long range fermion hopping", [FERMION] * 5, [(1.0, [("Cdag", 1), ("C", 5)])]),
    ("odd fermion parity", [FERMION] * 5, [(1.0, [("C", 3)])]),
    ("fermion density interaction", [FERMION] * 5,
     [(1.3, [("N", i)]) for i in range(1, 6)]
     + [(0.9, [("N", i), ("N", i + 1)]) for i in range(1, 5)]),
    ("spinful fermion hubbard", [ELECTRON] * 4,
     [(1.0, [("Cdagup", i), ("Cup", i + 1)]) for i in range(1, 4)]
     + [(1.0, [("Cdagdn", i), ("Cdn", i + 1)]) for i in range(1, 4)]
     + [(2.0, [("Nupdn", i)]) for i in range(1, 5)]),
    ("one site chain", [SPIN_HALF], [(1.0, [("Sz", 1)]), (0.5, [("Sx", 1)])]),
    ("two site chain", [SPIN_HALF] * 2,
     [(1.0, [("Sz", 1), ("Sz", 2)]), (0.2, [("Sx", 1)])]),
]


@pytest.mark.parametrize("label,codes,terms",
                         TERM_ZOO, ids=[z[0] for z in TERM_ZOO])
def test_to_mpo_is_exact_against_the_dense_reference(label, codes, terms):
    """to_mpo reproduces the operator itself, not merely something with the
    right spectrum -- checked entry by entry against a Kronecker product
    that never builds an MPO."""
    sites = SiteX(codes)
    ampo = AutoMPO.from_terms(sites, terms)
    got = mpo_dense(mb.to_mpo(ampo, cutoff=1e-14), sites)
    assert got == pytest.approx(ampo.dense_matrix(), abs=1e-12)


@pytest.mark.parametrize("label,codes,terms",
                         TERM_ZOO, ids=[z[0] for z in TERM_ZOO])
def test_automaton_is_exact_before_any_compression(label, codes, terms):
    """The machine itself is exact: to_mpo's compression sweep exists to
    honour cutoff/maxdim and to squeeze out redundancy the partial-prefix
    sharing cannot see, never to *reach* the right operator. Checked
    separately so a bug in the machine can't hide behind the sweep."""
    sites = SiteX(codes)
    ampo = AutoMPO.from_terms(sites, terms)
    got = mpo_dense(mb._automaton_mpo(ampo), sites)
    assert got == pytest.approx(ampo.dense_matrix(), abs=1e-12)


@pytest.mark.parametrize("n", [6, 14, 30])
def test_nearest_neighbour_heisenberg_reaches_bond_dimension_five(n):
    """The textbook constant, and the whole point of the rewrite: it is now
    reached *directly* rather than by compressing an intermediate of bond
    dimension 3(L-1). The bulk value is exactly 5 and does not grow with n;
    the two end bonds are smaller because the machine has no I (or no F)
    channel to carry there."""
    sites = SiteX([SPIN_HALF] * n)
    ampo = AutoMPO.from_terms(sites, heisenberg_terms(n))
    assert max(bond_dims(mb._automaton_mpo(ampo))) == 5
    assert max(bond_dims(mb.to_mpo(ampo, cutoff=1e-14))) == 5


def test_mpo_bond_dimension_does_not_grow_with_chain_length():
    """The scaling claim, stated as a test rather than as a timing: a
    nearest-neighbour Hamiltonian's MPO is the same size at L=100 as at
    L=10. (The old builder also ended here -- its cost was the O(L)
    intermediate it passed through, which is what no longer exists.)"""
    dims = []
    for n in (10, 100):
        sites = SiteX([SPIN_HALF] * n)
        ampo = AutoMPO.from_terms(sites, heisenberg_terms(n))
        dims.append(max(bond_dims(mb._automaton_mpo(ampo))))
    assert dims[0] == dims[1]


def test_a_one_site_chain_keeps_every_term():
    """Regression, found while rewriting this: on a chain with a single
    site the old builder silently kept only the *last* term. `to_mpo` had
    no bonds to sweep there, so `sum_many`'s concatenation was returned
    as-is, and 0.8*Sz + 0.6*Sx came back as 0.6*Sx alone. Through the
    public API that turned a 1-site Hamiltonian into a different operator
    (and, being non-Hermitian by accident, sent it to NH-DMRG)."""
    sites = SiteX([SPIN_HALF])
    ampo = AutoMPO.from_terms(sites, [(0.8, [("Sz", 1)]), (0.6, [("Sx", 1)])])
    assert mpo_dense(mb.to_mpo(ampo, cutoff=1e-14), sites) \
        == pytest.approx(ampo.dense_matrix(), abs=1e-12)


def test_empty_operator_is_still_the_zero_mpo():
    """An AutoMPO with no terms at all is a legitimate input (see
    _zero_mpo's docstring: arnolditk.py's orthogonalization produces one
    whenever a coefficient filters to zero) and must not go near the
    machine."""
    sites = SiteX([SPIN_HALF] * 4)
    mpo = mb.to_mpo(AutoMPO.from_terms(sites, []), cutoff=1e-14)
    assert mpo_dense(mpo, sites) == pytest.approx(np.zeros((16, 16)), abs=0.0)


@pytest.mark.parametrize("n", [4, 8])
def test_ground_state_still_matches_ed(n):
    """End to end, through the public API: the operator the machine builds
    is the one DMRG then solves."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    sc.set_hamiltonian(sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                           + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
                       + 0.3 * sc.Sz[0])
    assert sc.gs_energy(mode="DMRG") == pytest.approx(sc.gs_energy(mode="ED"),
                                                      abs=1e-6)
