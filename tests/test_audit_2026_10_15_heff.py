"""The effective exchange of the Hubbard dimer must be its singlet-triplet
splitting (2026-10-15).

get_heff fitted the spin operators inside the span of the four lowest
eigenstates. The singlet carries doubly occupied weight, where S0.S1 = 0,
so that fit returned J = (E_T-E_S)/(1/4+3w/4), with w the singly occupied
weight of the singlet: 2.1811 instead of 1.5616 at U=t, 0.9306 instead of
0.8284 at U=4t, 0.19949 instead of 0.19804 at U=20t. With
projector=single_occupancy_projector(fc) the eigenstates are projected on
the singly occupied model space and Lowdin orthonormalized (des
Cloizeaux), and the fitted H_eff has exactly the eigenvalues E_S, E_T.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy import effectivehamiltonian as eh


def dimer(U, t=1.0, lam=20.0):
    """Hubbard dimer, particle-hole symmetric U term, half filling pinned
    by lam*(N-2)^2 so that the N=1 and N=3 states stay above the triplet
    at small U."""
    fc = fermionchain.Spinful_Fermionic_Chain(2)
    h = t*fc.Cdagup[0]*fc.Cup[1] + t*fc.Cdagdn[0]*fc.Cdn[1]
    h = h + h.get_dagger()
    for i in range(2):
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
    ntot = fc.Ntot[0] + fc.Ntot[1]
    h = h + lam*(ntot-2)*(ntot-2)
    fc.set_hamiltonian(h)
    fc.nsweeps = 40
    fc.maxm = 120
    ops = [fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1], fc.Sz[0]*fc.Sz[1]]
    return fc, ops


def splitting(U, t=1.0):
    return 0.5*(np.sqrt(U**2 + 16*t**2) - U)


@pytest.mark.parametrize("mode,tol", [("ED", 1e-10), ("DMRG", 1e-6)])
@pytest.mark.parametrize("U", [1.0, 4.0, 20.0])
def test_dimer_exchange_equals_singlet_triplet_splitting(U, mode, tol):
    fc, ops = dimer(U)
    j = np.real(fc.get_heff(operators=ops, n=4, mode=mode,
                            projector=eh.single_occupancy_projector(fc)))
    assert j == pytest.approx([splitting(U)]*3, abs=tol)


@pytest.mark.parametrize("U", [1.0, 20.0])
def test_without_projector_the_fit_is_the_projected_operator_fit(U):
    """projector=None keeps the fit inside the eigenstates, which is the
    old answer, (E_T-E_S)/(1/4+3w/4)."""
    fc, ops = dimer(U)
    j = np.real(fc.get_heff(operators=ops, n=4, mode="ED"))
    e_s = 0.5*(U - np.sqrt(U**2 + 16))
    w = 1.0/(1.0 + (e_s/2.0)**2)
    assert j == pytest.approx([splitting(U)/(0.25 + 0.75*w)]*3, abs=1e-10)


def test_couplings_with_projector_reproduce_the_splitting():
    """The pairwise-products fit takes the same keyword."""
    fc, _ = dimer(1.0)
    coef = eh.get_effective_hamiltonian_couplings(fc, mode="ED", n=4,
            projector=eh.single_occupancy_projector(fc))
    for pair in [("S^x_1", "S^x_2"), ("S^y_1", "S^y_2"), ("S^z_1", "S^z_2")]:
        assert np.real(coef[pair]) == pytest.approx(splitting(1.0), abs=1e-10)


def test_manifold_larger_than_the_model_space_raises():
    """Five states cannot be represented in the four singly occupied
    states, so the fit refuses rather than returning couplings."""
    fc, ops = dimer(4.0)
    with pytest.raises(ValueError, match="model space"):
        fc.get_heff(operators=ops, n=5, mode="ED",
                    projector=eh.single_occupancy_projector(fc))
