"""Canonical form for MultiOperator (multioperatortk/canonical.py).

Two things are pinned here. First, that the rewrite is exact: putting an
operator in canonical form must not change the matrix it builds on the ED
backend, for every statistics dmrgpy supports, including the fermionic
signs that reordering two factors across sites picks up. Second, that the
Hermiticity proof is one-sided in the direction it claims -- a True is a
proof, a False only means not proven -- and in particular that it proves
the Heisenberg Hamiltonian, which is the case the sympy round trip this
replaced got wrong (it reported False for Sx[i]Sx[j]+Sy[i]Sy[j]+Sz[i]Sz[j]
and took 0.834 s at n=20 to do it, which is why mpsalgebra.exponential and
infinitechain.set_hamiltonian both carried a comment about not using it).
"""

import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain, bosonchain, parafermionchain
from dmrgpy.multioperator import MO2matrix
from dmrgpy.multioperatortk import canonical


def _spin_operators():
    sc = spinchain.Spin_Chain(["1/2"]*4)
    ops = {
            "heisenberg bond": sc.SS(0,1),
            "reversed product": sc.Sx[1]*sc.Sx[0],
            "three sites out of order": sc.Sz[2]*sc.Sx[0]*sc.Sy[1],
            "same site twice": sc.Sx[0]*sc.Sx[0],
            "full chain": sum([sc.SS(i,i+1) for i in range(3)],0),
            }
    return sc,ops


def _spinless_fermion_operators():
    fc = fermionchain.Fermionic_Chain(4)
    fc.set_hamiltonian(sum([fc.Cdag[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.C[i]
                            for i in range(3)],0))
    ops = {
            "hopping out of order": fc.Cdag[2]*fc.C[0],
            "reversed hopping": fc.C[0]*fc.Cdag[2],
            "four fermions": fc.Cdag[3]*fc.Cdag[1]*fc.C[2]*fc.C[0],
            "density times hopping": fc.N[2]*fc.Cdag[0]*fc.C[1],
            "occupation on one site": fc.Cdag[0]*fc.C[0],
            "anticommutator": fc.C[0]*fc.C[2]+fc.C[2]*fc.C[0],
            }
    return fc,ops


def _spinful_fermion_operators():
    fs = fermionchain.Spinful_Fermionic_Chain(3)
    fs.set_hamiltonian(sum([fs.Cdagup[i]*fs.Cup[i+1]+fs.Cdagup[i+1]*fs.Cup[i]
                            for i in range(2)],0))
    ops = {
            "up hopping out of order": fs.Cdagup[1]*fs.Cup[0],
            "down hopping and density": fs.Cdagdn[2]*fs.Cdn[0]*fs.Nup[1],
            }
    return fs,ops


def _boson_operators():
    bc = bosonchain.Bosonic_Chain(3,maxnb=[3,3,3])
    bc.set_hamiltonian(sum([bc.Adag[i]*bc.A[i+1]+bc.Adag[i+1]*bc.A[i]
                            for i in range(2)],0))
    ops = {
            "hopping out of order": bc.Adag[2]*bc.A[0],
            "density times hopping": bc.N[1]*bc.Adag[0]*bc.A[2],
            "occupation projector": bc.D[1][1]*bc.Adag[0],
            }
    return bc,ops


_BUILDERS = {
        "spin": _spin_operators,
        "spinless fermion": _spinless_fermion_operators,
        "spinful fermion": _spinful_fermion_operators,
        "boson": _boson_operators,
        }


@pytest.mark.parametrize("statistics",sorted(_BUILDERS))
def test_canonical_form_is_the_same_operator(statistics):
    """The rewrite sorts factors by site and signs the fermionic
    exchanges, so the matrix it builds has to be unchanged."""
    chain,ops = _BUILDERS[statistics]()
    obj = chain.get_ED_obj()
    for (label,op) in ops.items():
        m0 = np.array(MO2matrix(op,obj))
        m1 = np.array(MO2matrix(op.simplify(),obj))
        assert np.max(np.abs(m0-m1))<1e-10, label


def test_fermionic_sign_of_a_reordering():
    """Two annihilation operators on different sites anticommute, so
    their anticommutator has to cancel in the canonical form, while the
    same pair on one site is left untouched and does not."""
    fc = fermionchain.Fermionic_Chain(4)
    assert (fc.C[0]*fc.C[2]+fc.C[2]*fc.C[0]).is_zero()
    assert not (fc.C[0]*fc.C[2]-fc.C[2]*fc.C[0]).is_zero()
    # bosonic factors commute instead, so it is the commutator that goes
    sc = spinchain.Spin_Chain(["1/2"]*4)
    assert (sc.Sx[0]*sc.Sz[2]-sc.Sz[2]*sc.Sx[0]).is_zero()


def test_heisenberg_hamiltonian_is_proven_hermitian():
    """The case the sympy round trip got wrong, on every chain length,
    and the fermionic and bosonic analogues of it."""
    for n in (4,10,20):
        sc = spinchain.Spin_Chain(["1/2"]*n)
        h = sum([sc.SS(i,i+1) for i in range(n-1)],0)
        assert h.is_hermitian()
        assert (h + 0.3*sum([sc.Sz[i] for i in range(n)],0)).is_hermitian()
    fc = fermionchain.Fermionic_Chain(4)
    assert (fc.Cdag[0]*fc.C[1]+fc.Cdag[1]*fc.C[0]).is_hermitian()
    assert (fc.N[0]*fc.N[1]).is_hermitian()
    bc = bosonchain.Bosonic_Chain(3,maxnb=[3,3,3])
    assert (bc.Adag[0]*bc.A[1]+bc.Adag[1]*bc.A[0]).is_hermitian()


def test_a_non_hermitian_operator_is_not_proven_hermitian():
    sc = spinchain.Spin_Chain(["1/2"]*4)
    assert not (1j*sc.Sx[0]).is_hermitian()
    assert (1j*sc.Sx[0]).is_antihermitian()
    fc = fermionchain.Fermionic_Chain(4)
    assert not (fc.Cdag[0]*fc.C[1]).is_hermitian()
    assert not (fc.Cdag[0]*fc.C[1]+1j*fc.Cdag[1]*fc.C[0]).is_hermitian()


def test_the_proof_is_one_sided_and_the_chain_falls_back():
    """Sx Sx Sy is (1/4)Sy on a spin-1/2 site, so it is Hermitian, but
    saying so needs the local Hilbert space and not just the names. The
    symbolic test is allowed to miss it; Many_Body_Chain.is_hermitian,
    which probes numerically when the proof does not land, is not."""
    sc = spinchain.Spin_Chain(["1/2"]*4)
    op = sc.Sx[0]*sc.Sx[0]*sc.Sy[0]
    assert not op.is_hermitian() # not proven, which is allowed
    assert sc.is_hermitian(op) # the numerical probe sees through it
    # and the chain-level check still rejects a genuinely non-Hermitian one
    assert not sc.is_hermitian(1j*sc.Sx[0])


def test_an_unknown_name_is_never_reordered_or_proven():
    """Parafermionic operators reorder with a Z_n phase rather than a
    sign, so the canonical form leaves them exactly as written, and a
    term built out of them is never proven Hermitian."""
    pf = parafermionchain.Parafermionic_Chain(4,Z=3)
    h = pf.Sig[0]*pf.Sigd[1] + pf.Sig[1]*pf.Sigd[0]
    assert not h.is_hermitian()
    assert canonical.parity("Sig") is None
    # left spelled the way it was built, no reordering behind the caller
    terms = [tuple((o[0],o[1]) for o in t[1:]) for t in h.simplify().op]
    assert (("Sig",0),("SigDag",1)) in terms
    assert (("Sig",1),("SigDag",0)) in terms


def test_identities_and_zero_terms_are_dropped():
    """The "h = 0; h = h + term" idiom leaves a 0*identity placeholder,
    and an explicit identity factor multiplies nothing."""
    sc = spinchain.Spin_Chain(["1/2"]*4)
    h = sum([sc.SS(i,i+1) for i in range(3)],0)
    assert len(h.simplify().op)==len(h.op)-1 # the placeholder is gone
    idop = sc.get_operator("Id",0)
    assert len((idop*sc.Sx[1]).simplify().op[0])==2 # coefficient plus Sx
    assert (sc.Sx[0]-sc.Sx[0]).is_zero()


def test_equal_terms_are_collected():
    sc = spinchain.Spin_Chain(["1/2"]*4)
    h = sc.Sx[0]*sc.Sx[1] + sc.Sx[1]*sc.Sx[0] + 2.0*sc.Sx[0]*sc.Sx[1]
    out = h.simplify()
    assert len(out.op)==1
    assert abs(out.op[0][0]-4.0)<1e-12


@pytest.mark.parametrize("itensor_version",[2,3,"python"])
def test_ground_state_energy_is_unchanged(itensor_version):
    """The canonical form sits on the path every gs_energy() takes (the
    Hermiticity gate), so the energies it gates have to be untouched."""
    n = 6
    sc = spinchain.Spin_Chain(["1/2"]*n,itensor_version=itensor_version)
    h = sum([sc.SS(i,i+1) for i in range(n-1)],0)
    sc.set_hamiltonian(h)
    assert sc.gs_energy(mode="DMRG")==pytest.approx(sc.gs_energy(mode="ED"),
            abs=1e-6)
