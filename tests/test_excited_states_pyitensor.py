"""Excited states on the pure-Python backend (`itensor_version="python"`),
cross-checked against ED.

This file exists because that combination had no coverage at all --
`tests/_helpers.py`'s `versions=` kwarg spans the two compiled backends
only (see CLAUDE.md) -- and a real bug lived in the gap for a long time:
`pyitensor/dmrg.py::_bond_projections` conjugated the overlap-penalty
projector's environments once too few times, so for any complex-valued
state the penalty projected onto the wrong direction and barely
orthogonalized at all. Measured before the fix, on the 8-site spin-1/2
chain below: excited energies off by up to 0.55, and returned states with
mutual overlaps up to 0.81. The compiled backends were unaffected (they
use ITensor's own Weight overload), so every existing cross-backend test
passed throughout.

Two properties are asserted throughout, because the energies alone are not
enough to catch that class of bug -- an excited-state search that returns
the *same* state twice reports plausible energies:

- the energies match ED, index by index (a degenerate level appears once
  per member on both sides, see docs/ed_vs_dmrg_degenerate_multiplets.md);
- the returned states are mutually orthonormal.

The Hamiltonians deliberately include terms that make the state genuinely
complex (`Sy`, a Dzyaloshinskii-Moriya term), since the bug was exactly
invisible for real-valued tensors.
"""

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain


def _spin_chain(n=8, dm=0.0):
    sc = spinchain.Spin_Chain(["1/2"] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        if dm:   # Dzyaloshinskii-Moriya: makes the ground state complex
            h = h + dm*(sc.Sx[i]*sc.Sy[i+1] - sc.Sy[i]*sc.Sx[i+1])
    for i in range(n):
        h = h + 0.3*sc.Sz[i] + 0.13*(-1)**i*sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc


def _hubbard(cls, n=4):
    fc = cls(n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        for s in (0, 1):
            h = h + fc.Cdag[2*i+s]*fc.C[2*(i+1)+s] + fc.Cdag[2*(i+1)+s]*fc.C[2*i+s]
    for i in range(n):
        h = h + 2.0*fc.Cdag[2*i]*fc.C[2*i]*fc.Cdag[2*i+1]*fc.C[2*i+1]
        # Zeeman + chemical potential: lifts the spin/charge degeneracy, so
        # ED and DMRG can be compared index by index without the
        # degenerate-multiplet ambiguity.
        h = h + 0.35*(fc.Cdag[2*i]*fc.C[2*i] - fc.Cdag[2*i+1]*fc.C[2*i+1])
        h = h - 0.7*(fc.Cdag[2*i]*fc.C[2*i] + fc.Cdag[2*i+1]*fc.C[2*i+1])
    fc.set_hamiltonian(h)
    return fc


def _check(chain, nex, tol=1e-6):
    ref = np.array(chain.get_excited(n=nex, mode="ED")).real
    got = np.array(chain.get_excited(n=nex, mode="DMRG")).real
    assert got == pytest.approx(ref, abs=tol)
    return ref, got


@pytest.mark.parametrize("dm", [0.0, 0.25])
def test_spin_chain_excited_energies_match_ed(dm):
    _check(_spin_chain(dm=dm), nex=6)


@pytest.mark.parametrize("cls,label", [
    (fermionchain.Spinful_Fermionic_Chain, "interleaved"),
    (fermionchain.Spinful_Fermionic_Chain_Native, "native"),
])
def test_spinful_fermion_excited_energies_match_ed(cls, label):
    """Both spinful classes: the bug was found via the native one but hit
    the interleaved one harder (1.20 vs 0.54 before the fix)."""
    _check(_hubbard(cls), nex=8)


def test_returned_states_are_orthonormal():
    """The property the energies alone would not catch. Before the fix the
    largest off-diagonal overlap here was 0.81."""
    sc = _spin_chain(dm=0.25)
    n = 6
    _, wfs = sc.get_excited_states(n=n, purify=False)
    ov = np.array([[complex(wfs[a].dot(wfs[b])) for b in range(n)]
                   for a in range(n)])
    assert np.abs(np.diag(ov)) == pytest.approx(np.ones(n), abs=1e-6)
    off = np.abs(ov - np.diag(np.diag(ov)))
    assert off.max() < 1e-6


def test_penalty_projector_reproduces_the_exact_overlap():
    """Unit test of the fix itself, one level below the eigenvalues.

    `_bond_projections` returns, per already-found state, the vector
    `proj_k` for which `vdot(proj_k, v)` is the overlap `<wfs_k|Psi(v)>`
    of that state with the full MPS built by putting the two-site tensor
    `v` into psi's window. So contracting it against psi's OWN two-site
    tensor must return the plain full-MPS overlap `<wfs_k|psi>` -- and
    must return the same number at every bond, since that overlap is one
    scalar however the chain is sliced. Before the fix it was neither: off
    by 0.05-0.14 and different at every bond."""
    from dmrgpy.pyitensor import dmrg as pdmrg
    from dmrgpy.pyitensor.mpsalgebra import inner

    sc = _spin_chain(n=6, dm=0.25)
    sc.gs_energy()
    ch = sc._session
    psi = ch.wf0.copy(); psi.normalize()
    wk = ch._default_mps(); wk.normalize()
    H = ch._solve_H()

    exact = complex(inner(wk, psi))
    right_ov = [pdmrg._all_overlap_right(wk, psi)]
    left_ov = [{0: None}]
    left_env = {0: (None, None)}
    right_env = pdmrg._all_right_environments(H, psi)
    for i in range(1, psi.length()):
        Lh, Lbrah = left_env[i - 1]
        Rh, Rbrah = right_env[i + 2]
        _, order_in, _, x0 = pdmrg.two_site_heff(Lh, Lbrah, H, psi, i, Rh, Rbrah)
        proj = pdmrg._bond_projections([wk], left_ov, right_ov, i)[0]
        got = np.vdot(proj.transpose_to(order_in).reshape(-1), x0)
        assert got == pytest.approx(exact, abs=1e-10)
        left_env[i] = pdmrg._extend_left(Lh, Lbrah, H, psi, i)
        left_ov[0][i] = pdmrg._extend_overlap_left(left_ov[0][i - 1], wk, psi, i)
