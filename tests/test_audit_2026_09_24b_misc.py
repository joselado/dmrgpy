"""Regression tests for the `misc` cluster of the second 2026-09-24 hole
hunt (docs/audit_2026_09_24b_hole_hunt.md, findings 1, 6, 17 and 18).

   1. `mpsalgebra.disentangle_manifold` chose between `eigh` and `eig` on
      the bare, one-sided canonical-form proof, so a Hermitian operator the
      proof cannot see (`1j*Sx0*Sy0`, exactly `-Sz0/2`; the hopping
      `C0*Adag1 + Cdag0*A1`) went to `eig`, which does not orthogonalize
      inside a degenerate eigenspace: the output basis was non-orthonormal
      by 0.32 to 0.35 on the hunter's manifolds (0.165 on a Z3 chain). It
      now decides on the chain's own `is_hermitian` (proof, then probe) and
      takes `eigh` of the Hermitian part of the representation.
   6. `kpm_finite`'s `window_chain_kwargs` was a bare `setattr` loop, so a
      misspelled key was stored where nothing reads it and the default
      spectrum came back bit for bit; `itensor_version`/`mode` were accepted
      and ignored. Both now raise `TypeError`.
  17. Both correlator-based Kondo terms assumed an increasing `es`: on a
      coarse grid with a refined block appended, the second-order term was
      3.04 off a 6.97 peak and the potential term 0.654 off a 0.671 peak.
      Both helpers sort the grid now.
  18. The padding strip at the entry of the one-site TDVP route read the
      global `set_pad_bonds` flag, not the state, so a state padded earlier
      and evolved with the flag off ran one-site TDVP on the padded
      manifold, 0.19 off the unpadded trajectory on the 8-site Neel quench
      below. It is keyed on the state now, and an unpadded state is left
      bit for bit as it was.

Findings 1, 17 and the probe run on exact manifolds or exact Lehmann
correlators, so no convergence enters; 6 and 18 run on
itensor_version="python" with pinned seeds and sweep schedules.
"""

import warnings

import numpy as np
import pytest
import scipy.linalg as dlg

from dmrgpy import fermionchain, infinitechain, mpsalgebra, spinchain, timedependent
from dmrgpy.mpsalgebratk import disentangle
from dmrgpy.mpsalgebratk.disentangle import get_representation
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import chain as chainmod
from dmrgpy.pyitensor.mpscontainer import _link_at


# ------------------------------------------------------------ finding 1

@pytest.fixture(scope="module")
def spin_manifold():
    """The four eigenstates of a generic 2-site spin-1/2 Hamiltonian at full
    bond dimension (the hunter's script 08), orthonormal to 1e-15."""
    sc = spinchain.Spin_Chain([2]*2, itensor_version="python")
    h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[0]*sc.Sy[1] + 0.3*sc.Sz[0]*sc.Sz[1] \
        + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 8, 20
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        es, wfs = sc.get_excited_states(n=4)
    return sc, np.real(es), wfs


@pytest.fixture(scope="module")
def fermion_manifold():
    """The same for a 2-site spinless fermion chain with pairing."""
    fc = fermionchain.Fermionic_Chain(2, itensor_version="python")
    h = 0.8*fc.Cdag[0]*fc.C[1] + 0.8*fc.Cdag[1]*fc.C[0] \
        + 0.5*(fc.C[0]*fc.C[1] + fc.Cdag[1]*fc.Cdag[0]) \
        + 0.3*fc.N[0] - 0.4*fc.N[1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 8, 20
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        es, wfs = fc.get_excited_states(n=4)
    return fc, np.real(es), wfs


@pytest.fixture(scope="module")
def parafermion_manifold():
    """The nine eigenstates of a generic 2-site Z3 clock chain: every
    parafermionic name is off canonical.py's _PARITY table, so no operator
    on it is ever proven Hermitian (the reviewer's third family)."""
    from dmrgpy import parafermionchain
    pc = parafermionchain.Parafermionic_Chain(2, itensor_version="python")
    h = -(pc.Sig[0]*pc.Sigd[1] + pc.Sig[1]*pc.Sigd[0]) \
        - 0.7*(pc.Tau[0] + pc.Taud[0]) - 0.45*(pc.Tau[1] + pc.Taud[1]) \
        - 0.2*(pc.Sig[0] + pc.Sigd[0])
    pc.set_hamiltonian(h)
    pc.maxm, pc.nsweeps = 12, 20
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        es, wfs = pc.get_excited_states(n=9)
    return pc, np.real(es), wfs


def _manifold(request, which):
    return request.getfixturevalue(which + "_manifold")


def _gram(ws):
    return np.array([[a.dot(b) for b in ws] for a in ws])


def _coefficients(wfs, out):
    """C[i, j] = <w_i|out_j>, the output in the manifold's own coordinates"""
    return np.array([[a.dot(b) for b in out] for a in wfs])


def _unproven(chain, which):
    if which == "spin":
        return 1j*chain.Sx[0]*chain.Sy[0]            # exactly -Sz0/2
    if which == "parafermion":
        return chain.Sig[0] + chain.Sigd[0]          # spectrum 2,-1,-1 per site
    return chain.C[0]*chain.Adag[1] + chain.Cdag[0]*chain.A[1] # a hopping


@pytest.mark.parametrize("which", ["spin", "fermion", "parafermion"])
def test_unproven_hermitian_operator_gives_an_orthonormal_eigenbasis(request, which):
    chain, es, wfs = _manifold(request, which)
    A = _unproven(chain, which)
    assert not A.is_hermitian() # the proof cannot see it, which is the point
    np.random.seed(7)           # the probe draws its witness from np.random
    out = mpsalgebra.disentangle_manifold(wfs, A)
    assert np.max(np.abs(_gram(out) - np.eye(len(out)))) < 1e-12
    # the output spans the manifold, so H on it gives the levels back
    levels = np.sort(np.linalg.eigvalsh(get_representation(out, chain.hamiltonian)))
    assert np.max(np.abs(levels - np.sort(es))) < 1e-12
    # and A is diagonal on it, with its degenerate spectrum
    ma_out = get_representation(out, A)
    assert np.max(np.abs(ma_out - np.diag(np.diag(ma_out)))) < 1e-12


@pytest.mark.parametrize("which", ["spin", "fermion", "parafermion"])
def test_bare_proof_would_have_taken_eig(request, which, monkeypatch):
    """The pre-fix decision, rebuilt: the O(1) Gram error comes back."""
    chain, es, wfs = _manifold(request, which)
    monkeypatch.setattr(disentangle, "_is_hermitian", lambda w, A: A.is_hermitian())
    out = mpsalgebra.disentangle_manifold(wfs, _unproven(chain, which))
    # 0.35, 0.32 and 0.165 measured
    assert np.max(np.abs(_gram(out) - np.eye(len(out)))) > 1e-3


def _projectors(evals, V, tol=1e-8):
    """{eigenvalue: projector onto its eigenspace} from the columns of V"""
    out = {}
    for e in np.unique(np.round(np.real(evals)/tol)*tol):
        cols = V[:, np.abs(np.real(evals) - e) < 10*tol]
        out[e] = cols @ cols.conj().T
    return out


@pytest.mark.parametrize("which,spelling", [("spin", "non-degenerate"),
                                            ("spin", "degenerate"),
                                            ("fermion", "non-degenerate")])
def test_proven_operator_output_is_unchanged(request, which, spelling):
    """A proven operator takes eigh as before; only the matrix eigh reads
    moved, from ma's lower triangle to its Hermitian part, which differ at
    roundoff on an exact manifold. So the eigenvectors of a non-degenerate
    operator agree with the old ones up to a phase, and a degenerate
    operator's eigenspaces agree (the basis inside one is arbitrary, and
    a 1e-16 perturbation may rotate it)."""
    chain, es, wfs = _manifold(request, which)
    if which == "fermion":
        A = 0.3*chain.N[0] + 0.1*chain.N[1]
    elif spelling == "non-degenerate":
        A = 0.3*chain.Sz[0] + 0.1*chain.Sz[1]
    else:
        A = -0.5*chain.Sz[0]
    assert A.is_hermitian()
    ma = get_representation(wfs, A)
    old_e, old_V = dlg.eigh(ma) # the pre-fix branch, verbatim
    new_V = _coefficients(wfs, mpsalgebra.disentangle_manifold(wfs, A))
    if spelling == "non-degenerate":
        overlaps = np.abs(np.sum(old_V.conj()*new_V, axis=0))
        assert np.max(np.abs(overlaps - 1)) < 1e-12
    else:
        old_P, new_P = _projectors(old_e, old_V), _projectors(old_e, new_V)
        assert sorted(old_P) == sorted(new_P)
        for e in old_P:
            assert np.max(np.abs(old_P[e] - new_P[e])) < 1e-12


def test_non_hermitian_operator_still_takes_eig(spin_manifold):
    chain, es, wfs = spin_manifold
    A = chain.Sx[0] + 1j*chain.Sy[0] + 0.3*chain.Sz[0] + 0.1*chain.Sz[1]
    np.random.seed(7)
    new_V = _coefficients(wfs, mpsalgebra.disentangle_manifold(wfs, A))
    old_V = dlg.eig(get_representation(wfs, A))[1]
    assert np.max(np.abs(new_V - old_V)) < 1e-12


def test_states_without_a_chain_fall_back_to_the_bare_proof(spin_manifold):
    chain = spin_manifold[0]

    class Stateless:
        MBO = None

    assert disentangle._is_hermitian([Stateless()], chain.Sz[0])
    assert not disentangle._is_hermitian([Stateless()], 1j*chain.Sx[0]*chain.Sy[0])


# ------------------------------------------------------------ finding 6

@pytest.fixture(scope="module")
def heisenberg_ic():
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0])
    return ic


def _kpm_finite(ic, **wk):
    np.random.seed(3)
    x, y = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=4,
                         window_chain_kwargs=dict(maxm=10, nsweeps=4, **wk),
                         delta=0.2, es=np.linspace(-0.5, 3.0, 71))
    return np.real(np.asarray(y))


@pytest.mark.parametrize("typo", ["kpm_nscale", "kpmscale", "max_m", "nsweep"])
def test_misspelled_window_key_raises(heisenberg_ic, typo):
    with pytest.raises(TypeError, match=r"unknown window_chain_kwargs key\(s\): %s\." % typo):
        _kpm_finite(heisenberg_ic, **{typo: 3})


def test_every_unknown_window_key_is_named_sorted(heisenberg_ic):
    with pytest.raises(TypeError, match=r"key\(s\): kpm_nscale, max_m\."):
        _kpm_finite(heisenberg_ic, max_m=2, kpm_nscale=3)


@pytest.mark.parametrize("key,value", [("itensor_version", 3), ("mode", "ED")])
def test_backend_keys_are_rejected_by_name(heisenberg_ic, key, value):
    with pytest.raises(TypeError, match=r"cannot set %s: .*itensor_version=\"python\"" % key):
        _kpm_finite(heisenberg_ic, **{key: value})


@pytest.mark.parametrize("key", ["gs_energy", "_session"])
def test_method_or_private_window_key_raises(heisenberg_ic, key):
    # both pass hasattr(), and setting either would overwrite the chain's
    # own behaviour or state rather than a setting
    with pytest.raises(TypeError, match=r"unknown window_chain_kwargs key\(s\): %s" % key):
        _kpm_finite(heisenberg_ic, **{key: None})


def test_correct_window_key_still_moves_the_spectrum(heisenberg_ic):
    default = _kpm_finite(heisenberg_ic)
    tripled = _kpm_finite(heisenberg_ic, kpm_n_scale=3)
    # three times the moments at the same delta: the peak is 0.28 against
    # 0.80, measured
    assert np.max(np.abs(tripled - default)) > 0.5*np.max(default)


# ------------------------------------------------------------ finding 17

G, MUB = 2.0, 5.7883818066e-5 # eV/T


@pytest.fixture(scope="module")
def zeeman_spin():
    """The single S=1/2 at 10 T of tests/test_audit_2026_09_24_kondo.py"""
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
    return sc


def _both_terms(sc, es):
    from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
    from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
    eVs = np.linspace(-1e-3, 2e-3, 21)
    kw = dict(T0=1.0, mode="ED", submode="ED", delta=2e-6, es=es)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p = third_order_potential_dIdV_dc(sc, 0, eVs, 0.1, 0.3, **kw)
    s = second_order_dIdV_dc(sc, 0, eVs, U=0.3, **kw)
    return p, s, [w for w in caught if "third_order_potential_dIdV_dc" in str(w.message)]


_COARSE = np.linspace(-1e-3, 3e-3, 8000)
_FINE = np.linspace(1.0e-3, 1.4e-3, 4000) # refined around the 1.16 meV line


@pytest.mark.parametrize("label,grid", [
    ("coarse then fine", np.concatenate([_COARSE, _FINE])),
    ("fine then coarse", np.concatenate([_FINE, _COARSE])),
    ("descending", np.sort(np.concatenate([_COARSE, _FINE]))[::-1].copy()),
])
def test_kondo_terms_do_not_depend_on_the_order_of_es(zeeman_spin, label, grid):
    p_ref, s_ref, _ = _both_terms(zeeman_spin, np.sort(grid))
    p, s, caught = _both_terms(zeeman_spin, grid)
    # peaks 0.671 and 6.97; unsorted, these were 0.654 and 3.04 off
    assert np.max(np.abs(p - p_ref)) < 1e-12*np.max(np.abs(p_ref))
    assert np.max(np.abs(s - s_ref)) < 1e-12*np.max(np.abs(s_ref))
    # the sum-rule check sees the right weight, 0.7495 of 0.75
    assert caught == []


def test_kondo_terms_accept_es_as_a_list(zeeman_spin):
    grid = np.concatenate([_COARSE, _FINE])
    p_ref, s_ref, _ = _both_terms(zeeman_spin, np.sort(grid))
    p, s, _ = _both_terms(zeeman_spin, list(grid))
    assert np.max(np.abs(p - p_ref)) < 1e-12*np.max(np.abs(p_ref))
    assert np.max(np.abs(s - s_ref)) < 1e-12*np.max(np.abs(s_ref))


# ------------------------------------------------------------ finding 18

N18, K18, NT18, DT18 = 8, 4, 20, 0.05


def _prepared(start, pad_gs, K=K18):
    """(chain with the quench Hamiltonian set, its ground state before the
    quench), the ground state computed under set_pad_bonds(pad_gs). "neel"
    is the staggered-field product state (true bonds all 1), "dimer" the
    ground state of a dimerized Heisenberg chain in a staggered field."""
    np.random.seed(5)
    bk.set_pad_bonds(pad_gs)
    c = spinchain.Spin_Chain([2]*N18, itensor_version="python")
    c.maxm = K
    h0 = 0
    if start == "neel":
        c.nsweeps, c.cutoff = 10, 1e-12
        for i in range(N18):
            h0 = h0 + (-1)**i*c.Sz[i]
    else:
        c.nsweeps, c.cutoff = 12, 1e-10
        for i in range(N18 - 1):
            J = 1.0 if i % 2 == 0 else 0.2
            h0 = h0 + J*(c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1])
        for i in range(N18):
            h0 = h0 + 0.3*(-1)**i*c.Sz[i]
    c.set_hamiltonian(h0)
    wf = c.get_gs()
    h1 = 0
    for i in range(N18 - 1):
        h1 = h1 + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + 0.7*c.Sz[i]*c.Sz[i+1]
    for i in range(N18):
        h1 = h1 + 0.1*c.Sz[i]
    c.set_hamiltonian(h1)
    c.tevol_method = "TDVP_GSE"
    return c, wf


def _bonds(wf):
    return [_link_at(wf.cpp_handle, i, i + 1).dim for i in range(1, N18)]


def _evolve(start, pad_gs, sweeps, flag_off="clear", K=K18):
    """<Sz_0>(t) after the quench, with the flag off during the evolution
    (cleared, or on but suspended), plus the caller's wf bonds after it."""
    try:
        c, wf = _prepared(start, pad_gs, K)
        c.tdvp_gse_sweeps = sweeps
        if flag_off == "clear":
            bk.set_pad_bonds(None)
            out = timedependent.evolve_and_measure(c, operator=c.Sz[0], nt=NT18,
                                                   dt=DT18, wf=wf)
        else:
            with bk.pad_bonds_suspended():
                out = timedependent.evolve_and_measure(c, operator=c.Sz[0], nt=NT18,
                                                       dt=DT18, wf=wf)
        return np.real(np.asarray(out[1])), _bonds(wf)
    finally:
        bk.set_pad_bonds(None)


def _flag_gated_strip(psi):
    """The pre-fix _strip_bond_padding: a no-op whenever the flag is off."""
    if not bk.pad_bonds() or psi.length() < 2:
        return psi
    return chainmod._strip_sweep(psi)


@pytest.mark.parametrize("flag_off", ["clear", "suspend"])
def test_padded_state_evolved_with_the_flag_off_follows_the_unpadded_run(flag_off, monkeypatch):
    """tdvp_gse_sweeps=0: one-site TDVP from the Neel product state keeps
    bond dimension 1, so <Sz_0> stays at -1/2 unpadded; the padded state
    used to leave it, by 0.1925 over these 20 steps."""
    unpadded, _ = _evolve("neel", None, 0)
    padded, wf_bonds = _evolve("neel", K18, 0, flag_off)
    assert np.max(np.abs(padded - unpadded)) < 1e-12
    assert wf_bonds == [K18]*(N18 - 1) # the strip ran on a copy
    monkeypatch.setattr(chainmod, "_strip_bond_padding", _flag_gated_strip)
    before, _ = _evolve("neel", K18, 0, flag_off)
    assert np.max(np.abs(before - unpadded)) > 1e-2


def test_padded_entangled_state_evolved_with_the_flag_off_follows_the_unpadded_run():
    """The same on an entangled start (true bonds [2,4,8,5,8,4,2] at K=8),
    where the effect sat inside the method's own error: 7.0e-8 before."""
    unpadded, _ = _evolve("dimer", None, 0, K=8)
    padded, _ = _evolve("dimer", 8, 0, K=8)
    assert np.max(np.abs(padded - unpadded)) < 1e-12


@pytest.mark.parametrize("start,sweeps", [("neel", 0), ("neel", 3),
                                          ("dimer", 0), ("dimer", 3)])
def test_unpadded_runs_are_bit_identical(start, sweeps, monkeypatch):
    """With the flag off and nothing to strip, the state goes into the
    evolution exactly as it came in, as the flag-gated strip left it; an
    ungated sweep would re-gauge it (7e-16 to 4e-11, measured by the
    reviewer)."""
    K = K18 if start == "neel" else 8
    fixed, _ = _evolve(start, None, sweeps, K=K)
    monkeypatch.setattr(chainmod, "_strip_bond_padding", lambda psi: psi)
    today, _ = _evolve(start, None, sweeps, K=K)
    assert np.array_equal(fixed, today)
