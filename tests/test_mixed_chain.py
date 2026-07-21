"""Functional coverage for mixedchain.Mixed_Spin_Fermion_Chain, a chain
that mixes genuine spin sites and spinful-fermion locations (paired
spinless-fermion up/down sites, exactly like
fermionchain.Spinful_Fermionic_Chain) in the same chain.

Degenerate cases (all-spin / all-fermion sitesin) are cross-checked
against the existing, dedicated ED backends of spinchain.Spin_Chain /
fermionchain.Spinful_Fermionic_Chain. The genuinely mixed case has no
existing ED backend (see CLAUDE.md's mixed-chain notes), so it is
cross-checked against a small, from-scratch, dense NumPy Jordan-Wigner
reference built directly in this file (independent of
multioperatortk/jordanwigner.py) -- this is also what exercises the
case a Jordan-Wigner string has to cross a non-fermionic (spin) site,
which itensor_version=3 handles via ITensor's own SiteSet::op()
"F"-at-non-fermionic-site -> identity fallback (see mixedchain.py's
docstring).

Only itensor_version=3 (and, for one cross-check, "python") is
exercised: itensor_version=2 is not supported for this chain type (no
"F"-fallback) and is rejected at construction time, see mixedchain.py."""
import numpy as np
import pytest

from dmrgpy import mixedchain, spinchain, fermionchain

from _dense_jw_reference import SX, SY, SZ, embed, fermion_ops, gs as _gs

DMRG_TOL = 1e-6

# ---------------------------------------------------------------------
# Degenerate cases: mixed chain reduces to the existing pure classes
# ---------------------------------------------------------------------

def test_pure_spin_degenerate_case_matches_spin_chain_ed():
    """All-spin sitesin, non-uniform exchange/field: Mixed_Spin_Fermion_Chain
    DMRG (v3) must agree with spinchain.Spin_Chain's own ED backend."""
    labels = ["1/2", "1", "1/2"]

    sc = spinchain.Spin_Chain(labels)
    h_sc = (0.7*sc.SS(0, 1) + 1.3*sc.SS(1, 2)
            + 0.4*sc.Sz[0] - 0.9*sc.Sz[1] + 0.2*sc.Sz[2])
    sc.set_hamiltonian(h_sc)
    e_ed = sc.gs_energy(mode="ED")

    mc = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version=3)
    h_mc = (0.7*mc.SS(0, 1) + 1.3*mc.SS(1, 2)
            + 0.4*mc.Sz[0] - 0.9*mc.Sz[1] + 0.2*mc.Sz[2])
    mc.set_hamiltonian(h_mc)
    e_dmrg = mc.gs_energy(mode="DMRG")

    assert e_dmrg == pytest.approx(e_ed, abs=DMRG_TOL)


def test_pure_fermion_degenerate_case_matches_spinful_fermionic_chain_ed():
    """All-fermion sitesin, non-uniform hopping/on-site terms:
    Mixed_Spin_Fermion_Chain DMRG (v3) must agree with
    fermionchain.Spinful_Fermionic_Chain's own ED backend."""
    n = 3

    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h_fc = 0
    for i in range(n-1):
        t = 0.6 + 0.3*i
        h_fc = h_fc + t*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h_fc = h_fc + h_fc.get_dagger()
    for i in range(n):
        h_fc = h_fc + (0.2*i - 0.1)*fc.Ntot[i]
    fc.set_hamiltonian(h_fc)
    e_ed = fc.gs_energy(mode="ED")

    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F"]*n, itensor_version=3)
    mc.maxm = 60
    mc.nsweeps = 20
    h_mc = 0
    for i in range(n-1):
        t = 0.6 + 0.3*i
        h_mc = h_mc + t*(mc.Cdagup[i]*mc.Cup[i+1] + mc.Cdagdn[i]*mc.Cdn[i+1])
    h_mc = h_mc + h_mc.get_dagger()
    for i in range(n):
        h_mc = h_mc + (0.2*i - 0.1)*mc.Ntot[i]
    mc.set_hamiltonian(h_mc)
    e_dmrg = mc.gs_energy(mode="DMRG")

    assert e_dmrg == pytest.approx(e_ed, abs=DMRG_TOL)


# ---------------------------------------------------------------------
# Genuine mixed chain, cross-checked against the dense JW reference
# ---------------------------------------------------------------------

# Kondo-like model parameters shared by _build_kondo_mc() and
# _build_kondo_reference(), so the two Hamiltonians they build can never
# silently drift apart. `hf` is a small longitudinal field breaking the
# model's full SU(2) spin-rotation symmetry -- the ground state is
# otherwise part of a degenerate spin multiplet, so individual per-site
# Sz expectation values (unlike the energy or the SU(2)-invariant Ntot)
# would depend on which arbitrary member of the multiplet DMRG/eigh
# happens to converge to, making them impossible to cross-check
# reproducibly.
KONDO_PARAMS = dict(K1=0.5, K2=0.35, t=0.4, muA=0.15, muB=-0.2, hf=0.13)


def _build_kondo_reference(nsites, fermionic_mask, spin_pos, fermA, fermB):
    """Dense reference Hamiltonian/observables for the fermion-spin-fermion
    Kondo-like model used below. `fermA`=(up,dn) physical indices of the
    first fermion location, `fermB`=(up,dn) of the second, `spin_pos` the
    physical index of the spin location in between."""
    C, Cdag, N = fermion_ops(fermionic_mask)
    Sx_s = embed(SX, spin_pos, nsites)
    Sy_s = embed(SY, spin_pos, nsites)
    Sz_s = embed(SZ, spin_pos, nsites)

    def ferm_spin_ops(up, dn):
        sx = 0.5*(Cdag[up]@C[dn]) + 0.5*(Cdag[dn]@C[up])
        sy = -0.5j*(Cdag[up]@C[dn]) + 0.5j*(Cdag[dn]@C[up])
        sz = 0.5*N[up] - 0.5*N[dn]
        ntot = N[up] + N[dn]
        return sx, sy, sz, ntot

    sxA, syA, szA, ntotA = ferm_spin_ops(*fermA)
    sxB, syB, szB, ntotB = ferm_spin_ops(*fermB)

    p = KONDO_PARAMS
    H = (p["K1"]*(sxA@Sx_s + syA@Sy_s + szA@Sz_s)
         + p["K2"]*(Sx_s@sxB + Sy_s@syB + Sz_s@szB)
         + p["muA"]*ntotA + p["muB"]*ntotB
         + p["hf"]*(szA + Sz_s + szB))
    hop = p["t"]*(Cdag[fermA[0]]@C[fermB[0]] + Cdag[fermA[1]]@C[fermB[1]])
    H = H + hop + hop.conj().T

    observables = {"Sz_spin": Sz_s, "Ntot_A": ntotA, "Ntot_B": ntotB}
    return H, observables


def _build_kondo_mc(mc):
    """Same Kondo-like model, built with Mixed_Spin_Fermion_Chain's own
    operator lists (logical locations 0="F", 1="1/2", 2="F"), from the
    same KONDO_PARAMS as _build_kondo_reference()."""
    p = KONDO_PARAMS
    h = (p["K1"]*mc.SS(0, 1) + p["K2"]*mc.SS(1, 2)
            + p["muA"]*mc.Ntot[0] + p["muB"]*mc.Ntot[2])
    hop = p["t"]*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
    h = h + hop + hop.get_dagger()
    h = h + p["hf"]*(mc.Sz[0] + mc.Sz[1] + mc.Sz[2])
    return h


def test_mixed_chain_kondo_energy_and_vev_vs_dense_jw_reference():
    """fermion-spin-fermion chain (['F','1/2','F']): the hopping term
    between the two fermion locations has the spin location physically
    in between, so its Jordan-Wigner string must correctly pass through
    (no sign, since a spin site carries no fermionic parity) -- this is
    exactly the case flagged in mixedchain.py's docstring as relying on
    ITensor v3's "F"-at-non-fermionic-site fallback."""
    labels = ["F", "1/2", "F"]
    mc = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version=3)
    mc.maxm = 60
    mc.nsweeps = 20
    assert mc.sites == [0, 0, 2, 0, 0]
    assert mc.phys_index == [(0, 1), 2, (3, 4)]

    h_mc = _build_kondo_mc(mc)
    mc.set_hamiltonian(h_mc)
    e_dmrg = mc.gs_energy(mode="DMRG")
    sz_dmrg = mc.vev(mc.Sz[1]).real
    ntot_dmrg = mc.vev(mc.Ntot[0]).real

    fermionic_mask = [True, True, False, True, True]
    H_ref, obs = _build_kondo_reference(5, fermionic_mask, spin_pos=2,
            fermA=(0, 1), fermB=(3, 4))
    e_ref, gs = _gs(H_ref)
    sz_ref = np.vdot(gs, obs["Sz_spin"] @ gs).real
    ntot_ref = np.vdot(gs, obs["Ntot_A"] @ gs).real

    assert e_dmrg == pytest.approx(e_ref, abs=DMRG_TOL)
    assert sz_dmrg == pytest.approx(sz_ref, abs=1e-4)
    assert ntot_dmrg == pytest.approx(ntot_ref, abs=1e-4)


def test_mixed_chain_kondo_v3_vs_python_backend():
    """Same Kondo-like model, cross-checked between the ITensor v3 C++
    backend and the pure-Python pyitensor backend -- both already
    support heterogeneous per-site type lists (see mixedchain.py)."""
    labels = ["F", "1/2", "F"]

    mc3 = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version=3)
    mc3.maxm = 60
    mc3.nsweeps = 20
    mc3.set_hamiltonian(_build_kondo_mc(mc3))
    e3 = mc3.gs_energy(mode="DMRG")

    mcp = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version="python")
    mcp.maxm = 60
    mcp.nsweeps = 20
    mcp.set_hamiltonian(_build_kondo_mc(mcp))
    ep = mcp.gs_energy(mode="DMRG")

    assert e3 == pytest.approx(ep, abs=DMRG_TOL)


# ---------------------------------------------------------------------
# Robustness: itensor_version=2 guard, ED-not-implemented, masked-zero
# operators, and the shared generate_bilinear() indexing fix
# ---------------------------------------------------------------------

def test_itensor_version_2_rejected_at_construction():
    """itensor_version=2 has no fallback for a Jordan-Wigner string that
    crosses a spin site and hard-aborts the whole process; the
    constructor must reject it with a clean Python exception instead."""
    with pytest.raises(ValueError):
        mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"], itensor_version=2)


def test_itensor_version_2_rejected_by_setup_cpp():
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"], itensor_version=3)
    with pytest.raises(ValueError):
        mc.setup_cpp(2)


def test_get_ED_obj_raises_clear_error():
    """No ED backend exists yet for a mixed Hilbert space; gs_energy/vev
    with mode="ED" (or the automatic itensor_version==3, ns<3 fallback in
    mode.py) must fail with a clear NotImplementedError, not an
    AttributeError several call frames away."""
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"], itensor_version=3)
    with pytest.raises(NotImplementedError):
        mc.get_ED_obj()

    small = mixedchain.Mixed_Spin_Fermion_Chain(["F"], itensor_version=3)
    assert small.ns < 3 # triggers mode.py's automatic DMRG->ED fallback
    small.set_hamiltonian(0.3*small.Ntot[0])
    with pytest.raises(NotImplementedError):
        small.gs_energy(mode="DMRG")


def test_masked_operators_are_multioperators_not_plain_zero():
    """Fermion-only operators at spin locations must be real (zero)
    MultiOperator instances, not the literal Python int 0: multiplying
    two such masked entries (e.g. both operands at spin locations) used
    to degrade to a plain int (0*0==0), breaking any later
    MultiOperator-only method call."""
    # 3 spin locations (ns=3) so the chain is large enough to avoid
    # mode.py's separate itensor_version==3/ns<3 automatic ED fallback
    # (see test_get_ED_obj_raises_clear_error) and isolate this bug.
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["1/2", "1/2", "1/2"], itensor_version=3)
    term = mc.Cdagup[0]*mc.Cup[1] # both operands are masked (spin locations)
    dagger = term.get_dagger() # must not raise AttributeError
    mc.set_hamiltonian(mc.SS(0, 1) + term + dagger)
    e = mc.gs_energy(mode="DMRG") # the masked term contributes nothing
    e_ss_only = pytest.approx(-0.75, abs=DMRG_TOL) # Heisenberg dimer singlet, site 2 decoupled
    assert e == e_ss_only


def test_vev_on_masked_operator_returns_zero():
    """vev() on a masked (fermion-only) operator at a spin location must
    cleanly return 0.0 instead of raising a confusing RuntimeError from
    obj2MO silently dropping the requested MultiOperator name."""
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["1/2", "F"], itensor_version=3)
    mc.set_hamiltonian(mc.SS(0, 1))
    mc.gs_energy(mode="DMRG")
    assert mc.vev(mc.Ntot[0]) == pytest.approx(0.0, abs=1e-10) # location 0 is spin


def test_generate_bilinear_works_on_logical_location_lists():
    """generate_bilinear() used to index unconditionally up to self.ns
    (physical sites), crashing with IndexError whenever the operator
    lists passed in are shorter (logical-location-indexed, as for any
    chain with a fermion location)."""
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"], itensor_version=3)
    assert len(mc.Sx) < mc.ns
    h = mc.generate_bilinear(lambda i, j: 1.0 if i == j else 0.0, mc.Sx, mc.Sx)
    assert h is not None
