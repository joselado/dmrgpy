"""Coverage for Infinite_Many_Body_Chain's iMPS algebra: `apply_operator`,
`overlap` and `summps`, the public wrappers over pyitensor.idmrg's
apply_mpo/imps_overlap/imps_sum (gs_method="idmrg") and pyitensor.vumps's
apply_mpo/imps_sum (gs_method="vumps"). The module functions themselves are
covered by tests/test_infinite_chain.py, tests/test_vumps_apply_mpo.py and
tests/test_vumps_imps_sum.py; what is pinned here is the wrapper layer: the
MultiOperator-to-MPO construction (site order, same-site composition order,
the central-cell restriction), the normalization apply_operator adds (the
underlying apply_mpo refuses a state whose norm changed, so the operator is
divided by the measured per-cell norm and that norm is reported as `.eta`),
the derived chain's behaviour (vev/correlator read the new state, gs_energy
and the environment-based methods refuse it), the overlap as a fidelity per
site (the spectral radius of the mixed transfer matrix, real, so a state
and its spin-flipped copy, whose top mixed eigenvalues tie in modulus and
differ in sign, compare without error), and the backend/gs_method guards.

Cross-checks used: the exact D=1 field-polarized product state (Pauli-X
flips <Sz> exactly, Sp on the down state gives the up state, Sm*Sp leaves
it), unitary-operator invariants on the gapped TFIM at D>1 (<Sz> flips,
<Sz(0)Sz(r)> unchanged, eta 1), the exact squared norm of a non-unitary
operator (Sz on spin-1/2 has Sz^2 = 1/4, so eta = 1/4 and the state is the
Z-rotated one, with every Sz observable unchanged), the per-site fidelity's
three exact values (1, 1 up to gauge, 0 for orthogonal product states), and
a free-fermion two-site cell where the intra-cell hopping unitary commutes
with the intra-cell hopping, so that correlator is unchanged by it.
"""
import copy

import numpy as np
import pytest

from dmrgpy import infinitechain

GS_METHODS = ("idmrg", "vumps")
_FERMION_SITE_CODE = 0  # spinless fermion site (C/Cdag/F/N), as in test_idmrg_window_fermionic.py


def _field_chain(B, gs_method, spins=("1/2",), maxm=1):
    """H = -B*sum_i SzC[i]: an exact D=1 product state polarized along
    sign(B), on either solver."""
    ic = infinitechain.Infinite_Spin_Chain(list(spins), itensor_version="python")
    ic.gs_method = gs_method
    ic.maxm = maxm
    h = -B * ic.SzC[0]
    for i in range(1, len(spins)):
        h = h - B * ic.SzC[i]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _tfim_chain(g, gs_method, maxm):
    """H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] (the convention of
    test_vumps_apply_mpo.py), gapped at g=1.5 so both solvers converge at
    a small bond dimension."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.gs_method = gs_method
    ic.maxm = maxm
    ic.etol = 1e-11
    ic.maxiter = 200
    ic.set_hamiltonian(-4.0 * ic.SxC[0] * ic.SxR[0] - 2.0 * g * ic.SzC[0])
    ic.gs_energy()
    return ic


# -- apply_operator ---------------------------------------------------------

@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_apply_operator_pauli_x_flips_polarized_chain(gs_method):
    """D=1 field-polarized product state: Pauli-X (2*Sx) at every site
    flips <Sz> exactly, and being unitary gives eta exactly 1."""
    ic = _field_chain(5.0, gs_method)
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
    flipped = ic.apply_operator(2 * ic.SxC[0])
    assert flipped.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    assert flipped.eta == pytest.approx(1.0, abs=1e-8)
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)  # the original is untouched


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_apply_operator_conventions_on_down_state(gs_method):
    """The in/out layout and the same-site composition order, pinned on
    the exact down-polarized state: Sp turns it into the up state (a
    transposed layout would annihilate it), Sm*Sp (Sp applied first)
    leaves it unchanged (the reversed product annihilates it), and
    Sp*Sp, identically zero, raises ValueError rather than reaching the
    eigensolver with an all-zero transfer matrix."""
    down = _field_chain(-5.0, gs_method)
    assert down.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    up = down.apply_operator(down.get_operator("Sp", 0))
    assert up.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
    assert up.eta == pytest.approx(1.0, abs=1e-8)
    same = down.apply_operator(down.get_operator("Sm", 0) * down.get_operator("Sp", 0))
    assert same.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    assert same.eta == pytest.approx(1.0, abs=1e-8)
    with pytest.raises(ValueError):
        down.apply_operator(down.get_operator("Sp", 0) * down.get_operator("Sp", 0))


@pytest.mark.parametrize("gs_method", GS_METHODS)
@pytest.mark.parametrize("D", [2, 3])
def test_apply_operator_unitary_invariants_tfim(gs_method, D):
    """Genuinely entangled D>1 (gapped TFIM): Pauli-X at every site
    anticommutes with Sz, so <Sz> flips sign and <Sz(0)Sz(r)> is
    unchanged, and eta stays 1."""
    ic = _tfim_chain(1.5, gs_method, D)
    flipped = ic.apply_operator(2 * ic.SxC[0])
    assert flipped.eta == pytest.approx(1.0, abs=1e-6)
    sz = ic.vev("Sz", 0)
    assert abs(sz) > 0.3  # the paramagnet is polarized along the field
    assert flipped.vev("Sz", 0) == pytest.approx(-sz, abs=1e-6)
    for r in (1, 2, 3):
        c0 = ic.correlator("Sz", 0, "Sz", r)
        c1 = flipped.correlator("Sz", 0, "Sz", r)
        assert c1 == pytest.approx(c0, abs=1e-6)


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_apply_operator_non_unitary_is_normalized(gs_method):
    """A non-unitary operator: Sz at every site. Sz^2 = 1/4 on spin-1/2,
    so <psi|W^dag W|psi> per site is exactly 1/4, which is what `.eta`
    reports, while the stored state is the normalized one, Z|psi> with
    Z the Pauli matrix, so every Sz observable is unchanged. 0.9*Id gives
    eta 0.81 and the same physical state (fidelity 1)."""
    ic = _tfim_chain(1.5, gs_method, 2)
    w = ic.apply_operator(ic.SzC[0])
    assert w.eta == pytest.approx(0.25, abs=1e-8)
    assert w.vev("Sz", 0) == pytest.approx(ic.vev("Sz", 0), abs=1e-8)
    for r in (1, 2):
        assert w.correlator("Sz", 0, "Sz", r) == pytest.approx(
            ic.correlator("Sz", 0, "Sz", r), abs=1e-8)
    scaled = ic.apply_operator(0.9 * ic.get_operator("Id", 0))
    assert scaled.eta == pytest.approx(0.81, abs=1e-8)
    assert abs(ic.overlap(scaled)) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_apply_operator_two_site_cell(gs_method):
    """A ["1/2","1/2"] cell: one gate per cell, 2*SxC[0]*2*SxC[1], flips
    both sites; 2*SxC[0] alone flips only site 0, which pins the kron
    order of the cell operator (site 0 slowest)."""
    ic = _field_chain(5.0, gs_method, spins=("1/2", "1/2"))
    both = ic.apply_operator(2 * ic.SxC[0] * 2 * ic.SxC[1])
    assert both.eta == pytest.approx(1.0, abs=1e-8)
    assert both.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    assert both.vev("Sz", 1) == pytest.approx(-0.5, abs=1e-8)
    one = ic.apply_operator(2 * ic.SxC[0])
    assert one.eta == pytest.approx(1.0, abs=1e-8)
    assert one.vev("Sz", 0) == pytest.approx(-0.5, abs=1e-8)
    assert one.vev("Sz", 1) == pytest.approx(0.5, abs=1e-8)


def test_apply_operator_two_site_gate_grows_bond_dimension():
    """A sum of terms that does not factorize over the two sites, the
    unitary exp(-i*theta*S0.S1) written as alpha*Id + beta*S0.S1 (S.S has
    eigenvalues 1/4 and -3/4), tiled once per two-site cell on an
    entangled state: a chi_W>1 MPO, so the bond-growth path is exercised,
    and being unitary eta stays 1. S0.S1 commutes with the gate, so the
    intra-cell correlator is unchanged, while the inter-cell one moves.
    The gate grows the bond from maxm=4 to 16; truncating back to the
    default maxdim=maxm discards about 7.5e-4 of weight, which the gauge
    fixing cannot hold, and that is reported as a truncation error."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    ic.gs_method = "idmrg"
    h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
         + 0.4 * (ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0]))
    ic.set_hamiltonian(h)
    ic.maxm, ic.maxiter, ic.etol = 4, 100, 1e-12
    ic.gs_energy()
    theta = 0.37
    beta = np.exp(-1j * theta / 4) - np.exp(3j * theta / 4)
    alpha = np.exp(-1j * theta / 4) - beta / 4
    ss = ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
    gate = alpha * ic.get_operator("Id", 0) + beta * ss
    with pytest.raises(RuntimeError, match="discarded too much weight"):
        ic.apply_operator(gate, cutoff=1e-10)
    gated = ic.apply_operator(gate, cutoff=1e-10, maxdim=16)
    assert gated.eta == pytest.approx(1.0, abs=1e-6)
    intra = sum(ic.correlator(n, 0, n, 1) for n in ("Sx", "Sy", "Sz"))
    intra_gated = sum(gated.correlator(n, 0, n, 1) for n in ("Sx", "Sy", "Sz"))
    assert intra_gated == pytest.approx(intra, abs=1e-6)
    inter = ic.correlator("Sz", 1, "Sz", 1)
    inter_gated = gated.correlator("Sz", 1, "Sz", 1)
    assert abs(inter_gated - inter) > 1e-3


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_apply_operator_rejects_other_cells(gs_method):
    """A factor on the next or the previous cell is not a per-cell
    operator: ValueError, naming the offending site."""
    ic = _field_chain(5.0, gs_method)
    with pytest.raises(ValueError, match="outside the central cell"):
        ic.apply_operator(2 * ic.SxR[0])
    with pytest.raises(ValueError, match="outside the central cell"):
        ic.apply_operator(ic.SzC[0] * ic.SzL[0])


def test_apply_operator_fermionic_even_parity():
    """A spinless free-fermion two-site cell (gs_method="idmrg"): the
    intra-cell hopping unitary exp(-i*theta*h_cell), h_cell = Cdag0 C1 +
    Cdag1 C0, written as Id + (cos(theta)-1)*(N0 + N1 - 2 N0 N1) -
    i*sin(theta)*h_cell (h_cell^2 is the projector on the one-particle
    sector of the cell). It is parity-even, so it is accepted; unitary,
    so eta is 1; it commutes with h_cell, so <h_cell> is unchanged, and it
    conserves the cell's particle number. A parity-odd operator (Cdag
    alone) raises ValueError.

    maxm=4 deliberately: the gauge fixing inside apply_mpo is delicate
    once the grown bond dimension is far above the state's entanglement
    (its own docstring), and on this chain its isometry check fails at
    maxm=8 and 12 for every maxdim/cutoff tried, while at maxm=4 it holds
    for maxdim in (8, 16) with cutoff 1e-10 or 1e-8."""
    ic = infinitechain.Infinite_Many_Body_Chain(
        [_FERMION_SITE_CODE] * 2, itensor_version="python")
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter = 4, 60
    C = [ic.get_operator("C", i, "C") for i in range(2)]
    Cd = [ic.get_operator("Cdag", i, "C") for i in range(2)]
    CR = ic.get_operator("C", 0, "R")
    CdR = ic.get_operator("Cdag", 0, "R")
    N = [ic.get_operator("N", i, "C") for i in range(2)]
    h_cell = Cd[0] * C[1] + Cd[1] * C[0]
    ic.set_hamiltonian(-1.0 * h_cell - 0.6 * (Cd[1] * CR + CdR * C[1])
                       + 0.3 * (N[0] - N[1]))
    ic.gs_energy()
    theta = 0.41
    gate = (ic.get_operator("Id", 0)
            + (np.cos(theta) - 1.0) * (N[0] + N[1] - 2.0 * N[0] * N[1])
            - 1j * np.sin(theta) * h_cell)
    rotated = ic.apply_operator(gate, cutoff=1e-10, maxdim=2 * ic.maxm)
    assert rotated.eta == pytest.approx(1.0, abs=1e-6)

    def hop(chain):
        # <h_cell> = <Cdag_0 C_1> + <Cdag_1 C_0> = 2 Re <Cdag_0 C_1>
        return 2.0 * np.real(chain.correlator("Cdag", 0, "C", 1))

    def filling(chain):
        return chain.vev("N", 0) + chain.vev("N", 1)

    assert abs(hop(ic)) > 0.1
    assert hop(rotated) == pytest.approx(hop(ic), abs=1e-6)
    assert filling(rotated) == pytest.approx(filling(ic), abs=1e-6)
    with pytest.raises(ValueError, match="parity"):
        ic.apply_operator(Cd[0])


# -- the derived chain -------------------------------------------------------

@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_derived_chain_refuses_ground_state_methods(gs_method):
    """The chain apply_operator returns holds a state that is not an
    eigenstate: gs_energy and the methods reading the solver's own
    environments raise RuntimeError, vev/correlator/overlap/apply_operator
    keep working, and set_hamiltonian starts it over."""
    ic = _field_chain(5.0, gs_method)
    flipped = ic.apply_operator(2 * ic.SxC[0])
    assert flipped.e0 is None and flipped.converged is None
    with pytest.raises(RuntimeError, match="not a ground state"):
        flipped.gs_energy()
    with pytest.raises(RuntimeError, match="not a ground state"):
        flipped.kpm_finite("Sz", 0, "Sz", 0, 4)
    if gs_method == "vumps":
        with pytest.raises(RuntimeError, match="not a ground state"):
            flipped.excitation_energies(0.0)
    else:
        with pytest.raises(RuntimeError, match="not a ground state"):
            flipped.local_excitation_gap()
    twice = flipped.apply_operator(2 * flipped.SxC[0])
    assert twice.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
    assert abs(ic.overlap(twice)) == pytest.approx(1.0, abs=1e-8)
    flipped.set_hamiltonian(-5.0 * flipped.SzC[0])
    assert flipped.eta is None
    assert flipped.gs_energy() == pytest.approx(-2.5, abs=1e-8)


# -- overlap -------------------------------------------------------------------

@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_overlap_exact_values(gs_method):
    """Self-overlap 1; an identity-applied copy (re-gauged by the
    canonicalization) has fidelity of magnitude 1; the up- and
    down-polarized product states are orthogonal."""
    ic = _tfim_chain(1.5, gs_method, 2)
    assert ic.overlap(ic) == pytest.approx(1.0, abs=1e-8)
    same = ic.apply_operator(ic.get_operator("Id", 0))
    assert abs(ic.overlap(same)) == pytest.approx(1.0, abs=1e-6)
    up = _field_chain(10.0, gs_method)
    down = _field_chain(-10.0, gs_method)
    assert abs(up.overlap(down)) == pytest.approx(0.0, abs=1e-6)
    assert abs(up.overlap(up.apply_operator(2 * up.SxC[0]))) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_overlap_distinct_states_below_one(gs_method):
    """Two different TFIM ground states have a fidelity strictly between
    0 and 1, and it decreases as the two fields move apart (g=1.5 against
    g=2 and against g=3). The two solvers' D=3 states are different
    truncations of the exact state (iDMRG's fixed point is not the
    variational optimum), so their fidelities are not compared to each
    other: measured 0.9906 (idmrg) and 0.9953 (vumps) for g=1.5 against
    g=3."""
    a = _tfim_chain(1.5, gs_method, 3)
    near = abs(a.overlap(_tfim_chain(2.0, gs_method, 3)))
    far = abs(a.overlap(_tfim_chain(3.0, gs_method, 3)))
    assert 0.5 < far < near < 0.9999


def test_overlap_rejects_mismatched_chains():
    one = _field_chain(5.0, "vumps")
    two = _field_chain(5.0, "vumps", spins=("1/2", "1/2"))
    with pytest.raises(ValueError, match="unit cells"):
        one.overlap(two)
    other_method = _field_chain(5.0, "idmrg")
    with pytest.raises(ValueError, match="gs_method"):
        one.overlap(other_method)
    with pytest.raises(TypeError):
        one.overlap(1.0)


# -- summps --------------------------------------------------------------------

@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_summps_two_normalized_states_raise(gs_method):
    """Two normalized states tie in per-site norm, so the sum is a cat
    state the single-fixed-point construction cannot hold: RuntimeError,
    propagated from the underlying imps_sum. Every state the public API
    produces is normalized (a solved ground state, and apply_operator's
    output by construction), so there is no public route to the
    well-posed norm-mismatch case; the test below injects one."""
    up = _field_chain(10.0, gs_method)
    down = _field_chain(-10.0, gs_method)
    with pytest.raises(RuntimeError):
        up.summps(down)
    scaled = up.apply_operator(0.9 * up.get_operator("Id", 0))
    assert scaled.eta == pytest.approx(0.81, abs=1e-8)
    with pytest.raises(RuntimeError):
        down.summps(scaled)


@pytest.mark.parametrize("gs_method", GS_METHODS)
def test_summps_norm_mismatch_keeps_larger_branch(gs_method):
    """The well-posed case, reached by rescaling a chain's stored tensors
    by hand (the way the imps_sum examples do): summed with a 0.9x-scaled
    copy of a different TFIM ground state, the unscaled state survives
    alone, with its own <Sz>, <Sz(0)Sz(r)> and bond dimension, and the
    returned chain reports eta 1 and refuses gs_energy."""
    dom = _tfim_chain(1.5, gs_method, 2)
    other = _tfim_chain(2.0, gs_method, 2)
    scaled = copy.copy(other)
    if gs_method == "idmrg":
        from dmrgpy.pyitensor import idmrg
        from dmrgpy.pyitensor.tensor import ITensor
        res = other._result
        scaled._result = idmrg.PeriodicMPS(
            res.sites_uc, res.n_uc,
            [ITensor(T.inds, T.array * 0.9) for T in res.cell_list], eta=0.81)
    else:
        from dmrgpy.pyitensor import vumps
        res = other._vumps_result
        scaled._vumps_result = vumps.UniformMPS(
            res.sites_uc, res.n_uc, res.D, res.d_g, res.AL * 0.9, res.AR * 0.9,
            res.C * 0.9, res.AC * 0.81, eta=0.81)
    summed = dom.summps(scaled)
    assert summed.eta == pytest.approx(1.0, abs=1e-6)
    assert summed.vev("Sz", 0) == pytest.approx(dom.vev("Sz", 0), abs=1e-6)
    for r in (1, 2):
        assert summed.correlator("Sz", 0, "Sz", r) == pytest.approx(
            dom.correlator("Sz", 0, "Sz", r), abs=1e-6)
    with pytest.raises(RuntimeError, match="not a ground state"):
        summed.gs_energy()


# -- backend guards ------------------------------------------------------------

def test_itensor_version3_raises_not_implemented():
    """All three methods refuse itensor_version=3 before touching the C++
    backend, so this holds whether or not mpscpp3 is compiled."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(-5.0 * ic.SzC[0])
    other = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    other.set_hamiltonian(5.0 * other.SzC[0])
    with pytest.raises(NotImplementedError, match="python"):
        ic.apply_operator(2 * ic.SxC[0])
    with pytest.raises(NotImplementedError, match="python"):
        ic.overlap(other)
    with pytest.raises(NotImplementedError, match="python"):
        ic.summps(other)
    assert ic.e0 is None  # no solve was triggered


def test_multisite_vumps_state_raises_not_implemented():
    """A three-site cell runs the sequential multi-site VUMPS solver, whose
    per-site tensors the grouped iMPS algebra does not take."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * 3, itensor_version="python")
    ic.gs_method = "vumps"
    ic.maxm = 1
    ic.set_hamiltonian(-5.0 * (ic.SzC[0] + ic.SzC[1] + ic.SzC[2]))
    ic.gs_energy()
    with pytest.raises(NotImplementedError, match="multi-site"):
        ic.apply_operator(2 * ic.SxC[0])
    with pytest.raises(NotImplementedError, match="multi-site"):
        ic.overlap(ic)
