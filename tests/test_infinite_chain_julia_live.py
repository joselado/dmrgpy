"""Native Julia VUMPS on Infinite_Many_Body_Chain(..., itensor_version=
"julia_live") (mpsjulialive/vumps.jl), against the "python" VUMPS it ports
and against exact infinite-chain answers.

The two backends run the same algorithm, so on a model whose VUMPS fixed
point is unique (a gapped TFIM, a dimerized free-fermion chain) they land on
the same state: energy density, onsite values and correlators agree to the
convergence tolerance. The gapless Heisenberg chain is checked on its
energy only, since its state at finite D is one of a degenerate family and
the onsite values of two independent runs need not agree.

What the refusals of every other method look like is pinned in
tests/test_infinite_chain.py (they need no Julia).
"""

import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from _helpers import julia_live_param  # noqa: E402

from dmrgpy import infinitechain  # noqa: E402

BETHE = 0.25 - np.log(2)


def _tfim(itv, D, g=1.5):
    c = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itv)
    c.set_hamiltonian(-4 * c.SzC[0] * c.SzR[0] - 2 * g * c.SxC[0])
    c.maxm, c.maxiter, c.etol = D, 300, 1e-10
    return c


def _ssh(itv, D, t1=1.0, t2=0.4):
    c = infinitechain.Infinite_Many_Body_Chain([0, 0], itensor_version=itv)
    cdag = lambda g, i: c.get_operator("Cdag", i, group=g)
    a = lambda g, i: c.get_operator("C", i, group=g)
    c.set_hamiltonian(t1 * (cdag("C", 0) * a("C", 1) + cdag("C", 1) * a("C", 0))
                      + t2 * (cdag("C", 1) * a("R", 0) + cdag("R", 0) * a("C", 1)))
    c.maxm, c.maxiter, c.etol = D, 300, 1e-10
    return c


def _heisenberg(itv, D, n_uc=2):
    c = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc, itensor_version=itv)
    h = 0
    for i in range(n_uc - 1):
        h = h + c.SxC[i] * c.SxC[i + 1] + c.SyC[i] * c.SyC[i + 1] + c.SzC[i] * c.SzC[i + 1]
    j = n_uc - 1
    h = h + c.SxC[j] * c.SxR[0] + c.SyC[j] * c.SyR[0] + c.SzC[j] * c.SzR[0]
    c.set_hamiltonian(h)
    c.maxm, c.maxiter, c.etol = D, 300, 1e-10
    return c


@pytest.mark.parametrize("itv", [julia_live_param()])
@pytest.mark.parametrize("D", [1, 2, 4])
def test_tfim_matches_python_vumps(itv, D):
    py = _tfim("python", D)
    jl = _tfim(itv, D)
    assert jl.gs_energy() == pytest.approx(py.gs_energy(), abs=1e-8)
    assert jl.converged and py.converged
    assert jl._julia_vumps is not None and py._vumps_result is not None
    assert jl.vev("Sx", 0) == pytest.approx(py.vev("Sx", 0), abs=1e-6)
    # at D=1 the product state breaks the Ising symmetry (mean field orders
    # below g=2), and either sign of <Sz> is the same state up to the flip
    assert abs(jl.vev("Sz", 0)) == pytest.approx(abs(py.vev("Sz", 0)), abs=1e-6)
    for r in range(4):
        assert (jl.correlator("Sz", 0, "Sz", r)
                == pytest.approx(py.correlator("Sz", 0, "Sz", r), abs=1e-6))
        assert (jl.correlator("Sx", 0, "Sx", r)
                == pytest.approx(py.correlator("Sx", 0, "Sx", r), abs=1e-6))


@pytest.mark.parametrize("itv", [julia_live_param()])
def test_tfim_is_reproducible_over_random_starts(itv):
    """Each solve starts from its own random tensors; the converged energy
    does not depend on them."""
    es = [_tfim(itv, 4).gs_energy() for _ in range(3)]
    assert max(es) - min(es) < 1e-10


@pytest.mark.parametrize("itv", [julia_live_param()])
def test_dimerized_fermion_chain_matches_python_and_free_fermions(itv):
    """Two-site cell of spinless fermions with alternating hopping: the
    Jordan-Wigner string inside a cell and across cells enters both the
    Hamiltonian and <Cdag_i C_j>, and the energy density is known exactly,
    -(1/2) <|t1 + t2 e^{ik}|>_k."""
    py = _ssh("python", 4)
    jl = _ssh(itv, 4)
    e_jl = jl.gs_energy()
    assert e_jl == pytest.approx(py.gs_energy(), abs=1e-8)
    k = np.linspace(0, 2 * np.pi, 20001)[:-1]
    exact = -0.5 * np.mean(np.abs(1.0 + 0.4 * np.exp(1j * k)))
    assert e_jl == pytest.approx(exact, abs=1e-5)
    for p in (0, 1):
        assert jl.vev("N", p) == pytest.approx(py.vev("N", p), abs=1e-6)
        for r in range(5):
            assert (jl.correlator("Cdag", p, "C", r)
                    == pytest.approx(py.correlator("Cdag", p, "C", r), abs=1e-6))


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
@pytest.mark.parametrize("itv", [julia_live_param()])
def test_heisenberg_energy_matches_python_and_bethe_ansatz(itv):
    """Gapless, so the gauge mismatch stalls near 1e-4 at small D on both
    backends and `converged` is not compared; the energy is."""
    assert (_heisenberg(itv, 4).gs_energy()
            == pytest.approx(_heisenberg("python", 4).gs_energy(), abs=1e-6))
    assert _heisenberg(itv, 16).gs_energy() == pytest.approx(BETHE, abs=1e-4)


@pytest.mark.parametrize("itv", [julia_live_param()])
def test_matrix_free_branch_matches_dense_branch(itv):
    """Above a bond dimension of 24 (transfer maps) and a dimension of 256
    (H_AC, H_C) the Julia solver runs KrylovKit matrix-free; forcing that
    branch on small problems must give the same state as the dense one."""
    from dmrgpy.mpsjulialive.juliasession import Main
    dense = {}
    for name, mk, ops in (("tfim", lambda: _tfim(itv, 4), [("Sz", 0, "Sz"), ("Sx", 0, "Sx")]),
                          ("ssh", lambda: _ssh(itv, 4), [("Cdag", 0, "C"), ("Cdag", 1, "C")])):
        c = mk()
        dense[name] = (c.gs_energy(), [c.correlator(a, p, b, r) for a, p, b in ops for r in range(4)])
    old = Main.DmrgpyVumps.set_dense_thresholds(0, 0)
    try:
        for name, mk, ops in (("tfim", lambda: _tfim(itv, 4), [("Sz", 0, "Sz"), ("Sx", 0, "Sx")]),
                              ("ssh", lambda: _ssh(itv, 4), [("Cdag", 0, "C"), ("Cdag", 1, "C")])):
            c = mk()
            e = c.gs_energy()
            assert c.converged
            assert e == pytest.approx(dense[name][0], abs=1e-9)
            cors = [c.correlator(a, p, b, r) for a, p, b in ops for r in range(4)]
            assert np.allclose(cors, dense[name][1], atol=1e-6)
    finally:
        Main.DmrgpyVumps.set_dense_thresholds(old[0], old[1])


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
@pytest.mark.parametrize("itv", [julia_live_param()])
def test_matrix_free_branch_on_a_gapless_chain(itv):
    """Gapless Heisenberg on n_uc=2 at D=8, matrix-free everywhere: the
    energy must agree with the dense branch. Then GMRES is starved (Krylov
    dimension 2, one restart), so the environment solves miss and must
    fall back to the dense solve; the energy must still agree. At D=8 the
    state does not converge on any backend and a run lands on one of two
    nearby fixed points (-0.4429138 or -0.4429417, python included), so the
    tolerance is 1e-4 and not the solver precision."""
    from dmrgpy.mpsjulialive.juliasession import Main
    e_dense = _heisenberg(itv, 8).gs_energy()
    old = Main.DmrgpyVumps.set_dense_thresholds(0, 0)
    try:
        assert _heisenberg(itv, 8).gs_energy() == pytest.approx(e_dense, abs=1e-4)
        Main.seval("DmrgpyVumps.LINSOLVE_LIMITS[] = (2, 1)")
        Main.seval("DmrgpyVumps.LINSOLVE_FALLBACKS[] = 0")
        assert _heisenberg(itv, 8).gs_energy() == pytest.approx(e_dense, abs=1e-4)
        assert int(Main.seval("DmrgpyVumps.LINSOLVE_FALLBACKS[]")) > 0
    finally:
        Main.seval("DmrgpyVumps.LINSOLVE_LIMITS[] = (60, 200)")
        Main.DmrgpyVumps.set_dense_thresholds(old[0], old[1])


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
@pytest.mark.parametrize("itv", [julia_live_param()])
def test_redundant_bond_dimension_polarized_chain(itv):
    """The exact state is a product state, so at maxm=4 three directions
    carry no weight and the transfer matrix has a degenerate dominant
    eigenvalue; the bond candidate must take over and every run must land
    on the exact energy -1/2 - 1/8."""
    for _ in range(5):
        c = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itv)
        c.set_hamiltonian(-1.0 * c.SzC[0] - 0.5 * c.SzC[0] * c.SzR[0])
        c.maxm, c.maxiter, c.etol = 4, 300, 1e-10
        assert c.gs_energy() == pytest.approx(-0.625, abs=1e-10)
        assert c.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
