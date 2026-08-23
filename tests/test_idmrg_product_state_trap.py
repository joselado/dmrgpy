"""The growing algorithm must not fall into a product state it cannot leave.

Every particle-number-conserving Hamiltonian has the vacuum (and the filled
state) as an EXACT eigenstate, which makes it an absorbing fixed point of
iDMRG's growth loop: the local solve warm-started there returns it
immediately (a Krylov solve started on an exact eigenvector breaks down on
its first step), its Schmidt rank is 1, so every subsequent truncation has a
single nonzero singular value and nothing can grow the bond dimension back.
The loop then "converges" -- consecutive energy densities are identically
equal, so `etol` trips at the first opportunity -- and reports a product
state with `converged=True` in a fraction of a second.

Reported from a Kondo-chain calculation as a silently wrong ground state at
large |mu| (docs/known_issue_idmrg_product_state_collapse.md); the model
here is the cheapest one that reproduces it deterministically AND has an
exact closed-form answer to check against, so the test does not merely
assert "not the vacuum" but lands on the right number.

The model: a dimerized (SSH) spinless-fermion chain, `t1=1.0` intra-cell and
`t2=0.4` inter-cell, plus a uniform on-site `mu` on every site. Its lower
band is `-sqrt(t1^2 + t2^2 + 2*t1*t2*cos k)`, so with `+mu*N` added the
occupied states are those below zero, and both the filling and the energy
density follow from a one-dimensional integral (`_exact_band`). At `mu=1.25`
the true ground state has `n/site = 0.166092` at `e/site = -0.01648450` --
small enough that the vacuum is close in energy and the random start falls
into it, which before the fix it did on 6 runs out of 6, on BOTH backends.
"""
import numpy as np
import pytest

from dmrgpy import cppext
from dmrgpy import infinitechain

T1, T2, MU = 1.0, 0.4, 1.25


def _exact_band(mu, t1=T1, t2=T2):
    """(e/site, n/site) of the free chain from the band integral."""
    k = np.linspace(-np.pi, np.pi, 200001)
    e_k = -np.sqrt(t1 ** 2 + t2 ** 2 + 2 * t1 * t2 * np.cos(k)) + mu
    occ = e_k < 0
    e_cell = np.trapezoid(np.where(occ, e_k, 0.0), k) / (2 * np.pi)
    n_cell = np.trapezoid(occ.astype(float), k) / (2 * np.pi)
    return e_cell / 2, n_cell / 2


def _trapped_chain(itensor_version, mu=MU, maxm=16, maxiter=60):
    ic = infinitechain.Infinite_Many_Body_Chain([0, 0],
                                                 itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, 1e-11
    cdag = lambda g, i: ic.get_operator("Cdag", i, group=g)
    c = lambda g, i: ic.get_operator("C", i, group=g)
    n = lambda i: ic.get_operator("N", i, group="C")
    h = (T1 * (cdag("C", 0) * c("C", 1) + cdag("C", 1) * c("C", 0))
         + T2 * (cdag("C", 1) * c("R", 0) + cdag("R", 0) * c("C", 1))
         + mu * (n(0) + n(1)))
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


_BACKENDS = ["python"] + ([3] if cppext.available(3) else [])


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_does_not_collapse_into_the_vacuum(itensor_version):
    """The bare assertion the bug violates: the returned state must carry
    particles and must beat the vacuum's own energy (exactly 0 here)."""
    ic = _trapped_chain(itensor_version)
    n0 = float(np.real(ic.vev("N", 0)))
    assert n0 > 0.05, "collapsed into the vacuum (n={})".format(n0)
    assert float(np.real(ic.e0)) < -0.005


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_lands_on_the_exact_free_fermion_answer(itensor_version):
    """Not merely "not the vacuum": the escaped state is the right one.

    Tolerances are set by what `maxm=16` can actually resolve on a partially
    filled (gapless) band, not by what the fix can do: the "python" backend
    lands ~2e-4 from exact here and v3 within ~2e-5, and the filling
    converges more slowly than the energy. They are still two orders of
    magnitude tighter than the failure being guarded against, which returns
    the vacuum's own `e = 0` -- 1.6e-2 away."""
    e_exact, n_exact = _exact_band(MU)
    ic = _trapped_chain(itensor_version)
    assert float(np.real(ic.e0)) == pytest.approx(e_exact, abs=6e-4)
    assert float(np.real(ic.vev("N", 0))) == pytest.approx(n_exact, abs=0.03)


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_a_genuinely_product_ground_state_is_still_exact(itensor_version):
    """The case the noise could plausibly damage: a field-polarized XX chain
    whose true ground state IS a product state. The schedule keys on the
    state being a product state, so this model arms the noise and keeps
    re-arming it -- and must still come back exact once the cap is reached
    and the clean tail has run."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    ic.maxm, ic.maxiter, ic.etol = 4, 50, 1e-12
    ic.set_hamiltonian(-1.0 * (ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0])
                        - 6.0 * ic.SzC[0])
    ic.gs_energy()
    assert ic.vev("Sz", 0) == pytest.approx(0.5, abs=1e-8)
    assert ic.correlator("Sz", 0, "Sz", 3) == pytest.approx(0.25, abs=1e-8)
