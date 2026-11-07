"""The charge penalty of a conserved sector on itensor_version="python",
and the stop scale of every variational route that carries it.

On "python" a sector is a penalty lam*sum_k (Q_k-q_k)^2 added to the
operator the variational solves minimize (pyitensor/sector.py). Each local
Davidson solve stops once its Ritz residual is below a tolerance times the
scale of the operator, max|diag(H_eff)|. The penalty is large on every
configuration of a bond that is off the sector, so with it in that
diagonal the scale is that of the penalty, the solves stop after two or
three matvecs and the state keeps an admixture of the next state of the
sector. dmrg(penalty=) takes the penalty's diagonal out of the scale for
the ground-state route (pinned in test_sector_conservation_python.py);
the tests here pin the same for the excited-state route (dmrg_excited),
the generalized route (dmrg_generalized) and the one-site CBE route
(cbe.dmrg_cbe), and check the energies of each against ED restricted to
the sector.

The chain is the Heisenberg chain written as Sx Sx + Sy Sy + Sz Sz. At
Sz=0 its sector holds the global ground state, so the plain and the
sector solve see the same H and their median stop scales are comparable.
"""
import sys

import numpy as np
import pytest
import scipy.linalg

from dmrgpy import spinchain
from dmrgpy.multioperator import MO2matrix
from dmrgpy.pyitensor import cbe as cb
from dmrgpy.pyitensor import dmrg as dm

# Measured on d22: the three lowest sector energies of the 6-site chain
# reach sector ED to 4e-15 (run AB_after1; 1.2e-12 before the stop scale
# was fixed). The generalized eigenvalue of the 8-site chain in a sector
# was the quotient of the compressed penalized MPO, 3.8e-12 (Sz=0) and
# 4.7e-12 (Sz=2) off sector ED (run AB_after4); it is now the quotient
# under the plain H, 1.5e-14 (Sz=0) and 7.5e-15 (Sz=2) off (run F_fix).
EXCITED_TOL = 1e-8
GENERALIZED_TOL = 1e-13

# Plain (no sector) energies of the 10-site chain at maxm=40, nsweeps=8,
# with the global random stream seeded with 7 before each chain, recorded
# on d22 in a fresh process with the code before the penalty reached the
# excited, generalized and CBE routes (run AB_after4, where the new code
# gave the same numbers bit for bit). Outside a sector no penalty is
# passed. Inside pytest the earlier tests change how much of the random
# stream a solve draws, and another machine has another BLAS, so the
# comparison is at a tolerance well above that noise (1e-13 for the
# ground state, 7e-11 for the degenerate excited pair without a seed).
PLAIN_GS10 = -4.258035207277334
PLAIN_EX10 = [-4.258035207277337, -3.9306735895015885, -3.9306735895015565]
PLAIN_GEN10 = -6.283528222462229
PLAIN_CBE10 = -4.258035207277249
PLAIN_TOL = 1e-9


def heisenberg_chain(n, maxm=40, nsweeps=8):
    sc = spinchain.Spin_Chain([2] * n)
    sc.setup_python()
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    return sc, h, sum(sc.Sz)


def metric(sc, n):
    """A Hermitian positive-definite metric that conserves Sz: its
    spectrum lies within 1 +- 0.05*(n-1)."""
    a = 0
    for i in range(n - 1):
        a = a + 0.2 * sc.Sz[i] * sc.Sz[i+1]
    return 1 + a


def sector_block(sc, op, charge, target):
    ed = sc.get_ED_obj()
    M = np.array(MO2matrix(op, ed).todense())
    q = np.diag(np.array(MO2matrix(charge, ed).todense())).real
    sel = np.abs(q - target) < 1e-9
    return M[np.ix_(sel, sel)]


def spy_on(module, callers, record):
    """Wrap module._davidson_ground_state, recording the stop scale of each
    solve made from one of `callers` (the direct caller's function name)
    and the number of matvecs it took."""
    davidson = module._davidson_ground_state
    def spy(matvec, x0, diag, niter, rtol, hdiag=None, vtol=None):
        if sys._getframe(1).f_code.co_name not in callers:
            return davidson(matvec, x0, diag, niter, rtol, hdiag=hdiag, vtol=vtol)
        d = diag if hdiag is None else hdiag
        count = [0]
        def counted(v):
            count[0] += 1
            return matvec(v)
        out = davidson(counted, x0, diag, niter, rtol, hdiag=hdiag, vtol=vtol)
        record.append((float(np.max(np.abs(np.asarray(d)))), count[0]))
        return out
    return spy


def median_scale(record):
    return float(np.median([s for s, _ in record]))


def test_excited_states_stop_on_the_scale_of_h_not_of_the_penalty(monkeypatch):
    record = []
    monkeypatch.setattr(dm, "_davidson_ground_state",
                        spy_on(dm, ("_local_ground_state_penalized",), record))
    out = {}
    for sector in (False, True):
        sc, h, sz = heisenberg_chain(8)
        sc.set_hamiltonian(h)
        if sector:
            sc.set_conserved_sector(Sz=0)
        sc.gs_energy()
        record.clear()
        es, _ = sc.get_excited_states(n=2)
        out[sector] = (sorted(np.real(es)), median_scale(record))
    assert out[True][0][1] == pytest.approx(out[False][0][1], abs=1e-6)
    assert out[True][1] == pytest.approx(out[False][1], rel=0.5)


def test_generalized_solve_stops_on_the_scale_of_h_not_of_the_penalty(monkeypatch):
    record = []
    monkeypatch.setattr(dm, "_davidson_ground_state",
                        spy_on(dm, ("_local_ground_state",), record))
    out = {}
    for sector in (False, True):
        sc, h, sz = heisenberg_chain(8)
        sc.set_hamiltonian(h)
        if sector:
            sc.set_conserved_sector(Sz=0)
        record.clear()
        lam = sc.gs_energy_generalized(metric(sc, 8))
        out[sector] = (float(np.real(lam)), median_scale(record))
    assert out[True][0] == pytest.approx(out[False][0], abs=1e-6)
    assert out[True][1] == pytest.approx(out[False][1], rel=0.5)


def test_cbe_solve_stops_on_the_scale_of_h_not_of_the_penalty(monkeypatch):
    # cbe.py binds _davidson_ground_state by name at import, so it is
    # patched there and not in dmrg.py
    record = []
    monkeypatch.setattr(cb, "_davidson_ground_state",
                        spy_on(cb, ("_local_ground_state_one_site",), record))
    out = {}
    for sector in (False, True):
        sc, h, sz = heisenberg_chain(10)
        sc.dmrg_method = "CBE"
        sc.set_hamiltonian(h)
        if sector:
            sc.set_conserved_sector(Sz=0)
        record.clear()
        e = sc.gs_energy()
        out[sector] = (e, median_scale(record))
    assert out[True][0] == pytest.approx(out[False][0], abs=1e-6)
    assert out[True][1] == pytest.approx(out[False][1], rel=0.5)


@pytest.mark.parametrize("two_sz", [0, 2])
def test_three_lowest_sector_energies_match_sector_ed(two_sz):
    n = 6
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    es, _ = sc.get_excited_states(n=3)
    exact = np.linalg.eigvalsh(sector_block(sc, h, sz, two_sz / 2.0))[:3]
    assert sorted(np.real(es)) == pytest.approx(list(exact), abs=EXCITED_TOL)


@pytest.mark.parametrize("two_sz", [0, 2])
def test_generalized_eigenvalue_in_a_sector_matches_sector_ed(two_sz):
    """In a sector the reported lambda is <H>/<A> of the converged state
    under the unpenalized H, so it reaches the generalized eigenvalue of
    the sector block to the accuracy of the state."""
    n = 8
    sc, h, sz = heisenberg_chain(n)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=two_sz)
    A = metric(sc, n)
    lam = sc.gs_energy_generalized(A)
    Hb = sector_block(sc, h, sz, two_sz / 2.0)
    Ab = sector_block(sc, A, sz, two_sz / 2.0)
    exact = scipy.linalg.eigh(Hb, Ab, eigvals_only=True)[0]
    assert float(np.real(lam)) == pytest.approx(exact, abs=GENERALIZED_TOL)


def test_plain_solves_are_unchanged():
    """Outside a sector no penalty is passed anywhere, so the ground-state,
    excited-state, generalized and CBE solves of a chain without a sector
    give the numbers of the code before (see PLAIN_GS10)."""
    def seeded():
        np.random.seed(7)
        return heisenberg_chain(10)
    sc, h, sz = seeded()
    sc.set_hamiltonian(h)
    assert sc.gs_energy() == pytest.approx(PLAIN_GS10, abs=PLAIN_TOL)
    es, _ = sc.get_excited_states(n=3)
    assert sorted(np.real(es)) == pytest.approx(PLAIN_EX10, abs=PLAIN_TOL)
    sc, h, sz = seeded()
    sc.set_hamiltonian(h)
    lam = sc.gs_energy_generalized(metric(sc, 10))
    assert float(np.real(lam)) == pytest.approx(PLAIN_GEN10, abs=PLAIN_TOL)
    sc, h, sz = seeded()
    sc.dmrg_method = "CBE"
    sc.set_hamiltonian(h)
    assert sc.gs_energy() == pytest.approx(PLAIN_CBE10, abs=PLAIN_TOL)


def test_a_new_hamiltonian_sets_the_default_penalty_strength():
    """The default strength is terms_scale of the Hamiltonian being solved,
    not of the first one sent to the chain; an explicit
    set_sector_penalty() strength is kept across a new Hamiltonian."""
    from dmrgpy.pyitensor import sector as sectortk
    sc, h, sz = heisenberg_chain(6)
    sc.set_hamiltonian(h)
    sc.set_conserved_sector(Sz=0)
    e1 = sc.gs_energy()
    ses = sc._session
    p1, s1 = ses._penalty_mpo(), ses._penalty_strength
    assert s1 == sectortk.terms_scale(ses._h_terms)
    sc.set_hamiltonian(2.0 * h)
    assert sc.gs_energy() == pytest.approx(2.0 * e1, abs=1e-8)
    assert ses._penalty_mpo() is not p1
    assert ses._penalty_strength == pytest.approx(2.0 * s1)
    assert ses._penalty_strength == sectortk.terms_scale(ses._h_terms)
    sc.set_sector_penalty(100.0)
    sc.set_hamiltonian(h)
    sc.gs_energy()
    assert ses._penalty_strength == 100.0
