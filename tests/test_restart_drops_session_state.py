"""set_hamiltonian() (restart=True) starts the next solve from scratch on
every backend.

The v2/v3 session kept its own MPS across set_hamiltonian() and swept from
it, while Many_Body_Chain.restart() promised a cold recalculation. When
every term of H conserves total Sz, a two-site sweep, noise included,
cannot leave the Sz sector of its start, so a field sweep on one chain
object stayed in the sector of an earlier field: on a Kondo chain c S c S c
on v3, 5 of 12 sweeps came back 3.3e-2 too high at B=0.58 with the wrong
total Sz. "python" already dropped its state when the terms changed (2026-09
audit, finding 2, tests/test_audit_2026_09_pyitensor-evolution.py).
restart=False keeps the warm start on purpose and is not tested here.
"""

import pytest

from dmrgpy import spinchain, mixedchain
from dmrgpy import cppext

BACKENDS = [b for b in (2, 3, "python")
            if b == "python" or cppext.available(b)]


def _heisenberg(sc, n, B):
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i + 1] + sc.Sy[i]*sc.Sy[i + 1] \
              + sc.Sz[i]*sc.Sz[i + 1]
    return h + B*sum(sc.Sz[i] for i in range(n))


@pytest.mark.parametrize("v", BACKENDS)
def test_a_polarized_state_does_not_survive_a_new_hamiltonian(v):
    """The fully polarized state is an exact eigenstate of the Heisenberg
    chain at every field, and the only state of its Sz sector. On "python"
    the B=0 point used to come back as the ferromagnetic maximum (finding 2);
    on v2/v3 ITensor's Davidson escapes an exact eigenstate by randomizing a
    Krylov direction, so this case passed there even before."""
    n = 6
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm, sc.nsweeps = 40, 20
    for B in (3.0, 0.0):
        sc.set_hamiltonian(_heisenberg(sc, n, B))
        e = sc.gs_energy()
        assert e == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-6), \
            "B=%g on a reused chain" % B


@pytest.mark.parametrize("v", [b for b in BACKENDS if b != 2])
def test_a_field_sweep_on_a_kondo_chain_follows_the_ground_state(v):
    """The repro of the mixed chain: two spins 1/2 between three Hubbard
    sites, every term conserving total Sz, with a level crossing in total
    Sz between B=0.46 and B=0.58. On v3 a reused chain stayed in the old
    sector at B=0.58 in about 40% of the passes, so the field goes back and
    forth eight times (at 40% the old code passes all of them with
    probability 2e-2)."""
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F", "1/2", "F"],
                                             itensor_version=v)
    mc.maxm, mc.nsweeps = 60, 20
    K, t, U = 0.6, 0.5, 2.0
    h0 = 0
    for f in (0, 2, 4):
        h0 = h0 + U*mc.Nup[f]*mc.Ndn[f] - 0.5*U*mc.Ntot[f]
    for s in (1, 3):
        h0 = h0 + K*(mc.SS(s - 1, s) + mc.SS(s, s + 1))
    for f in (0, 2):
        hop = t*(mc.Cdagup[f]*mc.Cup[f + 2] + mc.Cdagdn[f]*mc.Cdn[f + 2])
        h0 = h0 + hop + hop.get_dagger()
    eed = {}
    for B in (0.46, 0.58):
        mc.set_hamiltonian(h0 - B*sum(mc.Sz))
        eed[B] = mc.gs_energy(mode="ED")
    for k in range(8):
        for B in (0.46, 0.58):
            mc.set_hamiltonian(h0 - B*sum(mc.Sz))
            e = mc.gs_energy()
            assert e == pytest.approx(eed[B], abs=1e-6), \
                "pass %d, B=%g on a reused chain" % (k, B)
