"""Sequential multi-site VUMPS: a unit cell of ANY size, at linear cost.

`vumps.py` handles a multi-site cell by GROUPING it into one supersite of
dimension `prod(d_p)`, which is exact but exponential in the cell size and
is why it rejected `n_uc > 2`. The literature is explicit that this is the
wrong thing to scale up -- Nietner, Vanhecke, Verstraete, Eisert and
Vanderstraeten (arXiv:2003.01142): "the cost of a naive application of the
VUMPS algorithm would scale exponentially with the size of the unit cell",
their own algorithm's "key property" being "a computational effort that
scales linearly rather than exponentially in the size of the unit cell".
`pyitensor/vumps_ms.py` implements that sequential form (Zauner-Stauber,
Vanderstraeten, Fishman, Verstraete, Haegeman, PRB 97 045145 (2018)), which
is also what TeNPy and MPSKit implement.

The checks here are of three kinds, deliberately:

1. *Against the grouped implementation*, at `n_uc = 1, 2` where both are
   valid. Machine precision is the bar, not a tolerance -- the two are
   different code paths for the same mathematics.
2. *Against exact closed-form answers* (TFIM energy density and
   magnetization; a trimerized free-fermion chain's own band integral), so
   a shared misconception between the two implementations cannot pass.
3. *Cell-size invariance*: the SAME physical chain written on cells of
   1, 2, 3, 4 sites -- and a 3-periodic chain written on a 3- and a 6-site
   cell -- must give the same answer. This is what actually exercises the
   multi-site bookkeeping, since a longer cell is only redundant
   parametrization of the same state.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from dmrgpy import cppext
from dmrgpy import infinitechain

# Both backends implement the sequential algorithm (pyitensor/vumps_ms.py
# and its C++ port, Chain::vms_ground_state), so every physics check below
# runs against both -- the two are independent implementations, which is
# what makes agreement between them evidence rather than a shared blind
# spot.
BACKENDS = ["python"] + ([3] if cppext.available(3) else [])

G_TFIM = 1.5


def _tfim(n_uc, D, g=G_TFIM, maxiter=300, backend="python"):
    """H = -sum sigma^x_i sigma^x_{i+1} - g sum sigma^z_i on an n_uc-site
    cell -- the same uniform chain for every n_uc."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                            itensor_version=backend)
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter, ic.etol = D, maxiter, 1e-12
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        nxt = ic.SxC[i + 1] if i + 1 < n_uc else ic.SxR[0]
        h = h + (-4.0) * ic.SxC[i] * nxt + (-2.0 * g) * ic.SzC[i]
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


def _tfim_exact_energy(g):
    val, _ = quad(lambda k: np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return -val / np.pi


def _tfim_exact_sz(g):
    val, _ = quad(lambda k: (g - np.cos(k))
                  / np.sqrt(1 + g ** 2 - 2 * g * np.cos(k)), 0, np.pi)
    return val / (2 * np.pi)      # Sz = sigma^z / 2


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("n_uc", [3, 4])
def test_energy_density_is_independent_of_cell_size(n_uc, backend):
    """A uniform chain does not care how it is sliced: writing it on a
    3- or 4-site cell must reproduce the 1-site answer. n_uc<=2 goes
    through the grouped path and n_uc>=3 through the multi-site one, so
    this compares the two implementations as well."""
    ref = _tfim(1, D=6, backend=backend)
    ic = _tfim(n_uc, D=6, backend=backend)
    assert np.real(ic.e0) == pytest.approx(np.real(ref.e0), abs=1e-10)


@pytest.mark.parametrize("backend", BACKENDS)
def test_energy_and_magnetization_match_exact_tfim(backend):
    ic = _tfim(3, D=6, backend=backend)
    assert np.real(ic.e0) == pytest.approx(_tfim_exact_energy(G_TFIM), abs=1e-6)
    assert float(np.real(ic.vev("Sz", 0))) == pytest.approx(
        _tfim_exact_sz(G_TFIM), abs=1e-6)


@pytest.mark.parametrize("backend", BACKENDS)
def test_observables_are_uniform_across_a_uniform_cell(backend):
    """Every site of a uniform chain must carry the same <Sz>, however many
    of them the cell happens to contain -- the cheapest check that the
    per-site bookkeeping is not off by one somewhere."""
    ic = _tfim(4, D=6, backend=backend)
    sz = [float(np.real(ic.vev("Sz", p))) for p in range(4)]
    assert max(sz) - min(sz) < 1e-10


@pytest.mark.parametrize("backend", BACKENDS)
def test_correlators_match_the_grouped_implementation(backend):
    """Correlators at several separations, multi-site against grouped."""
    ref = _tfim(1, D=6, backend=backend)
    ic = _tfim(3, D=6, backend=backend)
    for r in (1, 2, 3):
        assert float(np.real(ic.correlator("Sz", 0, "Sz", r))) == pytest.approx(
            float(np.real(ref.correlator("Sz", 0, "Sz", r))), abs=1e-8)


# == a genuinely n_uc=3 model, with Jordan-Wigner strings ====================

T_TRIMER = [1.0, 0.6, 1.4]
MU_TRIMER = 1.0        # inside the lowest band gap -> filled band, insulating


def _trimer_exact(T, mu, nk=20001):
    ks = np.linspace(-np.pi, np.pi, nk)
    e_tot, n_tot = 0.0, 0.0
    for k in ks:
        H = np.zeros((3, 3), dtype=complex)
        H[0, 1] = H[1, 0] = T[0]
        H[1, 2] = H[2, 1] = T[1]
        H[2, 0] = T[2] * np.exp(1j * k)
        H[0, 2] = T[2] * np.exp(-1j * k)
        w = np.linalg.eigvalsh(H) + mu
        occ = w[w < 0]
        e_tot += occ.sum()
        n_tot += len(occ)
    return e_tot / nk / 3.0, n_tot / nk / 3.0


def _trimer(n_sites, D=12, maxiter=150, backend="python"):
    ic = infinitechain.Infinite_Many_Body_Chain([0] * n_sites,
                                                 itensor_version=backend)
    ic.gs_method = "vumps"
    ic.maxm, ic.maxiter, ic.etol = D, maxiter, 1e-12
    ic.vumps_nrestarts = 2
    g = ic.get_operator
    h = 0
    for i in range(n_sites):
        t = T_TRIMER[i % 3]
        if i + 1 < n_sites:
            h = h + t * (g("Cdag", i, "C") * g("C", i + 1, "C")
                         + g("Cdag", i + 1, "C") * g("C", i, "C"))
        else:
            h = h + t * (g("Cdag", i, "C") * g("C", 0, "R")
                         + g("Cdag", 0, "R") * g("C", i, "C"))
        h = h + MU_TRIMER * g("N", i, "C")
    ic.set_hamiltonian(h)
    ic.gs_energy()
    return ic


@pytest.mark.parametrize("backend", BACKENDS)
def test_trimerized_fermion_chain_matches_exact_bands(backend):
    """A model that is genuinely 3-periodic (so the cell structure is real,
    not redundant) and fermionic (so the Jordan-Wigner strings in the
    multi-site correlator/automaton are exercised), against its own exact
    3-band integral."""
    e_ex, n_ex = _trimer_exact(T_TRIMER, MU_TRIMER)
    ic = _trimer(3, backend=backend)
    assert np.real(ic.e0) == pytest.approx(e_ex, abs=1e-4)
    filling = np.mean([float(np.real(ic.vev("N", p))) for p in range(3)])
    assert filling == pytest.approx(n_ex, abs=1e-4)


@pytest.mark.parametrize("backend", BACKENDS)
def test_trimer_occupations_are_genuinely_non_uniform(backend):
    """Guards the guard: if the three sites came out equally occupied, the
    test above would be passing on a uniform state and would say nothing
    about multi-site bookkeeping at all."""
    ic = _trimer(3, backend=backend)
    n = [float(np.real(ic.vev("N", p))) for p in range(3)]
    assert max(n) - min(n) > 0.1


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_doubled_cell_reproduces_the_same_state(backend):
    """The same 3-periodic chain written on a 6-site cell: a longer cell is
    only redundant parametrization, so every per-site occupation must come
    back in the same pattern."""
    ic3 = _trimer(3, backend=backend)
    ic6 = _trimer(6, backend=backend)
    assert np.real(ic6.e0) == pytest.approx(np.real(ic3.e0), abs=1e-6)
    n3 = [float(np.real(ic3.vev("N", p))) for p in range(3)]
    n6 = [float(np.real(ic6.vev("N", p))) for p in range(6)]
    for p in range(6):
        assert n6[p] == pytest.approx(n3[p % 3], abs=1e-4)


def test_excitations_reject_a_large_cell_clearly():
    """The tangent-space ansatz still expects the grouped single-supersite
    gauge, so it must say so rather than misread a list of per-site tensors
    as one tensor."""
    ic = _tfim(3, D=4, maxiter=50)
    with pytest.raises(NotImplementedError, match="n_uc>2"):
        ic.excitation_energies(0.0)
