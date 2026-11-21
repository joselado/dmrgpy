"""Findings 14 and 15 of docs/audit_2026_11_21_hole_hunt.md, on the thermal
Krylov route of the Kondo spectrum (kondospectrumtk/krylovkondo.py).

Finding 14: the opt-in tol= stop tested only the Ritz states already inside
the window of open final states, and counted an empty window as converged.
A Ritz value lies above the eigenvalue it approaches, so at the first check
the stop fired while a final state below max|eV| had no Ritz value in the
window yet, and its step was missing from dI/dV (14.3% of the peak on the
chain below). A Ritz state now counts when its energy minus its residual
reaches the window.

Finding 15: with an eigsh start the roundoff of the start in other S^z_tot
sectors kept the recursion going past the exhausted space (301 states where
13 suffice). With S^z_tot conserved, the initial state and the starts are
now put on their sectors.

The chain is the polarized 12-spin ferromagnetic chain in a field B along z.
S|GS> reaches only the 12 one-magnon states, and the lowest of them is the
k=0 magnon S^-_tot|GS>/sqrt(N), exactly B above the ground state, with
M2 = 1/2 |<f|S^-_0|GS>|^2 = 1/(2N) = 1/24: the analytic anchor."""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance, krylovkondo

from test_kondo_krylov import compare, heisenberg

meV = 1e-3
N = 12
J = -1.0*meV
B = 0.1*meV
# max|eV| between B and the lowest Ritz value of the space at the first
# check, B + 0.0112 meV: the window holds the k=0 magnon only
EMAX = 0.105*meV


@pytest.fixture(scope="module")
def ferromagnet():
    """J = -1 meV on every bond and the field -B Sz on every site: the
    ground state is the polarized one, and it is unique."""
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*N)
    h = 0
    for i in range(N-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(N):
        h = h - B*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def test_tol_stop_waits_for_a_ritz_state_whose_interval_reaches_the_window(ferromagnet):
    stop = krylovkondo.KrylovKondoSpectrum(ferromagnet, 0, 0., nkry=100, tol=1e-8)
    assert not stop.su2 and len(stop.weights) == 1
    # The regime of the finding. At the first check (10 blocks of one
    # column) the lowest Ritz value is above EMAX, so no Ritz value is in
    # the window and the old stop fired there, but its interval reaches
    # EMAX. A window that no interval reaches is still a stop, so a window
    # top of -1 eV returns the space of the first check.
    first = stop.state_spectrum(0, -1.)
    assert first.dim == 11
    e, r = first.e[1:], first.residual[1:]
    assert np.min(e) > EMAX # measured 0.11117 meV
    assert np.min(e - r) <= EMAX # measured 0.079 meV
    # with EMAX the stop waits, and the k=0 magnon is in the spectrum
    ss = stop.state_spectrum(0, EMAX)
    assert ss.dim > first.dim
    j = int(np.argmin(np.abs(ss.e - B)))
    assert abs(ss.e[j] - B) < 1e-10*B
    M2 = conductance._spin_matrix_elements_squared(ss)[:, 0]
    assert abs(M2[j] - 1./(2*N)) < 1e-12
    # and it is the only final state of the window
    assert int(np.sum(ss.e[1:] <= EMAX)) == 1


@pytest.mark.parametrize("order", [2, 3])
def test_tol_keeps_the_magnon_step(ferromagnet, order):
    """The public call with tol=1e-8 against the default one, which
    exhausts its space, and at order 2 against the analytic curve: the
    elastic 2 pi |<GS|Sz_0|GS>|^2 = 2 pi/4 at every bias, and the k=0
    magnon step 2 pi/24 above |eV| = B. No bias sits on B, where Theta0 is
    1/2."""
    eV = np.array([-EMAX, -0.1025*meV, -0.05*meV, 0.05*meV, 0.1025*meV, EMAX])
    kw = dict(site=0, Jrho_s=-0.05, T=0., order=order, mode="ED", submode="KRYLOV",
              nkry=100)
    _, a = ferromagnet.get_kondo_spectrum(eV, tol=1e-8, **kw)
    _, b = ferromagnet.get_kondo_spectrum(eV, **kw)
    # 0.143 (order 2) and 0.144 (order 3) of the peak before the fix
    assert np.max(np.abs(a - b)) < 1e-10*np.max(np.abs(b))
    if order == 2:
        exact = 2*np.pi*(0.25 + np.where(np.abs(eV) > B, 1./(2*N), 0.))
        assert np.max(np.abs(a - exact)) < 1e-10*np.max(exact)


def test_space_stops_at_the_one_magnon_states(ferromagnet):
    """The default call (no tol) stops when the one-magnon sector is
    spanned: the initial state and 12 Ritz states, the open-chain magnon
    band B + |J| (1 - cos(pi n/N)), n = 0..N-1, which S^-_0|GS> reaches
    entirely, with the weight sum rule 1/2 ||S^-_0|GS>||^2 = 1/2."""
    kks = krylovkondo.KrylovKondoSpectrum(ferromagnet, 0, 0., nkry=100)
    assert kks._sz_tot is not None
    ss = kks.state_spectrum(0)
    assert ss.dim <= 13 # 301 before the fix
    # all 12 one-magnon states: the Lanczos vectors from S^-_0|GS> are the
    # one-magnon states of sites 0, 1, 2, ..., so none is lost on the way
    assert ss.dim == 13
    band = B + abs(J)*(1. - np.cos(np.pi*np.arange(N)/N))
    assert np.max(np.abs(np.sort(ss.e[1:]) - band)) < 1e-10*abs(J)
    M2 = conductance._spin_matrix_elements_squared(ss)[:, 0]
    assert abs(np.sum(M2[1:]) - 0.5) < 1e-12


def test_sector_projection_moves_dIdV_at_roundoff(ferromagnet):
    """The start as eigsh returns it (the route before the fix, the
    projection switched off on the same object) against the projected one,
    at order 3 over the whole magnon band."""
    kks = krylovkondo.KrylovKondoSpectrum(ferromagnet, 0, 0., nkry=100)
    eV = np.linspace(-3*meV, 3*meV, 41) + 1.3e-5
    new = krylovkondo.krylov_kondo_dIdV(kks, eV, -0.05, order=3)
    kks._sz_tot = None
    old = krylovkondo.krylov_kondo_dIdV(kks, eV, -0.05, order=3)
    assert np.max(np.abs(new - old)) < 1e-12*np.max(np.abs(old)) # measured 1e-15


def dm_chain(n=6):
    """A spin-1/2 chain with uneven exchange, a Dzyaloshinskii-Moriya term
    along z and an uneven field along z: H is complex and conserves
    S^z_tot, so the projection runs on the complex route, with Sx|i> and
    Sy|i> put on the sectors m-1 and m+1."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        Jb = 1e-3*(1.0 + 0.3*np.cos(i))
        h = h + Jb*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
        h = h + 0.3e-3*(sc.Sx[i]*sc.Sy[i+1] - sc.Sy[i]*sc.Sx[i+1])
    for i in range(n):
        h = h + 1e-4*(1. + 0.5*np.cos(i))*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("T", [0., 1.0])
def test_projection_on_the_complex_route_matches_ed(T):
    sc = dm_chain()
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, T, nkry=40, p_cut=1e-14, su2=False)
    assert not kks.real and kks._sz_tot is not None
    sz = kks._sz_tot
    checked = 0
    for n in range(len(kks.weights)):
        psi = kks.V[:, n]/np.linalg.norm(kks.V[:, n])
        m = sz[int(np.argmax(np.abs(psi)))]
        if np.linalg.norm(psi[np.abs(sz - m) > 1e-9]) > krylovkondo.SECTOR_LEAK:
            continue # a mixture of sectors is not projected
        # the space stays in the sectors m-1, m and m+1
        assert kks.state_spectrum(n).dim <= int(np.sum(np.abs(sz - m) < 1. + 1e-9))
        checked += 1
    assert checked >= 1
    assert compare(sc, T) < 1e-11 # exhausted spaces: the full ED sums


def test_no_projection_without_sz_conservation():
    """A field in the xz plane or along y breaks S^z_tot, and nothing is
    projected; an SU(2) chain conserves it."""
    for field in ("xz", "y"):
        assert krylovkondo.KrylovKondoSpectrum(heisenberg(5, field), 1, 0.)._sz_tot is None
    assert krylovkondo.KrylovKondoSpectrum(heisenberg(5), 1, 0.)._sz_tot is not None
