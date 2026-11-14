"""The speedups of get_kondo_spectrum(mode="ED", submode="KRYLOV") from
docs/kondo_speedup_plan.md (item 11, the block Gram-Schmidt and the opt-in
residual stop of item 12), pinned by their numbers against the old form
kept here and by a count of what they save. The spaces here are truncated
on purpose (small nkry): every test of test_kondo_krylov.py exhausts its
space, and an exhausted space cannot see a change of the Krylov basis. One
test checks that the Krylov spectra of a real H take the real route of the
sums of conductance.py."""
import numpy as np
import pytest

from dmrgpy.kondospectrumtk import conductance, krylovkondo

from test_kondo_krylov import heisenberg

EV = np.linspace(-4e-3, 4e-3, 21) + 1.3e-5


def _old_state_spectrum(kks, k):
    """state_spectrum as of 2026-11-12: the Ritz vectors F = Q y formed,
    H@Q formed again, the residual ||H F - F theta|| explicit and the three
    matrices built on [psi, F] one by one."""
    psi = kks.V[:, k]
    psi = psi/np.linalg.norm(psi)
    H = kks.H
    Ei = float(np.real(psi.conj() @ (H @ psi)))
    if kks.real:
        Sp, Sm, Sz = kks._ladder
        starts = np.stack([Sp @ psi, Sm @ psi, Sz @ psi], axis=1)
    else:
        starts = np.stack([M @ psi for M in kks._S], axis=1)
    blocks = []

    def orth(V):
        scale = float(np.max(np.linalg.norm(V, axis=0))) if V.size else 0.
        if scale == 0.:
            return V[:, :0]
        for _ in range(2):
            V = V - np.outer(psi, psi.conj() @ V)
            for Q in blocks:
                V = V - Q @ (Q.conj().T @ V)
        U, s, _ = np.linalg.svd(V, full_matrices=False)
        return U[:, s > krylovkondo.DROP_TOL*scale]

    B = orth(starts)
    n = 0
    while B.shape[1] > 0 and n + B.shape[1] <= 3*kks.nkry:
        blocks.append(B)
        n += B.shape[1]
        B = orth(H @ B)
    Q = np.concatenate(blocks, axis=1)
    HQ = H @ Q
    Hk = Q.conj().T @ HQ
    theta, y = np.linalg.eigh((Hk + Hk.conj().T)/2.)
    F = Q @ y
    residual = np.linalg.norm(HQ @ y - F*theta[None, :], axis=0)
    basis = np.concatenate([psi[:, None], F], axis=1)
    if kks.real:
        Ap, Am, Az = (basis.T @ (M @ basis) for M in kks._ladder)
        mats = ((Ap + Am)/2., (Ap - Am)/2j, Az.astype(complex))
    else:
        mats = tuple(basis.conj().T @ (M @ basis) for M in kks._S)
    e = np.concatenate([[0.], theta - Ei])
    e[np.abs(e) <= kks.tol_degenerate] = 0.
    res0 = float(np.linalg.norm(H @ psi - Ei*psi))
    return krylovkondo.StateSpectrum(e, mats, kks.T, kks.kB,
                                     np.concatenate([[res0], residual]))


def _dIdV(ks):
    from dmrgpy.kondospectrumtk.stepfunctions import cached_fbuilder
    Fb = cached_fbuilder(ks.T, omega0=20e-3, Gamma0=5e-5, kB=ks.kB) if ks.T > 0 else None
    return (conductance.second_order_dIdV(ks, EV, U=0.3)
            + conductance.third_order_kondo_dIdV(ks, EV, -0.08, omega0=20e-3,
                                                 Gamma0=5e-5, Fb=Fb)
            + conductance.third_order_potential_dIdV(ks, EV, -0.08, 0.3, omega0=20e-3,
                                                     Gamma0=5e-5, Fb=Fb))


@pytest.mark.parametrize("field", ["xz", "y"])
def test_state_spectrum_keeps_the_numbers_on_a_truncated_space(field, monkeypatch):
    """Item 11 and the block Gram-Schmidt of block_krylov: the Ritz
    energies and dI/dV are the same to roundoff, and the matrices move by
    the roundoff of the basis, which the block Gram-Schmidt changes (3.6e-13
    of their largest element on this chain, against 7e-16 before it); in a
    real basis S- is the transpose of S+ and is not projected."""
    sc = heisenberg(6, field)
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, 1.0, nkry=6, p_cut=1e-6)
    assert kks.real == (field == "xz")
    calls = []
    project = krylovkondo._project
    monkeypatch.setattr(krylovkondo, "_project",
                        lambda *a: calls.append(1) or project(*a))
    for k in range(min(3, len(kks.weights))):
        old = _old_state_spectrum(kks, k)
        new = kks.state_spectrum(k)
        assert old.dim == new.dim == 19 # truncated: 63 states are reachable
        # the same space, so the Ritz energies agree to roundoff (8e-16)
        assert np.max(np.abs(old.e - new.e)) <= 1e-14*np.max(np.abs(old.e))
        for a, b in zip((old.Sx, old.Sy, old.Sz), (new.Sx, new.Sy, new.Sz)):
            # measured 3.6e-13 (field xz) and 2.9e-14 (y); two eigsh starts
            # of the same code differ by 9e-12 on 12 spins
            assert np.max(np.abs(a - b)) < 1e-12*np.max(np.abs(a))
        d0, d1 = _dIdV(old), _dIdV(new)
        assert np.max(np.abs(d1 - d0)) < 1e-13*np.max(np.abs(d0)) # measured 1.2e-15
        # the residual from the recurrence is the explicit one
        assert np.max(np.abs(new.residual - old.residual)) < 1e-12*krylovkondo.onenorm(kks.H) # 9e-15
    assert len(calls) == (2 if kks.real else 3)*min(3, len(kks.weights))


def test_the_krylov_spectrum_of_a_real_hamiltonian_is_of_the_real_form():
    """The third-order sums take conductance._spin_triplet's real route
    (Sx and Sz real, Sy i times a real matrix) for a real H, also in the
    SU(2) reduction, and the complex route for a field along y."""
    for sc, su2, real in ((heisenberg(5, "xz"), False, True),
                          (heisenberg(5), True, True),
                          (heisenberg(5, "y"), False, False)):
        kks = krylovkondo.KrylovKondoSpectrum(sc, 1, 1.0, nkry=6, su2=su2)
        for k in range(min(2, len(kks.weights))):
            assert conductance._spin_triplet(kks.state_spectrum(k))[1] == real




def test_residual_stop_is_the_fixed_size_it_stopped_at(monkeypatch):
    """Item 12, opt-in: in a narrow window (0.1 K, |eV| below 0.5 meV) the
    Ritz states that can be open final states converge long before the
    space is exhausted, the space stops there, and the spectrum is the one
    of a fixed space of that size."""
    monkeypatch.setattr(krylovkondo, "CHECK_EVERY", 2)
    sc = heisenberg(8, "xz")
    emax = 5e-4
    stop = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.1, nkry=60, p_cut=1e-6, tol=1e-8)
    full = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.1, nkry=60, p_cut=1e-6)
    assert full.tol is None # opt-in
    for k in range(len(stop.weights)):
        a = stop.state_spectrum(k, emax)
        b = full.state_spectrum(k, emax)
        assert b.dim == 181 # 3*nkry states and the initial one
        assert a.dim < b.dim # the stop fired
        assert a.dim % 3 == 1 # blocks of three columns, plus the initial state
        fixed = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.1, nkry=(a.dim - 1)//3,
                                                p_cut=1e-6).state_spectrum(k)
        # the same arithmetic (equal to the last bit here, 1e-14 for BLAS)
        assert np.max(np.abs(a.e - fixed.e)) <= 1e-14*np.max(np.abs(a.e))
        for x, y in zip((a.Sx, a.Sy, a.Sz), (fixed.Sx, fixed.Sy, fixed.Sz)):
            assert np.max(np.abs(x - y)) <= 1e-14*np.max(np.abs(x))
        # every Ritz state of the window is converged
        win = a.e[1:] <= stop.window(emax)
        assert np.all(a.residual[1:][win] < 1e-8*stop.scale)


def test_residual_stop_that_never_fires_changes_nothing(monkeypatch):
    """A tol no Ritz state reaches runs every check and never stops, and
    the spectrum is the one of the default call without checks."""
    monkeypatch.setattr(krylovkondo, "CHECK_EVERY", 2)
    sc = heisenberg(7, "y")
    never = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.5, nkry=30, p_cut=1e-6, tol=1e-300)
    plain = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.5, nkry=30, p_cut=1e-6)
    for k in range(min(3, len(plain.weights))):
        a, b = never.state_spectrum(k, 1e-3), plain.state_spectrum(k)
        assert a.dim == b.dim == 91
        assert np.max(np.abs(a.e - b.e)) <= 1e-14*np.max(np.abs(b.e))
        for x, y in zip((a.Sx, a.Sy, a.Sz), (b.Sx, b.Sy, b.Sz)):
            assert np.max(np.abs(x - y)) <= 1e-14*np.max(np.abs(y))


def test_residual_stop_moves_dIdV_by_the_tolerance():
    """The public call with tol=1e-8 against the default (no stop, 181
    states here), in narrow windows where the stop fires, at 0.1 K and at
    T=0, and in a wide one, where it does not, since every Ritz state there
    is a possible final state."""
    sc = heisenberg(8, "xz")
    for T, emax in ((0.1, 5e-4), (0., 3e-4), (1.0, 4e-3)):
        eV = np.linspace(-emax, emax, 21) + 1.3e-5
        kw = dict(site=1, Jrho_s=-0.08, U=0.3, T=T, omega0=20e-3, Gamma0=5e-5,
                  mode="ED", submode="KRYLOV", nkry=60, p_cut=1e-6)
        _, a = sc.get_kondo_spectrum(eV, tol=1e-8, **kw)
        _, b = sc.get_kondo_spectrum(eV, **kw)
        # measured 4e-12 (121 of 181 states), 1e-11 (31) and 0 (no stop)
        assert np.max(np.abs(a - b)) < 1e-8*np.max(np.abs(b))


def test_residual_stop_keywords():
    sc = heisenberg(4)
    for bad in (0., 1., -1e-8, True):
        with pytest.raises(ValueError, match="tol must be"):
            sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, T=0., submode="KRYLOV",
                                  tol=bad)
    with pytest.raises(TypeError, match=r"argument\(s\): tol\. "):
        sc.get_kondo_spectrum(EV, site=0, Jrho_s=-0.05, T=0., tol=1e-8)
    # opt-in: no stop unless tol is given
    import inspect
    from dmrgpy.kondospectrumtk.kondomethod import KondoSpectrumMixin
    sig = inspect.signature(KondoSpectrumMixin._get_kondo_spectrum_krylov)
    assert sig.parameters["tol"].default is None
    assert inspect.signature(krylovkondo.KrylovKondoSpectrum).parameters["tol"].default is None


class _CountingF():
    """An F that counts the entries it evaluates, keeping odd_part so that
    conductance._F_is_even reads it as the F it wraps."""
    def __init__(self, F):
        self.F, self.count = F, 0
        if hasattr(F, "odd_part"):
            self.odd_part = F.odd_part

    def __call__(self, x):
        self.count += np.size(x)
        return self.F(x)


@pytest.mark.parametrize("field", ["xz", "y"])
@pytest.mark.parametrize("T", [0., 1.0])
def test_transposed_exchange_table_keeps_the_numbers(field, T, monkeypatch):
    """Item 2 on the Krylov route: the exchange table of -eV on final x
    final states is the transpose of the one at +eV, so each bias needs
    n_f^2 fewer F values, and dI/dV moves at roundoff only."""
    from dmrgpy.kondospectrumtk.stepfunctions import cached_fbuilder
    sc = heisenberg(6, field)
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, T, nkry=6, p_cut=1e-6)
    ks = kks.state_spectrum(0)
    Fb = cached_fbuilder(T, omega0=20e-3, Gamma0=5e-5, kB=ks.kB) if T > 0 else None
    get_F = conductance._get_F
    counters = []
    monkeypatch.setattr(conductance, "_get_F",
                        lambda *a: counters.append(_CountingF(get_F(*a))) or counters[-1])
    d = [conductance.third_order_kondo_dIdV(ks, EV, -0.08, omega0=20e-3, Gamma0=5e-5,
                                            Fb=Fb, transpose_exchange=flag)
         for flag in (False, True)]
    nf = len(conductance._final_states(ks, conductance._occupied_states(ks), EV))
    assert counters[0].count - counters[1].count == nf*nf*len(EV)
    # measured 0 at T=0 (F0 is exactly even) and 2.3e-16 at 1 K
    assert np.max(np.abs(d[1] - d[0])) < 1e-13*np.max(np.abs(d[0]))


def test_only_the_krylov_route_transposes_the_exchange_table(monkeypatch):
    """The full ED route keeps its sum (n_f is far below dim there after
    the cut of _final_states, and item 2 gave 0 to 4%), the Krylov route
    takes the transposed table."""
    calls = []
    plain = conductance._kondo_one_block_transposed
    monkeypatch.setattr(conductance, "_kondo_one_block_transposed",
                        lambda *a: calls.append(1) or plain(*a))
    sc = heisenberg(5, "xz")
    kw = dict(site=1, Jrho_s=-0.08, U=0.3, T=1.0, omega0=20e-3, Gamma0=5e-5)
    sc.get_kondo_spectrum(EV, mode="ED", **kw)
    assert not calls
    sc.get_kondo_spectrum(EV, mode="ED", submode="KRYLOV", nkry=6, **kw)
    assert calls
