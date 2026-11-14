"""The speedups of get_kondo_spectrum(mode="ED", submode="KRYLOV") from
docs/kondo_speedup_plan.md (item 11), pinned by their numbers against the
old form kept here and by a count of what they save. The spaces here are
truncated on purpose (small nkry): every test of test_kondo_krylov.py
exhausts its space, and an exhausted space cannot see a change of the
Krylov basis. The last test checks that the Krylov spectra of a real H
take the real route of the sums of conductance.py."""
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
    """Item 11: the Krylov basis is the same to the last bit, the Ritz
    energies too, and the matrices move at roundoff only; in a real basis
    S- is the transpose of S+ and is not projected."""
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
        # the same basis: equal here to the last bit, but BLAS may round by
        # the alignment of the arrays, so the test asks for 1e-14
        assert np.max(np.abs(old.e - new.e)) <= 1e-14*np.max(np.abs(old.e))
        for a, b in zip((old.Sx, old.Sy, old.Sz), (new.Sx, new.Sy, new.Sz)):
            assert np.max(np.abs(a - b)) < 1e-14*np.max(np.abs(a)) # measured 7e-16
        d0, d1 = _dIdV(old), _dIdV(new)
        assert np.max(np.abs(d1 - d0)) < 1e-13*np.max(np.abs(d0)) # measured 6e-16
        # the residual from the recurrence is the explicit one
        assert np.max(np.abs(new.residual - old.residual)) < 1e-12*krylovkondo.onenorm(kks.H) # 9e-17
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


