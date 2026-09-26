"""The ED Kondo spectrum (kondospectrumtk/conductance.py) against a
brute-force second-order T-matrix of the paper's own model, Ternes, New
J. Phys. 17 063016 (2015): the tunnelling vertex
a^dag_s a_t (sigma.S/2 + U) and the sample exchange a^dag_s a_s sigma.S/2
are built as explicit Jordan-Wigner fermion x impurity operators, and the
second-order amplitude <F|V G0 V|I> is a matvec through every
intermediate Fock state. Which processes exist, their fermion signs and
their energy denominators are therefore read off the operators and H0,
never written by hand.

The fermions are a tip level k, the outgoing sample level k' and one more
sample level q, empty or doubly occupied. Every path through q is an
intermediate state whose q occupation differs from the initial one, with
E_I - E_M = c*eps_q + d; the band integral of 1/(c*eps_q + d) over the
empty (0,w0) or full (-w0,0) half is ln|d|-ln|d-w0| either way, whose
symmetric closed form is the paper's -F0(d). So each path contributes its
numerator times -F0(d_M), with d_M, where the log sits, taken from H0.
The dI/dV is then assembled as the paper assembles it: at T=0 the tip
electron sits at the tip's Fermi edge, eps_k = eV, the outgoing one at
eV - eps_if (t->s) or the incoming one at eV + eps_if (s->t), both
directions summed, the unpolarized electron trace taken twice ("SA
factor", conductance.py's docstring), and J_here = -Jrho_s (the paper
writes its sample exchange as a_s a_s^dag, the opposite sign).

This is what showed that the paper's eq. 25 puts the exchange diagram's
log at the wrong energy: until 2026-09-26 conductance.py used its
F(eV + eps_im) and missed this reference by 0.11 to 1.0 (5-20% of the
third-order term) on the models below; it now agrees to roundoff, and
_paper_eq25 keeps the old form here to show the reference can tell them
apart.
"""
import itertools
from types import SimpleNamespace

import numpy as np
import pytest

from dmrgpy.kondospectrumtk import conductance
from dmrgpy.kondospectrumtk.stepfunctions import F0, Theta0

MUB = 5.7883818060e-5 # eV/T
W0, G0 = 20e-3, 5e-6 # the paper's omega0, Gamma0
JRHO = -0.05
K, KP, Q = 0, 2, 4 # first mode of each spinful level: tip, outgoing, q
NMODES = 6
SIGMA = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], complex)]


def _fermions():
    a = np.array([[0., 1.], [0., 0.]]); z = np.diag([1., -1.]); e = np.eye(2)
    c = []
    for j in range(NMODES):
        out = np.array([[1.]])
        for n in range(NMODES):
            out = np.kron(out, z if n < j else (a if n == j else e))
        c.append(out)
    return c, [x.T.copy() for x in c]


C, CD = _fermions()
NUM = [CD[j] @ C[j] for j in range(NMODES)]
VAC = np.zeros(2**NMODES); VAC[0] = 1.


def _fock(ops):
    v = VAC.copy()
    for o in reversed(ops): v = o @ v
    return v


def _spin(s):
    m = np.arange(s, -s-1, -1); d = len(m)
    Sp = np.zeros((d, d), complex)
    for a in range(1, d): Sp[a-1, a] = np.sqrt(s*(s+1) - m[a]*(m[a]+1))
    return [(Sp + Sp.T)/2, (Sp - Sp.T)/(2j), np.diag(m).astype(complex)]


def _vertices(S, U):
    """tunnelling t->s and the sample exchange, with T = J = 1"""
    d = S[0].shape[0]
    tau = lambda lp, l: (0.5*sum(SIGMA[a][lp, l]*S[a] for a in range(3))
                         + (U*np.eye(d) if lp == l else 0.))
    exch = lambda mp, m: 0.5*sum(SIGMA[a][mp, m]*S[a] for a in range(3))
    Vts = sum(np.kron(CD[s+lp] @ C[K+l], tau(lp, l))
              for s in (KP, Q) for l, lp in itertools.product(range(2), repeat=2))
    Vss = sum(np.kron(CD[sp+mp] @ C[s+m], exch(mp, m))
              for s, sp in itertools.product((KP, Q), repeat=2)
              for m, mp in itertools.product(range(2), repeat=2))
    return Vts, Vts.conj().T, Vss


def _brute_dIdV(e, S, p, U, eVs):
    """second + third order dI/dV at T=0, both directions, in
    2*pi*e^2*T0^2/h units (conductance.py's)"""
    d = len(e)
    Vts, Vst, Vss = _vertices(S, U)
    Eimp = np.kron(np.ones(2**NMODES), e)
    nq = np.kron(NUM[Q] + NUM[Q+1], np.eye(d))
    out = np.zeros(len(eVs))
    for n, eV in enumerate(eVs):
        for i, f, direction in itertools.product(range(d), range(d), ("ts", "st")):
            if p[i] == 0.: continue
            eif = e[f] - e[i]
            th = Theta0(eV - eif) if direction == "ts" else Theta0(-eV - eif)
            if th == 0.: continue
            ekp = eV - eif if direction == "ts" else eV + eif
            H0 = np.kron(eV*(NUM[K] + NUM[K+1]) + ekp*(NUM[KP] + NUM[KP+1]),
                         np.eye(d)) + np.diag(Eimp)
            Vt = Vts if direction == "ts" else Vst
            acc = 0.
            for l, lp in itertools.product(range(2), repeat=2):
                src, dst = (K+l, KP+lp) if direction == "ts" else (KP+l, K+lp)
                A1 = A2 = 0.
                for qops in ([], [CD[Q], CD[Q+1]]):
                    vi = np.kron(_fock([CD[src]] + qops), np.eye(d)[i])
                    vf = np.kron(_fock([CD[dst]] + qops), np.eye(d)[f])
                    EI = np.real(vi @ H0 @ vi)
                    dM = EI - np.real(np.diag(H0))
                    through_q = np.abs(np.real(np.diag(nq)) - np.real(vi @ nq @ vi)) > 0.5
                    w = np.where(through_q, -F0(dM, omega0=W0, Gamma0=G0), 0.)
                    A1 = vf.conj() @ Vt @ vi
                    A2 = A2 + (vf.conj() @ Vss) @ (w*(Vt @ vi)) \
                            + (vf.conj() @ Vt) @ (w*(Vss @ vi))
                acc += abs(A1)**2 + 2*np.real(np.conj(A1)*(-JRHO)*A2)
            out[n] += p[i]*th*2*acc
    return 2*np.pi*out


def _ks(H, S, p=None):
    """a KondoSpectrum stand-in at T=0: eigenenergies from 0, spin
    operators in the eigenbasis, p on the ground state unless given"""
    e, v = np.linalg.eigh(H)
    e = e - e[0]
    Se = [v.conj().T @ x @ v for x in S]
    if p is None:
        p = (e <= 1e-9*e[-1]).astype(float); p = p/p.sum()
    return SimpleNamespace(e=e, p=np.asarray(p, float), dim=len(e), T=0.,
                           kB=8.617333262e-5, Sx=Se[0], Sy=Se[1], Sz=Se[2])


def _code_dIdV(ks, U, eVs):
    return (conductance.second_order_dIdV(ks, eVs, U=U)
            + conductance.third_order_kondo_dIdV(ks, eVs, JRHO, omega0=W0, Gamma0=G0)
            + conductance.third_order_potential_dIdV(ks, eVs, JRHO, U, omega0=W0,
                                                     Gamma0=G0))


def _paper_eq25(ks, U, eVs):
    """_code_dIdV with the exchange log at the paper's F(eV + eps_im)"""
    occ = conductance._occupied_states(ks)
    coeff = np.imag(conductance._triple_product_coefficients(ks, occ))/2.
    eps = ks.e[None, :] - ks.e[occ, None]
    def one(v):
        Th = Theta0(v[:, None, None] - eps[None])
        Fs = F0(v[:, None, None] - eps[None], W0, G0) + F0(v[:, None, None] + eps[None], W0, G0)
        return np.einsum('i,ifm,eif,eim->e', ks.p[occ], coeff, Th, Fs)
    return (conductance.second_order_dIdV(ks, eVs, U=U)
            + 4*np.pi*JRHO*(one(eVs) + one(-eVs))
            + conductance.third_order_potential_dIdV(ks, eVs, JRHO, U, omega0=W0,
                                                     Gamma0=G0))


S12, S1 = _spin(0.5), _spin(1.)
_I2 = np.eye(2)
_A = [np.kron(x, _I2) for x in S12]
_B = [np.kron(_I2, x) for x in S12]
MODELS = {
    # complex eigenvectors: the field is along x
    "S=1/2, Bx=10 T": (2*MUB*10*S12[0], S12),
    "S=1, D=-1 meV, B=(1,0,2) T": (-1e-3*S1[2]@S1[2] + 2*MUB*(S1[0] + 2*S1[2]), S1),
    "S=1/2 dimer, J=1 meV, Bz=5 T, tip on 0":
        (1e-3*sum(_A[a]@_B[a] for a in range(3)) + 2*MUB*5*(_A[2] + _B[2]), _A),
}
EVS = np.linspace(-4e-3, 4e-3, 21) + 1.3e-5 # off the exact step positions


@pytest.mark.parametrize("name", list(MODELS))
def test_ed_spectrum_matches_the_brute_force_t_matrix(name):
    H, S = MODELS[name]
    ks = _ks(H, S)
    U = 0.25 # covers the potential term too
    ref = _brute_dIdV(ks.e, [ks.Sx, ks.Sy, ks.Sz], ks.p, U, EVS)
    got = _code_dIdV(ks, U, EVS)
    assert np.max(np.abs(got - ref)) < 1e-10*np.max(np.abs(ref))
    # the reference tells the exchange diagram's two candidate logs apart
    assert np.max(np.abs(_paper_eq25(ks, U, EVS) - ref)) > 0.05


def test_every_occupied_initial_state_matches_the_brute_force_t_matrix():
    """The formula is per initial state i, and at T=0 only the ground
    state carries weight: give several states weight at once (not a
    Boltzmann distribution, just a test of the i dependence) so an
    exchange log measured from the wrong state would show."""
    H, S = MODELS["S=1, D=-1 meV, B=(1,0,2) T"]
    ks = _ks(H, S, p=[0.5, 0.3, 0.2])
    ref = _brute_dIdV(ks.e, [ks.Sx, ks.Sy, ks.Sz], ks.p, 0.1, EVS)
    got = _code_dIdV(ks, 0.1, EVS)
    assert np.max(np.abs(got - ref)) < 1e-10*np.max(np.abs(ref))


def test_zero_field_spin_half_is_where_the_paper_form_was_already_exact():
    """A free S=1/2 at B=0: every transition is elastic, so eq. 25's
    argument and the T-matrix's coincide, and so do the three."""
    ks = _ks(0*S12[2], S12)
    ref = _brute_dIdV(ks.e, [ks.Sx, ks.Sy, ks.Sz], ks.p, 0.25, EVS)
    assert np.allclose(_code_dIdV(ks, 0.25, EVS), ref, rtol=0., atol=1e-12)
    assert np.allclose(_paper_eq25(ks, 0.25, EVS), ref, rtol=0., atol=1e-12)
