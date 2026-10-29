"""Real-time evolution under a non-Hermitian Hamiltonian on "python".

pyitensor's TDVP exponentiated every local generator with a Lanczos
recursion, which keeps only the real part of <q|A|q> and discards, in its
full reorthogonalization, the projections a non-Hermitian A has on the
older Krylov vectors. On a 4-site Heisenberg chain with 0.7i*Sz_0 added,
<Sz_2(t)> came out 3.6e-3 off the exact normalized evolution while v3 was
9e-7 off. tdvp_step now decides from the MPO (mpsalgebra.mpo_is_hermitian)
and takes the Arnoldi exponentiator for a non-Hermitian one; see
docs/clik_mps_lindblad_plan.md, Phase 0.

Anchor: a dense expm of the 4-site and 8-site Hamiltonians. DMRG restores
the norm after every step, so the reference is the normalized
<psi|O|psi>/<psi|psi>.
"""

import numpy as np
import pytest
import scipy.linalg as dlg

from dmrgpy import cppext, spinchain, timedependent
from dmrgpy.pyitensor import tdvp as tdvptk
from dmrgpy.pyitensor.mpsalgebra import mpo_is_hermitian


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


_SX = np.array([[0, 1], [1, 0]]) / 2
_SY = np.array([[0, -1j], [1j, 0]]) / 2
_SZ = np.diag([0.5, -0.5])


def _site_op(o, i, n):
    m = np.array([[1.0]])
    for k in range(n):
        m = np.kron(m, o if k == i else np.eye(2))
    return m


def _fields(n):
    """A start Hamiltonian with a nondegenerate product ground state."""
    return [(0.3 + 0.1 * i, 0.4 - 0.07 * i) for i in range(n)]


def _start_h(ops, n):
    """sum_i hz_i Sz_i + hx_i Sx_i, on either side (ops gives Sz, Sx)."""
    return sum(hz * ops[0](i) + hx * ops[1](i)
               for i, (hz, hx) in enumerate(_fields(n)))


def _nh_h(ops, n, gamma, mul):
    """Heisenberg chain plus 0.5 Sx_0 + gamma Sz_0, with gamma complex;
    mul is the operator product on that side (@ dense, * MultiOperator)."""
    Sz, Sx, Sy = ops
    h = sum(mul(Sx(i), Sx(i + 1)) + mul(Sy(i), Sy(i + 1))
            + mul(Sz(i), Sz(i + 1)) for i in range(n - 1))
    return h + 0.5 * Sx(0) + gamma * Sz(0)


def _exact(n, gamma, nt, dt, site):
    ops = (lambda i: _site_op(_SZ, i, n), lambda i: _site_op(_SX, i, n),
           lambda i: _site_op(_SY, i, n))
    psi = np.linalg.eigh(_start_h(ops, n))[1][:, 0]
    U = dlg.expm(-1j * dt * _nh_h(ops, n, gamma, np.matmul))
    O = _site_op(_SZ, site, n)
    out = []
    for _ in range(nt):
        out.append(np.vdot(psi, O @ psi) / np.vdot(psi, psi))
        psi = U @ psi
    return np.array(out)


def _dmrg(version, n, gamma, nt, dt, site, tevol_method="TDVP"):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    sc.maxm = 64
    sc.tevol_method = tevol_method
    ops = (lambda i: sc.Sz[i], lambda i: sc.Sx[i], lambda i: sc.Sy[i])
    sc.set_hamiltonian(_start_h(ops, n))
    wf = sc.get_gs()
    h = _nh_h(ops, n, gamma, lambda a, b: a * b)
    _, c = timedependent.evolve_and_measure(sc, operator=sc.Sz[site], nt=nt,
                                            dt=dt, h=h, wf=wf)
    return np.array(c)


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
@pytest.mark.parametrize("gamma", [0.7j, 0.3 - 0.5j])
def test_nh_tdvp_follows_exact_normalized_evolution(version, gamma):
    """Was 3.6e-3 off on "python" at gamma=0.7i; v3 is 9e-7 off."""
    ex = _exact(4, gamma, 40, 0.05, 2)
    c = _dmrg(version, 4, gamma, 40, 0.05, 2)
    assert np.max(np.abs(c - ex)) < 3e-6


def test_nh_tdvp_python_matches_v3_beyond_one_block():
    """8 sites, so two-site blocks no longer span the whole space."""
    if not cppext.available(3):
        pytest.skip("needs the compiled itensor_version=3 extension")
    ex = _exact(8, 0.7j, 30, 0.05, 3)
    py = _dmrg("python", 8, 0.7j, 30, 0.05, 3)
    v3 = _dmrg(3, 8, 0.7j, 30, 0.05, 3)
    assert np.max(np.abs(py - v3)) < 1e-6
    assert np.max(np.abs(py - ex)) < 1e-5


def test_nh_tdvp_gse_follows_exact_normalized_evolution():
    """The one-site route (TDVP_GSE) takes the same exponentiator."""
    ex = _exact(4, 0.7j, 40, 0.05, 2)
    c = _dmrg("python", 4, 0.7j, 40, 0.05, 2, tevol_method="TDVP_GSE")
    assert np.max(np.abs(c - ex)) < 3e-6


def _mpo(sc, op):
    return sc._session._mpo(op.to_terms())


def test_mpo_is_hermitian_classifies():
    n = 12
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    heis = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
               + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    dmi = sum(sc.Sx[i] * sc.Sy[i + 1] - sc.Sy[i] * sc.Sx[i + 1]
              for i in range(n - 1))
    assert mpo_is_hermitian(_mpo(sc, heis))
    assert mpo_is_hermitian(_mpo(sc, heis + 0.3 * dmi))  # complex, Hermitian
    assert mpo_is_hermitian(_mpo(sc, 1e-9 * heis))       # relative, not absolute
    assert not mpo_is_hermitian(_mpo(sc, heis + 1j * sc.Sz[0]))
    assert not mpo_is_hermitian(_mpo(sc, heis + 1e-8j * sc.Sz[0]))
    assert not mpo_is_hermitian(_mpo(sc, 1e-9 * heis + 1e-15j * sc.Sz[0]))
    assert not mpo_is_hermitian(_mpo(sc, 1j * heis))


def _random_matrix(n, seed, hermitian):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (A + A.conj().T) / 2 if hermitian else A


@pytest.mark.parametrize("niter", [30, 4])
def test_arnoldi_exponentiator_matches_dense_expm(niter):
    """niter=4 forces the sub-stepping branch."""
    A = _random_matrix(40, 1, hermitian=False)
    v0 = np.random.default_rng(2).normal(size=40) + 0j
    for coeff in (-0.05j, 0.03 - 0.04j):
        got = tdvptk._arnoldi_expm_multiply(lambda x: A @ x, v0, coeff,
                                            niter=niter)
        ref = dlg.expm(coeff * A) @ v0
        assert np.linalg.norm(got - ref) < 1e-9 * np.linalg.norm(ref)


def test_arnoldi_agrees_with_lanczos_on_a_hermitian_operator():
    A = _random_matrix(40, 3, hermitian=True)
    v0 = np.random.default_rng(4).normal(size=40) + 0j
    lan = tdvptk._lanczos_expm_multiply(lambda x: A @ x, v0, -0.1j, niter=30)
    arn = tdvptk._arnoldi_expm_multiply(lambda x: A @ x, v0, -0.1j, niter=30)
    assert np.linalg.norm(lan - arn) < 1e-9


def test_hermitian_h_keeps_the_lanczos_route(monkeypatch):
    """hermitian=None on a Hermitian MPO never reaches Arnoldi, so every
    Hermitian result is unchanged by this fix."""
    def refuse(*args, **kwargs):
        raise AssertionError("Arnoldi used for a Hermitian H")
    monkeypatch.setattr(tdvptk, "_arnoldi_expm_multiply", refuse)
    c = _dmrg("python", 4, 0.0, 5, 0.05, 2)
    assert np.max(np.abs(c - _exact(4, 0.0, 5, 0.05, 2))) < 1e-6
