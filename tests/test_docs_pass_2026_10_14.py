"""Regressions for the bugs found while rewriting the user guide
(2026-10-14), one test per bug. Small systems, every one checked against
ED or an exact value."""

import contextlib
import io
import warnings

import numpy as np
import pytest

from dmrgpy import (cppext, fermionchain, fidelity, infinitechain, meanfield,
                    parafermionchain, spinchain)


def _quiet(f, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return f(*args, **kwargs)


def _cpp_versions():
    return [v for v in (2, 3) if cppext.available(v)]


def _heisenberg(n, **kwargs):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kwargs)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


def _interacting_fermions(**kwargs):
    fc = fermionchain.Fermionic_Chain(4, **kwargs)
    nn = lambda i, j: 1.0 if abs(i-j) == 1 else 0.0
    fc.set_hamiltonian(fc.generate_bilinear(nn, fc.Cdag, fc.C)
                       + 2.0*fc.generate_bilinear(nn, fc.N, fc.N))
    return fc


# 1. get_dm(inds=...) used C[i] for i in range(len(inds)), i.e. the first
# len(inds) sites whatever inds said.
@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_density_matrix_reads_the_requested_sites(mode):
    fc = fermionchain.Fermionic_Chain(6, itensor_version="python")
    hop = fc.generate_bilinear(lambda i, j: (1.0+0.3*i) if abs(i-j) == 1 else 0.0,
                               fc.Cdag, fc.C)
    fc.set_hamiltonian(hop + 0.4*fc.N[1] - 0.2*fc.N[3])
    inds = [1, 3, 4]
    wf = _quiet(fc.get_gs, mode=mode)
    dm = wf.get_dm(inds=inds)
    for a, i in enumerate(inds):
        for b, j in enumerate(inds):
            ref = fc.vev(fc.Cdag[j]*fc.C[i], mode="ED")
            assert dm[a, b] == pytest.approx(ref, abs=1e-6)


# 2. Z=2 parafermions (the transverse-field Ising chain in clock language)
# could not be built: v2/v3 aborted the process in SpinX, even for
# mode="ED", and "python" had no site of type code -1.
def _z2_energy(version, mode="DMRG"):
    n = 4
    pc = parafermionchain.Parafermionic_Chain(n, Z=2, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h - pc.Sig[i]*pc.Sigd[i+1]
    for i in range(n):
        h = h - 0.7*pc.Tau[i]
    pc.set_hamiltonian(h)
    return _quiet(pc.gs_energy, mode=mode)


def _tfim_exact(n=4, g=0.7):
    sx = np.array([[0, 1], [1, 0]])
    sz = np.diag([1, -1])
    def site(o, i):
        out = np.array([[1.0]])
        for k in range(n):
            out = np.kron(out, o if k == i else np.eye(2))
        return out
    h = sum(-site(sx, i) @ site(sx, i+1) for i in range(n-1))
    h = h + sum(-g*site(sz, i) for i in range(n))
    return np.linalg.eigvalsh(h)[0]


@pytest.mark.parametrize("version", ["python"] + _cpp_versions())
def test_z2_parafermion_chain_is_the_ising_chain(version):
    exact = _tfim_exact()
    assert _z2_energy(version, mode="ED") == pytest.approx(exact, abs=1e-10)
    assert _z2_energy(version) == pytest.approx(exact, abs=1e-8)


def test_unsupported_parafermion_order_raises():
    with pytest.raises(ValueError):
        parafermionchain.Parafermionic_Chain(3, Z=5)


# 3. v3 (and v2) gave ITensor's two Davidson steps per local solve, which
# left a chain whose Hamiltonian conserves a charge the sites do not short
# of its ground state at the default 15 sweeps (-1.6215 to -1.70156 over
# repeats, against ED -1.7015621187).
@pytest.mark.parametrize("version", _cpp_versions())
def test_default_sweeps_reach_ed_on_a_small_fermion_chain(version):
    ed = _quiet(_interacting_fermions(itensor_version=version).gs_energy,
                mode="ED")
    for _ in range(5):
        fc = _interacting_fermions(itensor_version=version)
        assert _quiet(fc.gs_energy) == pytest.approx(ed, abs=1e-8)


# 4. get_charge_gap(mode="ED") raised TypeError (sector_gap took no mode).
@pytest.mark.parametrize("mode", ["ED", "DMRG"])
def test_charge_gap_on_both_solvers(mode):
    fc = fermionchain.Fermionic_Chain(4, itensor_version="python")
    fc.set_hamiltonian(sum(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
                           for i in range(3)))
    # free fermions, levels 2cos(k pi/5): half filling, then one more
    exact = 2*np.cos(2*np.pi/5)
    assert _quiet(fc.get_charge_gap, d=1, mode=mode) == pytest.approx(exact, abs=1e-6)


# 5. get_fidelity(fmode="PT", n=...) dropped n, so the sum always ran over
# the default two states.
def test_perturbative_fidelity_uses_n_states():
    n = 3
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    h1 = 2*sum(sc.Sx[i] for i in range(n))
    lam, delta = 0.8, 1e-4
    sc.set_hamiltonian(h0 + lam*h1)
    H = sc.get_ED_obj().get_operator(h0 + lam*h1)
    H1 = sc.get_ED_obj().get_operator(h1)
    H = H.toarray() if hasattr(H, "toarray") else np.asarray(H)
    H1 = H1.toarray() if hasattr(H1, "toarray") else np.asarray(H1)
    es, vs = np.linalg.eigh(H)
    w = np.abs(vs[:, 0].conj() @ H1 @ vs)**2
    exact = np.sum(w[1:]/((es[1:]-es[0])**2 + delta**2))
    chi = _quiet(fidelity.get_fidelity, sc, h0, h1, lam, fmode="PT",
                 n=2**n, delta=delta, mode="ED")
    assert chi == pytest.approx(exact, rel=1e-6)


# 6. The mean-field chain's gs_energy() counted the decoupled exchange
# twice (-2.5 against -1.25 on a 6-site Neel state).
def test_mean_field_energy_is_the_energy_of_the_model():
    n = 6
    sc = _heisenberg(n)
    m0 = [[0, 0, 0.5*(-1)**i] for i in range(n)]
    mf = _quiet(meanfield.spinchain_meanfield, sc, p=0.0, m0=m0, mix=0.5,
                mode="ED")
    e_mf = _quiet(mf.gs_energy, mode="ED")
    assert e_mf == pytest.approx(-0.25*(n-1), abs=1e-6)
    assert e_mf == pytest.approx(
        _quiet(mf.vev, sc.hamiltonian, mode="ED").real, abs=1e-6)


# 7. get_full_SS_correlator(mode="ED") left the last state of the manifold
# set on the chain, and submode="ED" with dex= then raised TypeError, as it
# did after any set_gs().
def test_full_ss_correlator_leaves_the_chain_alone():
    sc = _heisenberg(3)  # a degenerate doublet
    kw = dict(mode="ED", submode="ED", name=(sc.Sz[0], sc.Sz[0]), dex=0.1)
    _, before = _quiet(sc.get_dynamical_correlator, **kw)
    _quiet(sc.get_full_SS_correlator, mode="ED")
    _, after = _quiet(sc.get_dynamical_correlator, **kw)
    assert np.max(np.abs(after-before)) < 1e-12
    sc.set_gs(_quiet(sc.get_gs, mode="ED"))
    _quiet(sc.get_dynamical_correlator, **kw)  # used to raise TypeError


# 8. A J1-J2 chain at J2=0.4 on a one-site cell (dimerized ground state,
# so no uniform one-site state is a fixed point) returned energies from
# -0.38 up to +0.15 at maxm=4 with only converged=False to say so.
def _j1j2(n_uc, maxm, maxiter=200):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"]*n_uc,
                                            itensor_version="python")
    h = 0
    for i in range(n_uc):
        for op in ("Sx", "Sy", "Sz"):
            o0 = ic.get_operator(op, i)
            for r, J in ((1, 1.0), (2, 0.4)):
                k = i + r
                h = h + J*o0*ic.get_operator(op, k % n_uc, group=k // n_uc)
    ic.set_hamiltonian(h)
    ic.maxm, ic.maxiter = maxm, maxiter
    return ic


def test_unconverged_infinite_chain_warns():
    ic = _j1j2(1, 4, maxiter=20)
    ic.vumps_nrestarts = 1
    with pytest.warns(RuntimeWarning, match="did not converge"):
        _quiet(ic.gs_energy)
    assert not ic.converged


def test_j1j2_converges_on_a_two_site_cell():
    ic = _j1j2(2, 4)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        e = _quiet(ic.gs_energy)
    assert ic.converged
    assert e == pytest.approx(-0.3801176696, abs=1e-6)
