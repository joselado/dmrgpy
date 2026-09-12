"""Regression tests for the 2026-09 audit's pyitensor time-evolution and
session-state cluster (findings 1, 2, 10 and 13 of
`docs/audit_2026_09_hole_hunt.md`).

All four are `itensor_version="python"` bugs, i.e. bugs of the backend a
bare `pip install dmrgpy` gets by default (`cppext.default_backend()`),
and all four are silent: a wrong number, never an exception. Each test
below therefore asserts the *right* number, with `itensor_version=3` or
ED as the reference wherever a reference is needed.

Chains are deliberately tiny -- at these sizes ED is exact, and (for
finding 1) a bond dimension above 2**(n/2) makes two-site TDVP exact, so
the assertions can be tight rather than tolerance-tuned.
"""

import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain, timedependent


# ---------------------------------------------------------------- helpers

def heisenberg(n, itensor_version, field=0.0, transverse=0.0):
    """Uniform S=1/2 Heisenberg chain, with an optional staggered Sz field
    and uniform Sx field (which break every symmetry that could otherwise
    make an evolution accidentally trivial)."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        if field: h = h + field * (-1) ** i * sc.Sz[i]
        if transverse: h = h + transverse * sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc, h


def _dense_mps(psi):
    """An MPS from pyitensor as a plain state vector, sites in order."""
    n = psi.length()
    T = psi.A(1)
    for i in range(2, n + 1):
        T = T * psi.A(i)
    sites = [next(ind for ind in psi.A(i).inds if ind.hastags("Site"))
             for i in range(1, n + 1)]
    order = [next(k for k, ind in enumerate(T.inds) if ind == s) for s in sites]
    return np.transpose(T.array, order).reshape(-1)


def _dense_mpo(H):
    """An MPO from pyitensor as a plain matrix (primed legs are the rows)."""
    n = H.length()
    T = H.A(1)
    for i in range(2, n + 1):
        T = T * H.A(i)
    ups, dns = [], []
    for i in range(1, n + 1):
        ss = [ind for ind in H.A(i).inds if ind.hastags("Site")]
        ups.append([x for x in ss if x.plev == 1][0])
        dns.append([x for x in ss if x.plev == 0][0])
    pos = lambda s: next(k for k, ind in enumerate(T.inds) if ind == s)
    arr = np.transpose(T.array, [pos(s) for s in ups] + [pos(s) for s in dns])
    d = int(np.prod(arr.shape[:n]))
    return arr.reshape(d, d)


# ------------------------------------------------- finding 1: TDVP order

def _tdvp_endpoint(itensor_version, dt, T=0.4, n=6):
    """<Sz0(T)> after evolving Sz[0]|gs> with TDVP at step dt."""
    sc, h = heisenberg(n, itensor_version, field=0.3, transverse=0.2)
    sc.maxm, sc.cutoff, sc.tevol_method = 64, 1e-14, "TDVP"
    wf = sc.applyoperator(sc.Sz[0], sc.get_gs())
    nt = int(round(T / dt)) + 1
    _ts, cs = timedependent.evolve_and_measure_dmrg(
            sc, operator=sc.Sz[0], wf=wf, nt=nt, dt=dt)
    return complex(cs[-1])


def test_python_tdvp_is_dt_independent_at_full_bond_dimension():
    """Two-site TDVP at a bond dimension large enough that no SVD ever
    truncates is *exact* (Lubich & Oseledets' projector-splitting
    exactness), so halving dt must not move the answer at all. On
    itensor_version="python" it used to move by exactly a factor of two --
    O(dt), one order worse even than the second-order integrator the
    module docstring describes -- because _half_sweep_lr was handed a
    *left*-canonical state (which is what applyMPO returns) where it needs
    a right-canonical one, making every right environment of the first
    step a non-orthogonal overlap matrix instead of an isometry. Only the
    first step of a trajectory was wrong, which is why the symptom at
    fixed total time looked like a first-order method.

    This pins the *order*, not a value: the three dt must agree with each
    other, and with v3, to near machine precision. Before the fix the
    dt=0.2 and dt=0.05 answers differed by 6.3e-03."""
    vals = {dt: _tdvp_endpoint("python", dt) for dt in (0.2, 0.1, 0.05)}
    ref = _tdvp_endpoint(3, 0.2)
    for dt, v in vals.items():
        assert v.real == pytest.approx(ref.real, abs=1e-8), \
            "python TDVP at dt=%g: %r vs v3 %r" % (dt, v, ref)
    spread = max(abs(a - b) for a in vals.values() for b in vals.values())
    assert spread < 1e-8


def test_python_tdvp_single_step_is_exact_at_full_bond_dimension():
    """The same defect seen at the primitive it lives in: one
    tdvp.tdvp_step() on a state that arrives *left*-canonical (which is
    what applyMPO returned when the bug was found, and what any
    position(n) leaves behind) must reproduce expm(-i*dt*H) exactly at
    full bond dimension. It used to be off by 2.9e-02/1.5e-02/7.4e-03 at
    dt=0.2/0.1/0.05 -- visibly first order. The gauge is set explicitly
    here rather than inherited from applyMPO, so that the invariant stays
    pinned however applyMPO's own output gauge changes."""
    from scipy.linalg import expm
    from dmrgpy.pyitensor.tdvp import tdvp_step
    n = 6
    sc, h = heisenberg(n, "python", field=0.3, transverse=0.2)
    sc.maxm, sc.cutoff = 128, 0.0
    sc.get_gs()
    session = sc._session
    psi0 = session.apply_pure_operator(session._mpo(sc.Sz[0].to_terms()),
                                       session.wf0)
    psi0.position(n)  # left-canonical: the gauge the bug needed
    assert psi0.center == n
    v0, Hd = _dense_mps(psi0), _dense_mpo(session.H)
    for dt in (0.2, 0.1, 0.05):
        out = tdvp_step(psi0.copy(), session.H, dt, cutoff=0.0, maxdim=128,
                        niter=50)
        err = np.linalg.norm(_dense_mps(out) - expm(-1j * dt * Hd) @ v0)
        assert err < 1e-9, "tdvp_step at dt=%g: error %.3e" % (dt, err)


# --------------------------------------- finding 2: stale DMRG start state

def test_set_hamiltonian_does_not_reuse_the_previous_ground_state():
    """A parameter sweep on one chain object must give the same answers as
    a fresh chain per point. The session used to keep the previous
    Hamiltonian's converged MPS as its DMRG start, and pyitensor has no
    noise term to escape with, so a start that happens to be an exact
    eigenstate of the new Hamiltonian left the solve stationary: the B=0
    point came back as the ferromagnetic *maximum* retained from the
    polarized end of the sweep."""
    n = 6
    reused, _ = heisenberg(n, "python")
    reused.maxm, reused.nsweeps = 40, 20
    for B in (3.0, 1.0, 0.5, 0.0):
        fresh, hf = heisenberg(n, "python")
        fresh.maxm, fresh.nsweeps = 40, 20
        hf = hf + B * sum(fresh.Sz[i] for i in range(n))
        fresh.set_hamiltonian(hf)
        e_ed = fresh.gs_energy(mode="ED")

        _, hr = heisenberg(n, "python")
        hr = 0
        for i in range(n - 1):
            hr = hr + reused.Sx[i] * reused.Sx[i + 1] \
                    + reused.Sy[i] * reused.Sy[i + 1] \
                    + reused.Sz[i] * reused.Sz[i + 1]
        hr = hr + B * sum(reused.Sz[i] for i in range(n))
        reused.set_hamiltonian(hr)
        assert reused.gs_energy() == pytest.approx(e_ed, abs=1e-6), \
            "B=%g on a reused chain" % B


def test_bandwidth_on_the_python_backend():
    """`Many_Body_Chain.bandwidth` is literally the two-set_hamiltonian-on-
    one-clone pattern of the test above (it solves for h and then for -h on
    the same clone), so the retained start state made it return 0.0 instead
    of the operator's actual spectral width."""
    n = 4
    sc, h = heisenberg(n, "python")
    sc.maxm, sc.nsweeps = 30, 20
    # exact max-min of the 4-site Heisenberg chain: 0.75 - (-1.6160254)
    assert sc.bandwidth(h) == pytest.approx(2.3660254, abs=1e-5)


# --------------------------------------------- finding 10: promote_to_dense

def _hubbard_nf3(itensor_version, U=1.7, n=3):
    fc = fermionchain.Spinful_Fermionic_Chain(n,
            itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagdn[i] * fc.Cdn[i + 1]
        h = h + fc.Cdagup[i + 1] * fc.Cup[i] + fc.Cdagdn[i + 1] * fc.Cdn[i]
    for i in range(n):
        h = h + U * fc.Nup[i] * fc.Ndn[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 64, 30
    return fc


def test_promote_to_dense_keeps_the_sector_energy():
    """promote_to_dense()'s own docstring guarantees that "the ground-state
    energy and wavefunction are kept, so a bare gs_energy() afterwards
    returns the sector's energy rather than re-solving unconstrained".
    itensor_version="python" used to re-solve and return the *global*
    ground state: the Hamiltonian re-send that promotion forces (the MPO
    was built on the QN indices) cleared the session's energy cache, and
    python confines a sector with a penalty on the variational solve
    rather than structurally, so nothing held the re-solve inside it.

    The n=3 Hubbard chain at Nf=3 is the case that can see this: unlike
    Sz=0 or spinless Nf=3, its sector energy (-1.94081402) and its global
    one (-2.33991308) genuinely differ."""
    for backend in (3, "python"):
        fc = _hubbard_nf3(backend)
        fc.set_conserved_sector(Nf=3)
        e_sector = fc.gs_energy()
        assert e_sector == pytest.approx(-1.9408140222, abs=1e-6)
        fc.promote_to_dense()
        assert fc.gs_energy() == pytest.approx(e_sector, abs=1e-6), \
            "%s re-solved unconstrained after promote_to_dense" % (backend,)
        assert fc.vev(sum(fc.N)).real == pytest.approx(3.0, abs=1e-5)


# ------------------------------------------------- finding 13: applyMPO

def test_applympo_is_exact_at_a_lossless_bond_dimension():
    """applyMPO used to be an unorthogonalized zip-up: it truncated each
    cut while building the product, in a gauge where the singular values
    are not the Schmidt values of the exact result, so discarding the
    smallest of them was not the optimal truncation. On an 8-site
    Heisenberg chain the exact H|psi> has Schmidt rank at most 16 at every
    bond -- maxdim=16 therefore discards literally nothing -- and the
    zip-up still lost 9.4e-05 in 2-norm."""
    n = 8
    sc, h = heisenberg(n, "python")
    sc.maxm, sc.cutoff = 64, 1e-14
    sc.get_gs()
    session = sc._session
    psi, H = session.wf0, session.H
    exact = _dense_mpo(H) @ _dense_mps(psi)

    # the premise: maxdim=16 is lossless here, verified rather than assumed
    for b in range(1, n):
        s = np.linalg.svd(exact.reshape(2 ** b, -1), compute_uv=False)
        assert (s > 1e-12 * s[0]).sum() <= 16

    from dmrgpy.pyitensor.mpsalgebra import applyMPO
    from dmrgpy.pyitensor.tensor import noPrime
    for maxdim in (16, 32, None):
        out = applyMPO(H, psi, maxdim=maxdim)
        for j in range(1, out.length() + 1):
            out.set_A(j, noPrime(out.A(j), "Site"))
        err = np.linalg.norm(_dense_mps(out) - exact)
        assert err < 1e-10, "applyMPO(maxdim=%s): error %.3e" % (maxdim, err)


def test_vev_npow_matches_the_explicitly_squared_operator():
    """The consumer that made finding 13 visible: `vev(h, npow=2)` applies
    the Hamiltonian MPO to the state, while `vev(h*h)` builds H^2 as an MPO
    and needs no application at all. They must agree; at stock defaults on
    the python backend they used to differ by 2.6e-09, and
    gs_energy_fluctuation() -- advertised as "a measure of how sharply the
    state is an eigenstate" -- read 5.1e-05 on a state whose true
    fluctuation is ~8e-06."""
    n = 10
    sc, h = heisenberg(n, "python")  # stock maxm=30, nsweeps=15
    assert sc.vev(h, npow=2).real == pytest.approx(sc.vev(h * h).real,
                                                   abs=1e-8)
    assert sc.gs_energy_fluctuation() < 1e-5
