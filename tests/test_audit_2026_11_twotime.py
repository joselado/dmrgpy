"""Pins for the two-time Kondo route (2026-11 audit, cluster twotime).

correlator-twotime-1: G(t2,tau) is built from three bra trajectories of
Sl|GS>, computed once, instead of six ket trajectories per t2 row. The old
construction is kept verbatim below as `_old_batches` and the new rows are
compared against it at saturation, with a real and a complex Hamiltonian.

correlator-twotime-4: an integer dt2 or dtau works like the float, and a bad
grid argument raises before the first TDVP step.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, cppext
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
from dmrgpy.kondospectrumtk.dmrgtwotime import _tdvp_trajectory, _EPS3, _AXES
from dmrgpy.kondospectrumtk.dmrgtwotime import _levi_civita_coeff_G_batches_dmrg

G = 2.0
MUB = 5.7883818066e-5

_BACKENDS = [
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3), reason="needs compiled v3"), id="v3"),
    pytest.param("python", id="pyitensor"),
]


# the construction before the fix, verbatim: six ket trajectories per row
def _old_batches(chain, site, dt2, n_t2_half, dtau,
                                       n_tau_half, t2_batch=1):
    """Generator of (t2_chunk, tau_row, G_chunk, Gx_chunk) rows -- see
    twotime.py's kondo_term_from_two_time for the contract on the
    (t2_chunk, G_chunk, Gx_chunk) it consumes -- built via real TDVP
    time evolution instead of edtwotimeref.py's eigenbasis shortcut.
    Gx_chunk is the exchange diagram's Gx(s,tau) at s = t2 (this module's
    step 6).
    t2_batch=1 here (each checkpoint's tau-trajectory is its own
    "chunk"): unlike the ED case, there is no cheap way to batch many t2
    checkpoints into one vectorized array, since each one requires its
    own independent tau-direction TDVP trajectory."""
    ops = {name: getattr(chain, name) for name in _AXES}
    E0 = chain.gs_energy()
    Hop = chain.toMPO(chain.hamiltonian - E0) # e0-shifted, matching the
                                               # e[0]=0 convention used
                                               # throughout kondospectrumtk
    gs = chain.get_gs()

    ref_wf = {l: ops[l][site]*gs for l in _AXES} # <GS|Sl, i.e. Sl|GS> (l fixed)

    t2_times, t2_wfs_by_j = {}, {}
    for j in _AXES:
        psi_j0 = ops[j][site]*gs
        times, wfs = _tdvp_trajectory(chain, Hop, psi_j0, dt2, n_t2_half)
        t2_times[j] = times
        t2_wfs_by_j[j] = wfs

    n_t2 = len(t2_times[_AXES[0]])
    for it2 in range(n_t2):
        t2 = t2_times[_AXES[0]][it2]
        G_row = np.zeros((1, 2*n_tau_half + 1), dtype=complex)
        tau_grid_row = None
        for jj, j in enumerate(_AXES):
            psi_j_t2 = t2_wfs_by_j[j][it2]
            for kk, k in enumerate(_AXES):
                if jj == kk: continue # eps_jkl==0 for every l here -- skip
                                       # the expensive tau trajectory itself,
                                       # not just the (already-guarded) sum
                branch = ops[k][site]*psi_j_t2
                tau_times, tau_wfs = _tdvp_trajectory(chain, Hop, branch,
                                                       dtau, n_tau_half)
                if tau_grid_row is None: tau_grid_row = tau_times
                for ll, l in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    overlap = np.array([ref_wf[l].dot(wf) for wf in tau_wfs])
                    G_row[0, :] += c*overlap
        # exchange diagram (step 6): Gx(s,tau) = <psi_l(-s-tau)|Sk|psi_j(-s)>
        # at s = t2; the t2 grid is symmetric, so -s is the mirrored index
        mirror = n_t2 - 1 - it2
        kets = {(j, k): ops[k][site]*t2_wfs_by_j[j][mirror]
                for j in _AXES for k in _AXES if j != k}
        Gx_row = np.zeros((1, 2*n_tau_half + 1), dtype=complex)
        for ll, l in enumerate(_AXES):
            # checkpoints psi_l(-s+tau') for tau' on the tau grid; reversed,
            # entry n is psi_l(-s-tau_n) (the grid is symmetric)
            _, bras = _tdvp_trajectory(chain, Hop, t2_wfs_by_j[l][mirror],
                                       dtau, n_tau_half)
            bras = bras[::-1]
            for jj, j in enumerate(_AXES):
                for kk, k in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    Gx_row[0, :] += c*np.array([b.dot(kets[(j, k)])
                                                for b in bras])
        yield np.array([t2]), tau_grid_row, G_row, Gx_row


def _chain(version, n=4, complex_field=False):
    sc = spinchain.Spin_Chain(["1/2"]*n, itensor_version=version)
    h = G*MUB*10.0*sc.Sz[0]
    if complex_field:
        h = h + 3e-5*(0.6*sc.Sx[1] + 0.8*sc.Sy[1])
    for i in range(n - 1):
        h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    sc.maxm = 40
    sc.nsweeps = 12
    return sc


@pytest.mark.parametrize("version", _BACKENDS)
@pytest.mark.parametrize("complex_field", [False, True])
def test_G_from_bra_trajectories_equals_ket_trajectories(version, complex_field,
                                                       monkeypatch):
    sc = _chain(version, complex_field=complex_field)
    sc.get_gs()
    args = (sc, 0, 500.0, 2, 300.0, 3)
    count = {"old": 0, "new": 0}
    def counted(key, f):
        def g(*a, **k):
            count[key] += 1
            return f(*a, **k)
        return g
    monkeypatch.setitem(globals(), "_tdvp_trajectory", counted("old", _tdvp_trajectory))
    monkeypatch.setattr(dt, "_tdvp_trajectory", counted("new", dt._tdvp_trajectory))
    rows_old = list(_old_batches(*args))
    rows_new = list(_levi_civita_coeff_G_batches_dmrg(*args))
    # the saving of finding 11: three bra trajectories per row, not six ket
    # ones (48 before)
    assert count == {"old": 48, "new": 21}
    assert len(rows_old) == len(rows_new) == 5
    scale = max(np.abs(r[2]).max() for r in rows_old)
    assert scale > 0.1
    for a, b in zip(rows_old, rows_new):
        assert np.allclose(a[0], b[0])
        assert np.allclose(a[1], b[1], atol=1e-9)
        assert np.abs(a[2] - b[2]).max() < 1e-8*scale
        assert np.abs(a[3] - b[3]).max() < 1e-12


def _kw(**k):
    d = dict(omega0=2e-3, Gamma0=2e-4, dt2=1000., n_t2_half=1, dtau=1000.,
             n_tau_half=3)
    d.update(k)
    return d


@pytest.mark.parametrize("name", ["dt2", "dtau"])
def test_integer_step_equals_float_step(name):
    sc = _chain("python", n=3)
    sc.get_gs()
    eV = np.linspace(-2e-3, 2e-3, 3)
    a = dt.two_time_kondo_term_dmrg(sc, 0, eV, **_kw(**{name: 1000.}))
    b = dt.two_time_kondo_term_dmrg(sc, 0, eV, **_kw(**{name: 1000}))
    assert np.allclose(a, b, rtol=0, atol=1e-14)


@pytest.mark.parametrize("bad", [dict(dt2=0.), dict(dt2=-1.), dict(dtau=0),
                                 dict(n_t2_half=0), dict(n_tau_half=0),
                                 dict(n_t2_half=1.5), dict(n_tau_half=-2),
                                 dict(dt2=float("nan"))])
def test_bad_grid_raises_before_any_tdvp_step(bad, monkeypatch):
    sc = _chain("python", n=3)
    sc.get_gs()

    def boom(*a, **k):
        raise AssertionError("a TDVP trajectory started")
    monkeypatch.setattr(dt, "_tdvp_trajectory", boom)
    with pytest.raises(ValueError):
        dt.two_time_kondo_term_dmrg(sc, 0, np.array([1e-3]), **_kw(**bad))
