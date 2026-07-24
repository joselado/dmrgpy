"""Coverage for the DMRG two-time-correlator construction of the
third-order Kondo term (kondospectrumtk/dmrgtwotime.py), following
Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.

This module was written against the verified DMRG API (toMPO, tdvp_step,
MPS operator application/overlap) but could not be exercised until a
compiled itensor_version=3 backend was available. Once compiled, running
this test surfaced three real, now-fixed bugs, none of which showed up
in earlier ED-only testing:

  - tdvp_step silently renormalizes its output to unit norm on every
    single call. Sj|GS> (the state each two-time trajectory starts from)
    is generally NOT unit-norm (Sj isn't norm-preserving), so letting
    that renormalization stand silently discarded the true amplitude at
    every step -- confirmed directly, it produced overlaps exactly a
    factor of 2 too large end to end for a state of norm 0.5. Fixed in
    _tdvp_trajectory by evolving a normalized copy and rescaling every
    returned checkpoint by the true original norm.
  - The forward/backward time-stepping loop shared one running `times`
    list, so the "backward" branch kept counting down from the forward
    branch's endpoint instead of from t=0 -- it never actually reached
    negative times at all (e.g. n_half=3, dt=1000 produced
    times=[0,0,1000,1000,2000,2000,3000], not the intended symmetric
    [-3000,...,3000]). Fixed by walking forward and backward from t=0
    independently and merging afterward.
  - twotime.py's t2-integral used np.trapz per chunk, which is exactly 0
    for a single-point chunk (no width to integrate over) -- and DMRG
    chunks are always single-point (each t2 checkpoint is its own real
    trajectory), so this silently zeroed the entire third-order Kondo
    term end to end. Fixed in kondo_term_from_two_time with a uniform
    Riemann-sum accumulation (trapezoidal endpoint correction only at
    the true global first/last t2 points) that is well defined for any
    chunk size.

After these fixes, the two-time G(t2,tau) construction matches the ED
reference (edtwotimeref.py) to ~1e-9-1e-10 pointwise, and the full
eV-swept third-order Kondo term matches to ~1e-10 when compared against
a grid-consistent ED reference (same t2/tau grid convention on both
sides -- an earlier attempt using edtwotimeref.py's own
width/npts-based grid builder directly showed an O(0.01-0.1) spurious
mismatch purely from the two grid builders using different endpoint
conventions at nominally "matching" resolution, not a physics bug; see
the grid-consistent construction below, which sidesteps that by
building both G(t2,tau) datasets on the literal same array).

Every test below is parametrized over itensor_version in (3, "python"):
dmrgtwotime.py never hardcodes which compiled/pure-Python backend it's
driving (see that module's own docstring), so it was tried directly
against a pyitensor (itensor_version="python", no compiled extension
needed) chain too and found to work unmodified, matching the ED
reference even more tightly (~1e-14) than the compiled v3 backend does.
The itensor_version=3 parametrization is still skipped when no compiled
backend is available; itensor_version="python" always runs.
"""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy import cppext
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.twotime import kondo_term_from_two_time
from dmrgpy.kondospectrumtk.edtwotimeref import _levi_civita_coeff_G_chunk

G = 2.0
MUB = 5.7883818066e-5 # eV/T

_BACKENDS = [
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="needs a compiled itensor_version=3 backend"), id="v3"),
    pytest.param("python", id="pyitensor"),
]


def _build_chain(itensor_version):
    """A 3-site (not 1-site: ITensor v3's two-site DMRG/TDVP can't handle
    chains shorter than 3 sites -- confirmed directly, a 1-site chain hit
    an internal C++ index error building the Hamiltonian MPO; pyitensor
    has no such limitation, but the same chain size is used for both
    backends here for a fair/uniform comparison) S=1/2 chain: a
    Zeeman-split impurity at site 0 (the tip-coupled site) with a weak
    Heisenberg coupling to two spectator sites, giving DMRG a
    properly-sized system to operate on without changing the impurity
    physics much."""
    sc = spinchain.Spin_Chain(["1/2", "1/2", "1/2"], itensor_version=itensor_version)
    h = G*MUB*10.0*sc.Sz[0]
    for i in range(2):
        h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_tdvp_trajectory_matches_ed_time_evolution(itensor_version):
    """Direct check of the fixed _tdvp_trajectory: <ref|exp(-iHt)Sj|GS>
    at every checkpoint of a forward+backward trajectory, for an
    operator combination with a genuinely nonzero overlap (some
    combinations vanish by a symmetry of this particular Hamiltonian --
    confirmed directly; Sx-Sx does not)."""
    from dmrgpy.kondospectrumtk.dmrgtwotime import _tdvp_trajectory

    sc = _build_chain(itensor_version)
    gs = sc.get_gs()
    E0 = sc.gs_energy()
    Hop = sc.toMPO(sc.hamiltonian - E0)
    ks = KondoSpectrum(sc, site=0, T=0.0)

    psi0 = sc.Sx[0]*gs
    dt, n_half = 1000.0, 3
    times, wfs = _tdvp_trajectory(sc, Hop, psi0, dt, n_half)
    assert np.allclose(times, dt*np.arange(-n_half, n_half+1))

    ref = sc.Sx[0]*gs
    ref_vec = ks.Sx[:, 0]
    v0 = ks.Sx[:, 0]
    for t, wf in zip(times, wfs):
        got = ref.dot(wf)
        expected = np.conjugate(ref_vec).dot(v0*np.exp(-1j*ks.e*t))
        assert got == pytest.approx(expected, abs=1e-6)


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_two_time_G_matches_ed_reference_pointwise(itensor_version):
    """The full Levi-Civita-contracted G(t2,tau) construction
    (_levi_civita_coeff_G_batches_dmrg), compared point by point against
    the ED reference on a small grid."""
    from dmrgpy.kondospectrumtk.dmrgtwotime import _levi_civita_coeff_G_batches_dmrg

    sc = _build_chain(itensor_version)
    sc.get_gs()
    ks = KondoSpectrum(sc, site=0, T=0.0)

    dt2, n_t2_half = 500.0, 2
    dtau, n_tau_half = 300.0, 2
    dmrg_rows = list(_levi_civita_coeff_G_batches_dmrg(
            sc, 0, dt2, n_t2_half, dtau, n_tau_half))

    e = ks.e
    axes = ["Sx", "Sy", "Sz"]
    Sops = {"Sx": ks.Sx, "Sy": ks.Sy, "Sz": ks.Sz}
    eps3 = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]: eps3[a, b, c] = 1.
    for a, b, c in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]: eps3[a, b, c] = -1.

    def coeffG_ed(t2, tau):
        total = 0j
        for jj, j in enumerate(axes):
            v_t2 = Sops[j][:, 0]*np.exp(-1j*e*t2)
            for kk, k in enumerate(axes):
                phi_tau = (Sops[k]@v_t2)*np.exp(-1j*e*tau)
                for ll, l in enumerate(axes):
                    c = eps3[jj, kk, ll]
                    if c == 0.: continue
                    total += c*np.conjugate(Sops[l][:, 0]).dot(phi_tau)
        return total

    maxdiff = 0.
    for t2_chunk, tau_row, G_chunk in dmrg_rows:
        t2 = t2_chunk[0]
        for it, tau in enumerate(tau_row):
            diff = abs(G_chunk[0, it] - coeffG_ed(t2, tau))
            maxdiff = max(maxdiff, diff)
    assert maxdiff < 1e-6


@pytest.mark.parametrize("itensor_version", _BACKENDS)
def test_two_time_kondo_term_dmrg_matches_grid_consistent_ed_reference(itensor_version):
    """End-to-end: the full eV-swept third-order Kondo term via
    two_time_kondo_term_dmrg, against an ED reference built on the
    literal same (t2,tau) grid (not edtwotimeref.py's own width/npts
    grid builder, whose endpoint convention differs from
    dmrgtwotime.py's at nominally "matching" resolution -- see this
    module's docstring)."""
    from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg

    sc = _build_chain(itensor_version)
    sc.get_gs()
    ks = KondoSpectrum(sc, site=0, T=0.0)

    omega0, Gamma0 = 2e-3, 5e-6
    n_t2_half, n_tau_half = 10, 15
    dt2 = 25./Gamma0/n_t2_half
    dtau = (2*np.pi/2e-5)/n_tau_half
    eVs = np.linspace(-2e-3, 2e-3, 5)

    dmrg_term = two_time_kondo_term_dmrg(
            sc, 0, eVs, omega0=omega0, Gamma0=Gamma0,
            dt2=dt2, n_t2_half=n_t2_half, dtau=dtau, n_tau_half=n_tau_half)

    t2_grid = dt2*np.arange(-n_t2_half, n_t2_half+1)
    tau_grid = dtau*np.arange(-n_tau_half, n_tau_half+1)
    def batches():
        yield t2_grid, _levi_civita_coeff_G_chunk(ks, t2_grid, tau_grid)
    ed_term = kondo_term_from_two_time(t2_grid, tau_grid, batches(), eVs,
                                        omega0, Gamma0)

    assert np.max(np.abs(dmrg_term - ed_term)) < 1e-4
