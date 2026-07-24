import numpy as np
from .twotime import kondo_term_from_two_time

# DMRG (itensor_version=3) construction of the two-time third-order Kondo
# term (twotime.py), replacing edtwotimeref.py's exact eigenbasis time
# evolution with real TDVP time evolution -- avoiding excited-state
# enumeration/diagonalization entirely, as requested. Reuses the exact
# same kernel machinery (theta0_filter, K_W) already validated against
# the excited-state-sum third_order_kondo_dIdV via edtwotimeref.py, since
# that machinery only ever consumes a discretely-sampled G(t2,tau) array
# -- it does not care how G was computed.
#
# IMPORTANT CAVEAT: this module was written against the verified,
# already-used-elsewhere DMRG API (Many_Body_Chain.toMPO, the C++
# tdvp_step single-step primitive already driven manually by tdz.py's
# _advance_complex_time_step, MPS operator application/overlap) but could
# NOT be executed or numerically validated in the environment this was
# developed in (no compiled C++ backend available: cppext.available(3) is
# False there). Treat it as a best-effort implementation that follows the
# same patterns already proven correct elsewhere in this codebase, not as
# independently verified the way edtwotimeref.py/twotime.py's ED path is
# (which matches the excited-state-sum reference to ~0.01-1%, depending
# on grid resolution). Run test_kondo_spectrum_dmrgtwotime.py (skipped
# automatically when no C++ backend is compiled) in an environment with
# itensor_version=3 compiled to actually check this.
#
# Algorithm (the "checkpoint-and-branch" construction from the design
# discussion): for each j in {Sx,Sy,Sz} (the eq. "3rd-normal" vertex
# applied at t=0):
#   1. build |psi_j(0)> = Sj|GS>
#   2. step |psi_j(0)> forward AND backward in time (tdvp_step with
#      +dt2/-dt2) to get a trajectory of checkpoints |psi_j(t2)> across
#      the full t2_grid (both signs, since G(t2,tau) is needed for
#      negative t2 too -- see twotime.py's K_W kernel, defined for all
#      real t2)
#   3. at each checkpoint, apply each k in {Sx,Sy,Sz}: |phi_jk(t2)> =
#      Sk|psi_j(t2)>
#   4. step |phi_jk(t2)> forward AND backward across the full tau_grid,
#      taking the overlap with each of the three fixed reference states
#      <GS|Sl (l in {Sx,Sy,Sz}) at every tau checkpoint -- this directly
#      gives G_jkl(t2,tau) for that (j,k,l) triple
#   5. combine via the Levi-Civita contraction (as edtwotimeref.py does)
#      into the single coeffG(t2,tau) twotime.py's kernels consume.
#
# This costs, per t2 checkpoint (of which there are len(t2_grid)), one
# tau-direction trajectory per (j,k) pair (9 combinations, several
# vanishing by antisymmetry but computed generically here) -- i.e.
# O(len(t2_grid)) separate short TDVP trajectories, each
# O(len(tau_grid)) steps. This is the "N_t2 separate DMRG runs" cost
# flagged from the start of this feature's design discussion.

_AXES = ("Sx", "Sy", "Sz")

_EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    _EPS3[_i, _j, _k] = 1.
for _i, _j, _k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    _EPS3[_i, _j, _k] = -1.


def _tdvp_trajectory(chain, Hop, wf0, dt, n_half):
    """(times, wfs): 2*n_half+1 checkpoints of wf0 evolved under Hop,
    spanning t in [-n_half*dt, +n_half*dt] in steps of dt, via repeated
    single-step tdvp_step calls (the same primitive tdz.py's
    _advance_complex_time_step drives manually; here with a plain real
    dt, forward for positive steps and backward -- i.e. dt<0, which
    tdvp_step already supports as a "possibly complex/signed" time step,
    per its own docstring -- for negative ones)."""
    from .. import mps as mpsmod
    times = [0.0]
    wfs = [wf0]
    wf = wf0
    for _ in range(n_half):
        handle = chain._session.tdvp_step(Hop.cpp_handle, wf.cpp_handle, dt)
        wf = mpsmod.MPS(chain, cpp_handle=handle).copy()
        times.append(times[-1] + dt)
        wfs.append(wf)
    wf = wf0
    for _ in range(n_half):
        handle = chain._session.tdvp_step(Hop.cpp_handle, wf.cpp_handle, -dt)
        wf = mpsmod.MPS(chain, cpp_handle=handle).copy()
        times.append(times[-1] - dt)
        wfs.append(wf)
    order = np.argsort(times)
    times = np.array(times)[order]
    wfs = [wfs[i] for i in order]
    return times, wfs


def _levi_civita_coeff_G_batches_dmrg(chain, site, dt2, n_t2_half, dtau,
                                       n_tau_half, t2_batch=1):
    """Generator of (t2_chunk, G_chunk) pairs -- see twotime.py's
    kondo_term_from_two_time for the contract -- built via real TDVP
    time evolution instead of edtwotimeref.py's eigenbasis shortcut.
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
                branch = ops[k][site]*psi_j_t2
                tau_times, tau_wfs = _tdvp_trajectory(chain, Hop, branch,
                                                       dtau, n_tau_half)
                if tau_grid_row is None: tau_grid_row = tau_times
                for ll, l in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    overlap = np.array([ref_wf[l].dot(wf) for wf in tau_wfs])
                    G_row[0, :] += c*overlap
        yield np.array([t2]), tau_grid_row, G_row


def two_time_kondo_term_dmrg(chain, site, eVs, omega0=20e-3, Gamma0=5e-6,
                              dt2=1.0, n_t2_half=200, dtau=1.0,
                              n_tau_half=200):
    """DMRG counterpart of edtwotimeref.two_time_kondo_term_ed -- see
    this module's docstring for the algorithm and its important
    untested-in-development caveat. dt2/n_t2_half and dtau/n_tau_half
    set the (uniform) time grids in each leg; see twotime.py's module
    docstring for the resolution/range requirements (K_W needs t2
    spacing finer than 1/omega0 and a range wider than several/Gamma0;
    the Hilbert-transform-based Theta0 filter is comparatively
    forgiving)."""
    def batches():
        for t2_chunk, _tau_row, G_chunk in _levi_civita_coeff_G_batches_dmrg(
                chain, site, dt2, n_t2_half, dtau, n_tau_half):
            yield t2_chunk, G_chunk

    # kondo_term_from_two_time only needs t2_grid/tau_grid for their
    # spacing (dt2/dtau) and midpoint check, both fully determined by the
    # step sizes and half-widths -- no need to wait for the actual
    # trajectory output to build them
    t2_grid = dt2*np.arange(-n_t2_half, n_t2_half + 1)
    tau_grid = dtau*np.arange(-n_tau_half, n_tau_half + 1)
    return kondo_term_from_two_time(t2_grid, tau_grid, batches(), eVs,
                                     omega0, Gamma0)
