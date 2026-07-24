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
# VALIDATED (see test_kondo_spectrum_dmrgtwotime.py): after a compiled
# itensor_version=3 backend became available, this module's G(t2,tau)
# matched the ED reference (edtwotimeref.py) to ~1e-9-1e-10 pointwise,
# and the full eV-swept third-order Kondo term matched a grid-consistent
# ED reference to ~1e-10. Getting there surfaced three real bugs, none of
# which had shown up in the ED-only testing this module was originally
# written against -- see _tdvp_trajectory's and
# twotime.kondo_term_from_two_time's docstrings for the details (in
# short: tdvp_step silently renormalizes every step to unit norm, which
# discards Sj|GS>'s true amplitude unless corrected for explicitly; a
# forward/backward time-stepping bug meant "backward" checkpoints never
# actually reached negative times; and a per-chunk np.trapz integral is
# exactly 0 for the single-t2-point chunks this module's real-time
# evolution necessarily produces). A 1-site test chain also hit an
# internal ITensor v3 error unrelated to any of the above -- see
# _build_chain in the test file for why chains need >=3 sites.
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
    per its own docstring -- for negative ones).

    tdvp_step renormalizes its output to unit norm at every single call
    (confirmed directly: even a near-zero dt=0.01 step takes an
    input norm^2=0.25 state straight to norm^2=1) -- correct/harmless for
    evolving an already-normalized physical state, but wf0 here is
    Sj|GS> or Sk|psi(t2)>, which is generally NOT unit-norm (Sj isn't
    norm-preserving). Silently letting that renormalization stand would
    throw away wf0's true amplitude at every step -- confirmed directly,
    it produced overlaps exactly a factor of 2 too large end to end for a
    wf0 of norm 0.5. Fixed by normalizing wf0 before the loop and
    rescaling every returned checkpoint by the true original norm
    afterward, so the trajectory is evolved on a unit-norm state
    throughout (where tdvp_step's renormalization is a no-op) but
    reported at wf0's actual scale."""
    from .. import mps as mpsmod
    norm0 = np.sqrt(wf0.dot(wf0).real)
    wf0_unit = wf0*(1./norm0) if norm0 > 1e-12 else wf0

    def _walk(step):
        # forward (step=dt) or backward (step=-dt) checkpoints from t=0,
        # NOT including t=0 itself -- each direction must restart its own
        # time counter at 0 rather than continuing the other direction's
        # (confirmed directly: sharing one running `times` list across
        # both loops made the "backward" branch keep counting down from
        # the forward branch's endpoint instead of from 0, so it never
        # actually reached negative times at all -- e.g. n_half=3,
        # dt=1000 produced times=[0,0,1000,1000,2000,2000,3000], not the
        # intended [-3000,...,0,...,3000])
        t, wf, out_t, out_wf = 0.0, wf0_unit, [], []
        for _ in range(n_half):
            handle = chain._session.tdvp_step(Hop.cpp_handle, wf.cpp_handle, step)
            wf = mpsmod.MPS(chain, cpp_handle=handle).copy()
            t = t + step
            out_t.append(t); out_wf.append(wf)
        return out_t, out_wf

    fwd_t, fwd_wf = _walk(dt)
    bwd_t, bwd_wf = _walk(-dt)
    times = np.array(bwd_t[::-1] + [0.0] + fwd_t)
    wfs = bwd_wf[::-1] + [wf0_unit] + fwd_wf
    wfs = [norm0*wf for wf in wfs]
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
