import numpy as np
from .twotime import kondo_term_from_two_time

# Pure-ED reference implementation of the two-time third-order Kondo term
# (twotime.py), used to validate that construction independently of the
# excited-state-sum third_order_kondo_dIdV (conductance.py) -- and,
# eventually, as the small-system ground truth for the DMRG two-leg
# time-evolution construction (kondospectrumtk/dmrgtwotime.py). Computes
# G(t2,tau)=<GS|Sl(t2+tau)Sk(t2)Sj(0)|GS> via exact time evolution in the
# Hamiltonian eigenbasis already available on a KondoSpectrum (ks.e,
# ks.Sx/Sy/Sz), rather than by any actual real-time-evolution algorithm --
# appropriate here since the point is to test the two-time/kernel
# machinery itself, not ED's own (already independently validated)
# ability to get transition matrix elements.

_EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    _EPS3[_i, _j, _k] = 1.
for _i, _j, _k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    _EPS3[_i, _j, _k] = -1.

_AXES = ("Sx", "Sy", "Sz")


def _levi_civita_coeff_G_chunk(ks, t2_chunk, tau_grid):
    """coeffG(t2,tau) = sum_jkl eps_jkl <GS|Sl(t2+tau)Sk(t2)Sj(0)|GS>,
    the full Levi-Civita-contracted two-time correlator (eq. "3rd-normal"'s
    triple product, now as a function of two real times instead of three
    discrete states), on one t2-chunk x the full tau_grid."""
    e = ks.e
    ops = {"Sx": ks.Sx, "Sy": ks.Sy, "Sz": ks.Sz} # ops[name][a,b] = <a|S|b>
    phase_t2 = np.exp(-1j*np.outer(t2_chunk, e)) # (nb, dim)
    phase_tau = np.exp(-1j*np.outer(tau_grid, e)) # (ntau, dim)
    out = np.zeros((len(t2_chunk), len(tau_grid)), dtype=complex)
    for jj, j in enumerate(_AXES):
        v0 = ops[j][:, 0] # Sj|GS> in the eigenbasis
        v_t2 = phase_t2*v0[None, :] # Sj|GS> evolved to each t2
        for kk, k in enumerate(_AXES):
            phi_t2 = v_t2 @ ops[k].T # Sk applied at each t2
            for ll, l in enumerate(_AXES):
                c = _EPS3[jj, kk, ll]
                if c == 0.: continue
                ref_l = ops[l][:, 0] # <GS|Sl, for the final overlap
                weighted = phi_t2*np.conjugate(ref_l)[None, :]
                out += c*(weighted @ phase_tau.T)
    return out


def two_time_kondo_term_ed(ks, eVs, omega0=20e-3, Gamma0=5e-6,
                            t2_width=None, t2_npts=200_000,
                            tau_width=None, tau_npts=4_000,
                            t2_batch=2_000):
    """Term(eV) for every eV in `eVs` (see
    twotime.kondo_term_from_two_time's return-value convention) computed
    via exact-diagonalization two-time evolution. Requires ks.T==0 (only
    the ground state is a meaningful initial state for this single-
    three-point-function construction; see kondospectrumtk/dmrgkondo.py
    docs for why T=0 is the natural regime for this method in general).

    t2_width/tau_width default to scales tied to Gamma0/omega0 (see
    module docstring in twotime.py for why K_W needs t2 resolution finer
    than 1/omega0 and range wider than ~1/Gamma0). G(t2,tau) is built
    once per t2-chunk and reused for the whole eVs sweep."""
    if ks.T != 0.: raise ValueError("two_time_kondo_term_ed requires T=0")
    if t2_width is None: t2_width = 40./Gamma0
    if tau_width is None: tau_width = 2*np.pi/1e-5
    t2_grid = np.linspace(-t2_width, t2_width, t2_npts)
    tau_grid = np.linspace(-tau_width, tau_width, tau_npts, endpoint=False)

    def batches():
        for start in range(0, t2_npts, t2_batch):
            t2_chunk = t2_grid[start:start+t2_batch]
            yield t2_chunk, _levi_civita_coeff_G_chunk(ks, t2_chunk, tau_grid)

    return kondo_term_from_two_time(t2_grid, tau_grid, batches(), eVs,
                                     omega0, Gamma0)
