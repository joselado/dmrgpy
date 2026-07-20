"""
Complex-time-evolution dynamical correlator (submode "TDZ"), following
Cao, Lu, Stoudenmire & Parcollet, "Dynamical correlation functions from
complex time evolution", arXiv:2311.10909.

Idea: instead of evolving |psi(t)> along the real-time axis (submode
"TD", timedependent.py), evolve along a complex-time contour
z(t,alpha0) = Integral_0^t exp(-i*alpha0*f(t')) dt', f(t)=exp(-t*omega0)
(the paper's "decreasing-angle" contour, its Appendix A). Since
Im(z(t,alpha0))<0, this evolution damps high-energy content as it goes,
so the MPS bond dimension needed for a given accuracy grows far more
slowly than under real-time evolution (the paper reports chi~20-30 vs
chi~500-700 for comparable accuracy on the Anderson impurity model). The
true real-time (alpha0=0) correlator is then recovered order by order via
a perturbative Taylor expansion in alpha0 around this simulated contour
(Eq. 6), see _reconstruct_real_axis below.

Each order n of that expansion needs two ingredients (Eq. 7 and
Appendix B):
  - phi^(n)(t,alpha0) = <GS|O1 H^n|psi(t,alpha0)> -- since H is
    Hermitian this is just <H^n(O1^dagger GS)|psi(t,alpha0)>, so
    {H^n(O1^dagger GS)}_{n=0..n_max} can be precomputed ONCE (repeated
    MPO application to a fixed vector, independent of t) and reused as a
    fixed bra at every time step -- the per-step cost is then just
    n_max+1 overlaps against the evolving state, not n_max+1 fresh MPO
    applications per step.
  - J^(n)(t,alpha0) = -i*d^n/dalpha0^n[z(t,alpha0)] -- a pure scalar
    contour integral, no tensor-network work at all.

This module implements the algorithm as plain Python using only the
already-exposed toMPO()/StaticOperator/.dot() primitives (the same
pattern cvm.py's cvm_correction_vector uses for its own hand-rolled
Python-level CG algorithm), plus one new backend primitive per
itensor_version: a single complex-time-step propagator
(self._session.tdvp_step(H_handle, wf_handle, dz), see
pyitensor/chain.py and, once added, mpscpp3/mpscpp2's own pybind11
bindings) -- distinct from quench_tdvp/evolve_and_measure_tdvp
(timedependent.py), which assume one constant *real* dt for their whole
internal loop and thus can't be reused here, since the per-step contour
increment dz_k varies with t_k (f(t) is not constant).
"""
import numpy as np

from . import operatornames
from . import multioperator
from .timedependent import _fourier_transform_correlator

_MAX_SUPPORTED_ORDER = 4  # explicit Appendix-B g^(n) formulas only go this far


def _cumtrapz0(y, dt):
    """Cumulative trapezoidal integral of y (uniform grid, spacing dt),
    with out[0]=0 and len(out)==len(y)."""
    out = np.zeros(len(y), dtype=complex)
    out[1:] = np.cumsum(0.5*(y[:-1]+y[1:])*dt)
    return out


def _reconstruct_real_axis(alpha0, n_max, Jn, phi):
    """
    Combine phi^(n)(t,alpha0) (n=0..n_max, overlaps of H^n(B|GS>) against
    the complex-time-evolved state) and J^(n)(t,alpha0) (n=1..n_max, pure
    contour-integral scalars) into the real-time (alpha0->0) correlator,
    via the explicit order-by-order reconstruction of arXiv:2311.10909's
    Eq. 6, with the g^(n) terms given explicitly in its Appendix B
    (reproduced here verbatim for n=1..4; the paper finds n_max<=4-5
    always suffices for alpha0<~0.3, so this is a hardcoded cap rather
    than a fully general Faa-di-Bruno recursion for arbitrary n).

    phi[0](t,alpha0) is the raw complex-time correlator itself (n=0, no
    Hamiltonian insertion): G>(t,alpha0) = -i*phi[0](t,alpha0) is the
    base point of the Taylor series in Eq. 6, i.e.
    G>(t,0) ~= -i*(phi[0](t,alpha0) + sum_n g^(n)(t,alpha0)) -- the
    literal g^(n) formulas hardcoded below (Appendix B) are consistent
    with the paper's own Eq. 3/4/7/B1 definitions only once this overall
    "-i" on the *whole* bracket is included (confirmed independently by
    an exact single-eigenstate toy-model derivation: expanding
    G>(alpha0)=-i*exp(-i*z(alpha0)*E1) directly in alpha0 and comparing
    to the phi^(n)/J^(n) machinery pins down this placement uniquely).
    This function instead returns the *bare* correlator
    C(t,0) := i*G>(t,0) = phi[0] + sum_n g^(n), matching this codebase's
    own evolution_dmrg_DC/quench_tdvp convention (no -i prefactor)
    rather than the paper's own G> convention directly, so the result
    can feed the same downstream FFT tail unchanged (see
    dynamical_correlator_tdz) -- note this is *not* i*phi[0]+i*sum_n
    g^(n): the two i's from C:=i*G> and G>=-i*(...) cancel exactly.
    """
    if n_max > _MAX_SUPPORTED_ORDER:
        raise NotImplementedError(
                "TDZ only implements the explicit n<=%d Taylor "
                "reconstruction terms given in arXiv:2311.10909's "
                "Appendix B; the paper itself finds this always "
                "suffices for alpha0<~0.3" % _MAX_SUPPORTED_ORDER)
    gsum = 0.0
    if n_max >= 1:
        gsum = gsum - alpha0*Jn[1]*phi[1]
    if n_max >= 2:
        gsum = gsum + (alpha0**2/2.)*(Jn[2]*phi[1] + Jn[1]**2*phi[2])
    if n_max >= 3:
        gsum = gsum + (-alpha0**3/6.)*(Jn[3]*phi[1] + 3*Jn[1]*Jn[2]*phi[2]
                + Jn[1]**3*phi[3])
    if n_max >= 4:
        gsum = gsum + (alpha0**4/24.)*(Jn[4]*phi[1]
                + (4*Jn[1]*Jn[3] + 3*Jn[2]**2)*phi[2]
                + 6*Jn[1]**2*Jn[2]*phi[3] + Jn[1]**4*phi[4])
    return phi[0] + gsum


def _advance_complex_time_step(self, Hop, wf, dz, do_gse=False):
    """
    One complex-time evolution step of size dz (a possibly-complex
    scalar), advancing wf under exp(-i*dz*Hop). Uses two-site TDVP
    (whose forward/backward coefficients are already generic to complex
    time steps -- see pyitensor/tdvp.py's module docstring and
    mpscpp3/TDVP/README.md, which documents its own "t" argument as
    "real, imaginary, or complex") when available, matching exactly the
    same TDVP-vs-Taylor-MPO choice self.tevol_method/timedependent.py
    already makes for real-time evolution. Falls back to the MPO-Taylor
    evolve_taylor_step() (also generic to complex z -- see
    mpscpp2/mpscpp3's own evoloperator()/pyitensor's _evoloperator())
    otherwise -- mpscpp2's only route to TDZ, since it has no TDVP at
    all.

    self.tevol_method="TDVP_GSE" instead runs one-site TDVP; when
    do_gse=True (the caller gates this to the first self.tdvp_gse_sweeps
    steps, mirroring quench_tdvp_gse()/evolve_and_measure_tdvp_gse()) a
    global_subspace_expand() call precedes the one-site step, since dz's
    per-step contour increment varies with t here (unlike quench_tdvp_gse's
    fixed real dt), which is exactly why this stays a Python-level driver
    loop instead of a C++ one like those two. Same itensor_version support
    as plain "TDVP" (3 or "python" only) -- see evolution_dmrg_DC's
    docstring for why a v2 port isn't available here.
    """
    if self.itensor_version in (3, "python") and self.tevol_method == "TDVP":
        handle = self._session.tdvp_step(Hop.cpp_handle, wf.cpp_handle, dz)
    elif self.itensor_version in (3, "python") and self.tevol_method == "TDVP_GSE":
        wf_handle = wf.cpp_handle
        if do_gse:
            wf_handle = self._session.global_subspace_expand(Hop.cpp_handle,
                    wf_handle, self.tdvp_gse_krylov_order, self.tdvp_gse_cutoff, 0)
        handle = self._session.tdvp_step(Hop.cpp_handle, wf_handle, dz, 1)
    else:
        handle = self._session.evolve_taylor_step(Hop.cpp_handle, wf.cpp_handle, dz)
    from . import mps as mpsmod
    return mpsmod.MPS(self, cpp_handle=handle).copy()


def _complex_time_correlator(self, A, B, alpha0, n_max, dt, nt, omega0):
    """
    Single-direction complex-time-evolution correlator
    C(t,0) ~ <B GS|exp(-iHt)A GS> (bare convention, matching
    evolution_dmrg_DC/quench_tdvp -- see this module's docstring),
    reconstructed from an alpha0!=0 complex-time simulation. A,B are
    MultiOperators; the caller is responsible for any dagger convention
    (dynamical_correlator_tdz passes A=O1.get_dagger(),B=O2, exactly
    mirroring evolution_dmrg_DC's own A/B convention).

    Returns (ts,cs) with ts=dt*arange(nt), matching evolution_DC's own
    (pre-conjugation) return shape so the two submodes can share the same
    downstream windowing/FFT tail.
    """
    self.get_gs()  # ensures self.wf0 (ground state) and self.e0 are set
    Hop = self.toMPO(self.hamiltonian - self.e0)  # build once, reused every step
    wf_g = self.wf0

    wf_fwd = self.toMPO(A)*wf_g  # |psi(0,alpha0)> = A|GS>, forward-evolved
    bras = [self.toMPO(B)*wf_g]  # bras[0] = B|GS> (n=0, no H inserted)
    for _n in range(n_max):
        bras.append(Hop*bras[-1])  # bras[n] = H^n(B|GS>), built once, not per step

    ts = dt*np.arange(nt+1)
    f_t = np.exp(-ts*omega0)  # decreasing-angle contour f_D(t), Appendix A
    dz_dt = np.exp(-1j*alpha0*f_t)  # dz(t,alpha0)/dt
    z = _cumtrapz0(dz_dt, dt)
    dz = np.diff(z)  # per-step complex contour increment, nt steps

    Jn = {}
    for n in range(1, n_max+1):
        integrand = -1j*((-1j*f_t)**n)*dz_dt
        Jn[n] = _cumtrapz0(integrand, dt)

    phi = {n: np.zeros(nt+1, dtype=complex) for n in range(n_max+1)}
    for k in range(nt+1):
        for n in range(n_max+1):
            phi[n][k] = bras[n].dot(wf_fwd)
        if k < nt:
            wf_fwd = _advance_complex_time_step(self, Hop, wf_fwd, dz[k],
                    do_gse=(k < self.tdvp_gse_sweeps))

    cs = _reconstruct_real_axis(alpha0, n_max, Jn, phi)
    return ts[:nt], cs[:nt]


def dynamical_correlator_tdz(self, name="XX", es=None, alpha0=0.1, n_max=4,
        dt=0.1, tmax=None, nt=None, delta=5e-2, damping_periods=6,
        window=[-1, 10], factor=1, **kwargs):
    """
    Dynamical correlator via complex-time evolution + perturbative
    real-axis reconstruction (submode "TDZ"; Cao, Lu, Stoudenmire &
    Parcollet, "Dynamical correlation functions from complex time
    evolution", arXiv:2311.10909). See this module's own docstring for
    the algorithm.

    alpha0: complex-time contour angle parameter (paper default range
        0.1-0.3). Larger alpha0 reduces the bond dimension needed to
        reach a given accuracy further, but needs a higher n_max to
        reconstruct the real axis accurately (see the paper's Fig. 4).
    n_max: order of the Taylor-in-alpha0 reconstruction, <=4 (see
        _reconstruct_real_axis's docstring for why 4 is the current cap).
    dt,tmax/nt: as in submode "TD" (timedependent.py) -- dt is the
        real-time step of the underlying TDVP propagator (the complex
        contour is layered on top of it, not a replacement for it); nt
        defaults from damping_periods/delta/dt exactly as "TD" does, if
        tmax is not given either.
    es,delta,damping_periods,window,factor: passed straight through to
        the same windowing/FFT tail "TD" uses
        (timedependent._fourier_transform_correlator) -- delta here is
        the *final* Lorentzian broadening applied to the reconstructed
        real-time correlator before the FFT, a separate knob from alpha0
        (which only controls entanglement growth during the simulation
        itself, not the output broadening).

    Note (current scope): only the "greater" branch
    G>_{O1,O2}(t) ~ <B GS|exp(-iHt)A GS> is computed (A=O1^dagger, B=O2,
    exactly evolution_dmrg_DC's own convention) and treated as the full
    time-domain correlator fed to the FFT tail -- the same simplification
    timedependent.py's own "TD" submode already makes (it computes but
    never combines in a ga=conj(gr) "lesser" branch, see its
    dynamical_correlator()). A genuine G< second run (the paper's Eq. 5,
    alpha0->-alpha0 with O1<->O2 swapped) is a possible follow-up, not
    required to match "TD"'s existing fidelity.
    """
    if nt is None:
        if tmax is None: nt = int(damping_periods/delta/dt)
        else: nt = int(tmax/dt)
    tmax_eff = nt*dt
    omega0 = 2.*np.pi/tmax_eff  # paper's own convention: lowest resolvable frequency

    name = operatornames.str2MO(self, name, **kwargs)
    O1, O2 = name[0], name[1]
    A = O1.get_dagger()  # matches evolution_dmrg_DC's own A=O1.get_dagger()
    B = O2

    ts, cs = _complex_time_correlator(self, A, B, alpha0, n_max, dt, nt, omega0)
    cs = cs.real - 1j*cs.imag  # match evolution_DC's own conjugation convention
    return _fourier_transform_correlator(ts, cs, dt, es=es, window=window,
            delta=delta, factor=factor)
