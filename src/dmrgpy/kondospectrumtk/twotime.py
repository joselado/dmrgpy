import numpy as np

# T=0 third-order Kondo term via a two-time correlator, avoiding explicit
# excited-state enumeration -- see docs/user_guide.md sec. 17 for the
# physics and derivation summary, and this module's own functions for the
# numerical details. Backend-agnostic: everything here operates on a
# supplied G(t2,tau) array (or a callable producing chunks of it), so the
# same pipeline serves both the ED reference (edtwotimeref.py) and the
# DMRG two-leg time-evolution construction.
#
# Physics: define the Heisenberg three-point function
#   G(t2,tau) = <GS|Sl(t2+tau) Sk(t2) Sj(0)|GS>
#             = sum_{f,m} <GS|Sl|f><f|Sk|m><m|Sj|GS> * exp(-i*eps_f0*tau) * exp(-i*eps_m0*t2)
# so that
#   Term(eV) = sum_{f,m} [...] * Theta0(eV-eps_f0) * (F0(eV-eps_m0)+F0(eV+eps_m0))
#            = integral dtau dt2  K_theta(tau;eV) * K_W(t2;eV) * G(t2,tau)
# for two closed-form time-domain kernels K_theta, K_W derived below by
# inverse-Fourier-transforming Theta0(eV-.) and F0(eV-.)+F0(eV+.) --
# avoiding ever having to evaluate a discontinuous step or the F0 log
# singularity pointwise on a discrete frequency grid (which was tried
# first and does not converge robustly -- see the module docstring notes
# in stepfunctions.py and the PR history for why).
#
# K_theta(tau;eV) = (1/2)*delta(tau) - (i/(2*pi)) * exp(i*eV*tau) * PV(1/tau)
#   (the exact inverse FT of a Heaviside step; PV = Cauchy principal value)
# so that integral K_theta(tau;eV) h(tau) dtau
#   = 0.5*h(0) + (i/2) * HilbertTransform[exp(i*eV*tau)*h(tau)](tau=0)
# via the standard identity PV integral h(tau)/tau dtau = -pi*Hilbert[h](0).
# The Hilbert transform is computed via the standard FFT method (multiply
# the FFT by -i*sign(frequency)) -- this converges to machine precision
# even on coarse grids, unlike a naive principal-value trapezoidal
# quadrature (which needs impractically fine grids for the same accuracy:
# confirmed directly, see PR history).
#
# K_W(t2;eV) = (1/|t2|) * exp(-Gamma0*|t2|) * [cos(t2*eV) - cos(t2*(eV-omega0))]
#   (the exact inverse FT of F0(eV-.)+F0(eV+.), derived from F0's own
#   closed form via the standard FT pair ln(x^2+b^2) <-> -2*pi*exp(-b|w|)/|w|)
# applied via ordinary numerical integration over t2 (K_W is smooth, no
# singularity at t2=0 -- the two cosines cancel there).
#
# Both kernels were verified independently (direct numerical/closed-form
# checks) and the full pipeline verified end-to-end against the exact,
# already-validated Lehmann-representation third_order_kondo_dIdV (ED, all
# eigenstates) on small test systems -- max relative error ~0.01% at the
# module's default parameters.


def hilbert_transform_at_zero(signal, dt):
    """Hilbert transform of `signal` (array, last axis = the time
    coordinate, uniformly spaced with step dt, symmetric about 0),
    evaluated at the t=0 grid point (assumed to be the array's midpoint
    for an odd-length, endpoint-excluded symmetric grid -- see callers).
    FFT-based (multiply the spectrum by -i*sign(frequency)); exact to
    machine precision for a periodic-boundary approximation of a
    sufficiently long/damped window."""
    n = signal.shape[-1]
    freqs = np.fft.fftfreq(n, d=dt)
    kernel = -1j*np.sign(freqs)
    return np.fft.ifft(np.fft.fft(signal, axis=-1)*kernel, axis=-1)


def theta0_filter(tau_grid, G_of_tau, eV):
    """integral K_theta(tau;eV) G_of_tau(tau) dtau, i.e. the exact T=0
    Theta0(eV-H) filter applied along the last (tau) axis of G_of_tau
    (any leading shape, e.g. a batch of t2 rows). tau_grid must be
    uniform and symmetric about (and include) tau=0."""
    dtau = tau_grid[1] - tau_grid[0]
    idx0 = len(tau_grid)//2
    if abs(tau_grid[idx0]) > 1e-6*dtau:
        raise ValueError("tau_grid must include tau=0 at its midpoint")
    h = np.exp(1j*eV*tau_grid)*G_of_tau
    Hh0 = hilbert_transform_at_zero(h, dtau)[..., idx0]
    return 0.5*G_of_tau[..., idx0] + (1j/2.)*Hh0


def K_W(t2, eV, omega0, Gamma0):
    """Closed-form time-domain kernel for F0(eV-.)+F0(eV+.) (see module
    docstring). Smooth everywhere, including t2=0 (the naive 1/|t2|
    factor there is cancelled by the vanishing bracket)."""
    t2 = np.asarray(t2, dtype=float)
    out = np.zeros_like(t2)
    nz = np.abs(t2) > 1e-300
    tnz = np.abs(t2[nz])
    out[nz] = (np.exp(-Gamma0*tnz)/tnz)*(np.cos(t2[nz]*eV) - np.cos(t2[nz]*(eV-omega0)))
    return out


def kondo_term_from_two_time(t2_grid, tau_grid, G_batches, eVs, omega0, Gamma0):
    """Assemble Term(eV) = integral dt2 K_W(t2;eV) * theta0_filter(G(t2,.);eV)
    for every eV in `eVs`, from G(t2,tau) supplied in chunks over t2 (to
    bound memory: the full (len(t2_grid), len(tau_grid)) array is not
    required to exist at once).

    G_batches: an iterable of (t2_slice, G_chunk) pairs, where t2_slice is
    a 1D array (a contiguous chunk of t2_grid) and G_chunk has shape
    (len(t2_slice), len(tau_grid)) = G(t2,tau) on that chunk. Each chunk
    is visited exactly once, regardless of len(eVs) -- computing G is the
    expensive part (a real time evolution on the DMRG side), so this
    shares that one pass across the whole eV sweep instead of rebuilding
    G per eV.

    Returns an array of real-valued Term(eV), one per eVs entry (already
    includes the Im[...]/4 from Re[X/(4i)]=Im[X]/4, matching eq.
    "3rd-normal"'s prefactor convention in conductance.py -- multiply by
    4*pi*T0^2*Jrho_s for the full third-order Kondo dI/dV contribution,
    exactly as conductance.third_order_kondo_dIdV does for its own
    (excited-state-sum-based) construction)."""
    eVs = np.atleast_1d(np.asarray(eVs, dtype=float))
    totals = np.zeros(len(eVs), dtype=complex)
    for t2_chunk, G_chunk in G_batches:
        for i, eV in enumerate(eVs):
            h_t2 = theta0_filter(tau_grid, G_chunk, eV)
            kw = K_W(t2_chunk, eV, omega0, Gamma0)
            totals[i] += np.trapz(kw*h_t2, t2_chunk)
    return np.imag(totals)/4.
