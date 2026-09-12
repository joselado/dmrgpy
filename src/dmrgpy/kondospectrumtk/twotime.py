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
#   g(eV) = sum_{f,m} [...] * Theta0(eV-eps_f0) * (F0(eV-eps_m0)+F0(eV+eps_m0))
#         = integral dtau dt2  K_theta(tau;eV) * K_W(t2;eV) * G(t2,tau)
# with the measured term summed over both tunneling directions,
#   Term(eV) = g(eV) + g(-eV)
# (see conductance.py's module docstring for why the Kondo term's two
# directions add rather than cancel, and for the Im[...]/2 normalization)
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
# K_W(t2;eV) = cos(eV*t2) * { exp(-Gamma0*|t2|)/|t2|
#                 - (2/(pi*|t2|)) [sin(z) Ci(z) - cos(z) si(z)] },  z = omega0*|t2|
#   (si(z) = Si(z) - pi/2; the exact inverse FT of F0(eV-.)+F0(eV+.) for
#   F0(x) = ln(omega0+|x|) - 1/2 ln(x^2+Gamma0^2): the first term is the
#   standard pair ln(x^2+b^2) <-> -2*pi*exp(-b|t|)/|t|, the second is
#   2*int_0^inf ln(omega0+x) cos(xt) dx = -(2/t) int_0^inf sin(xt)/(omega0+x) dx
#   by parts, Gradshteyn-Ryzhik 3.722.1). The two 1/|t2| pieces cancel at
#   t2->0, leaving an integrable log, K_W ~ -(2*omega0/pi) ln(omega0|t2|):
#   F0 decays only as omega0/|x| beyond the band, so its transform is not
#   finite at t2=0. On the uniform t2 grid, which contains t2=0 exactly,
#   that one point takes the average of K_W over its own cell (K_W's `dt`
#   argument), which is what a Riemann sum wants of an integrable
#   singularity; every other point is evaluated directly.
#   Until 2026-09-12 F0 was the electron-like ln|omega0-x| - ln|x| (see
#   stepfunctions.py's module docstring), whose transform
#   (1/|t2|) exp(-Gamma0|t2|) [cos(t2 eV) - cos(t2 (eV-omega0))] is
#   smooth at t2=0 -- the band-edge cosine there is the sharp-cutoff
#   singularity at x=omega0 that the current F0 does not have.
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


def K_W(t2, eV, omega0, Gamma0, dt=None, n_avg=8):
    """Closed-form time-domain kernel for F0(eV-.)+F0(eV+.) (see module
    docstring). Log-singular (integrably) at t2=0, so on a uniform grid
    of spacing `dt` the points within n_avg cells of it (t2=0 included)
    are returned as cell averages, (1/dt) int_{t2-dt/2}^{t2+dt/2} K_W --
    what a Riemann sum wants of a kernel that is not smooth on the cell
    scale there (measured on a pure cosine: the plain midpoint values
    left a 0.9% error at omega0*dt=0.5, the cell averages 0.02%). A t2
    array containing an exact 0 requires dt."""
    from scipy.special import sici
    from scipy.integrate import quad
    t2 = np.asarray(t2, dtype=float)
    out = np.zeros_like(t2)
    if dt is not None:
        near = np.abs(t2) < (n_avg + 0.5)*dt
        if near.any():
            far = ~near
            out[far] = K_W(t2[far], eV, omega0, Gamma0)
            for i in np.nonzero(near)[0]:
                a, b = t2[i] - dt/2., t2[i] + dt/2.
                pts = [0.] if a < 0. < b else None
                cell, _ = quad(lambda t: K_W(np.array([t]), eV, omega0, Gamma0)[0],
                               a, b, points=pts, limit=200)
                out[i] = cell/dt
            return out
    nz = np.abs(t2) > 1e-300
    ta = np.abs(t2[nz])
    z = omega0*ta
    small = z < 1e-4
    band = np.empty_like(ta)
    si, ci = sici(z[~small])
    band[~small] = (2./(np.pi*ta[~small]))*(np.sin(z[~small])*ci
                                           - np.cos(z[~small])*(si - np.pi/2))
    # z->0: the 1/t of the band term cancels the 1/t of the Gamma0 term
    # exactly; evaluate the difference by its series instead of by
    # cancellation (gamma = Euler's constant)
    gamma = 0.5772156649015329
    ts = ta[small]
    band[small] = 1./ts + (2*omega0/np.pi)*(gamma - 1. + np.log(omega0*ts)) - omega0**2*ts/2.
    main = np.exp(-Gamma0*ta)/ta
    main[small] = 1./ts - Gamma0 + Gamma0**2*ts/2.
    out[nz] = np.cos(t2[nz]*eV)*(main - band)
    if not nz.all():
        raise ValueError("K_W diverges (logarithmically) at t2=0: pass "
                         "dt, the grid spacing, to get its cell average")
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
    G per eV. Chunks may be as small as a single t2 point (as they are on
    the DMRG side, where every t2 checkpoint is its own real trajectory) --
    the t2 integral is done as a uniform Riemann sum weighted by t2_grid's
    own spacing (with the standard trapezoidal half-weight only at the two
    true endpoints of the *full* t2_grid, not per chunk), which is well
    defined for any chunk size; a per-chunk np.trapz (tried first) silently
    returns exactly 0 for every single-point chunk -- confirmed directly,
    it produced an all-zero result end to end -- since a trapezoidal rule
    needs at least two points to have any width to integrate over.

    Returns an array of real-valued Term(eV), one per eVs entry (already
    includes the Im[...]/2 spin-average normalization and the sum over
    both tunneling directions, matching conductance.py's own conventions
    exactly -- multiply by 4*pi*T0^2*Jrho_s for the full third-order
    Kondo dI/dV contribution, as
    conductance.third_order_kondo_dIdV does for its own
    (excited-state-sum-based) construction).

    Both tunneling directions are handled here rather than by the callers
    so that edtwotimeref.py, dmrgtwotime.py and Spin_Chain's mode="DMRG"
    path all inherit them from one place. Only the two closed-form
    kernels depend on eV, so the s->t direction is obtained by evaluating
    them at -eV as well (eq. "sym_z": the Kondo term's two directions add,
    giving a result even in eV) -- G(t2,tau), the expensive part, is still
    built in a single pass over G_batches regardless of how many eV
    points are swept."""
    eVs = np.atleast_1d(np.asarray(eVs, dtype=float))
    nev = len(eVs)
    both = np.concatenate([eVs, -eVs]) # t->s, then s->t
    totals = np.zeros(len(both), dtype=complex)
    dt2 = t2_grid[1] - t2_grid[0]
    t2_first, t2_last = t2_grid[0], t2_grid[-1]
    for t2_chunk, G_chunk in G_batches:
        weights = np.full(len(t2_chunk), dt2)
        weights[np.isclose(t2_chunk, t2_first)] *= 0.5
        weights[np.isclose(t2_chunk, t2_last)] *= 0.5
        for i, eV in enumerate(both):
            h_t2 = theta0_filter(tau_grid, G_chunk, eV)
            kw = K_W(t2_chunk, eV, omega0, Gamma0, dt=dt2)
            totals[i] += np.sum(kw*h_t2*weights)
    out = np.imag(totals)/2. # SA factor 2 -- see conductance.py's docstring
    return out[:nev] + out[nev:]
