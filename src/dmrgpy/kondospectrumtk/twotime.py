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
#             = sum_{f,m} X_fm * exp(-i*eps_f0*tau) * exp(-i*eps_m0*t2),
#   X_fm = <GS|Sl|f><f|Sk|m><m|Sj|GS>,
# and the same function on a sheared grid,
#   Gx(s,tau) = G(-s, tau+s)
#             = sum_{f,m} X_fm * exp(-i*eps_f0*tau) * exp(-i*(e_f-e_m)*s),
# so that
#   g(eV) = sum_{f,m} X_fm * Theta0(eV-eps_f0) * (F0(eV-eps_m0) + F0(eV-(e_f-e_m)))
#         = integral dtau dt  K_theta(tau;eV) * K_F(t;eV) * [G(t,tau) + Gx(t,tau)],
# the direct diagram coming from G and the exchange diagram from Gx (see
# conductance.py's module docstring for why the exchange log sits at
# eV-(e_f-e_m), which is where this departs from the paper's eq. 25),
# with the measured term summed over both tunneling directions,
#   Term(eV) = g(eV) + g(-eV)
# (see conductance.py's module docstring for why the Kondo term's two
# directions add rather than cancel, and for the Im[...]/2 normalization)
# for two closed-form time-domain kernels K_theta, K_F derived below by
# inverse-Fourier-transforming Theta0(eV-.) and F0(eV-.) --
# avoiding ever having to evaluate a discontinuous step or the F0 log
# singularity pointwise on a discrete frequency grid (which was tried
# first and does not converge robustly -- see the module docstring notes
# in stepfunctions.py and the PR history for why).
#
# Until 2026-09-26 only G was built, and the kernel along t2 was
# K_W = K_F(eV) + K_F(-eV), the transform of F0(eV-eps_m0)+F0(eV+eps_m0):
# the paper's eq. 25 for the exchange diagram, whose log differs from
# eV-(e_f-e_m) whenever f and GS are not degenerate. A kernel along t2
# cannot express the new argument at all, since e_f-e_m is not a frequency
# of G's t2 axis; it is the frequency of Gx's s axis.
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
# K_F(t;eV) = (1/2) exp(i*eV*t) * { exp(-Gamma0*|t|)/|t|
#                 - (2/(pi*|t|)) [sin(z) Ci(z) - cos(z) si(z)] },  z = omega0*|t|
#   (si(z) = Si(z) - pi/2; the exact inverse FT of F0(eV-.), i.e.
#   integral K_F(t;eV) exp(-i*b*t) dt = F0(eV-b), for
#   F0(x) = ln(omega0+|x|) - 1/2 ln(x^2+Gamma0^2): the braces are the
#   transform of 2*F0, the first term by the standard pair
#   ln(x^2+b^2) <-> -2*pi*exp(-b|t|)/|t|, the second by
#   2*int_0^inf ln(omega0+x) cos(xt) dx = -(2/t) int_0^inf sin(xt)/(omega0+x) dx
#   by parts, Gradshteyn-Ryzhik 3.722.1, and the phase shifts it to eV).
#   The two 1/|t| pieces cancel at t->0, leaving an integrable log,
#   K_F ~ -(omega0/pi) ln(omega0|t|): F0 decays only as omega0/|x| beyond
#   the band, so its transform is not finite at t=0. On the uniform t
#   grid, which contains t=0 exactly, the points next to it take the
#   average of K_F over their own cell (K_F's `dt` argument), which is
#   what a Riemann sum wants of an integrable singularity; every other
#   point is evaluated directly.
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


def K_F(t, eV, omega0, Gamma0, dt=None, n_avg=8):
    """Closed-form time-domain kernel of F0(eV-.), complex (see module
    docstring): integral K_F(t;eV) exp(-i*b*t) dt = F0(eV-b). Applied to
    G's t2 axis it gives the direct diagram's F0(eV-eps_m0), and to Gx's s
    axis the exchange diagram's F0(eV-(e_f-e_m)). Log-singular
    (integrably) at t=0, so on a uniform grid of spacing `dt` the points
    within n_avg cells of it (t=0 included) are returned as cell
    averages, (1/dt) int_{t-dt/2}^{t+dt/2} K_F -- what a Riemann sum wants
    of a kernel that is not smooth on the cell scale there (measured on a
    pure cosine: the plain midpoint values left a 0.9% error at
    omega0*dt=0.5, the cell averages 0.02%). A t array containing an
    exact 0 requires dt."""
    from scipy.special import sici
    from scipy.integrate import quad
    t = np.asarray(t, dtype=float)
    out = np.zeros(t.shape, dtype=complex)
    if dt is not None:
        near = np.abs(t) < (n_avg + 0.5)*dt
        if near.any():
            far = ~near
            out[far] = K_F(t[far], eV, omega0, Gamma0)
            for i in np.nonzero(near)[0]:
                a, b = t[i] - dt/2., t[i] + dt/2.
                pts = [0.] if a < 0. < b else None
                re, _ = quad(lambda u: K_F(np.array([u]), eV, omega0, Gamma0)[0].real,
                             a, b, points=pts, limit=200)
                im, _ = quad(lambda u: K_F(np.array([u]), eV, omega0, Gamma0)[0].imag,
                             a, b, points=pts, limit=200)
                out[i] = (re + 1j*im)/dt
            return out
    nz = np.abs(t) > 1e-300
    ta = np.abs(t[nz])
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
    out[nz] = 0.5*np.exp(1j*eV*t[nz])*(main - band)
    if not nz.all():
        raise ValueError("K_F diverges (logarithmically) at t=0: pass "
                         "dt, the grid spacing, to get its cell average")
    return out


def kondo_term_from_two_time(t2_grid, tau_grid, G_batches, eVs, omega0, Gamma0):
    """Assemble
        Term(eV) = integral dt K_F(t;eV) * [theta0_filter(G(t,.);eV)
                                            + theta0_filter(Gx(t,.);eV)]
    for every eV in `eVs`, from G(t2,tau) and Gx(s,tau) supplied in
    chunks over the time t (t2 for G, s for Gx, on the one grid t2_grid;
    chunked to bound memory: the full (len(t2_grid), len(tau_grid))
    arrays are not required to exist at once).

    G_batches: an iterable of (t_slice, G_chunk, Gx_chunk) triples, where
    t_slice is a 1D array (a contiguous chunk of t2_grid) and G_chunk,
    Gx_chunk have shape (len(t_slice), len(tau_grid)): G(t2,tau) and
    Gx(s,tau) = G(-s,tau+s) on that chunk (see the module docstring; G
    carries the direct diagram, Gx the exchange one). A plain (t_slice,
    G_chunk) pair is refused, since G alone cannot produce the exchange
    diagram. Each chunk
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
    totals = np.zeros(2*nev, dtype=complex) # t->s, then s->t
    dt2 = t2_grid[1] - t2_grid[0]
    t2_first, t2_last = t2_grid[0], t2_grid[-1]
    for batch in G_batches:
        if len(batch) != 3:
            raise ValueError(
                "G_batches must yield (t_slice, G_chunk, Gx_chunk) triples: "
                "the exchange diagram's log sits at eV-(e_f-e_m), which "
                "needs Gx(s,tau)=G(-s,tau+s) (see twotime.py's module "
                "docstring); G alone gives only the paper's eq. 25 form")
        t2_chunk, G_chunk, Gx_chunk = batch
        weights = np.full(len(t2_chunk), dt2)
        weights[np.isclose(t2_chunk, t2_first)] *= 0.5
        weights[np.isclose(t2_chunk, t2_last)] *= 0.5
        # the Theta0 filter is linear and both diagrams share the t kernel,
        # so filter the sum once rather than each diagram separately
        GG = np.asarray(G_chunk) + np.asarray(Gx_chunk)
        for i, eV in enumerate(eVs):
            # K_F(t;-eV) = conj(K_F(t;eV)) (the braces are real), so the
            # s->t direction reuses the t->s kernel: its cell averages next
            # to t=0 are adaptive quadratures of an oscillating log and
            # were most of the cost
            kf = K_F(t2_chunk, eV, omega0, Gamma0, dt=dt2)
            totals[i] += np.sum(kf*theta0_filter(tau_grid, GG, eV)*weights)
            totals[nev+i] += np.sum(np.conj(kf)*theta0_filter(tau_grid, GG, -eV)
                                    *weights)
    out = np.imag(totals)/2. # SA factor 2 -- see conductance.py's docstring
    return out[:nev] + out[nev:]
