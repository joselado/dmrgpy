"""The leading-log flow of leadinglog._joint_flow on a device, through JAX.

leadinglog._joint_flow takes this route when the array backend of pyitensor
is JAX (pyitensor.backend.set_backend("jax")), and the flow then runs on a
GPU when JAX has one. The equations, the layout of the state and the
integrator are those of the NumPy route: the same right-hand side written
in jax.numpy, and the same DOP853 steps (_Dop853 is SciPy's RungeKutta and
DOP853 with its coefficients, its initial step, its error norm and its step
control, with the state on the device). So the two routes give the same
numbers to roundoff, and the same steps unless one sits at the edge of
acceptance.

What runs where. Every array of the size of the state stays on the device.
The host keeps the scalars of the integrator, builds the tables of the
thermal kernel at each cutoff (thermalkernel.ThermalKernel.table, a few
tens of kB), and receives one scalar per stage attempt (the error norm) and
one per step (the sample coupling). The tip part of the right-hand side is
computed per chunk of biases (and of occupied columns), with each chunk
below CHUNK_BYTES, in one of two forms (TIP_FORM, _DeviceFlow._tip): as
leadinglog._tip_bracket does for a few occupied columns, and as one matrix
product per bias for many.

Memory. The integrator holds the 13 stages of a step, the state at its
start and the new state, and one stage argument, about 16 copies of the
state where SciPy's DOP853 holds about 24 (leadinglog.memory_model). At a
stop of the flow three more stages are computed for the dense output.

What it is for. The d^3 part of the flow is four 2d x 2d products per
right-hand side, which a GPU with full double precision does far faster
than a CPU. On a card with weak double precision the products are not
faster than a few cores, and what the device gains there is the rest of the
right-hand side and the integrator, which run on one core in the NumPy
route. README, section 7, has the measurements.
"""
import warnings

import numpy as np
from scipy.integrate._ivp import dop853_coefficients as _dc
from scipy.optimize import brentq

from ..pyitensor import backend

_EPS = np.finfo(float).eps

# bytes of one chunk of the tip part (the operands of its products and
# their kernel)
CHUNK_BYTES = 512e6

# how the tip part is computed (_DeviceFlow._tip): "columns" sums over the
# electron index first, as leadinglog._tip_bracket does, with an operand of
# 4 d^2 numbers per occupied column and the exchange kernel built again for
# every chunk of columns; "block" applies the kernel to the vertex and does
# one product per bias for all the columns, with twice the arithmetic.
# "auto" takes "columns" up to BLOCK_FROM - 1 occupied columns.
TIP_FORM = "auto"
BLOCK_FROM = 4

# scipy.integrate._ivp.rk
SAFETY, MIN_FACTOR, MAX_FACTOR = 0.9, 0.2, 10.
_NS = _dc.N_STAGES # 12
_ERROR_EXPONENT = -1./(7 + 1) # error_estimator_order = 7


def _nonzero(c, K):
    """The pairs (c_j, K_j) with c_j != 0, as an array and a tuple; one
    pair at least (the dense output at x=0 has every c_j zero)."""
    idx = [j for j in range(len(K)) if c[j] != 0.] or [0]
    return np.array([c[j] for j in idx]), tuple(K[j] for j in idx)


class _Kernels():
    """The jitted linear algebra of the integrator, shared by every flow."""
    _built = None

    @classmethod
    def get(cls):
        if cls._built is None:
            cls._built = cls()
        return cls._built

    def __init__(self):
        jax, jnp = backend._JAX, backend.xp()

        def lincomb(c, K):
            acc = c[0]*K[0]
            for j in range(1, len(K)):
                acc = acc + c[j]*K[j]
            return acc

        # y + (sum_j c_j K_j)*h, as np.dot(K.T, c)*h in rk_step
        self.comb = jax.jit(lambda y, h, c, K: y + lincomb(c, K)*h)

        def errsq(y, y_new, rtol, atol, c5, K5, c3, K3):
            scale = atol + jnp.maximum(jnp.abs(y), jnp.abs(y_new))*rtol
            e5 = lincomb(c5, K5)/scale
            e3 = lincomb(c3, K3)/scale
            return jnp.sum(e5*e5), jnp.sum(e3*e3)
        self.errsq = jax.jit(errsq)

        def scaled_sq(v, y0, rtol, atol):
            x = v/(atol + jnp.abs(y0)*rtol)
            return jnp.sum(x*x)
        self.scaled_sq = jax.jit(scaled_sq)
        self.diff_sq = jax.jit(lambda a, b, y0, rtol, atol: scaled_sq(a - b, y0, rtol, atol))
        self.dense = jax.jit(lambda y_old, cy, dy, c, K: y_old + cy*dy + lincomb(c, K))
        self.sub = jax.jit(lambda a, b: a - b)


class _Dop853():
    """scipy.integrate.DOP853 with the state on the device. fun(t, y)
    returns a device array; the scalars stay on the host."""

    def __init__(self, fun, t0, y0, t_bound, rtol, atol):
        if rtol < 100*_EPS:
            warnings.warn("At least one element of `rtol` is too small. Setting "
                          "`rtol = np.maximum(rtol, %g)`." % (100*_EPS))
            rtol = 100*_EPS
        self.kern = _Kernels.get()
        self.nfev = 0
        self._fun = fun
        self.t, self.y, self.t_bound = float(t0), y0, float(t_bound)
        self.rtol, self.atol = float(rtol), float(atol)
        self.direction = float(np.sign(t_bound - t0)) if t_bound != t0 else 1.
        self.n = int(y0.size)
        self.f = self.fun(self.t, y0)
        self.h_abs = self._initial_step()
        self.t_old = self.y_old = self.h_previous = self.K = None
        self.status = "running"

    def fun(self, t, y):
        self.nfev += 1
        return self._fun(t, y)

    def _norm(self, sq):
        return float(np.sqrt(sq))/self.n**0.5

    def _initial_step(self):
        """scipy.integrate._ivp.common.select_initial_step, max_step=inf"""
        k = self.kern
        interval = abs(self.t_bound - self.t)
        if interval == 0.:
            return 0.
        d0 = self._norm(k.scaled_sq(self.y, self.y, self.rtol, self.atol))
        d1 = self._norm(k.scaled_sq(self.f, self.y, self.rtol, self.atol))
        h0 = 1e-6 if (d0 < 1e-5 or d1 < 1e-5) else 0.01*d0/d1
        h0 = min(h0, interval)
        y1 = k.comb(self.y, h0*self.direction, np.ones(1), (self.f,))
        f1 = self.fun(self.t + h0*self.direction, y1)
        d2 = self._norm(k.diff_sq(f1, self.f, self.y, self.rtol, self.atol))/h0
        del y1, f1
        if d1 <= 1e-15 and d2 <= 1e-15:
            h1 = max(1e-6, h0*1e-3)
        else:
            h1 = (0.01/max(d1, d2))**(1./(7 + 1))
        return min(100*h0, h1, interval)

    def step(self):
        """RungeKutta._step_impl and OdeSolver.step"""
        k = self.kern
        t, y = self.t, self.y
        # drop what the last step kept for its dense output
        self.y_old = self.K = None
        min_step = 10*np.abs(np.nextafter(t, self.direction*np.inf) - t)
        h_abs = max(self.h_abs, min_step)
        accepted = rejected = False
        while not accepted:
            if h_abs < min_step:
                self.status = "failed"
                return "Required step size is less than spacing between numbers."
            h = h_abs*self.direction
            t_new = t + h
            if self.direction*(t_new - self.t_bound) > 0:
                t_new = self.t_bound
            h = t_new - t
            h_abs = np.abs(h)
            K = [self.f]
            for s in range(1, _NS):
                c, Ks = _nonzero(_dc.A[s, :s], K)
                K.append(self.fun(t + _dc.C[s]*h, k.comb(y, h, c, Ks)))
            c, Ks = _nonzero(_dc.B, K)
            y_new = k.comb(y, h, c, Ks)
            K.append(self.fun(t + h, y_new))
            c5, K5 = _nonzero(_dc.E5, K)
            c3, K3 = _nonzero(_dc.E3, K)
            e5, e3 = (float(v) for v in k.errsq(y, y_new, self.rtol, self.atol,
                                                 c5, K5, c3, K3))
            if e5 == 0. and e3 == 0.:
                error_norm = 0.
            else:
                error_norm = np.abs(h)*e5/np.sqrt((e5 + 0.01*e3)*self.n)
            if error_norm < 1:
                if error_norm == 0:
                    factor = MAX_FACTOR
                else:
                    factor = min(MAX_FACTOR, SAFETY*error_norm**_ERROR_EXPONENT)
                if rejected:
                    factor = min(1, factor)
                h_abs *= factor
                accepted = True
            else:
                h_abs *= max(MIN_FACTOR, SAFETY*error_norm**_ERROR_EXPONENT)
                rejected = True
        self.h_previous = h
        self.y_old, self.K = y, K
        self.t_old, self.t, self.y = t, t_new, y_new
        self.h_abs = h_abs
        self.f = K[-1]
        if self.direction*(self.t - self.t_bound) >= 0:
            self.status = "finished"
        return None

    def dense_output(self):
        """DOP853._dense_output_impl and Dop853DenseOutput as one function
        of t. The interpolant is written as a combination of y_old, y and
        the 16 stages, so that F is never stored."""
        k = self.kern
        K = list(self.K)
        h = self.h_previous
        for s in range(_NS + 1, _dc.N_STAGES_EXTENDED):
            c, Ks = _nonzero(_dc.A[s, :s], K)
            K.append(self.fun(self.t_old + _dc.C[s]*h, k.comb(self.y_old, h, c, Ks)))
        dy = k.sub(self.y, self.y_old)
        t_old, y_old = self.t_old, self.y_old

        def at(t):
            x = (t - t_old)/h
            # Dop853DenseOutput: y_old + sum_r p_r F_r, with p_r below
            p = np.array([x, x*(1 - x), x**2*(1 - x), x**2*(1 - x)**2,
                          x**3*(1 - x)**2, x**3*(1 - x)**3, x**4*(1 - x)**3])
            # F0 = dy, F1 = h K0 - dy, F2 = 2 dy - h (f + K0), F3.. = h D K
            cK = h*(p[3:] @ _dc.D)
            cK[0] += h*(p[1] - p[2])
            cK[_NS] -= h*p[2] # f = K[12]
            c, Ks = _nonzero(cK, K)
            return k.dense(y_old, p[0] - p[1] + 2*p[2], dy, c, Ks)
        return at


def _integrate(solver, sample_coupling, ref):
    """leadinglog._integrate on the device solver."""
    g_old = sample_coupling(solver.y)
    max_coupling = g_old
    y, D_stop = solver.y, None
    while solver.status == "running":
        message = solver.step()
        if solver.status == "failed":
            raise RuntimeError("leading-log flow failed: "+str(message))
        y = solver.y
        g_new = sample_coupling(y)
        if g_old <= ref and g_new > ref:
            dense = solver.dense_output()
            s_stop = brentq(lambda s: ref - sample_coupling(dense(s)),
                            solver.t_old, solver.t, xtol=4*_EPS, rtol=4*_EPS)
            y = dense(s_stop)
            D_stop = float(np.exp(s_stop))
            max_coupling = max(max_coupling, sample_coupling(y))
            break
        max_coupling = max(max_coupling, g_new)
        g_old = g_new
    return y, D_stop, max_coupling


def _thermal_static(k):
    """The scalars of a ThermalKernel that its tables are read with."""
    return (float(k.h), int(k.nn), float(k.omega0), float(k.Gamma0)**2,
            float(k.t0), float(k.iht), int(k.nf))


def _tables(k, D, Dlow=-1.):
    """The tables of a ThermalKernel at a cutoff, on the device; () at T=0."""
    if k is None:
        return ()
    cn, cf = k.table(D, Dlow)
    jnp = backend.xp()
    return (jnp.asarray(cn), jnp.asarray(cf))


class _DeviceFlow():
    """The right-hand side of leadinglog._joint_flow in jax.numpy, for one
    flow: the constants on the device and the jitted functions."""

    def __init__(self, e, w, cols, omega0, Gamma0, real, kernel, onshell_kernel):
        jax, jnp = backend._JAX, backend.xp()
        self.jnp = jnp
        d, n, m = len(e), len(w), len(cols)
        self.d, self.n, self.m, self.real = d, n, m, real
        self.omega0, self.Gamma0 = float(omega0), float(Gamma0)
        self.kernel, self.onshell_kernel = kernel, onshell_kernel
        self.st_tip = None if kernel is None else _thermal_static(kernel)
        self.st_on = None if onshell_kernel is None else _thermal_static(onshell_kernel)
        item = 8 if real else 16
        form = TIP_FORM if TIP_FORM != "auto" else ("block" if m >= BLOCK_FROM else "columns")
        if form not in ("block", "columns"):
            raise ValueError("TIP_FORM must be \"auto\", \"block\" or \"columns\", got %r"
                             % (TIP_FORM,))
        self.form = form

        def chunks(k):
            nk = -(-n//k)
            return jnp.asarray(np.concatenate([w, np.full(nk*k - n, w[-1])]).reshape(nk, k))
        e = np.asarray(e, dtype=float)
        cols = np.asarray(cols, dtype=np.int64)
        base = (jnp.asarray(e), jnp.asarray(w), jnp.asarray(e[:, None] - e[None, :]))
        if form == "block":
            # biases per chunk: of the direct sum, where the operand and the
            # result are 2d x 2m per bias, and of the exchange sum, where
            # the operand is 2d x 2d per bias and its kernel d x d
            self.kd = int(max(1, min(n, CHUNK_BYTES//(3*4*d*m*item))))
            self.kc = int(max(1, min(n, CHUNK_BYTES//((4*item + 8)*d*d))))
            # the columns (n,b) of the block form, n = 0, 1, for b in cols
            colsel = np.concatenate([cols, d + cols])
            self.consts = base + (jnp.asarray(e[cols]), jnp.asarray(colsel), chunks(self.kd),
                                  chunks(self.kc))
        else:
            # occupied columns per chunk: P and Q are 4 mc d^2 numbers each;
            # biases per chunk of the exchange sum: a d x kc x d kernel
            self.mc = int(max(1, min(m, CHUNK_BYTES//(2*4*d*d*item))))
            nch = -(-m//self.mc)
            colp = np.concatenate([cols, np.full(nch*self.mc - m, cols[-1])])
            self.kc = int(max(1, min(n, CHUNK_BYTES//(8*d*d))))
            self.consts = base + (jnp.asarray(colp.reshape(nch, self.mc)), chunks(self.kc))
        nb = 4*d*d
        self.nb = nb

        def rhs(D, y, tab_on, tab_tip, consts):
            return self._rhs(D, y, tab_on, tab_tip, consts)
        self.rhs = jax.jit(rhs)
        self.tip = jax.jit(lambda Gb, Ab, D, Dlow, tab, consts:
                           self._tip(Gb, Ab, D, Dlow, tab, consts, True))
        if real:
            def coupling(y):
                B = y[:nb].reshape(2*d, 2*d)
                return jnp.maximum(jnp.maximum(
                    jnp.max(jnp.abs((B[:d, d:] + B[d:, :d])/2.)),
                    jnp.max(jnp.abs((B[:d, d:] - B[d:, :d])/2.))),
                    jnp.max(jnp.abs((B[:d, :d] - B[d:, d:])/2.)))
        else:
            def coupling(y):
                X = self._complex(y[:2*nb]).reshape(4, d, d)
                return jnp.max(jnp.abs(X[1:]))
        self.coupling = jax.jit(coupling)

    # the kernels, line for line those of leadinglog.py and thermalkernel.py
    def _k_two(self, x, D):
        """leadinglog._kernel"""
        jnp = self.jnp
        ax = jnp.abs(x)
        return 1./(D + jnp.sqrt(ax*ax + self.Gamma0*self.Gamma0)) - 1./(D + self.omega0 + ax)

    def _k_int(self, x, Dlow, Dhigh):
        """leadinglog._kernel_integral"""
        jnp = self.jnp
        ax = jnp.abs(x)
        r = jnp.sqrt(ax*ax + self.Gamma0*self.Gamma0)
        return (jnp.log((Dhigh + r)/(Dlow + r))
                - jnp.log((Dhigh + self.omega0 + ax)/(Dlow + self.omega0 + ax)))

    def _k_one(self, ax, D, Dlow, omega0, G2, integral):
        """the kernel of leadinglog._hole_sum and thermalkernel._k0"""
        jnp = self.jnp
        r = jnp.sqrt(ax*ax + G2)
        if integral:
            return jnp.log((D + r)/(Dlow + r)) - jnp.log((D + omega0 + ax)/(Dlow + omega0 + ax))
        return (omega0 + ax - r)/((D + r)*(D + omega0 + ax))

    def _kval(self, ax, D, Dlow, tab, st, integral):
        """thermalkernel._kval on an array"""
        jnp = self.jnp
        cn, cf = tab
        h, nn, omega0, G2, t0, iht, nf = st
        near = ax < h*nn
        s = jnp.where(near, ax/h, 0.)
        j = jnp.minimum(s.astype(jnp.int64), nn - 1)
        ds = s - j
        vn = ((cn[j, 0]*ds + cn[j, 1])*ds + cn[j, 2])*ds + cn[j, 3]
        axf = jnp.where(near, h*nn, ax)
        k = self._k_one(axf, D, Dlow, omega0, G2, integral)
        s2 = (jnp.log(axf) - t0)*iht
        j2 = s2.astype(jnp.int64)
        inside = j2 < nf
        j2 = jnp.minimum(j2, nf - 1)
        ds2 = s2 - j2
        vf = jnp.where(inside,
                       k*(1. + ((cf[j2, 0]*ds2 + cf[j2, 1])*ds2 + cf[j2, 2])*ds2 + cf[j2, 3]),
                       k)
        return jnp.where(near, vn, vf)

    def _kdir(self, x, D, Dlow, tab, integral):
        """the kernel of the direct sum of the tip vertex"""
        if self.st_tip is not None:
            return self._kval(self.jnp.abs(x), D, Dlow, tab, self.st_tip, integral)
        return self._k_int(x, Dlow, D) if integral else self._k_two(x, D)

    def _khole(self, ax, D, Dlow, tab, integral):
        """the kernel of the exchange sum of the tip vertex"""
        if self.st_tip is not None:
            return self._kval(ax, D, Dlow, tab, self.st_tip, integral)
        return self._k_one(ax, D, Dlow, self.omega0, self.Gamma0*self.Gamma0, integral)

    def _kon(self, dE, D, tab):
        """the kernel of the on-shell flow at the energy differences"""
        if self.st_on is not None:
            return self._kval(self.jnp.abs(dE), D, -1., tab, self.st_on, False)
        return self._k_two(dE, D)

    # block form, as leadinglog._block, _unblock, _swap and _bracket
    def _block(self, A):
        jnp = self.jnp
        return jnp.concatenate([
            jnp.concatenate([A[0] + A[3], A[1] - 1j*A[2]], axis=1),
            jnp.concatenate([A[1] + 1j*A[2], A[0] - A[3]], axis=1)], axis=0)

    def _unblock(self, B):
        d = self.d
        return self.jnp.stack([(B[:d, :d] + B[d:, d:])/2., (B[:d, d:] + B[d:, :d])/2.,
                               0.5j*(B[:d, d:] - B[d:, :d]), (B[:d, :d] - B[d:, d:])/2.])

    def _swap(self, B):
        jnp, d = self.jnp, self.d
        return jnp.concatenate([jnp.concatenate([B[:d, :d], B[d:, :d]], axis=1),
                                jnp.concatenate([B[:d, d:], B[d:, d:]], axis=1)], axis=0)

    def _bracket(self, X, Y):
        return X @ Y - self._swap(self._swap(Y) @ self._swap(X))

    def _complex(self, y):
        """the complex numbers of a float view, as ndarray.view(complex)"""
        y2 = y.reshape(-1, 2)
        return backend._JAX.lax.complex(y2[:, 0], y2[:, 1])

    def _floats(self, z):
        """the float view of complex numbers, as ndarray.view(float)"""
        return self.jnp.stack([z.real, z.imag], axis=-1).reshape(-1)

    def _tip(self, Gb, Ab, D, Dlow, tab, consts, integral):
        """leadinglog._tip_bracket: (n,2,d,2,m), the tip vertex at each
        bias with its columns cut, in the form self.form (TIP_FORM)."""
        if self.form == "block":
            return self._tip_block(Gb, Ab, D, Dlow, tab, consts, integral)
        return self._tip_columns(Gb, Ab, D, Dlow, tab, consts, integral)

    def _tip_block(self, Gb, Ab, D, Dlow, tab, consts, integral):
        """The kernel of either sum depends on the impurity indices only, so
        it multiplies a vertex elementwise, the same for the four electron
        blocks. The direct sum at a bias w is G (A_cols * K) with
        K[c,b] = K(w-(e_c-e_b)), and the exchange sum is
        swap((swap(A) * K') swap(G)_cols) with K'[a,c] = K(w-(e_a-e_c)),
        which is _bracket's hole channel with the kernel on A. Both are
        matrix products per bias, the exchange kernel is built once per
        bias for every occupied column, and nothing grows with their
        number but the products."""
        jnp, lax = self.jnp, backend._JAX.lax
        d, n, m = self.d, self.n, self.m
        e, ecol, colsel, wdir, wexc = consts[0], consts[3], consts[4], consts[5], consts[6]
        Ac = Ab[:, colsel] # (2d, 2m)
        sG = self._swap(Gb)[:, colsel]
        sA = self._swap(Ab)

        def direct(wk):
            # x[k,c,j] = w_k - (e_c - e_cols_j), as in _tip_bracket
            x = wk[:, None, None] - (e[None, :, None] - ecol[None, None, :])
            K = jnp.tile(self._kdir(x, D, Dlow, tab, integral), (1, 2, 2))
            Y = jnp.moveaxis(Ac[None]*K, 0, 1).reshape(2*d, -1) # (2d, kd 2m)
            return jnp.moveaxis((Gb @ Y).reshape(2*d, -1, 2*m), 1, 0)

        def exchange(wk):
            # x[k,a,c] = (w_k - e_a) + e_c, as in _hole_sum
            ax = jnp.abs((wk[:, None, None] - e[None, :, None]) + e[None, None, :])
            K = jnp.tile(self._khole(ax, D, Dlow, tab, integral), (1, 2, 2))
            X = (sA[None]*K).reshape(-1, 2*d) # (kc 2d, 2d)
            return (X @ sG).reshape(-1, 2*d, 2*m) # rows (n,a), columns (i,j)

        def over(f, ws):
            out = f(ws[0]) if ws.shape[0] == 1 else lax.map(f, ws)
            return out.reshape(-1, 2*d, 2*m)[:n]
        part = over(direct, wdir).reshape(n, 2, d, 2, m)
        hole = over(exchange, wexc).reshape(n, 2, d, 2, m)
        return part - hole.transpose(0, 3, 2, 1, 4)

    def _hole(self, Q, D, Dlow, tab, integral, consts):
        """out[k,a,l] = sum_c Q[a,l,c] K(w[k]-(e[a]-e[c])), for a real Q,
        over the chunks of biases: _hole_sum of leadinglog.py."""
        jnp, lax = self.jnp, backend._JAX.lax
        e, wchunks = consts[0], consts[4]

        def chunk(wk):
            # [a,k,c], with x = (w_k - e_a) + e_c as in _hole_sum
            ax = jnp.abs((wk[None, :, None] - e[:, None, None]) + e[None, None, :])
            K = self._khole(ax, D, Dlow, tab, integral)
            return jnp.einsum('akc,alc->kal', K, Q)
        out = chunk(wchunks[0]) if wchunks.shape[0] == 1 else lax.map(chunk, wchunks)
        return out.reshape((-1,) + out.shape[-2:])[:self.n]

    def _tip_columns(self, Gb, Ab, D, Dlow, tab, consts, integral):
        """leadinglog._tip_bracket line for line, per chunk of occupied
        columns: the electron index is summed before the kernel is applied,
        so the products are half those of _tip_block, but each chunk of
        columns has its operands P and Q of 4 mc d^2 numbers and builds the
        exchange kernel again."""
        jnp, lax = self.jnp, backend._JAX.lax
        d, n, mc = self.d, self.n, self.mc
        e, w, colchunks = consts[0], consts[1], consts[3]
        G4 = Gb.reshape(2, d, 2, d)
        A4 = Ab.reshape(2, d, 2, d)
        cplx = not self.real

        def chunk(cj):
            x = w[None, None, :] - (e[None, :, None] - e[cj][:, None, None]) # [j,c,k]
            Kp = self._kdir(x, D, Dlow, tab, integral)
            P = jnp.einsum('iakc,kcnj->jianc', G4, A4[:, :, :, cj]).reshape(mc, 4*d, d)
            Q = jnp.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cj]).reshape(d, 4*mc, d)
            if not cplx:
                part = jnp.matmul(P, Kp).reshape(mc, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
                hole = self._hole(Q, D, Dlow, tab, integral, consts)
                return part - hole.reshape(n, d, 2, 2, mc).transpose(0, 2, 1, 3, 4)
            part = jnp.matmul(P.real, Kp) + 1j*jnp.matmul(P.imag, Kp)
            part = part.reshape(mc, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
            hole = self._hole(jnp.concatenate([Q.real, Q.imag], axis=1), D, Dlow, tab,
                              integral, consts)
            hole = (hole[:, :, :4*mc] + 1j*hole[:, :, 4*mc:]).reshape(n, d, 2, 2, mc)
            return part - hole.transpose(0, 2, 1, 3, 4)
        if colchunks.shape[0] == 1:
            return chunk(colchunks[0])
        out = lax.map(chunk, colchunks) # (nch, n, 2, d, 2, mc)
        out = jnp.moveaxis(out, 0, 4).reshape(n, 2, d, 2, -1)
        return out[..., :self.m]

    def _rhs(self, D, y, tab_on, tab_tip, consts):
        """the right-hand side of _joint_flow, both routes"""
        jnp, d, nb = self.jnp, self.d, self.nb
        dE = consts[2]
        Kf = self._kon(dE, D, tab_on)
        if self.real:
            X = y[:2*nb].reshape(2, 2*d, 2*d)
            Gb, Ab = X[0], X[1]
            GKb = Gb*jnp.tile(Kf, (2, 2))
            out = jnp.concatenate([self._bracket(Gb, GKb).ravel(),
                                   self._bracket(GKb, Ab).ravel(),
                                   self._tip(Gb, Ab, D, -1., tab_tip, consts, False).ravel()])
            return out*(-D)
        z = self._complex(y)
        X = z[:2*nb].reshape(2, 4, d, d)
        Gb, GKb, Ab = self._block(X[0]), self._block(X[0]*Kf), self._block(X[1])
        out = jnp.concatenate([self._unblock(self._bracket(Gb, GKb)).ravel(),
                               self._unblock(self._bracket(GKb, Ab)).ravel(),
                               self._tip(Gb, Ab, D, -1., tab_tip, consts, False).ravel()])
        return self._floats(out*(-D))

    def fun(self, s, y):
        """fun(t, y) of the integrator, with t = ln D"""
        D = float(np.exp(s))
        return self.rhs(D, y, _tables(self.onshell_kernel, D), _tables(self.kernel, D),
                        self.consts)

    def couplings(self, y):
        """G and A_on in block form on the device"""
        d, nb = self.d, self.nb
        if self.real:
            X = y[:2*nb].reshape(2, 2*d, 2*d)
            return X[0], X[1]
        X = self._complex(y[:4*nb]).reshape(2, 4, d, d)
        return self._block(X[0]), self._block(X[1])


def memory_model(d, m, real=False):
    """Bytes one flow on the device holds for an impurity of d states and m
    occupied columns, as (fixed, per incoming energy), the device form of
    leadinglog.memory_model: about 20 copies of the state (13 stages, the
    state at the start and at the end of a step, a stage argument, and at a
    stop three more stages and the difference for the dense output), whose
    couplings are 8 d^2 numbers and whose tip vertex is 4 d m per energy;
    about ten 2d x 2d matrices of the right-hand side; three chunks of the
    tip part. Real numbers when the flow takes the real route. On the RTX
    A2000 of d22 it gives 6.1 GB for 8 spins 1/2 with 223 occupied states and
    122 energies, where 6.15 GB were measured; for one occupied state it is
    above the peak (3.3 GB against 2.4 GB at 10 spins), since the chunks of a
    small impurity are smaller than CHUNK_BYTES."""
    item = 8 if real else 16
    return (20*8 + 10*4)*d*d*item + 3*CHUNK_BYTES, 20*4*d*m*item


def joint_flow(ks, w, cols, Jrho_s, omega0, Gamma0, rtol, atol, pole_threshold,
               kernel, onshell_kernel):
    """leadinglog._joint_flow on the device, called from there with its
    arguments checked. Returns the tip vertex (n,4,d,m) on the host and the
    info of the flow."""
    from . import leadinglog as ll
    jnp = backend.xp()
    S = ll._spin_ops(ks)
    d = S.shape[-1]
    e = np.ascontiguousarray(ks.e, dtype=float)
    w = np.ascontiguousarray(w, dtype=float)
    cols = np.asarray(cols, dtype=int)
    n, m = len(w), len(cols)
    Dmax, Dmin = ll._energy_range(ks, omega0, w, Gamma0)
    smax = np.max(np.abs(S))/2.
    ref = pole_threshold*smax
    real = ll._block_is_real(S)
    A_bare = ll._bare_vertex(S, 1., 0.)
    Ab0 = ll._block(A_bare)
    T0 = Ab0.reshape(2, d, 2, d)[:, :, :, cols]
    if real:
        y0 = np.concatenate([ll._block(ll._bare_vertex(S, ll.SIGN*Jrho_s, 0.)).ravel(),
                             Ab0.ravel(),
                             np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()]).real.copy()
    else:
        y0 = np.concatenate([ll._bare_vertex(S, ll.SIGN*Jrho_s, 0.).ravel(), A_bare.ravel(),
                             np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()]).view(float)
    flow = _DeviceFlow(e, w, cols, omega0, Gamma0, real, kernel, onshell_kernel)
    solver = _Dop853(flow.fun, np.log(Dmax), jnp.asarray(y0), np.log(Dmin), rtol, atol)
    del y0
    y, D_stop, max_coupling = _integrate(solver, lambda y: float(flow.coupling(y)), ref)
    nfev = solver.nfev
    del solver
    tip = y[2*flow.nb:] if real else y[4*flow.nb:]
    if D_stop is not None and D_stop > Dmin:
        Gb, Ab = flow.couplings(y)
        extra = flow.tip(Gb, Ab, D_stop, Dmin, _tables(kernel, D_stop, Dmin), flow.consts)
        if not real:
            extra = flow._floats(extra.ravel())
        tip = tip + extra.ravel()
    tip = np.asarray(tip)
    Tb = (tip if real else tip.view(complex)).reshape(n, 2, d, 2, m)
    return ll._cut_columns(Tb), {"D_stop": D_stop,
                                 "max_coupling": float(max_coupling/(smax if smax > 0. else 1.)),
                                 "nfev": nfev}
