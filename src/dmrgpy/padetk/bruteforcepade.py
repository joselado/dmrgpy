"""N-point Pade approximant through given complex points.

The Thiele continued fraction (the Vidberg-Serene construction) passes
exactly through every point (z_i, u_i), and is then evaluated anywhere in
the complex plane, typically just above the real axis. The continuation is
ill conditioned: the farther the data points are from where it is evaluated,
the more it smooths, so a spectrum continued from far above the axis loses
its sharp features. The coefficients are computed in double precision.
"""
import numpy as np


def _coefficients(zs, us):
    n = len(zs)
    g = np.array(us, dtype=complex)
    a = np.zeros(n, dtype=complex)
    a[0] = g[0]
    for p in range(1, n):
        # g_p(z_i) = (g_{p-1}(z_{p-1}) - g_{p-1}(z_i)) / ((z_i - z_{p-1}) g_{p-1}(z_i))
        gp = g.copy()
        for i in range(p, n):
            gp[i] = (g[p-1] - g[i])/((zs[i] - zs[p-1])*g[i])
        g = gp
        a[p] = g[p]
    return a


def fit(zs, us):
    """Return p(z), the Pade approximant through the points (zs, us)."""
    zs = np.asarray(zs, dtype=complex)
    us = np.asarray(us, dtype=complex)
    if len(zs) != len(us) or len(zs) == 0:
        raise ValueError("pade fit: zs and us must have the same nonzero length")
    if np.any(us == 0):
        raise ValueError("pade fit: a data value is exactly zero")
    a = _coefficients(zs, us)
    n = len(a)
    def p(z):
        z = complex(z)
        out = 1.
        for k in range(n-1, 0, -1): # evaluate the continued fraction from the bottom
            out = 1. + a[k]*(z - zs[k-1])/out
        return a[0]/out
    return p
