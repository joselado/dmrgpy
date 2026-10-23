"""Maximum-entropy density with prescribed power moments on an interval.

Of all densities p(x) on [a,b] with int x^k p dx = mu_k for k=0..n-1, the
one of largest entropy -int p log p is p(x) = exp(sum_k lambda_k x^k). The
multipliers minimize the convex dual F(lambda) = int exp(sum lambda_k x^k)
dx - sum lambda_k mu_k, whose gradient is the moment mismatch. We map
[a,b] to [-1,1] and use Legendre polynomials instead of powers, which is the
same family of densities but keeps the Newton steps well conditioned, and
evaluate every integral on a Gauss-Legendre grid.

Same interface as PyMaxEnt (Saad and Ruai), which dmrgpy called before:
reconstruct(mu, bnds) returns (sol, lambdas), sol a vectorized callable
that is zero outside bnds.
"""
import warnings
from math import comb

import numpy as np
from numpy.polynomial import legendre


def _shifted_moments(mu, a, b):
    """Power moments of t = (2x-a-b)/(b-a) from those of x."""
    al, be = 2./(b-a), -(a+b)/(b-a)
    n = len(mu)
    return np.array([sum(comb(k, j)*al**j*be**(k-j)*mu[j]
                         for j in range(k+1)) for k in range(n)])


def reconstruct(mu, bnds=(-1., 1.), nquad=400, tol=1e-10, maxiter=200):
    """Maximum-entropy density with moments mu[0..n-1] on bnds=[a,b].

    mu[0] is the total weight. Returns (sol, lambdas): sol(x) is the
    density, lambdas the multipliers in the Legendre basis of the mapped
    interval. Warns when Newton does not reach tol, which happens when the
    moments belong to no smooth density (a few isolated points, or noise)."""
    mu = np.asarray(mu, dtype=float)
    a, b = float(bnds[0]), float(bnds[1])
    if not b > a: raise ValueError("reconstruct: bnds must have b > a")
    n = len(mu)
    if n < 1: raise ValueError("reconstruct: needs at least one moment")
    if not mu[0] > 0: raise ValueError("reconstruct: mu[0] must be positive")
    mt = _shifted_moments(mu, a, b)
    # Legendre moments L_k = int P_k(t) q(t) dt, from the power moments
    P = [legendre.leg2poly(np.eye(n)[k]) for k in range(n)]
    L = np.array([np.dot(P[k], mt[:len(P[k])]) for k in range(n)])
    t, w = legendre.leggauss(nquad)
    V = np.array([legendre.legval(t, np.eye(n)[k]) for k in range(n)])
    lam = np.zeros(n)
    lam[0] = np.log(mu[0]/2.) # uniform start
    def dual(l):
        q = np.exp(np.clip(l @ V, -700., 700.))
        return np.sum(w*q) - l @ L, q
    F, q = dual(lam)
    for it in range(maxiter):
        g = V @ (w*q) - L
        if np.max(np.abs(g)) < tol*mu[0]: break
        H = (V*(w*q)) @ V.T
        step = np.linalg.solve(H + 1e-14*np.trace(H)*np.eye(n), g)
        s = 1.
        while True: # backtracking on the convex dual
            Fn, qn = dual(lam - s*step)
            if Fn <= F - 1e-4*s*(g @ step) or s < 1e-10: break
            s *= 0.5
        lam, F, q = lam - s*step, Fn, qn
    else:
        warnings.warn("maxent reconstruct: Newton did not converge in %d "
                      "iterations (moment mismatch %.2e of the weight); the "
                      "moments may belong to no smooth density on %s"
                      % (maxiter, np.max(np.abs(g))/mu[0], [a, b]),
                      RuntimeWarning, stacklevel=2)
    al = 2./(b-a)
    def sol(x):
        x = np.asarray(x, dtype=float)
        tx = al*x - (a+b)/(b-a)
        out = al*np.exp(np.clip(legendre.legval(tx, lam), -700., 700.))
        return np.where((x >= a) & (x <= b), out, 0.)
    return sol, lam
