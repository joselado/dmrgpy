# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# METTS (Minimally Entangled Typical Thermal States) on the live Julia
# backend (itensor_version="julia_live") vs exact ED thermal averages --
# the julia_live counterpart of ../metts_VS_exact (itensor_version in
# ("python",3)). See src/dmrgpy/mpsjulialive/metts.jl for the algorithm
# itself: imaginary-time TDVP evolution of a classical product state by
# e^{-beta H/2} (reusing tdvp.jl's tdvp_step, no new evolution primitive),
# followed by a sequential ("perfect sampling") collapse back down to a
# new product state, alternating the collapse basis between Sz and Sx from
# one sample to the next for ergodicity -- a value-level port of
# pyitensor/metts.py's algorithm (White & Stoudenmire, New J. Phys. 12,
# 055026 (2010), arXiv:1002.1305), not an independent implementation.
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain

n = 4
J, B = 1.0, 0.3

sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
sc.setup_julia()
h = 0
for i in range(n - 1):
    h = h + J*sc.Sx[i]*sc.Sx[i+1] + J*sc.Sy[i]*sc.Sy[i+1] + J*sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + B*sc.Sz[i]
sc.set_hamiltonian(h)

Ts = [0.5, 1.0, 2.0, 3.0]
ed_szs, metts_szs, metts_sz_errs = [], [], []
ed_es, metts_es, metts_e_errs = [], [], []

print("%6s  %12s  %20s  %12s  %20s" % ("T", "ED <Sz0>", "METTS <Sz0>", "ED <E>", "METTS <E>"))
for T in Ts:
    ed_sz = sc.vev(sc.Sz[0], mode="ED", T=T).real
    ed_e = sc.vev(h, mode="ED", T=T).real
    # Sz and h are measured together on one shared METTS sample chain
    # (see vevtk/mettsvev.py's docstring) rather than resampling once per
    # operator -- the nwarmup+nsamples imaginary-time evolutions are the
    # expensive part, so this halves the wall-clock time for the same
    # result.
    (metts_sz, err_sz), (metts_e, err_e) = sc.metts_vev(
        [sc.Sz[0], h], T, nsamples=300, nwarmup=30,
        dbeta_half_step=0.08, seed=42)
    print("%6.2f  %12.5f  %8.5f +/- %8.5f  %12.5f  %8.5f +/- %8.5f"
          % (T, ed_sz, metts_sz.real, err_sz, ed_e, metts_e.real, err_e))

    ed_szs.append(ed_sz) ; metts_szs.append(metts_sz.real) ; metts_sz_errs.append(err_sz)
    ed_es.append(ed_e) ; metts_es.append(metts_e.real) ; metts_e_errs.append(err_e)

    # generous tolerance: METTS is a Monte Carlo method, so exact
    # agreement isn't expected, only agreement within a comfortable
    # multiple of its own (Markov-correlated, likely optimistic) reported
    # standard error -- see tests/test_metts_vev.py for tighter,
    # seed-pinned regression bounds.
    assert abs(metts_sz.real - ed_sz) < 0.08, \
        "METTS <Sz0> disagrees with ED at T=%g (itensor_version=julia_live)" % T
    assert abs(metts_e.real - ed_e) < 0.08, \
        "METTS <E> disagrees with ED at T=%g (itensor_version=julia_live)" % T

print("\nTEST PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
ax1.plot(Ts, ed_szs, "k-o", label="ED")
ax1.errorbar(Ts, metts_szs, yerr=metts_sz_errs, fmt="s", capsize=3, label="METTS (julia_live)")
ax1.set_xlabel("T") ; ax1.set_ylabel(r"$\langle S_0^z\rangle$") ; ax1.legend() ; ax1.grid(alpha=0.3)

ax2.plot(Ts, ed_es, "k-o", label="ED")
ax2.errorbar(Ts, metts_es, yerr=metts_e_errs, fmt="s", capsize=3, label="METTS (julia_live)")
ax2.set_xlabel("T") ; ax2.set_ylabel(r"$\langle H\rangle$") ; ax2.legend() ; ax2.grid(alpha=0.3)

fig.suptitle("METTS (Julia/ITensors.jl backend) vs exact ED, n=%d Heisenberg+field chain" % n)
fig.tight_layout()
plt.show()
