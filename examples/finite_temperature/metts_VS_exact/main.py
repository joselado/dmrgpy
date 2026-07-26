# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# METTS (Minimally Entangled Typical Thermal States) vs exact ED thermal
# averages, following E.M. Stoudenmire and Steven R. White, "Minimally
# entangled typical thermal state algorithms", New J. Phys. 12, 055026
# (2010), arXiv:1002.1305. See src/dmrgpy/pyitensor/metts.py for the
# algorithm itself: imaginary-time TDVP evolution of a classical product
# state by e^{-beta H/2}, followed by a sequential ("perfect sampling")
# collapse back down to a new product state, alternating the collapse
# basis between Sz and Sx from one sample to the next for ergodicity.
#
# itensor_version="python" and 3 both implement METTS (see
# vevtk/mettsvev.py's docstring; mpscpp3/chain_session.h's Chain::metts_vev
# is a direct port of the same algorithm onto real ITensor v3) -- this
# example cross-checks both against the exact ED Boltzmann average for a
# small Heisenberg chain + field, where the full spectrum is cheap to
# diagonalize directly.
import numpy as np
from dmrgpy import spinchain, cppext

n = 4
J, B = 1.0, 0.3

versions = ["python"]
if cppext.available(3):
    versions.append(3)

Ts = [0.5, 1.0, 3.0]
for version in versions:
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    if version == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(version=version)
    h = 0
    for i in range(n - 1):
        h = h + J*sc.Sx[i]*sc.Sx[i+1] + J*sc.Sy[i]*sc.Sy[i+1] + J*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n):
        h = h + B*sc.Sz[i]
    sc.set_hamiltonian(h)

    print("\n=== itensor_version=%s ===" % (version,))
    print("%6s  %12s  %20s  %12s  %20s" % ("T", "ED <Sz0>", "METTS <Sz0>", "ED <E>", "METTS <E>"))
    for T in Ts:
        ed_sz = sc.vev(sc.Sz[0], mode="ED", T=T).real
        ed_e = sc.vev(h, mode="ED", T=T).real
        metts_sz, err_sz = sc.metts_vev(sc.Sz[0], T, nsamples=300, nwarmup=30,
                                         dbeta_half_step=0.08, seed=42)
        metts_e, err_e = sc.metts_vev(h, T, nsamples=300, nwarmup=30,
                                       dbeta_half_step=0.08, seed=42)
        print("%6.2f  %12.5f  %8.5f +/- %8.5f  %12.5f  %8.5f +/- %8.5f"
              % (T, ed_sz, metts_sz.real, err_sz, ed_e, metts_e.real, err_e))

        # generous tolerance: METTS is a Monte Carlo method, so exact
        # agreement isn't expected, only agreement within a comfortable
        # multiple of its own (Markov-correlated, likely optimistic)
        # reported standard error -- see tests/test_metts_vev.py for
        # tighter, seed-pinned regression bounds.
        assert abs(metts_sz.real - ed_sz) < 0.08, \
            "METTS <Sz0> disagrees with ED at T=%g (itensor_version=%s)" % (T, version)
        assert abs(metts_e.real - ed_e) < 0.08, \
            "METTS <E> disagrees with ED at T=%g (itensor_version=%s)" % (T, version)

print("\nTEST PASSED")
