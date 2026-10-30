import sys, time; sys.path.insert(0, "/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import spinchain
from dmrgpy.lindblad import Lindblad_Chain
n = 6
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.5*sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
sp = sc.Sx[0] + 1j*sc.Sy[0]
sm = sc.Sx[n-1] - 1j*sc.Sy[n-1]
jumps = [0.6*sp, 0.6*sm] + [0.4*sc.Sz[i] for i in range(n)]
lc = Lindblad_Chain(sc, jumps)
t = time.time()
lam, res, obs = lc.get_spectrum()
gap = -lam[1].real
mz = [lc.steady_state_vev(sc.Sz[i]) for i in range(n)]
print("time %.0f s maxm %d" % (time.time() - t, sc.maxm))
lamE, _, obsE = lc.get_spectrum(mode="ED")
print("gap %.4f ED %.4f" % (gap, -lamE[1].real), "lam0", lam[0])
print("mz", np.round(np.real(mz), 4), "err %.1e" % max(abs(np.array(mz) - np.array([obsE(sc.Sz[i], 0) for i in range(n)]))))
