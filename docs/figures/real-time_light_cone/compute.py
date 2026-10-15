# A domain wall of the spin-1/2 XX chain released at t=0, evolved with TDVP,
# against the exact result: the XX chain is free fermions with hopping 1/2,
# so <Sz_i>(t) = C_ii(t) - 1/2 with C(t) = U C(0) U^dagger, U = exp(-i h t).
import numpy as np
import scipy.linalg
from dmrgpy import spinchain, timedependent

n = 40
dt, nt, nchunks = 0.05, 5, 64  # the profile is stored every nt*dt

sc = spinchain.Spin_Chain(["S=1/2"]*n)
sc.maxm, sc.tevol_method = 100, "TDVP"
# initial state: all up on the left half, all down on the right half
eps = np.where(np.arange(n) < n//2, 1.0, -1.0)
sc.set_hamiltonian(sum(-eps[i]*sc.Sz[i] for i in range(n)))
wf = sc.get_gs()
sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                       for i in range(n-1)))


def profile(wf):
    return [sc.aMb(wf, sc.Sz[i], wf).real for i in range(n)]


sz = [profile(wf)]
for _ in range(nchunks):
    _, _, wf = timedependent.evolve_and_measure(sc, operator=sc.Sz[n//2], nt=nt,
                                                dt=dt, wf=wf, return_wf=True)
    sz.append(profile(wf))
ts = nt*dt*np.arange(nchunks + 1)

h = 0.5*(np.eye(n, k=1) + np.eye(n, k=-1))
c0 = np.diag((eps > 0).astype(float))
U = [scipy.linalg.expm(-1j*h*t) for t in ts]
sz_exact = [np.real(np.diag(u @ c0 @ u.conj().T)) - 0.5 for u in U]

np.savez("data.npz", n=n, ts=ts, sz=sz, sz_exact=sz_exact)
