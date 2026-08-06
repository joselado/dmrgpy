# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical METTS (real-time finite-temperature dynamical correlators) on
# the live Julia backend (itensor_version="julia_live") vs an exact ED
# reference -- the julia_live counterpart of ../dynamical_metts_VS_ED
# (itensor_version in ("python",3)). See src/dmrgpy/mpsjulialive/metts.jl's
# metts_dynamical_correlator for the algorithm itself: for every METTS
# sample |phi> sampled by the same Markov chain ../metts_julia_VS_ED
# already uses, set |v(0)>=B|phi>, |w(0)>=|phi>, real-time-evolve both
# under H (reusing tdvp.jl's tdvp_step, no new evolution primitive), and
# average C(t)=<w(t)|A|v(t)> over samples to estimate
# C_AB(t) = <A(t)B>_beta -- a value-level port of
# pyitensor/metts.py's "Dynamical METTS" section (Wang, McClarty, Dankova,
# Honecker & Wietek, "Spectroscopy and complex-time correlations using
# minimally entangled typical thermal states", arXiv:2405.18484, Sec. II /
# Algorithm 1), not an independent implementation.
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain
from dmrgpy.edtk.edchain import EDOperator

n = 4
J, B_field = 1.0, 0.4
T = 1.0
nt, dt = 8, 0.2


def ed_time_domain_reference(sc, name, T, ts):
    """Exact C_AB(t)=<A(t)B>_T via a full Boltzmann-weighted Lehmann sum
    over every ED eigenstate -- A,B used exactly as given (no dagger),
    same convention edtk/dynamics.py's dynamical_correlator_finite_T
    uses."""
    edobj = sc.get_ED_obj()
    emu, vs = edobj.get_diagonalized_hamiltonian()
    a0 = EDOperator(name[0], edobj).SO
    b0 = EDOperator(name[1], edobj).SO
    e0 = np.min(emu)
    ex = emu - e0
    beta = 1.0/T
    w = np.exp(-beta*ex) ; w = w/w.sum()
    U = np.array(vs) ; Uh = np.conjugate(U.T)
    A = Uh@a0@U ; B = Uh@b0@U
    out = np.zeros(len(ts), dtype=complex)
    for k, t in enumerate(ts):
        phase = np.exp(1j*(ex[:, None]-ex[None, :])*t)
        out[k] = np.sum(w[:, None]*phase*A*B.T)
    return out


sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
sc.setup_julia()
h = 0
for i in range(n-1):
    h = h + J*sc.Sx[i]*sc.Sx[i+1] + J*sc.Sy[i]*sc.Sy[i+1] + J*sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + B_field*sc.Sz[i]
sc.set_hamiltonian(h)

name = (sc.Sz[0], sc.Sz[0])
ts, means, stderrs = sc.metts_dynamical_correlator(
    name, T, nt=nt, dt=dt, nsamples=150, nwarmup=30,
    dbeta_half_step=0.08, seed=42)
ref = ed_time_domain_reference(sc, name, T, ts)

print("\n=== itensor_version=julia_live, T=%g ===" % T)
print("%6s  %22s  %22s  %10s" % ("t", "METTS C(t)", "ED C(t)", "stderr"))
for k in range(nt):
    print("%6.2f  %9.5f%+9.5fj  %9.5f%+9.5fj  %10.5f"
          % (ts[k], means[k].real, means[k].imag, ref[k].real, ref[k].imag, stderrs[k]))

diff = np.abs(means-ref)
tol = np.maximum(5*stderrs, 0.03)  # 5-sigma headroom, floored for the near-zero-stderr t=0 point
assert np.all(diff <= tol), \
    "dynamical METTS disagrees with the exact ED reference beyond " \
    "the expected statistical tolerance (itensor_version=julia_live)"

print("\nTEST PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
ax1.plot(ts, ref.real, "k-o", label="ED")
ax1.errorbar(ts, means.real, yerr=stderrs, fmt="s", capsize=3, label="METTS (julia_live)")
ax1.set_xlabel("t") ; ax1.set_ylabel(r"Re $C_{S_0^zS_0^z}(t)$") ; ax1.legend() ; ax1.grid(alpha=0.3)

ax2.plot(ts, ref.imag, "k-o", label="ED")
ax2.errorbar(ts, means.imag, yerr=stderrs, fmt="s", capsize=3, label="METTS (julia_live)")
ax2.set_xlabel("t") ; ax2.set_ylabel(r"Im $C_{S_0^zS_0^z}(t)$") ; ax2.legend() ; ax2.grid(alpha=0.3)

fig.suptitle("Dynamical METTS (Julia/ITensors.jl backend) vs exact ED, n=%d Heisenberg+field chain, T=%g" % (n, T))
fig.tight_layout()
plt.show()
