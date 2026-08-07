# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical METTS (real-time finite-temperature dynamical correlators) vs
# an exact ED reference, following Wang, McClarty, Dankova, Honecker &
# Wietek, "Spectroscopy and complex-time correlations using minimally
# entangled typical thermal states", arXiv:2405.18484, Sec. II /
# Algorithm 1. See src/dmrgpy/pyitensor/metts.py's
# "Dynamical METTS" section for the algorithm itself: for every METTS
# sample |psi_i> sampled by the same Markov chain the static METTS
# algorithm (metts_VS_exact/main.py) already uses, set |v_i(0)>=B|psi_i>,
# |w_i(0)>=|psi_i>, real-time-evolve both under H, and average
# C^i(t)=<w_i(t)|A|v_i(t)> over samples to estimate
# C_AB(t) = <A(t)B>_beta.
#
# itensor_version="python" and 3 both implement this (mpscpp3/
# chain_session.h's Chain::metts_dynamical_correlator is a direct port of
# the same algorithm onto real ITensor v3) -- this example cross-checks
# both against an exact ED reference: a full Boltzmann-weighted Lehmann
# sum over every eigenstate (edtk/dynamics.py's
# dynamical_correlator_finite_T is the frequency-domain version of this
# same sum; the time-domain sum used directly here avoids an FFT
# round-trip for this comparison), for a small Heisenberg chain + field
# where the full spectrum is cheap to diagonalize directly.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain, cppext
from dmrgpy.edtk.edchain import EDOperator

n = 4
J, B_field = 1.0, 0.4
T = 1.0
nt, dt = 8, 0.2

versions = ["python"]
if cppext.available(3):
    versions.append(3)


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


results = {}
for version in versions:
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    if version == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(version=version)
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

    print("\n=== itensor_version=%s, T=%g ===" % (version, T))
    print("%6s  %22s  %22s  %10s" % ("t", "METTS C(t)", "ED C(t)", "stderr"))
    for k in range(nt):
        print("%6.2f  %9.5f%+9.5fj  %9.5f%+9.5fj  %10.5f"
              % (ts[k], means[k].real, means[k].imag, ref[k].real, ref[k].imag, stderrs[k]))

    diff = np.abs(means-ref)
    tol = np.maximum(5*stderrs, 0.03)  # 5-sigma headroom, floored for the near-zero-stderr t=0 point
    assert np.all(diff <= tol), \
        "dynamical METTS disagrees with the exact ED reference beyond " \
        "the expected statistical tolerance (itensor_version=%s)" % (version,)

    results[version] = (ts, means, stderrs, ref)

print("\nTEST PASSED")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
plotted_ref = False
for version, (ts, means, stderrs, ref) in results.items():
    if not plotted_ref:
        ax1.plot(ts, ref.real, "k-o", label="ED")
        ax2.plot(ts, ref.imag, "k-o", label="ED")
        plotted_ref = True
    ax1.errorbar(ts, means.real, yerr=stderrs, fmt="s", capsize=3,
                 label="METTS (%s)" % (version,))
    ax2.errorbar(ts, means.imag, yerr=stderrs, fmt="s", capsize=3,
                 label="METTS (%s)" % (version,))
ax1.set_xlabel("t") ; ax1.set_ylabel(r"Re $C_{S_0^zS_0^z}(t)$") ; ax1.legend() ; ax1.grid(alpha=0.3)
ax2.set_xlabel("t") ; ax2.set_ylabel(r"Im $C_{S_0^zS_0^z}(t)$") ; ax2.legend() ; ax2.grid(alpha=0.3)
fig.suptitle("Dynamical METTS vs exact ED, n=%d Heisenberg+field chain, T=%g" % (n, T))
fig.tight_layout()
plt.show()
