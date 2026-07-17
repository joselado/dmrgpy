# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Benchmark: speed and accuracy of the TDVP real-time-evolution engine
# (mpscpp3/TDVP/, sc.tevol_method="TDVP", the itensor_version=3 default)
# versus the legacy MPO-Taylor method it replaced as default
# (sc.tevol_method="MPO", still the backup -- see CLAUDE.md's "real-time
# MPS evolution defaults to TDVP" and timedependent.py's
# evolve_and_measure_dmrg()), for chains of n=10, 20 and 30 sites.
#
# Accuracy is assessed two different ways depending on system size:
#  - n=10: exact diagonalization is still cheap (2**10=1024-dim Hilbert
#    space), so the DMRG trajectory is compared directly against the ED
#    one (same pattern as main.py in this folder).
#  - n=20, n=30: ED is no longer feasible. Instead, each method evolves
#    the state forward for nt steps and then *backward* (dt -> -dt) for
#    nt more steps, and the result is overlapped against the original
#    state. Since the true time-evolution operator is unitary, a forward
#    +backward round trip is the identity in exact arithmetic; how far
#    the overlap |<psi0|psi_roundtrip>| deviates from 1 measures each
#    method's own numerical error, without needing an independent
#    reference. This also happens to expose an existing difference
#    between the two methods: TDVP normalizes at every step
#    (chain_session.h's tdvp_step(), "DoNormalize"=true) while the
#    MPO-Taylor evolve_and_measure() never renormalizes the state, so its
#    norm drifts over many steps and the round-trip overlap can end up
#    *larger* than 1 in magnitude, not just smaller -- hence measuring
#    the error as abs(1-overlap) rather than 1-overlap.
import time
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

ED_MAX_N = 10 # largest system size for which ED is still used as reference
NT = 30 # number of time steps (each direction, for the round trip)
DT = 0.05 # time step
MAXM = 60 # bond dimension, kept fixed across sizes for a fair comparison


def make_chain(n):
    sc = spinchain.Spin_Chain([2 for i in range(n)]) # S=1/2 chain
    sc.maxm = MAXM
    return sc


def hamiltonians(sc,n):
    h0 = 0 # Neel-favoring staggered field, used to prepare the initial state
    for i in range(n): h0 = h0+(-1)**i*sc.Sz[i]
    h1 = 0 # Heisenberg Hamiltonian, quenched into
    for i in range(n-1):
        h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h0,h1


def benchmark_size(n):
    print("="*70)
    print("System size n =",n)
    sc = make_chain(n)
    h0,h1 = hamiltonians(sc,n)
    sc.set_hamiltonian(h0)
    wf = sc.get_gs() # DMRG ground state, used as the common starting point
    wfED = sc.get_gs(mode="ED") if n<=ED_MAX_N else None
    sc.set_hamiltonian(h1) # quench Hamiltonian
    op = sc.Sz[0]

    out = {}
    for method in ["TDVP","MPO"]:
        sc.tevol_method = method
        t0 = time.time()
        ts,sz,wf_fwd = timedependent.evolve_and_measure(sc,operator=op,
                nt=NT,dt=DT,wf=wf,return_wf=True)
        t_forward = time.time()-t0

        if wfED is not None:
            _tsED,szED = timedependent.evolve_and_measure(sc,operator=op,
                    nt=NT,dt=DT,wf=wfED,mode="ED")
            error = np.max(np.abs(sz-szED))
            label = "max|DMRG-ED|"
        else:
            _ts2,_sz2,wf_back = timedependent.evolve_and_measure(sc,operator=op,
                    nt=NT,dt=-DT,wf=wf_fwd,return_wf=True)
            overlap = wf.dot(wf_back)
            error = abs(1.0-overlap)
            label = "|1-round_trip_overlap|"

        out[method] = dict(t_forward=t_forward,error=error,label=label)
        print("  %-4s: forward time = %7.3fs   %s = %.3e"%(
            method,t_forward,label,error))

    speedup = out["MPO"]["t_forward"]/out["TDVP"]["t_forward"]
    print("  Speedup (MPO forward time / TDVP forward time) = %.2fx"%speedup)
    return out


if __name__=="__main__":
    results = {}
    for n in [10,20,30]:
        results[n] = benchmark_size(n)

    print("="*70)
    print("Summary")
    print("%4s  %10s  %10s  %8s  %14s  %14s"%(
        "n","TDVP t[s]","MPO t[s]","speedup","TDVP error","MPO error"))
    for n,res in results.items():
        tt,tm = res["TDVP"]["t_forward"],res["MPO"]["t_forward"]
        print("%4d  %10.3f  %10.3f  %7.2fx  %14.3e  %14.3e"%(
            n,tt,tm,tm/tt,res["TDVP"]["error"],res["MPO"]["error"]))

    # Regression guard: TDVP should stay accurate (either against ED, or
    # via the forward/backward round trip) at every tested size -- this is
    # a much tighter tolerance than the MPO-Taylor backup can meet (see
    # the printed comparison above), by design.
    tol = 1e-4
    for n,res in results.items():
        err = res["TDVP"]["error"]
        assert err<tol, "TDVP error %.3e too large at n=%d (tol=%.1e)"%(err,n,tol)

    print("TEST PASSED")
