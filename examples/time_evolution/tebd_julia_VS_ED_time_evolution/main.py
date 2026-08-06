# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test: real-time quench evolution on the live Julia backend
# (itensor_version="julia_live") with tevol_method="TEBD" -- 2nd-order
# Trotter (odd/even bond) gate evolution, mpsjulialive/tebd.jl -- must
# agree with exact diagonalization, mirroring
# ../tebd_VS_ED_time_evolution's own itensor_version="python" check and
# ../tdvp_gse_julia_VS_ED_time_evolution's julia_live structure.
#
# Unlike the other two backends' TEBD (which consume
# MultiOperator.to_terms()'s Jordan-Wigner-*predressed* "A"/"Adag"/"F"
# term list), tebd.jl resolves the Jordan-Wigner string itself, from
# scratch, off the *raw* "C"/"Cdag" term list mpo.py's text_mpo() already
# serializes for the MPO path -- real ITensors.jl's builtin "Fermion" site
# type only defines op("C",..)/op("Cdag",..)/op("F",..), not dmrgpy's own
# "A"/"Adag" names (see mpsjulialive/tebd.jl's header comment). This
# script sweeps both a spin chain and a fermion hopping chain to exercise
# that path (the fermion case is the one that would break if the from-
# scratch Jordan-Wigner threading were wrong).

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain, fermionchain
from dmrgpy import timedependent

n = 4  # small enough for exact diagonalization
nt, dt = 50, 0.05


def spin_quench():
    """Neel-favoring field -> nearest-neighbor Heisenberg, same setup as
    ../tebd_VS_ED_time_evolution."""
    sc = spinchain.Spin_Chain([2 for _ in range(n)])  # S=1/2
    sc.setup_julia()
    sc.tevol_method = "TEBD"
    h0, h1 = 0, 0
    for i in range(n):
        h0 = h0 + (-1)**i*sc.Sz[i]
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h0)
    wf, wfED = sc.get_gs(), sc.get_gs(mode="ED")
    sc.set_hamiltonian(h1)
    ts, sz = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
    tsED, szED = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt,
            wf=wfED, mode="ED")
    return "spin chain", ts, sz, tsED, szED


def fermion_quench():
    """Staggered onsite field -> hopping + staggered onsite Hamiltonian --
    exercises tebd.jl's Jordan-Wigner threading (see this file's own
    header comment). Finer dt than spin_quench()'s: confirmed directly
    that dt=0.05 accumulates just over 1e-4 of 2nd-order Trotter error by
    t=2.0 for this Hamiltonian's larger coupling scale (hopping=1 across
    3 bonds plus a 0.3 onsite term); halving dt (same total evolved time)
    brings it comfortably inside tolerance."""
    nt_f, dt_f = 100, 0.02
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_julia()
    fc.tevol_method = "TEBD"
    h0, h1 = 0, 0
    for i in range(n):
        h0 = h0 + (-1)**i*fc.N[i]
    for i in range(n - 1):
        h1 = h1 + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h1 = h1 + 0.3*sum(fc.N[i] for i in range(n))
    fc.set_hamiltonian(h0)
    wf, wfED = fc.get_gs(), fc.get_gs(mode="ED")
    fc.set_hamiltonian(h1)
    ts, n2 = timedependent.evolve_and_measure(fc, operator=fc.N[2], nt=nt_f, dt=dt_f, wf=wf)
    # edtk's wf= unwrapping only matches an exact edchain.State, which
    # Fermionic_State (its own subclass) fails -- a pre-existing,
    # unrelated ED-backend gap; wfED.v sidesteps it (see
    # tests/test_julia_live.py::test_tebd_matches_ed_fermion_chain).
    tsED, n2ED = timedependent.evolve_and_measure(fc, operator=fc.N[2], nt=nt_f, dt=dt_f,
            wf=wfED.v, mode="ED")
    return "fermion chain", ts, n2, tsED, n2ED


results = [spin_quench(), fermion_quench()]

tol = 1e-4
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, (label, ts, y, tsED, yED) in zip(axes, results):
    diff = np.max(np.abs(y - yED))
    print("%-14s max |diff| (julia TEBD vs ED) = %.3e" % (label, diff))
    assert diff < tol, "%s: TEBD (julia_live) disagrees with ED by %g (tol=%g)" \
            % (label, diff, tol)
    ax.plot(ts, y.real, label="DMRG (TEBD, julia_live)", c="blue")
    ax.scatter(tsED, yED.real, label="ED", c="red", s=14, zorder=3)
    ax.set_xlabel("time") ; ax.set_title(label)
    ax.legend(fontsize=8) ; ax.grid(alpha=0.3)
axes[0].set_ylabel(r"$\langle S_0^z\rangle(t)$")
axes[1].set_ylabel(r"$\langle N_2\rangle(t)$")
fig.suptitle("julia_live TEBD real-time quench vs exact diagonalization (n=%d)" % n)
fig.tight_layout()
print("TEST PASSED")
plt.show()
