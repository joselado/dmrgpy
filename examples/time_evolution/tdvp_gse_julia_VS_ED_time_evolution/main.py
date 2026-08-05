# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test: real-time quench evolution on the live Julia backend
# (itensor_version="julia_live") with tevol_method="TDVP_GSE" -- one-site
# TDVP plus Krylov global subspace expansion (arXiv:2005.06104) -- must
# agree with exact diagonalization, same as ../tdvp_gse_VS_ED_time_evolution
# already checks for itensor_version=3 (this script mirrors it).
#
# On this backend neither piece is a dmrgpy port: ITensorMPS.jl ships the
# expansion as expand(psi, H; alg="global_krylov") (citing the same
# Yang-White paper) and its tdvp() takes nsite=1 directly, so
# mpsjulialive/tdvp.jl only wires the two together in the same
# "expand-then-step for the first tdvp_gse_sweeps steps" structure the
# other backends use.
#
# The starting state here is a genuine PRODUCT state (bond dimension 1),
# deliberately: one-site TDVP conserves bond dimension exactly, so without
# the expansion it is *structurally* unable to represent the entangled
# state the quench produces, and the run is visibly wrong rather than
# merely inaccurate. The third curve below (tdvp_gse_sweeps=0, i.e. the
# same one-site integrator with the expansion switched off) shows exactly
# that, which is what makes this a real check that the expansion runs at
# all rather than silently falling through to the two-site path.
#
# Note this is also where itensor_version=3 has a known, isolated
# bond-dimension-1 failure (see ../tdvp_gse_VS_ED_time_evolution's own
# long comment, which is why *that* script has to mix in an XX+YY coupling
# to give its starting state bond dimension >1). The Julia route handles
# the product-state start fine, so no such workaround is needed here.

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain
from dmrgpy import timedependent

n = 6            # small enough for exact diagonalization
nt = 40 ; dt = 0.05


def build(tevol_method, gse_sweeps):
    sc = spinchain.Spin_Chain([2 for i in range(n)]) # S=1/2
    sc.setup_julia()
    sc.tevol_method = tevol_method
    sc.tdvp_gse_sweeps = gse_sweeps
    sc.tdvp_gse_krylov_order = 3
    sc.tdvp_gse_cutoff = 1e-10
    h0 = 0 ; h1 = 0
    for i in range(n): # pure staggered field -> exact product ground state
        h0 = h0 + (-1)**i*sc.Sz[i]
    for i in range(n-1): # quench to Heisenberg
        h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()             # DMRG ground state of h0 (bond dim 1)
    wfED = sc.get_gs(mode="ED")  # ED ground state, same Hamiltonian
    sc.set_hamiltonian(h1)       # quench
    return sc,wf,wfED


runs = [("TDVP",0,          "two-site TDVP"),
        ("TDVP_GSE",nt,     "one-site TDVP + GSE"),
        ("TDVP_GSE",0,      "one-site TDVP, GSE off")]

results = []
for method,gse_sweeps,label in runs:
    sc,wf,wfED = build(method,gse_sweeps)
    op = sc.Sz[0]
    ts,sz = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
    tsED,szED = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,
            wf=wfED,mode="ED")
    diff = np.max(np.abs(sz-szED))
    print("%-24s (tdvp_gse_sweeps=%2d)  max |diff| vs ED = %.3e"
          %(label,gse_sweeps,diff))
    results.append((label,ts,sz,tsED,szED,diff))

tol = 1e-4
# the two real integrators must track ED...
assert results[0][5] < tol, "julia two-site TDVP disagrees with ED"
assert results[1][5] < tol, "julia TDVP_GSE disagrees with ED"
# ...and one-site TDVP with the expansion switched off must NOT, or the
# expansion isn't doing anything and the test above proves nothing
assert results[2][5] > 1e-2, ("one-site TDVP without GSE tracked ED anyway "
        "-- the expansion is not what made the TDVP_GSE run work")

print("TEST PASSED")

# visualize the comparison behind the asserts above
plt.figure(figsize=(7,4.5))
plt.scatter(results[0][3],results[0][4].real,label="ED",c="black",zorder=3,s=14)
for label,ts,sz,_tsED,_szED,diff in results:
    plt.plot(ts,sz.real,label="%s (max diff %.1e)"%(label,diff))
plt.xlabel("time") ; plt.ylabel(r"$\langle S_0^z\rangle(t)$")
plt.title("julia_live real-time evolution from a product state (n=%d)"%n)
plt.legend(fontsize=8) ; plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
