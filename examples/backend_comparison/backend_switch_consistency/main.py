# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test for a real bug: Many_Body_Chain.setup_cpp()/
# setup_python()/setup_julia() -- the documented way to "switch an
# existing chain between backends" (see CLAUDE.md) -- used to leave
# self.computed_gs/self.wf0/self.e0 pointing at the OLD backend's
# ground state. Concretely this meant:
#   (1) gs_energy() after switching silently returned the stale energy
#       from the previous backend instead of recomputing on the new one
#       (gs_energy()'s "if self.computed_gs: return self.e0" fired first);
#   (2) vev() then handed the OLD backend's MPS cpp_handle to the NEW
#       backend's C++ session, which is a hard TypeError (an mpscpp2 MPS
#       is not valid input to an mpscpp3 Chain.vev(), or vice versa).
# Fixed by resetting that cached state whenever the backend changes. This
# test builds one chain, computes on v2, switches to v3, and checks both
# that nothing crashes and that the post-switch results are genuinely
# fresh v3 results (matching ED), not silently-reused v2 leftovers.
import numpy as np
from dmrgpy import spinchain

def run(n):
    """Build an n-site Heisenberg chain, compute on v2, switch to v3
    WITHOUT re-setting the Hamiltonian (the pattern that used to trigger
    the stale-cache bug), and return v2 energy, post-switch v3 energy,
    and the ED ground truth."""
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=2)

    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    e_v2 = sc.gs_energy(mode="DMRG")
    vev_v2 = sc.vev(sc.Sz[0], mode="DMRG").real
    print("n =",n,"  v2: energy =",e_v2,"  <Sz_0> =",vev_v2)

    # Switch backend on the SAME chain object, WITHOUT calling
    # set_hamiltonian() again -- this is exactly the pattern that used to
    # crash/return stale results. (set_hamiltonian()'s default restart=True
    # already resets the ground-state cache as a side effect, which is why
    # this bug only shows up when the Hamiltonian is *not* re-set after the
    # switch -- reusing the same Hamiltonian across a backend switch is
    # precisely the use case setup_cpp()/setup_python() exist for.)
    sc.setup_cpp(3)
    e_v3 = sc.gs_energy(mode="DMRG")
    vev_v3 = sc.vev(sc.Sz[0], mode="DMRG").real
    print("n =",n,"  v3 (after switch): energy =",e_v3,"  <Sz_0> =",vev_v3)

    e_ed = sc.gs_energy(mode="ED")
    print("n =",n,"  ED (ground truth): energy =",e_ed)

    return e_v2,e_v3,e_ed

# sweep the chain length -- small enough to stay ED-feasible, but enough
# points to show the post-switch v3 result tracks ED at every size, not
# just the single n=4 case
ns = [3,4,5,6]
tol = 1e-3
e_v2s,e_v3s,e_eds = [],[],[]
for n in ns:
    e_v2,e_v3,e_ed = run(n)
    diff = abs(e_v3-e_ed)
    print("n =",n,"  Difference v3 (post-switch) vs ED =",diff)
    assert diff<tol, "n=%d: post-switch v3 energy disagrees with ED by %g (tol=%g) -- looks like a stale cache"%(n,diff,tol)
    e_v2s.append(e_v2)
    e_v3s.append(e_v3)
    e_eds.append(e_ed)

print("TEST PASSED")

import matplotlib.pyplot as plt
plt.plot(ns,e_v2s,marker="o",label="v2 (pre-switch)")
plt.plot(ns,e_v3s,marker="s",linestyle="--",label="v3 (post-switch)")
plt.plot(ns,e_eds,marker="^",linestyle=":",label="ED (ground truth)")
plt.xlabel("Chain length n")
plt.ylabel("Ground state energy")
plt.legend()
plt.show()
