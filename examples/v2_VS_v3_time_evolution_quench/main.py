# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Compare ITensor v2 vs v3 vs the pure-Python backend (itensor_version=
# "python") on real-time evolution (quench dynamics), mirroring
# examples/time_evolution/time_evolution_quench/main.py: prepare the
# ground state of one Hamiltonian, then time-evolve it under a different
# Hamiltonian and measure a local operator (Chain::quench()/
# Chain::evolve_and_measure() in chain_session.h; pyitensor's two-site
# Krylov TDVP, see tdvp.py). v3 and the python backend both default to
# TDVP (Many_Body_Chain.tevol_method, see CLAUDE.md); v2 has no TDVP and
# always uses the MPO-Taylor-expansion path -- so some v2-vs-others
# difference beyond pure numerical noise is expected here, unlike the
# other v2_VS_v3_* examples.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import timedependent

n = 6

def evolve(itensor_version):
    sc = spinchain.Spin_Chain([2 for i in range(n)],itensor_version=itensor_version)
    h0 = 0
    for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i] # Neel-favoring field
    h1 = 0
    for i in range(n-1):
        h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs() # start from the Neel-like ground state
    sc.set_hamiltonian(h1) # quench to the Heisenberg Hamiltonian
    op = sc.Sz[0]
    nt = 30 ; dt = 1e-2 # small system/few steps, just to compare backends
    (ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=nt,dt=dt,wf=wf)
    return np.array(ts),np.array(sz)

ts3,sz3 = evolve(3)
tspy,szpy = evolve("python")

print("<Sz_0>(t) sample (ITensor v3)  :",sz3[:5].real)
print("<Sz_0>(t) sample (pure Python) :",szpy[:5].real)
print("Max abs difference v3 vs pure Python (both TDVP) =",np.max(np.abs(sz3-szpy)))

# v2's ED fallback (when the mpscpp2 extension isn't compiled -- see
# mode.py) doesn't support evolve_and_measure() at all (self._session is
# None in ED mode, since real-time evolution was never wired up on that
# path); not related to the python backend, so just skip it gracefully
# instead of crashing the whole script.
try:
    ts2,sz2 = evolve(2)
    print("<Sz_0>(t) sample (ITensor v2)  :",sz2[:5].real)
    print("Max abs difference v2 vs v3          =",np.max(np.abs(sz2-sz3)))
except AttributeError as e:
    print("ITensor v2 comparison skipped ({}) -- likely v2's ED fallback, "
          "see mode.py".format(e))
