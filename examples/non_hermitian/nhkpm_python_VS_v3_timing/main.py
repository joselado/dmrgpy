# Non-Hermitian KPM (NH-KPM) dynamical correlator: pure-Python (pyitensor)
# backend vs ITensor v3 (mpscpp3), on a non-uniform (site-dependent
# hopping) interacting fermionic chain with a non-uniform imaginary
# onsite (gain/loss) potential.
#
# This is the same algorithm as examples/non_hermitian/nhkpm_v3_VS_ED
# (src/dmrgpy/nonhermitian/kpm.py, src/dmrgpy/pyitensor/chain.py's
# Chain.nhkpm_moments, mpscpp3/chain_session.h's Chain::nhkpm_moments),
# now exercised on the pure-Python backend as well, and timed against v3
# to gauge the cost of running NH-KPM with zero compiler/pybind11
# dependency (see docs/documentation.md's "Backend performance" section
# for the general itensor_version=3-vs-"python" comparison; NH-KPM's own
# repeated-per-frequency moment recursion makes this comparison
# particularly relevant, since it is far more matvec-heavy per correlator
# than the Hermitian KPM path).

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
from dmrgpy import fermionchain

n = 6

def build_chain(version):
    fc = fermionchain.Fermionic_Chain(n)
    if version=="python": fc.setup_python()
    else: fc.setup_cpp(version=version)
    h = 0
    # non-uniform (site-dependent) real hopping amplitudes
    ts = [1.0+0.15*i for i in range(n-1)]
    for i in range(n-1):
        h = h + ts[i]*(fc.Cdag[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.C[i])
    # non-uniform imaginary onsite (gain/loss) potential: alternating
    # sign with a site-dependent magnitude, not a uniform staggered
    # pattern (which would instead risk a Re-degenerate ground state,
    # see nhdmrg_VS_ED_VS_arnoldi and the nh_dmrg_v3_port memory notes)
    gammas = [0.5*(-1)**i*(1.0+0.1*i) for i in range(n)]
    for i in range(n):
        h = h + 1j*gammas[i]*fc.Cdag[i]*fc.C[i]
    # weak interaction: keeps the ground state genuinely many-body
    # (entangled) without growing the bond dimension needed too much
    V = 0.3
    for i in range(n-1):
        h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    # small real symmetry-breaking field: without it this particular
    # model has an EXACT real-part degeneracy at the ground state
    # (checked via ED: two lowest eigenvalues share Re(E) to 1e-15,
    # differing only in Im(E)) -- nhdmrg()'s random-start Arnoldi search
    # then has no way to consistently prefer one member of that manifold
    # over the other, so independent v3/python runs can silently converge
    # to different-but-equally-valid members and disagree with each
    # other (see the nh_dmrg_v3_port memory notes' Re-degeneracy
    # pitfall). This field opens a real-part gap of ~0.06, well above
    # nhdmrg()'s convergence tolerance.
    h = h + 0.05*fc.N[0]
    fc.set_hamiltonian(h)
    return fc

es = np.linspace(0.,4.,10) # frequency mesh
kwargs = dict(E_max=12,n=60,delta=0.3,es=es)

fc3 = build_chain(3)
assert not fc3.is_hermitian(fc3.hamiltonian)
name3 = [fc3.N[0],fc3.N[0]] # <N_0(z) N_0(0)>
t0 = time.time()
es3,ys3 = fc3.get_dynamical_correlator(mode="DMRG",submode="KPM",
        name=name3,**kwargs)
t_v3 = time.time()-t0

fcpy = build_chain("python")
namepy = [fcpy.N[0],fcpy.N[0]]
t0 = time.time()
espy,yspy = fcpy.get_dynamical_correlator(mode="DMRG",submode="KPM",
        name=namepy,**kwargs)
t_py = time.time()-t0

print(f"n={n} sites, non-uniform hopping + non-uniform imaginary onsite potential")
print(f"itensor_version=3 : {t_v3:.2f}s")
print(f"itensor_version=\"python\": {t_py:.2f}s  ({t_py/t_v3:.1f}x slower than v3)")
print()
print(f"{'e':>6} {'v3':>28} {'python':>28}")
for e,y3,ypy in zip(es,ys3,yspy):
    print(f"{e:6.3f} {y3:28.6f} {ypy:28.6f}")

assert np.all(np.isfinite(ys3.real)) and np.all(np.isfinite(ys3.imag))
assert np.all(np.isfinite(yspy.real)) and np.all(np.isfinite(yspy.imag))

reldiff = np.max(np.abs(ys3-yspy))/np.max(np.abs(ys3))
print()
print("max relative difference (v3 vs python):",reldiff)
assert reldiff<1e-6, "v3 and python NH-KPM dynamical correlators disagree"
print("NH-KPM (v3 vs python) regression test PASSED")
