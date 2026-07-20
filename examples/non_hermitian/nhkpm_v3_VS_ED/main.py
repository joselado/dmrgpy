# Non-Hermitian Kernel Polynomial Method (NH-KPM) dynamical correlator,
# cross-checked between the ED backend and the ITensor v3 (mpscpp3)
# backend.
#
# NH-KPM (src/dmrgpy/nonhermitian/kpm.py, src/dmrgpy/algebra/kpm.py's
# get_mu_n_nh/get_spec_kpm_nh, and Chain::nhkpm_moments in
# mpscpp3/chain_session.h) is a port of the biorthogonal Chebyshev-moment
# algorithm from NHKPM.jl (https://github.com/GUANGZECHEN/NHKPM.jl, itself
# implementing Phys. Rev. Lett. 130, 100401) to dmrgpy. It computes the
# dynamical correlator <psi_L|A(z) B|psi_R> at complex frequency
# z=e0+e+1j*delta for a non-Hermitian Hamiltonian, with (psi_L,psi_R) the
# biorthogonal ground state pair from the existing NH-DMRG solver
# (nhdmrg()) on the DMRG side, or from algebra.biorthogonal_ground_state
# on the ED side.
#
# Unlike the ordinary (Hermitian) KPM dynamical correlator, the expansion
# operator (z*Id-H)/E_max depends on z itself, so the moments are
# recomputed from scratch at every frequency point rather than amortized
# once over the whole spectrum - this mirrors the reference algorithm's
# own cost profile, and is why this example keeps both the chain and the
# frequency mesh small.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

n = 4

def build_chain(version):
    fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
    fc.setup_cpp(version=version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
    for i in range(n):
        h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i] # staggered imaginary potential
    for i in range(n-1):
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    # tiny symmetry-breaking field: without it the ground state is part of
    # a Re-degenerate manifold (the staggered +-i potential's conjugate
    # partners share the same real part), and independently diagonalizing
    # H and H^dagger on the ED side (or independently converging psil/psir
    # on the DMRG side) can then land on different members of that
    # manifold, so <psi_L|psi_R> collapses towards 0 instead of the
    # biorthogonal pair actually agreeing - see nh_dmrg_v3_port memory notes
    # about the same Re-degeneracy pitfall in nhdmrg() itself.
    h = h + 0.03*fc.N[0]
    fc.set_hamiltonian(h)
    return fc

es = np.linspace(0.,4.,12) # frequency mesh
name = None # set below, same pair of operators for both backends
kwargs = dict(E_max=10,n=80,delta=0.3,es=es)

fc_ed = build_chain(3) # itensor_version only matters for the DMRG call below
assert not fc_ed.is_hermitian(fc_ed.hamiltonian)
name = [fc_ed.N[0],fc_ed.N[0]] # <N_0(z) N_0(0)>
es_ed,ys_ed = fc_ed.get_dynamical_correlator(mode="ED",submode="KPM",
        name=name,**kwargs)

fc_v3 = build_chain(3)
name3 = [fc_v3.N[0],fc_v3.N[0]]
es_v3,ys_v3 = fc_v3.get_dynamical_correlator(mode="DMRG",submode="KPM",
        name=name3,**kwargs)

assert np.all(np.isfinite(ys_ed.real)) and np.all(np.isfinite(ys_ed.imag))
assert np.all(np.isfinite(ys_v3.real)) and np.all(np.isfinite(ys_v3.imag))

reldiff = np.max(np.abs(ys_v3-ys_ed))/np.max(np.abs(ys_ed))
print("NH-KPM dynamical correlator, ED vs itensor_version=3:")
for e,yed,yv3 in zip(es,ys_ed,ys_v3):
    print(f"  e={e:.3f}  ED={yed:.6f}  v3={yv3:.6f}")
print("max relative difference:",reldiff)
assert reldiff<1e-6, "ED and v3 NH-KPM dynamical correlators disagree"
print("NH-KPM (ED vs v3) regression test PASSED")
