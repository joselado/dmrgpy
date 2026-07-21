# Add the root path of the dmrgpy library, and of tests/ for the shared
# dense Jordan-Wigner reference oracle (tests/_dense_jw_reference.py) --
# kept in one place rather than duplicated here, since it's the
# physics-critical part of the cross-check.
import os ; import sys
sys.path.append(os.getcwd()+'/../../src')
sys.path.append(os.getcwd()+'/../../tests')

# Regression test / demo for mixedchain.Mixed_Spin_Fermion_Chain: a chain
# mixing a genuine spin-1/2 site and spinful-fermion locations (paired
# spinless-fermion up/down sites, see mixedchain.py), used here for a
# small Kondo-lattice-like model -- a local moment (real spin site)
# exchange-coupled to two conduction-electron sites, with a direct
# hopping between the two fermion sites that has the spin site physically
# in between (exercising the Jordan-Wigner string crossing a
# non-fermionic site, see mixedchain.py's docstring).
#
# Ground-state energy and a couple of expectation values are checked
# against the same dense Jordan-Wigner reference tests/test_mixed_chain.py
# uses, since there is no dedicated ED backend for this chain type yet.
import numpy as np
from dmrgpy import mixedchain
from _dense_jw_reference import SX, SY, SZ, fermion_ops, embed

labels = ["F", "1/2", "F"] # fermion - spin - fermion
mc = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version=3)
mc.maxm = 60
mc.nsweeps = 20
assert mc.sites == [0, 0, 2, 0, 0] # 2 fermion sites, 1 spin site, 2 fermion sites
assert mc.phys_index == [(0, 1), 2, (3, 4)]

# Same parameters as tests/test_mixed_chain.py's KONDO_PARAMS. `hf` is a
# small longitudinal field: lifts the SU(2) ground-multiplet degeneracy
# so per-site Sz has a single, reproducible value (see
# tests/test_mixed_chain.py's KONDO_PARAMS docstring for why).
K1, K2, t, muA, muB, hf = 0.5, 0.35, 0.4, 0.15, -0.2, 0.13
h = K1*mc.SS(0, 1) + K2*mc.SS(1, 2) + muA*mc.Ntot[0] + muB*mc.Ntot[2]
hop = t*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2]) # crosses the spin site
h = h + hop + hop.get_dagger()
h = h + hf*(mc.Sz[0] + mc.Sz[1] + mc.Sz[2])

mc.set_hamiltonian(h)
e_dmrg = mc.gs_energy(mode="DMRG")
sz_spin_dmrg = mc.vev(mc.Sz[1]).real
ntot_fermA_dmrg = mc.vev(mc.Ntot[0]).real

print("DMRG ground-state energy", e_dmrg)
print("<Sz> at the spin location", sz_spin_dmrg)
print("<Ntot> at the first fermion location", ntot_fermA_dmrg)

# --- independent dense Jordan-Wigner reference ---
nsites = 5 # fAup, fAdn, spin, fBup, fBdn
fermionic_mask = [True, True, False, True, True]
C, Cdag, N = fermion_ops(fermionic_mask)
Sx_s, Sy_s, Sz_s = embed(SX,2,nsites), embed(SY,2,nsites), embed(SZ,2,nsites)

def ferm_spin_ops(up, dn):
    sx = 0.5*(Cdag[up]@C[dn]) + 0.5*(Cdag[dn]@C[up])
    sy = -0.5j*(Cdag[up]@C[dn]) + 0.5j*(Cdag[dn]@C[up])
    sz = 0.5*N[up] - 0.5*N[dn]
    ntot = N[up] + N[dn]
    return sx, sy, sz, ntot

sxA, syA, szA, ntotA = ferm_spin_ops(0, 1)
sxB, syB, szB, ntotB = ferm_spin_ops(3, 4)

H = (K1*(sxA@Sx_s + syA@Sy_s + szA@Sz_s)
     + K2*(Sx_s@sxB + Sy_s@syB + Sz_s@szB)
     + muA*ntotA + muB*ntotB + hf*(szA + Sz_s + szB))
hopref = t*(Cdag[0]@C[3] + Cdag[1]@C[4])
H = H + hopref + hopref.conj().T

es, vs = np.linalg.eigh(H)
e_ref = es[0]
gs = vs[:, 0]
sz_ref = np.vdot(gs, Sz_s@gs).real
ntot_ref = np.vdot(gs, ntotA@gs).real

print("Reference ground-state energy", e_ref)
print("Reference <Sz> at the spin location", sz_ref)
print("Reference <Ntot> at the first fermion location", ntot_ref)

assert abs(e_dmrg-e_ref) < 1e-5
assert abs(sz_spin_dmrg-sz_ref) < 1e-4
assert abs(ntot_fermA_dmrg-ntot_ref) < 1e-4
print("OK")
