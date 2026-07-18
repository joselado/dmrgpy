# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: entanglement entropy (single-site, pair, and mutual
# information, via the generic reduced-density-matrix machinery in
# entropy.py/densitymatrix.py's reduced_dm_projective -- works on any
# wavefunction type, DMRG or ED, since it only uses vev-style operator
# algebra, not a C++ session) must agree between DMRG (v2, v3,
# pure-Python) and ED. None of the existing entanglement_entropy/*
# examples check against ED -- they only plot a single DMRG backend's
# result. This directly exercises entropy.py, whose formula was one of
# the real bugs fixed in the "Many_Body_Chain call-graph review" pass
# (that specific fix was in entanglement.py's Peschel's-formula path,
# a different function -- but the same class of bug, a wrong entropy
# formula, could just as easily hit entropy_dm() here).
#
# Important: entropy_dm() is the SAME function regardless of whether the
# reduced density matrix came from a DMRG or an ED wavefunction, so a
# DMRG-vs-ED comparison alone cannot catch a bug in that shared formula
# (confirmed directly: injecting a uniform scaling bug into entropy_dm()
# and rerunning this test, DMRG and ED still agreed with each other --
# both went through the same broken formula). The known-value check
# below (an exact 2-site AFM singlet has entanglement entropy = log(2),
# independent of dmrgpy's own machinery) is what actually catches a
# formula-level bug; the DMRG-vs-ED loop below it only catches DMRG
# under-converging or diverging from ED given a correct formula.
import numpy as np
from dmrgpy import spinchain

# --- known-value sanity check: 2-site AFM Heisenberg singlet ---
# ground state is exactly (|ud>-|du>)/sqrt(2) -- maximally entangled,
# so the entropy of either site (equivalently, of the pair with the
# rest of a 2-site system) is exactly log(2), by hand, independent of
# any of dmrgpy's own reduced-density-matrix code.
sc2 = spinchain.Spin_Chain(["S=1/2","S=1/2"])
h2 = sc2.Sx[0]*sc2.Sx[1] + sc2.Sy[0]*sc2.Sy[1] + sc2.Sz[0]*sc2.Sz[1]
sc2.set_hamiltonian(h2)
wf2 = sc2.get_gs(mode="ED")
s_singlet = sc2.get_site_entropy(wf2,0)
print("2-site AFM singlet: site(0) entropy = %.6f (exact: log(2) = %.6f)"%(s_singlet,np.log(2)))
tol_exact = 1e-4
diff_exact = abs(s_singlet-np.log(2))
assert diff_exact<tol_exact, "singlet entropy %.6f != log(2) (diff=%g, tol=%g) -- entropy_dm() formula looks wrong"%(s_singlet,diff_exact,tol_exact)

n = 4 # small enough for exact diagonalization
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + 0.3*sc.Sz[0] # small field, breaks the fully-symmetric special
                      # case so entanglement isn't trivially uniform
sc.set_hamiltonian(h)

wf_ed = sc.get_gs(mode="ED")
s_site_ed = sc.get_site_entropy(wf_ed,0)
s_pair_ed = sc.get_pair_entropy(wf_ed,0,1)
s_mutual_ed = sc.get_mutual_information(wf_ed,0,2)
print("ED: site(0)=%.6f  pair(0,1)=%.6f  mutual(0,2)=%.6f"%(s_site_ed,s_pair_ed,s_mutual_ed))

tol = 1e-3
for v in [2,3,"python"]:
    sc.set_hamiltonian(h) # re-set: switching backend invalidates the
                           # cached ground state (see restart(), called
                           # by set_hamiltonian) so each backend gets a
                           # fresh DMRG run
    if v!="python": sc.setup_cpp(v)
    else: sc.setup_python()
    wf = sc.get_gs(mode="DMRG")
    s_site = sc.get_site_entropy(wf,0)
    s_pair = sc.get_pair_entropy(wf,0,1)
    s_mutual = sc.get_mutual_information(wf,0,2)
    print("v%s: site(0)=%.6f  pair(0,1)=%.6f  mutual(0,2)=%.6f"%(v,s_site,s_pair,s_mutual))

    d_site = abs(s_site-s_site_ed)
    d_pair = abs(s_pair-s_pair_ed)
    d_mutual = abs(s_mutual-s_mutual_ed)
    assert d_site<tol, "v%s site entropy disagrees with ED by %g (tol=%g)"%(v,d_site,tol)
    assert d_pair<tol, "v%s pair entropy disagrees with ED by %g (tol=%g)"%(v,d_pair,tol)
    assert d_mutual<tol, "v%s mutual information disagrees with ED by %g (tol=%g)"%(v,d_mutual,tol)

print("TEST PASSED")
