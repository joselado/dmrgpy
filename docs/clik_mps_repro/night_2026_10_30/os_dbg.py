import sys; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.tdvp import tdvp_step
from dmrgpy.pyitensor.gse import global_subspace_expand
from dmrgpy.pyitensor.mpsalgebra import inner
from dmrgpy.pyitensor.metts import product_state
from dmrgpy.pyitensor.mpscontainer import _link_at
L, N = 3, 2
herm = sys.argv[1] == "herm"
bc, lc = bh(L, N)
dc = lc.doubled; S = dc._session
G = S._mpo(lc.L.to_terms())
if herm:   # a Hermitian generator with the same charges: the doubled-chain H x 1 + 1 x H^T
    from dmrgpy.lindblad import on_sublattice, transpose
    G = S._mpo((on_sublattice(bc.hamiltonian, 0) + on_sublattice(transpose(bc.hamiltonian), 1)).to_terms())
nk = sum(dc.N[2*i] for i in range(L)); nb = sum(dc.N[2*i+1] for i in range(L))
Q = S._mpo((((nk - N)*(nk - N)) + ((nb - N)*(nb - N))).to_terms())
p, u = lc._seeds(*bh_start(L, N))
psi = product_state(S.sites, u); psi.normalize()
for m in range(3):
    psi = global_subspace_expand(G, psi, 3, 1e-10, maxdim=64, bond_maxdim=64)
    psi, _ = clik._step(psi, G, 0.5, 1e-10, 64, 50, not herm and False or herm)
bonds = lambda x: [_link_at(x, i, i+1).dim for i in range(1, x.length())]
print("start bonds", bonds(psi), "dev %.1e" % abs(inner(psi, Q, psi)))
g = global_subspace_expand(G, psi, 3, 1e-10, maxdim=64, bond_maxdim=64)
print("GSE: bonds", bonds(g), "1-|<psi|gse>| %.1e" % (1-abs(inner(psi, g))/np.sqrt(abs(inner(g,g)))), "dev %.1e" % abs(inner(g, Q, g)))
nrm = lambda x: x * (1/np.sqrt(abs(inner(x, x))))
for dt in (0.5, 0.1, 0.02):
    z = 1j*dt
    A = nrm(tdvp_step(psi.copy(), G, z, 1e-10, 64, hermitian=herm))
    B = nrm(tdvp_step(psi.copy(), G, z, 1e-10, 64, num_center=1, hermitian=herm))
    C = nrm(tdvp_step(g.copy(), G, z, 1e-10, 64, num_center=1, hermitian=herm))
    print("dt %.2f  two-one: %.1e  two-(gse+one): %.1e  dev(gse+one) %.1e" % (dt, 1-abs(inner(A,B)), 1-abs(inner(A,C)), abs(inner(C,Q,C))))
