# Bond entanglement entropies S(l) of open critical chains (XX, Heisenberg,
# critical Ising) and of the gapped spin-1 Heisenberg chain, with the exact
# XX entropies from free fermions as the reference.
import numpy as np
from dmrgpy import spinchain

L = 100
names = ["xx", "heisenberg", "ising", "spin1"]
maxm = {"xx": 200, "heisenberg": 200, "ising": 120, "spin1": 100}
ell = np.arange(1, L)


def entropies(name):
    sc = spinchain.Spin_Chain(["S=1" if name == "spin1" else "S=1/2"]*L)
    xy = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] for i in range(L-1))
    zz = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1))
    if name == "xx":
        h = xy
    elif name == "ising":
        h = zz + 0.5*sum(sc.Sx[i] for i in range(L))  # critical field
    else:
        h = xy + zz
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm[name], 20, 1e-12
    wf = sc.get_gs()
    return [wf.get_bond_entropy(i, i+1) for i in range(L-1)]


entropy = [entropies(name) for name in names]

# exact XX chain: free fermions with hopping 1/2 at half filling, S(l) from
# the eigenvalues nu of the correlation matrix restricted to the block
hop = 0.5*(np.eye(L, k=1) + np.eye(L, k=-1))
occ = np.linalg.eigh(hop)[1][:, :L//2]
cmat = occ @ occ.T
s_exact = []
for l in ell:
    nu = np.clip(np.linalg.eigvalsh(cmat[:l, :l]), 1e-15, 1-1e-15)
    s_exact.append(-np.sum(nu*np.log(nu) + (1-nu)*np.log(1-nu)))

np.savez("data.npz", L=L, names=names, ell=ell, entropy=entropy,
         s_exact=s_exact)
