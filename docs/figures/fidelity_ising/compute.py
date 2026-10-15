# Fidelity susceptibility of the open transverse-field Ising chain
# H = sum sz_i sz_i+1 + lambda sum sx_i (Pauli matrices) across lambda = 1,
# with ED on the shortest chain as the reference.
import numpy as np
from dmrgpy import spinchain, fidelity

ns = [12, 16, 24, 32]
# lambda starts at 0.8: below about 0.75 the Z2 doublet of the ordered
# chain splits by less than tol at n=32, and chi mixes the two states
lams = np.round(np.linspace(0.8, 1.4, 49), 8)
dl, tol = 0.02, 1e-4


def chain(n):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    sc.maxm, sc.nsweeps = 60, 20
    h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    h1 = 2*sum(sc.Sx[i] for i in range(n))
    sc.set_hamiltonian(h0 + h1)
    return sc, h0, h1


def chi(n, **kw):  # chi, total over the chain, at every lambda
    sc, h0, h1 = chain(n)
    return [fidelity.get_fidelity(sc, h0, h1, l, dl=dl, delta=tol,
                                  fmode="derivative", **kw) for l in lams]


chi_dmrg = [chi(n) for n in ns]
chi_ed = chi(ns[0], mode="ED")

np.savez("data.npz", lam=lams, ns=ns, chi=chi_dmrg, chi_ed=chi_ed)
