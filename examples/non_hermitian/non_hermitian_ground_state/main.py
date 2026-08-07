# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

from dmrgpy import mpsalgebra
# now define the total spin operator
hi = sum(sc.Sx)*sum(sc.Sx) + sum(sc.Sy)*sum(sc.Sy) + sum(sc.Sz)*sum(sc.Sz)
sc.maxm = 10

def get(Starget):
    """Return S(S+1) of the state picked by get_gs() for a "dummy"
    non-Hermitian Hamiltonian H + i*(S^2-Starget*(Starget+1)), which adds
    a purely imaginary shift to every total-spin sector except S=Starget
    (whose eigenvalues stay real, Im(e)=0). The Re(e) part of every
    eigenvalue is untouched by this shift, so this only demonstrates a
    genuine selection of the S=Starget sector if the underlying NH-DMRG
    solver targets Im(e)=0 rather than simply smallest Re(e)."""
    ss1 = Starget*(Starget+1.) # S(S+1) for the targeted spin
    Hn  = h + 1j*(hi-ss1) # non-Hermitian Hamiltonian
    # (use set_hamiltonian, not a bare attribute assignment, so that
    # get_gs() below is forced to recompute for each target instead of
    # returning the cached ground state from a previous call)
    sc.set_hamiltonian(Hn)
    wf = sc.get_gs()
    return wf.dot(hi*wf).real

# compute S(S+1) of the selected eigenvector, targeting S=1
print("S(S+1)",get(1.))

# sweep over the total-spin sectors that are actually reachable for
# n=4 spin-1/2 sites (S=0,1,2) and see which sector get_gs() actually
# lands on for each target -- nhdmrg.py's docstring documents that this
# backend's non-Hermitian solver targets the eigenvalue of smallest
# Re(e) (not Im(e)=0), which this Hamiltonian's construction does not
# change across sectors, so all three targets are expected to converge
# to the same true (S=0) ground state rather than tracking the diagonal
Stargets = [0.,1.,2.]
achieved = [get(S) for S in Stargets]
for S,a in zip(Stargets,achieved):
    print("Target S =",S,"achieved S(S+1) =",a)

import matplotlib.pyplot as plt
expected = [S*(S+1) for S in Stargets]
plt.plot(expected,expected,"k--",label="Ideal (S(S+1)=target)")
plt.scatter(expected,achieved,c="red",s=80,label="get_gs() result")
plt.xlabel("Target S(S+1)")
plt.ylabel("Achieved S(S+1)")
plt.legend()
plt.show()









