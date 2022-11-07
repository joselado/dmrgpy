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
# we are goign to define a "dummy" non-Hermitian Hamiltonian, just the
# original Hamiltonian plus S^2
ss1 = 1.*(1.+1) # S(S+1) for S=1
Hn  = h + 1j*(hi-ss1) # non-Hermitian Hamiltonian
# the previous Hamiltonian has complex eigenvalues for S!=1 sectors
# now perform the non-Hermitian diagonalization
# We will target the most negative eigenvalues with Im(e)=0
# delta is the tolerancy for the imaginary part
sc.maxm = 10
sc.hamiltonian = Hn
wf = sc.get_gs()

# now compute S(S+1) to check we get the "right" eigenvector
print("S(S+1)",wf.dot(hi*wf).real)









