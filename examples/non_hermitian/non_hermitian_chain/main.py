# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4
fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
mh = np.zeros((n,n),dtype=np.complex) # TB matrix
for i in range(n-1):
    mh[i,i+1] = 1.0
    mh[i+1,i] = 1.0
for i in range(n): mh[i,i] = 1*1j*(-1)**i
h = 0 # initialize Hamiltonian
for i in range(n):
    for j in range(n):
        h = h + mh[i,j]*fc.Cdag[i]*fc.C[j]
for i in range(n-1): h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
from dmrgpy import mpsalgebra
# the GS mode targets the state with minimum Re(E)

fc.set_hamiltonian(h) ; print(fc.get_excited(mode="ED",n=3)) 
#exit()
#fc.set_hamiltonian(h) 
#e,wf = fc.get_excited_states(mode="ED",n=1)
#wf=wf[0]
#print(wf.dot(sum(fc.N)*wf)) ; exit()
#exit()

#e,wf = mpsalgebra.mpsarnoldi(fc,h,mode="GS",n=10,verbose=2,nwf=3,maxit=3)
fc.maxm = 20
e,wf = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,verbose=2,n=3,
maxit=3)
e,wf = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,verbose=2,n=3,
maxit=3,wfs=wf)

print("MPS Energy",e)
#import scipy.linalg as lg
#es = lg.eigvals(mh) # eigenvalues
#print("Exact TB energy",np.sum(es[es.real<0.0]))









