## DMRGPY ##

# Summary #

This is a Python library to compute quasi-one-dimensional
spin chains and fermionic systems using matrix product states
with the matrix renormalization group as implemented in ITensor. Most
of the computations can be performed both with DMRG and exact
diagonalization for small systems, which allows to benchmark the
results.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development.

# How to install #

The script install.sh will compile both ITensor and a C++ program
that uses it. Afterwards, it is only required to add to the .bashrc
the following line

export DMRGROOT=PATH_TO_DMRGPY"/src"

After this, you can write in your Python script

import os
import sys
sys.path.append(os.environ["DMRGROOT"])

And import the sublibrary that you want, for example

from dmrgpy import spinchain

# Capabilities #
- Ground state energy
- Excitation gap
- Excited states
- Static correlation functions
- Time evolution and measurements
- Dynamical correlation functions computed with the Kernel polynomial method
- Dynamical correlation functions with time dependent DMRG


# Examples

## Ground state energy of an S=1/2 spin chain
```python
from dmrgpy import spinchain
spins = [2 for i in range(30)] # 2*S+1=2 for S=1/2
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
print("Ground state energy",sc.gs_energy())
```

## Static correlator of an S=1 spin chain
```python
from dmrgpy import spinchain
spins = [3 for i in range(30)] # 2*S+1=3 for S=1
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
pairs = [(0,i) for i in range(30)] # between the edge and the rest
cs = sc.get_correlator(pairs)
```

## Ground state energy of a bilinear-biquadratic Hamiltonian
```python
from dmrgpy import spinchain
ns = 6 # number of sites in the spin chain
spins = [3 for i in range(ns)] # S=1 chain
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
h = 0 # initialize Hamiltonian
Si = [sc.Sx,sc.Sy,sc.Sz] # store the three components
for i in range(ns-1): # loop 
    for S in Si: h = h + S[i]*S[i+1]  # bilinear
    for S in Si: h = h + 1./3.*S[i]*S[i+1]*S[i]*S[i+1]  # biquadratic
sc.set_hamiltonian(h) # create the Hamiltonian
print("Energy with DMRG",sc.gs_energy(mode="DMRG"))
print("Energy with ED",sc.gs_energy(mode="ED"))
```

## Magnetization an S=1 spin chain with an edge magnetic field
```python
from dmrgpy import spinchain
spins = [3 for i in range(40)] # 2*S+1=3 for S=1
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
sc.set_fields(lambda i: [0,0,(i==0)*0.01]) # only in the first site
mx,mx,mz = sc.get_magnetization()
print("Mz",mz)
```

## Bond dimension energy convergence for an S=1/2 Heisenberg chain
```python
from dmrgpy import spinchain
spins = [2 for i in range(30)] # 2*S+1=2 for S=1/2
for maxm in [1,2,5,10,20,30,40]: # loop over bond dimension
  sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
  sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
  sc.maxm = maxm # set the bond dimension
  e = sc.gs_energy() # get the ground state energy
  print("Energy",e,"for bond dimension",maxm)
```


## Excited states with DMRG and ED 
```python
from dmrgpy import spinchain
spins = [2 for i in range(12)] # 2*S+1=2 for S=1/2
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
es1 = sc.get_excited(n=6,mode="DMRG")
es2 = sc.get_excited(n=6,mode="ED")
print("Excited states with DMRG",es1)
print("Excited states with ED",es2)
```

## Singlet-triplet gap of the Haldane Heisenberg S=1 spin chain
```python
from dmrgpy import spinchain
# Haldane chain with S=1/2 on the edge to remove the topological modes
spins = [2]+[3 for i in range(40)]+[2]
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
es = sc.get_excited(n=2,mode="DMRG")
gap = es[1]-es[0] # compute gap
print("Gap of the Haldane chain",gap)
```

## Edge dynamical correlator of a Haldane chain
```python
from dmrgpy import spinchain
spins = [3 for i in range(40)] # 2*S+1=3 for S=1
sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
sc.get_dynamical_correlator(i=0,j=0,name="ZZ")
```


## Spin and charge correlator of the 1D Hubbard model
```python
from dmrgpy import fermionchain
n = 20 # number of sites
fc = fermionchain.Spinful_Fermionic_Hamiltonian(n)
# first neighbor hopping
fc.set_hoppings_spinful(lambda i,j: (abs(i-j)==1)*1.0) 
# Hubbard term
fc.set_hubbard_spinful(lambda i,j: ((i-j)==0)*1.0) 
pairs = [(0,i) for i in range(n)]
# compute the two correlators
zz = fc.get_correlator(pairs=pairs,name="ZZ")
dd = fc.get_correlator(pairs=pairs,name="densitydensity")
print("Spin correlators",zz)
print("Density correlators",dd)
```


## Generic interacting fermionic Hamiltonian
```python
import numpy as np
from dmrgpy import fermionchain
n = 6 # number of different spinless fermionic orbitals
# fc is an object that contains the information of the many body system
fc = fermionchain.Fermionic_Hamiltonian(n) # create the object
# create a random Hermitian hopping matrix
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # make it Hermitian
fc.set_hoppings(lambda i,j: m[i,j])
def vijkl(i,j,k,l):
    """Function defining the many body interaction"""
    if i==j and k==l and abs(i-k)==1: return 1.0
    else: return 0.0
fc.set_vijkl(vijkl) # add interaction term
print("GS energy with ED",fc.gs_energy(mode="ED")) # energy with exact diag
print("GS energy with DMRG",fc.gs_energy(mode="DMRG")) # energy with DMRG
```

