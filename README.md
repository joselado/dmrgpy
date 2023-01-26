## DMRGPY ##

# Summary #

This is a Python library to compute quasi-one-dimensional
spin chains and fermionic systems using matrix product states
with the density matrix renormalization group as implemented in ITensor
(C++ and Julia versions). Most
of the computations can be performed both with DMRG and exact
diagonalization for small systems, which allows to benchmark the
results.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development.

# How to install #

## Linux and Mac ##

Execute the script 

```bash
python install.py 
```

and it will compile both ITensor and a C++ program
that uses it. 

If your default C++ compiler is not g++ (version 6 or higher), execute the installation script providing the specific compiler to use (g++-6 for example)

```bash
python install.py gpp=g++-6 
```

Alternatively, in case you just want to use the Julia version,
execute the script 

```bash
python install_julia.py
```

The installation scripts will
also add dmrgpy to the PYTHONPATH of the python interpreter you used
to execute them.

Afterwards you can import the dmrgpy sublibrary that you want, for example

```python
from dmrgpy import spinchain
```

## Windows ##

For using this program in Windows, the easiest solution is to create a virtual
machine using [Virtual Box](https://www.virtualbox.org/), installing
a version of [Ubuntu](https://releases.ubuntu.com/20.04/)
in that virtual machine, and following the previous
instructions.

# Capabilities #
- Possible models include spinless fermions, spinful fermions, spins, parafermions and bosons
- Ground state energy
- Ground state wavefunction
- Excitation energies
- Excited wavefunctions
- Arbitrary expectation values, including static correlation functions
- Time evolution of arbitrary states
- MPS algebra: sum of MPS, application of operators, exponential and inverse
- MPO algebra: sums, products, trace, trace of inverse for generic operators
- Dynamical correlation functions computed with the Kernel polynomial method
- Dynamical correlation functions with time dependent DMRG
- Generic operator distributions computed with the Kernel polynomial method
- Iterative MPS Hermitian and non-Hermitian diagonalization solvers 
- Hermitian and non-Hermitian degeneracy detection

# Examples

## Ground state energy of an S=1/2 spin chain
```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
print("Ground state energy",sc.gs_energy())
```


## Static correlator of an S=1/2 spin chain
```python
from dmrgpy import spinchain
n = 30
spins = ["S=1/2" for i in range(n)] # S=1 in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
cs = [sc.vev(sc.Sz[0]*sc.Sz[i]).real for i in range(n)]
```
![Alt text](images/S12chain_correlator.png?raw=true "Static correlator in an S=1/2 chain, showing power-law decay of correlations")


## Static correlator of an S=1 spin chain
```python
from dmrgpy import spinchain
n = 30
spins = ["S=1" for i in range(n)] # S=1 in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian
cs = [sc.vev(sc.Sz[0]*sc.Sz[i]).real for i in range(n)]
```
![Alt text](images/S1chain_correlator.png?raw=true "Static correlator in an S=1 chain, showing coupling of the emergent edge excitations")


## Conformal field theory central charge of a critical Ising model
```python
from dmrgpy import spinchain
n = 100 # number of sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # initialize
for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1] # Ising coupling
for i in range(n): h = h + 0.5*sc.Sx[i] # transverse field
sc.set_hamiltonian(h) # set the Hamiltonian
sc.maxm = 200 # increase bond dimension for a critical system
wf = sc.get_gs() # compute ground state
print("Central charge",wf.get_CFT_central_charge()) # compute central charge
```

## Ground state energy of a bilinear-biquadratic Hamiltonian
```python
from dmrgpy import spinchain
ns = 6 # number of sites in the spin chain
spins = ["S=1" for i in range(ns)] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
Si = [sc.Sx,sc.Sy,sc.Sz] # store the three components
for i in range(ns-1): # loop 
    for S in Si: h = h + S[i]*S[i+1]  # bilinear
    for S in Si: h = h + 1./3.*S[i]*S[i+1]*S[i]*S[i+1]  # biquadratic
sc.set_hamiltonian(h) # create the Hamiltonian
print("Energy with DMRG",sc.gs_energy(mode="DMRG"))
print("Energy with ED",sc.gs_energy(mode="ED"))
```

## Magnetization of an S=1 spin chain with an edge magnetic field
```python
from dmrgpy import spinchain
n = 40
spins = ["S=1" for i in range(n)] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
h = h + sc.Sz[0]*0.1 # edge magnetic field
sc.set_hamiltonian(h) # create the Hamiltonian
mz = [sc.vev(sc.Sz[i]).real for i in range(n)]
print("Mz",mz)
```

## Bond dimension energy convergence for an S=1/2 Heisenberg chain
```python
from dmrgpy import spinchain
import numpy as np
n= 30 # size of the chain
spins = ["S=1/2" for i in range(n)] # S=1/2 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
bds = range(3,20,2) # bond dimension
es,des = [],[] # storage of energies and fluctuations
for maxm in bds: # loop over bond dimension
  sc.set_hamiltonian(h) # create the Hamiltonian
  sc.maxm = maxm # set the bond dimension
  e = sc.gs_energy() # get the ground state energy
  wf = sc.get_gs() ; de = wf.dot(h*(h*wf)) # Energy square
  de = np.sqrt(np.abs(de-e**2)) # energy fluctuation
  es.append(e/n) # store energy
  des.append(de/n) # energy fluctuation
```

![Alt text](images/bond_dimension.png?raw=true "Convergence of the energy as a function of the bond dimension for an S=1/2 chain")

## Excited states with DMRG and ED 
```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(12)] # 2*S+1=2 for S=1/2
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
es1 = sc.get_excited(n=6,mode="DMRG")
es2 = sc.get_excited(n=6,mode="ED")
print("Excited states with DMRG",es1)
print("Excited states with ED",es2)
```

## Singlet-triplet gap of the Haldane Heisenberg S=1 spin chain
```python
from dmrgpy import spinchain
# Haldane chain with S=1/2 on the edge to remove the topological modes
spins = ["S=1/2"]+["S=1" for i in range(40)]+["S=1/2"]
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1): 
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
es = sc.get_excited(n=2,mode="DMRG")
gap = es[1]-es[0] # compute gap
print("Gap of the Haldane chain",gap)
```



## Local dynamical spin correlator of an S=1/2 chain
```python
import numpy as np
from dmrgpy import spinchain
n = 40
# create an S=1/2 spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
zs = [] # empty list
for i in range(n): # loop over sites
  name = (sc.Sz[i],sc.Sz[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s) # store
```

![Alt text](images/dyn_corr_spatial_long.png?raw=true "Dynamical spin correlator for different sites of an S=1/2 chain")

## Local dynamical spin correlator of an S=1 chain
```python
import numpy as np
from dmrgpy import spinchain
n = 40
# create an S=1/2 spin chain
spins = ["S=1" for i in range(n)] # spin 1/2 heisenberg chain
# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
zs = [] # empty list
for i in range(n): # loop over sites
  name = (sc.Sz[i],sc.Sz[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s) # store
```

![Alt text](images/dyn_corr_spatial_long_S1.png?raw=true "Dynamical spin correlator for different sites of an S=1 chain")



## Local dynamical spin correlator of an S=1/2 chain with a S=1 impurity
```python
import numpy as np
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(14)] # spin 1/2 heisenberg chain
spins = spins + ["S=1"] + spins # put S=1 in the middle
n = len(spins) # total number of spins
# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
zs = [] # empty list
for i in range(n): # loop over sites
  name = (sc.Sz[i],sc.Sz[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s.real) # store
```

![Alt text](images/dyn_corr_spatial_impurity.png?raw=true "Dynamical spin correlator for an S=1/2 with an S=1 impurity in the middle")



## Non-local dynamical spin correlator of an S=1/2 chain
```python
import numpy as np
from dmrgpy import spinchain
n = 10
# create an S=1/2 spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
for i in range(n): # loop over sites
  name = (sc.Sz[0],sc.Sz[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s) # store
```

![Alt text](images/dyn_corr_nonlocal_spatial.png?raw=true "Dynamical non-local spin correlator for different sites of an S=1/2 chain")


## Bulk and edge dynamical correlator of a Haldane chain
```python
from dmrgpy import spinchain
n = 20 ; spins = ["S=1" for i in range(n)] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
(e0,d0) = sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]))
(eb,db) = sc.get_dynamical_correlator(name=(sc.Sz[n//2],sc.Sz[n//2]))
```

![Alt text](images/dyn_corr_bulk_edge.png?raw=true "Dynamical spin correlator for different sites of an S=1 chain")

## Spin and charge correlator of the 1D Hubbard model
```python
from dmrgpy import fermionchain
n = 20 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
# first neighbor hopping
h = 0
for i in range(n-1):
  h = h + fc.Cdagup[i]*fc.Cup[i+1]
  h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # Make Hermitian
# Hubbard term
for i in range(n):
  h = h + 4.*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
fc.set_hamiltonian(h) # initialize the Hamiltonian
# compute the two correlators
zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)]
cc = [fc.vev(fc.Cdagup[0]*fc.Cup[i]).real for i in range(n)]
```

![Alt text](images/hubbard_correlator.png?raw=true "Spin and charge correlator in the Hubbard model")



## Spin correlator in the Hubbard model as function of the interaction
```python
from dmrgpy import fermionchain
import numpy as np
n = 14 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
# first neighbor hopping
h = 0
for i in range(n-1):
  h = h + fc.Cdagup[i]*fc.Cup[i+1]
  h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # Make Hermitian
# Hubbard term
hU = 0
for i in range(n):
  hU = hU + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)

zzs = [] # storage for correlators
Us = np.linspace(0.,4.,6) # Hubbard Us 
for U in Us:
  fc.set_hamiltonian(h+U*hU) # initialize the Hamiltonian
  zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)]
  zzs.append(zz) # store zz correlator
```

![Alt text](images/hubbard_correlator_VS_U.png?raw=true "Spin correlator in the Hubbard model for different interactions U")




## Generic interacting fermionic Hamiltonian
```python
import numpy as np
from dmrgpy import fermionchain
n = 6 # number of different spinless fermionic orbitals
# fc is an object that contains the information of the many body system
fc = fermionchain.Fermionic_Chain(n) # create the object
h = 0
# create random hoppings
for i in range(n):
  for j in range(i):
    h = h + fc.Cdag[i]*fc.C[j]*np.random.random()
# create random density interactions
for i in range(n):
  for j in range(i):
    h = h + fc.N[i]*fc.N[j]*np.random.random()
h = h + h.get_dagger() # make the Hamiltonian Hermitian
fc.set_hamiltonian(h) # set the Hamiltonian in the object
print("GS energy with ED",fc.gs_energy(mode="ED")) # energy with exact diag
print("GS energy with DMRG",fc.gs_energy(mode="DMRG")) # energy with DMRG
```


# Choosing between the C++ and Julia backend #
The library uses ITensor in the background. Currently dmrgpy allows to 
choose between
ITensor2 (C++), or ITensors (Julia). The default version executed is
the the C++ v2 version, if you want to instead use the Julia version
execute the method ".setup_julia()", for example

```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
sc.setup_julia()
```

and all the subsequent computations will be performed with Julia

