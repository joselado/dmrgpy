## DMRGPY ##

# Summary #

DMRGPY is a Python library to study quasi-one-dimensional quantum
many-body systems with matrix product states: spin chains, spinless and
spinful fermions, bosons and parafermions. Ground states and excitations
are computed with the density matrix renormalization group (DMRG), based
on ITensor. Small systems can also be solved with exact diagonalization
(ED), so you can check the DMRG results. Chains can be finite, or
infinite with iDMRG and VUMPS.

Many examples are in the `examples` folder.

# Disclaimer #

This library is still under heavy development.

# How to install #

## With pip ##

```bash
pip install dmrgpy
```

This installs the Python part of the library. It does not need a
compiler. It contains exact diagonalization and a DMRG/TDVP code written
in pure Python, which is used by default when the C++ part is not
compiled. You can also ask for it explicitly:

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
```

The C++ version of ITensor is not in the pip package, because it is large
and must be compiled on your machine. To get it, clone this repository
and follow the instructions below. Optional extras:
`pip install "dmrgpy[julia]"` for the Julia backend, `[jax]` to run the
pure-Python backend on a GPU, `[stats]`, `[parallel]`, or `[full]` for
all of them.

## Linux and Mac (compiled C++ backend) ##

Run the script

```bash
python install.py 
```

It first checks that everything needed is present (a C++ compiler,
LAPACK/BLAS, `pybind11`, `make`). It finds the right compiler for your
Python, also the conda compiler if you use Anaconda or Miniconda. Then it
compiles ITensor and the C++ extension.

By default it compiles ITensor v3. ITensor v2 is still supported, and you
can choose:

```bash
python install.py --itensor-version=2     # only ITensor v2
python install.py --itensor-version=3     # only ITensor v3 (the default)
python install.py --itensor-version=both  # both, one after the other
```

To use a specific C++ compiler (g++ 6 or newer, or a matching Clang):

```bash
python install.py --gpp=g++-6 
```

If LAPACK/BLAS is not found, give the path to OpenBLAS:

```bash
python install.py --openblas --openblas_libdir=/path/to/openblas/lib --openblas_includedir=/path/to/openblas/include
```

To only check the compiler and LAPACK/BLAS, without compiling ITensor
(which takes several minutes):

```bash
python install.py --doctor
```

To use the Julia version instead, run

```bash
python install_julia.py
```

Both scripts add dmrgpy to the PYTHONPATH of the Python you used. After
that you can import dmrgpy, for example

```python
from dmrgpy import spinchain
```

### HPC clusters with environment modules (e.g. Aalto's Triton) ###

On clusters with environment modules (Lmod) there may be no C++ compiler
until you load a module. On Triton, for example:

```bash
module load scicomp-python-env gcc
python3 install.py
```

The compiler that comes with a module Python (for example Triton's
`scicomp-python-env`) can be broken, with its `cc1plus` missing. In that
case `install.py` uses the compiler of the system or of the loaded
module, and prints a note. You must load the same modules (or activate
the same conda environment) in every later session and job script that
imports `dmrgpy`, because the compiled extension depends on them.

If you do not need the C++ backend on the cluster, `pip install dmrgpy`
works without loading any module.

## Windows ##

`pip install dmrgpy` works on Windows, because ED and the pure-Python
DMRG/TDVP (`itensor_version="python"`) do not need a compiler.

The C++ backends do not work on Windows. If you need them, install
[Ubuntu](https://releases.ubuntu.com/20.04/) in a virtual machine (for
example with [Virtual Box](https://www.virtualbox.org/)) and follow the
Linux instructions.

# Tutorials #
You can find several tutorials [here](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024), on the following topics

- [Many-body quantum magnets](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/quantum_magnetism.ipynb)
- [Many-body correlated fermionic systems](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/quantum_interacting_fermions.ipynb)
- [Tensor networks for many-body quantum magnets](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/mps_quantum_magnets.ipynb)
- [Tensor networks for many-body correlated fermionic systems](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/mps_many_body_fermionic.ipynb)


# Capabilities #
- Models: spinless fermions, spinful fermions, spins, parafermions and bosons
- The same script runs on several backends (ITensor C++ v2 and v3, ITensors.jl in Julia, a pure-Python DMRG, and ED), so results can be compared between them
- The pure-Python backend can run on a GPU with JAX
- Calculations in a fixed sector of particle number and/or total Sz
- Ground state energy and wavefunction
- Excitation energies and excited wavefunctions
- Expectation values of any operator, including static correlation functions
- Time evolution of any state, with TDVP, TEBD or a Taylor expansion of the evolution operator
- MPS algebra: sums of MPS, application of operators, exponential and inverse
- MPO algebra: sums, products, trace, and trace of the inverse
- Dynamical correlation functions, with several methods chosen with `submode=`: kernel polynomial method (`"KPM"`), real-time (`"TD"`) and complex-time (`"TDZ"`) evolution, correction vector (`"CVM"`, `"CVM_explicit"`, `"CVMimag"`), root-N Krylov correction vector (`"ROOTN"`), sum over excited states (`"EX"`) and maximum entropy (`"maxent"`)
- Energy truncation of the KPM Chebyshev vectors (`kpm_energy_truncate`, only with `itensor_version=3` and `"python"`), for a higher energy resolution
- Distribution of any operator with the kernel polynomial method
- Hermitian and non-Hermitian iterative MPS eigensolvers, and detection of degeneracies
- Infinite chains with iDMRG and VUMPS: ground state, excitation gaps, dynamical correlators and the dynamical structure factor S(k,w), directly in the thermodynamic limit
- Finite temperature with METTS and thermal purification
- Entanglement entropy, reduced density matrices, mutual information and the central charge of a critical chain
- Topological diagnostics: symmetry sectors, edge and zero modes, Zak phase
- Non-Hermitian systems: ground and excited states, dynamical correlators, skin effect and degeneracies

# Examples

The snippets below are also in `examples/readme_examples/` as scripts
you can run. The `examples/` folder has more than 100 other scripts, one
per model or feature.

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

## Fixing the particle number of a fermionic chain
```python
from dmrgpy import fermionchain
n = 8 # number of spinless fermionic sites
# quantum-number sectors need itensor_version=3, "python", or mode="ED"
fc = fermionchain.Fermionic_Chain(n,itensor_version="python")
h = 0 # initialize Hamiltonian
for i in range(n-1): h = h + fc.Cdag[i]*fc.C[i+1] # first neighbor hopping
h = h + h.get_dagger() # make the Hamiltonian Hermitian
fc.set_hamiltonian(h) # set the Hamiltonian
for nf in range(n+1): # loop over particle-number sectors
  fc.set_conserved_sector(Nf=nf) # confine the calculation to Nf particles
  print("Energy with",nf,"particles",fc.gs_energy())
fc.set_conserved_sector() # back to the full Hilbert space
print("Global ground state",fc.gs_energy())
```

`set_conserved_sector` takes `Nf` (number of fermions, for spinless and
spinful fermions), `Sz` (total spin, given as the integer `2*Sz` like in
ITensor, for spins and spinful fermions) and `Nb` (number of bosons). You
can combine them, for example `set_conserved_sector(Nf=8,Sz=0)` for a
Hubbard chain. Every operator you use (the Hamiltonian, and the operators
in `vev` or in correlators) must conserve these quantities, otherwise you
get a `ValueError` with the name of the operator.

Exact diagonalization also works in a sector, so the same script runs
without any compiled backend:

```python
fc = fermionchain.Fermionic_Chain(8)
fc.mode = "ED" # this chain is solved by exact diagonalization
fc.set_hamiltonian(h)
fc.set_conserved_sector(Nf=3)
print(fc.gs_energy(mode="ED")) # the three-particle ground state
```

For `Spinful_Fermionic_Chain`, fixing the total `Sz` only works with ED.
The DMRG version of this chain is made of spinless fermion sites, which
only carry `Nf`.

## Ground state of an infinite chain with iDMRG
```python
from dmrgpy import infinitechain
import numpy as np
# a one-site unit cell = a uniform, infinite Heisenberg chain
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg" # growing algorithm ("vumps" is the default)
# C = the central cell, R = the next cell to the right
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h) # set the Hamiltonian of the infinite chain
ic.maxm = 30 # bond dimension
print("Energy per site",ic.gs_energy())
print("Exact (Bethe ansatz)",0.25-np.log(2))
print("<Sz(0)Sz(r)>",[ic.correlator("Sz",0,"Sz",r).real for r in range(1,5)])
```

An infinite chain is defined by its unit cell, not by its length, so all
results are in the thermodynamic limit (the energy is per site). Infinite
chains work with `itensor_version="python"` (the default for them) and
`itensor_version=3`. More examples are in `examples/idmrg/`.


# Choosing a backend #

All backends run the same Python code. You choose one with
`itensor_version`:

- `3`: ITensor v3 in C++ (`python install.py`). The default when it is compiled.
- `2`: ITensor v2 in C++ (`python install.py --itensor-version=2`). No conserved sectors, and time evolution only with the Taylor-expanded evolution operator.
- `"python"`: DMRG and TDVP written in Python. No compiler needed, but slower than C++. The default when the C++ part is not compiled. It can run on a GPU with JAX.
- `"julia_live"`: ITensors.jl through `juliacall` (`python install_julia.py`). Only part of the features are implemented, the others raise an error.

Exact diagonalization (`mode="ED"`) is used for small systems and to check
DMRG. If a C++ backend is not compiled, or the chain is too short for it,
dmrgpy uses ED instead and prints a note.

To change the backend of an existing chain:

```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
sc.setup_cpp(version=2) # switch to ITensor v2 (C++)
# sc.setup_python()     # or: pure Python, no compiler needed
# sc.setup_julia()      # or: ITensors.jl (Julia)
```

Most methods also accept `mode="ED"`, to check a single result with
exact diagonalization, for example `sc.gs_energy(mode="ED")`.
