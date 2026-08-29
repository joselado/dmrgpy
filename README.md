## DMRGPY ##

# Summary #

This is a Python library to compute quasi-one-dimensional
spin chains and fermionic systems using matrix product states
with the density matrix renormalization group as implemented in ITensor
(C++ v2/v3, Julia, or a pure-Python/NumPy reimplementation that needs
no compiler; see "Choosing a backend" below). Chains can be finite or,
through iDMRG/VUMPS, infinite (see "Infinite chains" below). Most
of the computations can be performed both with DMRG and exact
diagonalization for small systems, which allows to benchmark the
results.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development.

# How to install #

## With pip ##

```bash
pip install dmrgpy
```

This installs the pure-Python part of the library, which needs no compiler
and no C++ toolchain. Both the exact-diagonalization backend and the
pure-Python DMRG/TDVP backend work out of the box; select the latter with

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
```

The compiled ITensor (C++) backends are deliberately *not* shipped on PyPI
-- they are built against a large vendored copy of ITensor and are compiled
in place from a clone of this repository, as described below. Optional
extras are available for the non-default backends and tools:
`pip install "dmrgpy[julia]"` (Julia/ITensors.jl backend), `[jax]` (the
GPU array backend of the pure-Python backend, see "Running the
pure-Python backend on a GPU" below), `[stats]`, `[parallel]`, or
`[full]` for all of them.

## Linux and Mac (compiled C++ backend) ##

Execute the script 

```bash
python install.py 
```

It first checks that everything needed to build is present (a C++
compiler, LAPACK/BLAS, `pybind11`, `make`), auto-detecting the right C++
compiler for the Python interpreter you used to run it (including a
conda-provided compiler if you're using Anaconda/Miniconda), and only then
compiles ITensor and the in-process C++ extension that uses it.

By default this builds the ITensor v3 backend (`itensor_version=3`). The
older ITensor v2 backend is still supported, and either or both can be
selected explicitly

```bash
python install.py --itensor-version=2     # only ITensor v2
python install.py --itensor-version=3     # only ITensor v3 (the default)
python install.py --itensor-version=both  # both, one after the other
```

If you need to override the auto-detected C++ compiler (version 6 or higher, or a matching Clang), pass it explicitly

```bash
python install.py --gpp=g++-6 
```

If LAPACK/BLAS can't be found automatically, point the script at OpenBLAS explicitly

```bash
python install.py --openblas --openblas_libdir=/path/to/openblas/lib --openblas_includedir=/path/to/openblas/include
```

To check what the script would use (compiler, BLAS/LAPACK) without
actually running the several-minutes-long ITensor build, use

```bash
python install.py --doctor
```

Alternatively, in case you just want to use the Julia version,
execute the script

```bash
python install_julia.py
```

Either installation script also adds dmrgpy to the PYTHONPATH of the
python interpreter you used to execute it. Afterwards you can import the
dmrgpy sublibrary that you want, for example

```python
from dmrgpy import spinchain
```

### HPC clusters with environment modules (e.g. Aalto's Triton) ###

On clusters using Lmod/environment modules, no C++ compiler may be on
PATH at all until a module is loaded. On Triton, for example:

```bash
module load scicomp-python-env gcc
python3 install.py
```

If you're using a module-provided Python distribution (e.g. Triton's
`scicomp-python-env`), its bundled conda-style compiler wrapper can
sometimes be present but broken (missing its `cc1plus` backend). If so,
`install.py` automatically falls back to the system/module-loaded
compiler and prints a note about it. Whichever compiler ends up being
used, and whichever modules/conda environment were loaded at install
time, must be loaded/activated again in every later session (including
job scripts) that imports `dmrgpy`, since the compiled extension depends
on them.

If you don't need the compiled C++ backend on the cluster, `pip install
dmrgpy` (see above) works directly in your Triton Python environment
without loading any modules at all, using the pure-Python/ED backends
instead.

## Windows ##

`pip install dmrgpy` (see above) works out of the box on Windows: the
package is pure Python, and the ED and pure-Python (`itensor_version=
"python"`) DMRG/TDVP backends need no compiler.

The compiled C++ backends (`install.py`) are POSIX-only and not supported
on Windows. If you need one of them, the easiest solution is to create a
virtual machine using [Virtual Box](https://www.virtualbox.org/),
installing a version of [Ubuntu](https://releases.ubuntu.com/20.04/) in
that virtual machine, and following the Linux instructions above.

# Tutorials #
You can find several tutorials [here](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024), in particular organized around the following topics

- [Many-body quantum magnets](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/quantum_magnetism.ipynb)
- [Many-body correlated fermionic systems](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/quantum_interacting_fermions.ipynb)
- [Tensor networks for many-body quantum magnets](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/mps_quantum_magnets.ipynb)
- [Tensor netowrks for many-body correlated fermionic systems](https://github.com/joselado/Advanced_Computational_Methods_Physics_2024/blob/main/jupyter-notebooks/mps_many_body_fermionic.ipynb)


# Capabilities #
- Possible models include spinless fermions, spinful fermions, spins, parafermions and bosons
- Interchangeable backends for the same API: compiled ITensor C++ (v2 and v3), a live Julia/ITensors.jl session, a pure-Python/NumPy reimplementation needing no compiler, and exact diagonalization; results can be cross-checked between any of them
- Optional GPU execution of the pure-Python backend, by putting its tensors on a JAX device instead of in host memory
- Optional conserved quantum-number sectors: confine an entire calculation to a fixed particle number and/or total Sz
- Ground state energy
- Ground state wavefunction
- Excitation energies
- Excited wavefunctions
- Arbitrary expectation values, including static correlation functions
- Time evolution of arbitrary states, with TDVP, TEBD or a Taylor-expanded evolution operator
- MPS algebra: sum of MPS, application of operators, exponential and inverse
- MPO algebra: sums, products, trace, trace of inverse for generic operators
- Dynamical correlation functions computed with the Kernel polynomial method
- Dynamical correlation functions with time dependent DMRG
- Several interchangeable algorithms for the same dynamical correlator, selected with `submode=`: the Kernel polynomial method (`"KPM"`), real-time (`"TD"`) and complex-time (`"TDZ"`) evolution, correction-vector methods (`"CVM"`, `"CVM_explicit"`, `"CVMimag"`), a root-N Krylov correction vector (`"ROOTN"`), an explicit sum over DMRG excited states (`"EX"`), and maximum-entropy reconstruction (`"maxent"`)
- Optional energy truncation of the Chebyshev vectors in the KPM correlator (`kpm_energy_truncate`, on `itensor_version=3` and `"python"` only), so the expansion window can be narrowed for higher spectral resolution
- Generic operator distributions computed with the Kernel polynomial method
- Iterative MPS Hermitian and non-Hermitian diagonalization solvers 
- Hermitian and non-Hermitian degeneracy detection
- Infinite chains with iDMRG/VUMPS: ground states, excitation gaps and dynamical correlators directly in the thermodynamic limit, including a momentum-resolved dynamical structure factor S(k,w) built from the quasiparticle branches' own exact spectral weights
- Finite-temperature calculations via METTS and thermal purification
- Entanglement entropy, reduced density matrices, mutual information and CFT central charge extraction
- Topological diagnostics: parity/symmetry sectors, edge/zero modes and Zak phase
- Non-Hermitian systems: DMRG ground and excited states, dynamical correlators, and skin-effect/degeneracy analysis

# Examples

The snippets below are mirrored as self-contained scripts under
`examples/readme_examples/`; the full `examples/` directory has 100+
further scripts, one per physical model or feature.

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
# quantum-number sectors need itensor_version=3 or "python"
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

`set_conserved_sector` accepts `Nf` (particle number, spinless and
spinful fermions), `Sz` (total spin projection, in ITensor's integer
`2*Sz` units, for spin chains and spinful fermions) and `Nb` (boson
number), alone or in combination, e.g.
`set_conserved_sector(Nf=8,Sz=0)` for a Hubbard chain. Every operator
used on the chain (the Hamiltonian, and anything passed to
`vev`/correlators) must itself conserve the requested quantities, or a
`ValueError` names the offending operator.

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

An `Infinite_Many_Body_Chain` is defined by a single unit cell rather
than by a total length, so everything it returns is a thermodynamic-limit
quantity (the energy is *per site*). See "Infinite chains" below for what
each backend and each `gs_method` supports.


# Choosing a backend #
The library uses ITensor in the background, and exposes the same
public API regardless of which backend actually runs the calculation.
Available backends are

- **ITensor v3 (C++)**, the default (`itensor_version=3`), compiled by
  `python install.py`. An in-process pybind11 extension. Real-time
  evolution defaults to TDVP (`tevol_method="TDVP"`), with TEBD
  (`"TEBD"`, nearest-neighbor Hamiltonians only), `"AUTO"` (TEBD when
  the Hamiltonian allows it, TDVP otherwise), one-site TDVP with global
  subspace expansion (`"TDVP_GSE"`) and the Taylor-expanded
  evolution-operator MPO (`"MPO"`) also available. This is the only
  compiled backend with quantum-number sectors
  (`set_conserved_sector`).
- **ITensor v2 (C++)** (`itensor_version=2`), compiled with `python
  install.py --itensor-version=2` (or `--itensor-version=both` to build
  v2 and v3 together). Real-time evolution always uses the
  Taylor-expanded evolution-operator MPO (no TDVP/TEBD support), and it
  has no quantum-number sectors.
- **Pure Python** (`pyitensor`, `itensor_version="python"`), a
  from-scratch NumPy/SciPy reimplementation of the ITensor v3 API subset
  dmrgpy needs: ground/excited-state DMRG, the same
  TDVP/TEBD/AUTO/TDVP_GSE/MPO time-evolution methods, KPM dynamical
  correlators, METTS, quantum-number sectors, iDMRG/VUMPS, and more.
  Needs no C++ compiler or compiled extension at all, at the cost of
  being slower than the compiled backends. It can also put its tensors
  on a GPU, see below.
- **Julia (ITensors.jl)** (`itensor_version="julia_live"`), a live
  in-process Julia session driven through `juliacall`, set up with
  `python install_julia.py` instead of `install.py`. It implements a
  subset of the API rather than all of it. For instance,
  `get_dynamical_correlator` supports `submode` `"KPM"`, `"CVM"`,
  `"TDZ"`, `"EX"` and `"maxent"` there, and raises `NotImplementedError`
  for the others. It also never falls back to another backend: a
  feature that isn't implemented there raises instead.
- **Exact diagonalization (ED)**, a pure Python/NumPy/SciPy fallback.
  Used automatically for small systems, and automatically in place of a
  C++ backend whose extension wasn't compiled (and for `itensor_version=3`
  on chains of fewer than 3 sites, which ITensor v3's two-site DMRG
  cannot handle), so a script keeps running, just slower, even
  without a full build. This automatic fallback covers the C++ backends
  only: the Julia backend never falls back, and neither does a chain with
  a conserved sector set, since ED would silently answer with the
  *global* ground state instead of the requested sector's.

The default backend is ITensor v3 (C++). Switch an existing chain to
another one with `.setup_cpp(version=2)`, `.setup_python()` or
`.setup_julia()`, for example

```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
sc.setup_cpp(version=2) # switch to ITensor v2 (C++)
# sc.setup_python()     # or: pure Python, no compiler needed
# sc.setup_julia()      # or: ITensors.jl (Julia)
```

and all subsequent computations on `sc` will be performed with that
backend. Most methods also accept a `mode="DMRG"|"ED"` kwarg to
cross-check a single call against ED without switching the chain's
backend, for example `sc.gs_energy(mode="ED")`.


# Infinite chains #

Besides the finite-chain classes above, `infinitechain.py` provides
`Infinite_Many_Body_Chain` (and its `Infinite_Spin_Chain` subclass), a
translationally-invariant chain defined by a single `n_uc`-site unit
cell instead of a fixed total length, solved directly in the
thermodynamic limit. It is a separate object rather than a
`Many_Body_Chain` subclass, since a fixed number of sites, an ED
cross-check and the finite-chain dispatch machinery have no meaning for
an infinite system. Hamiltonians are written with the L/C/R operator
convention (`SzC[i]` for the central cell, `SzR[i]` for the next cell
to the right, `SzL[i]` for the previous one); see the iDMRG example
above.

Two ground-state algorithms are available, selected with
`ic.gs_method`:

- `"vumps"` (the default), a direct fixed-bond-dimension solver
- `"idmrg"`, the growing infinite-size DMRG algorithm

Only `itensor_version="python"` (the default here) and
`itensor_version=3` are supported; anything else raises
`NotImplementedError`. What each combination supports:

| method | `gs_method` | `"python"` | `3` |
|---|---|---|---|
| `gs_energy` (energy per site) | either | yes | yes |
| `vev`, `correlator` | either | yes | yes |
| `excitation_energies`, `excitation_gap` | `"vumps"` | yes | yes |
| `spectral_weights`, `dynamical_structure_factor` | `"vumps"` | yes | no |
| `local_excitation_gap` | `"idmrg"` | yes | `window=0` only |
| `td_dynamical_correlator` | `"idmrg"` | yes | yes |
| `kpm_finite` | n/a | yes | yes |

`kpm_finite` is a finite-window approximation: it builds a window of
`n_window` unit cells with open boundaries and runs the ordinary
finite-chain KPM machinery on it, so it depends on neither the chain's
`itensor_version` nor its `gs_method`. `td_dynamical_correlator`, by
contrast, uses infinite boundary conditions and real-time evolution,
while `spectral_weights`/`dynamical_structure_factor` give exact
thermodynamic-limit delta peaks from the quasiparticle branches
themselves. Unit cells with `n_uc>2` need `gs_method="vumps"`; the
growing `"idmrg"` algorithm rejects them.

Runnable examples live in `examples/idmrg/`.

# Running the pure-Python backend on a GPU #

`itensor_version="python"` has a second, orthogonal axis: *which array
library* its tensors are made of. `pyitensor/backend.py` selects it
process-wide, before a chain is built:

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")   # or "numpy", the default
print(backend.device_info()) # e.g. "jax: cuda:0"
backend.set_pad_bonds(20)    # freeze every bond dimension (see below)

from dmrgpy import spinchain          # nothing else changes
spins = ["S=1/2" for i in range(10)]
sc = spinchain.Spin_Chain(spins,itensor_version="python")
h = 0
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20                 # must match set_pad_bonds above
print("Ground state energy",sc.gs_energy())
```

This is a device axis, not a new backend: it changes nothing about
backend dispatch, only where the arrays live. It needs `jax` with its
CUDA plugin (`pip install "dmrgpy[jax]"` installs `jax` itself). The
compiled C++ backends and the ED path are unaffected by it.

**Whether it helps depends on bond dimension, not on chain length.**
Each device operation carries a fixed dispatch cost regardless of how
small the tensor is, so many small tensors lose and few large ones win.
Measured on an H200: below a bond dimension of ~120 the GPU is *slower*
than NumPy, ~120-160 is break-even, and above that it pays off (~5x at
240, and ~20x at 480 on a calculation that genuinely needs that bond
dimension; a uniform Heisenberg ground state converges well below it
and cannot show a speedup at any `maxm`). A naive design that converted arrays per call was
5-11x slower than the host path, which is why there is exactly one
conversion point in the engine and everything downstream of it stays on
the device.

Two knobs matter, and they work together rather than separately:
`backend.set_pad_bonds(K)` freezes every bond at `K` so the engine stops
minting a new tensor shape (and hence a new compiled kernel) every time
a bond dimension changes, and `backend.set_jit()` fuses the hot inner
kernels. The default `set_jit("auto")` turns jitting on exactly when
`set_pad_bonds` has been set. On a cold run (a script that starts,
computes and exits) the two together were measured at 6.4-12.1x, where
either alone was worth only 2.4-3.8x; inside a long-running process the
trade reverses and `set_jit(True)` *without* padding is the better
choice.

The full measured comparison, including which models are and are not
meaningful GPU benchmarks, is in `docs/gpu_cpu_performance.md`, and the
design in `docs/documentation.md`.

# Tests and benchmarks #

The `tests/` directory holds a pytest suite that cross-checks the DMRG
backends against ED (and against each other) on small systems:

```bash
python run_tests.py     # or: pytest tests
```

`benchmarks/run_benchmarks.py` answers a different question: not "is
this backend correct" but "which backend is fastest, and by how much".
It sweeps a uniform S=1/2 Heisenberg chain over a range of chain
lengths, timing a ground-state energy, a static correlator and a KPM
dynamical correlator on every backend available in the current
environment, cross-checking each against ED, and writes a LaTeX report
(compiled to PDF if `pdflatex` is available) with tables and plots:

```bash
cd benchmarks && python run_benchmarks.py --help
```

Note that DMRG is made of a great many *small* dense linear-algebra
calls, so more BLAS threads is usually slower rather than faster, and
under oversubscription (a shared node, or several dmrgpy runs at once)
dramatically so. Pin threads before importing numpy when timing
anything:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 script.py
```

After running examples, `python clean.py` removes the working
directories and stray output files they leave behind.
