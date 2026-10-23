# Performance

On this page we make a calculation faster without changing its result:
one BLAS thread per process, several points of a scan in parallel, and a
GPU for the pure-Python backend.

## BLAS threads

DMRG makes many small matrix operations, and on small matrices several
BLAS threads are slower than one. It matters most when the machine is
shared, for example on a cluster node or with several runs in parallel,
where a calculation can become 3 to 17 times slower. Pin one thread
before Python imports numpy:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 myscript.py
```

Inside a script, `blasthreads.limit(1)` does the same for one block. It
needs the package `threadpoolctl`, and without it does nothing and gives a
warning:

```python
from dmrgpy import blasthreads
with blasthreads.limit(1):
    e0 = sc.gs_energy()
```

dmrgpy never changes the number of threads by itself.

## Parallel scans

The points of a parameter scan are independent, so they can run in
parallel processes. `parallel.pcall(f, args)` applies `f` to every entry
of `args` and returns the list of results, in order:

```python
from dmrgpy import parallel
parallel.cores = 4                 # processes; 1 runs in serial

def energy(b):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    sc.set_hamiltonian(h0 + b*h1)
    return sc.gs_energy()

es = parallel.pcall(energy, bs)    # one process per point
```

`f` builds its own chain, since a chain holds the state of one solve and
is not shared between processes. Pin one BLAS thread per process as
above, or the processes compete for the cores. `parallel.cores` is 1 by
default, and the processes need the package `multiprocess`; without it
`pcall` runs in serial. See `examples/utilities/parallel`.

## GPU

The `itensor_version="python"` backend can keep its tensors on a GPU. It
needs `jax` with CUDA support and a GPU with fast double precision, as in data-centre cards. A consumer card with slow double precision does not give a speedup. The switch is global and is set before
building a chain:

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")       # "numpy", the default, goes back to the CPU
print(backend.device_info())     # for example "jax: cuda:0"
```

The rest of the script does not change, and the results agree with the
CPU ones to about $10^{-11}$. A GPU helps only when the bond dimension is
large: below about 120 it is slower than the CPU, and at 240 it is about
5 times faster. So it pays for real-time evolution, where the bond
dimension grows, and rarely for a 1D ground state. Infinite chains and
non-Hermitian DMRG run but do not get faster. For a script that does one
calculation and exits, `backend.set_pad_bonds(maxm)` fixes all tensor
shapes so that each kernel is compiled only once. With padding on, the
default `backend.set_jit("auto")` also fuses the kernels. For a KPM
correlator, set `kpmmaxm` equal to `maxm`. See `docs/gpu_cpu_performance.md`
for the measurements.
