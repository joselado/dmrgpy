# BLAS threads

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
