# GPU

The `itensor_version="python"` backend can keep its tensors on a GPU. It
needs `jax` with CUDA support. The switch is global and is set before
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
correlator,
set `kpmmaxm` equal to `maxm`. See `docs/gpu_cpu_performance.md` for the
measurements.
