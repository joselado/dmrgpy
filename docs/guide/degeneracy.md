# Ground-state degeneracy

A degenerate ground state signals a broken symmetry or a topological
phase. Examples are the two ordered states of an Ising chain and the four
edge states of an open Haldane chain, whose splitting is shown in
[Excited states](excited-states.md#low-lying-spectrum). We count the levels close to $E_0$,

$$g\approx\sum_i\exp\!\left[-\left(\frac{(E_i-E_0)^2}{\delta}\right)^2\right],$$

adding excited states $E_i$ until the count stops growing:

```python
g = sc.get_gs_degeneracy()             # close to an integer
g = sc.get_gs_degeneracy(delta=1e-4)   # a stricter window
```

A level counts when $|E_i-E_0|$ is smaller than about $\sqrt{\delta}$,
which is 0.1 at the default `delta=1e-2`. A result that is not close to
an integer means that a level sits at the edge of the window. On a finite
chain a degeneracy is often split by an amount that decays with the
length, so choose `delta` relative to that splitting. `mode="ED"` uses
exact diagonalization.
