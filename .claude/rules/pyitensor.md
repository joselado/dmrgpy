---
paths:
  - "src/dmrgpy/pyitensor/**"
  - "tests/test_mpo_automaton_builder.py"
  - "tests/test_audit_*pyitensor*.py"
  - "tests/test_audit_*mpobuilder*.py"
  - "tests/test_*gpu*.py"
---

# The pure-Python backend (`pyitensor/`)

These are the notes CLAUDE.md used to carry inside its "DMRG, pure Python"
bullet; they load when a file under `pyitensor/` is read, written or
edited. What `pyitensor/` is, and where it sits in the dispatch, is still
in CLAUDE.md.

## A new method also runs on the GPU

New `pyitensor` methods are written to run on the GPU as well, whenever
that is possible, locked by the author on 2026-10-24 ("whenever we
implement a new method, if possible, we should write it so that it also
works on the GPU", for pyitensor). The way we do this is by taking the
array namespace from `pyitensor/backend.py` (`xp = backend.xp()`, once
per function) instead of `np`, moving host data to the backend at one
entry point, bringing back only scalars and small matrices
(`backend.scalar`, `backend.to_host`), and never writing into an array
in place without `backend.setblock`/`backend.setrow`; a new method is
then checked under `backend.set_backend("jax")` against its NumPy
result. `pyitensor/task.py` is a worked example. When a step cannot run
on the device, say so in the method's docstring.

## The MPO builder

MPO construction is not one of the simplifications `"python"` takes:
`pyitensor/mpobuilder.py`'s `to_mpo` assembles the finite-state machine over the terms' partial products directly (the shape ITensor's
own `toMPO(...,{"Exact",true})` produces), O(L^2) instead of the old
O(L^4) sum-and-compress, keeping two sweeps afterwards purely to honour
`cutoff`/`maxdim`, the first exact and only the return one truncating (a
truncating first sweep compressed any operator whose coefficients were
all below about 2e-7 into a different one), and an operator whose largest
coefficient is below 1 is built at unit scale and scaled back, like
v2/v3's `to_mpo_unit()`. Two rules in there are load-bearing and easy to
get subtly wrong: the coefficient goes on the transition into the final
state (so terms differing only by a coefficient still share partial
states), and transitions into that state accumulate while structural
ones are assigned (two syntactically identical terms trace the same path,
so their coefficients must sum while their shared structure must not be
written twice). The old construction is kept as `_sum_of_term_mpos`, the
reference the machine is tested against in
`tests/test_mpo_automaton_builder.py`.
