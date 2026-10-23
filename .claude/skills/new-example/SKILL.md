---
name: new-example
description: Add or revise a script under examples/ in dmrgpy. Use this whenever the user asks for an example, a demo, a script showing how to use a feature, a backend-versus-backend comparison, or a regression script in the spirit of tests/, and also when a new feature has just landed and an example for it is the natural next step. The examples tree has rules that are easy to miss (it must plot, it must pin its own sweep schedule, it must import this worktree and not site-packages) and an audit found 105 scripts that had quietly broken the first of them.
---

# Adding an example

`examples/` is the second regression surface, and the thing that separates it
from `tests/` is that a human is expected to look at the result. That is the
premise behind most of the rules below.

## Where it goes

`examples/<theme>/<name>/main.py`, one self-contained script per folder. The
themes are `groundstate`, `staticcorrelators`, `dynamical_correlator`,
`spin_models`, `fermion_models`, `boson_models`, `time_evolution`,
`excited_states`, `entanglement`, `topological`, `kondo`, `non_hermitian`,
`finite_temperature`, `magnetization`, `idmrg`, `backend_comparison`,
`algebra`, `parity`, `pyitensor`, `utilities` and `readme_examples`. Pick the
one that matches the physics, not the backend, except for a script whose whole
point is comparing backends. A comparison of the two compiled backends is named
`v2_VS_v3_<thing>` and lives in its theme folder. A script that mirrors a
snippet from `README.md` goes in `readme_examples/` and has to stay in step with
the README.

## The script

Open with the path line the tree already uses, which is depth three:

```python
# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
```

**Keep it short and simple** (the author's standing rule, 2026-10-23: "from
now on, i want all the examples like this"). An example is 40 lines, 50 at
most, blank lines and comments included. It is written the way a condensed
matter physicist who is not a professional programmer and not a native English
speaker writes it: plain English, short sentences, simple code. It focuses on
the physics and on the one capability it demonstrates, and leaves out
everything else: history, audit notes, implementation explanations, timing
code, loops over backends, helper functions, try/except, long docstrings. The
two approved models are `staticcorrelators/spin_correlator_short/main.py` and
`dynamical_correlator/spin_spectrum_short/main.py`; read them before writing.

Then a header comment of three to six lines saying the physics: what is
computed, the formula if it is short, what it is used for or what measures it.
Build the Hamiltonian with explicit loops, with an inline comment only where a
line is not obvious (`sc.maxm = 40 # bond dimension`). Build the smallest
system that shows the effect.

**Pin the sweep schedule.** Set `nsweeps`, `maxm` and any tolerance explicitly
rather than inheriting the library defaults, and assert at a tolerance the
pinned schedule actually reaches. `boson_models/v2_VS_v3_boson` inherited the
defaults and its assert then fired non-deterministically on a clean tree, in
two runs of five for one reader and five of seven for another: the two backends
were converging to the same answer, just not far enough to be compared where it
compared them. It now pins `nsweeps=80`, `maxm=100` and asserts at 1e-5.

**Plot what you computed.** Any script that produces a sequence of values ends
with a `matplotlib` figure of it: a quantity against time, site, frequency or a
coupling; the reference as a thick line and the other method dashed on top of
it. One simple figure, `plt.plot`, axis labels, a legend and `plt.show()`, no
`savefig` and no styling.

A script whose entire output is one scalar is the single exception, and it is a
last resort rather than a default: before taking it, look for the cheap axis
that is almost always already there, such as chain length, a field or coupling,
time, bond dimension, site index, or ED against DMRG. The audit that found 105
print-only scripts found that most of them had already computed the array and
simply never drew it.

**Assert where the script is a regression guard.** Several examples carry real
asserts and are meant as tests in the same spirit as `tests/`, for instance
`time_evolution/tdvp_VS_ED_time_evolution`,
`backend_comparison/backend_switch_consistency`,
`staticcorrelators/static_correlator_VS_ED` and
`entanglement/entanglement_entropy_VS_ED`. Those are the templates to copy. The
assert stays as the pass or fail guard and the plot goes in after it, so the one
script serves as an automated check and as a visual sanity check when run by
hand. In the short style this is one printed "Max difference" number and one
`assert` on it.

## Run it before declaring it done

From the example's own directory, with threads pinned, this worktree's `src`
forced onto the path, and a non-interactive backend so `plt.show()` does not
block:

```bash
cd examples/<theme>/<name>
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg \
  PYTHONPATH=<this worktree>/src python3 main.py
```

The `PYTHONPATH` is not belt and braces. The script's own line *appends* to
`sys.path`, and `site-packages/dmrgpy` is typically a symlink into one specific
checkout, so a plain `python3 main.py` from another worktree silently exercises
the wrong code with no error at all.

Afterwards run `python clean.py` from the repo root, which removes the generated
working directories and stray `ERROR`/`*.OUT` files the run leaves behind.

Say what the numbers came out as when you report back, not just that it ran.
