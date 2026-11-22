---
paths:
  - "pyproject.toml"
  - "release_pypi.py"
  - ".gitignore"
  - "src/dmrgpy/juliapkg.json"
  - "src/dmrgpy/mpsjulialive/ITensorNHDMRG/**"
  - "src/dmrgpy/mpsjulialive/juliasession.py"
  - "src/dmrgpy/algebra/kpmextrapolate.py"
---

# Packaging for PyPI (`pyproject.toml`)

These are the notes CLAUDE.md used to carry under "Packaging / PyPI"; they
load when `pyproject.toml`, `release_pypi.py`, `.gitignore` or
`juliapkg.json` is read or edited.

`pyproject.toml` builds a **pure-Python** distribution: a `py3-none-any`
wheel (~450KB) containing the ED backend, the `itensor_version="python"`
pyitensor DMRG/TDVP backend, and the Julia bridge. It deliberately ships no
C++ at all. `[tool.setuptools.packages.find]` excludes `dmrgpy.mpscpp2` and
`dmrgpy.mpscpp3`, which are ~400MB of vendored ITensor plus a build that
needs a compiler, LAPACK/BLAS and pybind11, none of which belongs in a pip
install. Users who want the compiled backends clone the repo and run
`install.py`; without them everything still works via `"python"`/ED
(`mode.py`), so the wheel is fully usable as installed.

```bash
python -m build                  # -> dist/*.whl + dist/*.tar.gz
python -m twine check dist/*     # metadata/README validation
python -m twine upload --repository testpypi dist/*   # dry run on TestPyPI first
python -m twine upload dist/*    # the real thing
```

Things in here that are load-bearing and easy to break:

- **`numba` is a core dependency, not an extra.** `algebra/kpmextrapolate.py`
  imports it at module level and sits on the unconditional
  `manybodychain -> dynamics -> kpmdmrg -> algebra.kpm` import path, so
  without it even `import dmrgpy.spinchain` raises `ModuleNotFoundError`.
  This is unrelated to pyitensor's own opt-in numba/JAX matvec kernel
  (`pyitensor/kernels.py`), which is properly `try/except`-guarded and
  genuinely optional.
- **`requires-python = ">=3.10"` is set by numba**, which declares
  `>=3.10`. Claiming a lower floor does not make dmrgpy work on 3.9; pip
  just backtracks onto an ancient, untested numba/numpy pair.
- **`setuptools>=77` in `[build-system]`**, not the more common 61: the
  `[project]` table uses PEP 639 licensing (bare SPDX `license` string plus
  `license-files`), which older setuptools rejects outright.
- **The `julia` extra is `juliacall`, not `julia`.** The live backend
  migrated from PyJulia to juliacall/PythonCall.jl. The only surviving
  `import julia` is `juliarun.py`, part of the legacy subprocess path that
  is not reachable from the public API, so PyJulia is intentionally not
  pulled in.
- **`src/dmrgpy/juliapkg.json` must stay in `package-data`.** juliacall/
  juliapkg reads it at import time to resolve and precompile ITensors,
  ITensorMPS 0.4, KrylovKit (used directly by `mpsjulialive/vumps.jl` and
  `cbe.jl`) and VectorInterface; if it is missing from the wheel,
  `pip install dmrgpy[julia]` yields a Julia backend with no ITensors.jl.
  `.gitignore`'s unanchored `julia*` pattern would also swallow it, hence
  the explicit `!src/dmrgpy/juliapkg.json` un-ignore. Do not remove either
  one.
- **ITensorNHDMRG.jl is a copy, not a dependency.** Its compat entry
  refuses ITensorMPS 0.4, so its `src/` lives unchanged in
  `src/dmrgpy/mpsjulialive/ITensorNHDMRG/` (LICENSE and commit hash in
  `UPSTREAM.md`), `juliasession.py` includes it as a module before any
  other Julia file, and `get_gs.jl`/`nhdmrg.jl` load it with
  `using .ITensorNHDMRG`. Its files are in `package-data`; its
  VectorInterface dependency is in `juliapkg.json`. The folder name must
  not start with `julia`, because of the same `.gitignore` pattern.
- **The classifiers say `OS Independent`**: the pip-installed package (ED
  and `itensor_version="python"` backends) works on Windows. Nothing in the
  Python layer shells out any more (`sites.py::initialize()` used to run
  `mkdir -p` via `subprocess`, removed in `1f7b7a9` along with the
  `.mpsfolder` it created and `Many_Body_Chain.clean()`, whose only job was
  removing that folder). The compiled C++ backends (`install.py`) remain
  POSIX-only and are not part of the wheel regardless.
- Every module in the wheel must at least byte-compile, or pip prints a
  `SyntaxError` on install; `python -m compileall` over the installed
  package is the cheap check.

Verify a release candidate by installing the built wheel into a clean
virtualenv (not from `src/`, so import bugs and missing package-data
actually surface) and cross-checking the two shipped backends against each
other: `itensor_version="python"` DMRG against `mode="ED"` on a small chain
should agree to ~1e-15. `dist/`/`build/` are gitignored; build artifacts are
not committed. A given version can never be re-uploaded to PyPI (or
TestPyPI) once published, so bump the version rather than trying to replace
a release.
