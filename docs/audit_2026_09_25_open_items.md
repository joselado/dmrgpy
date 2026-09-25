# Audit, 2026-09-25: ten open items of the five hole hunts

This record is the sixth pass over the Python layer, and the first that hunted
nothing new: it took the ten most urgent of the items the five previous records
(`audit_2026_08_hole_hunt.md`, `audit_2026_09_hole_hunt.md` and the three
`audit_2026_09_24*_hole_hunt.md`) had left open in a Status line or recorded as
an unreviewed lead, ranked by how silently wrong each one is and how many callers
it reaches. It ran on 2026-09-25 on `8dd2198` (clean tree, both compiled
extensions current, `make -n pybind` reporting nothing to do in `mpscpp2` and
`mpscpp3`), as a workflow of eleven agents: one fix agent per cluster, which
reproduced each item on `8dd2198` before touching it and kept a repro whose
before and after outputs are spliced below verbatim, then one reviewer per
cluster briefed to refute both that each defect was real as stated and that
each fix holds, then a repair pass wherever the reviewer found a fix incomplete.
Three clusters ran in parallel in the shared tree, and the C++ half of the
small-units item ran last and alone, since it rebuilt both extensions.

This file is the evidence, not a task list, the same convention as the five
earlier records. Every item is fixed or narrowed and fixed, each entry carries a
`**Status**` line, and the regressions live in
`tests/test_audit_2026_09_25_<cluster>.py` for `construction`, `session` and
`scale`. None of the ten was refuted; the reviewers narrowed four (`maxde`
lost its example and user-guide halves, which `8dd2198` had already fixed;
`init-kwargs` lost the claim that `mode=` behaves as an assignment in the first
pass, which the repair then made true; `generalized-cache` gained a
non-Hermitian half the first pass had missed; `small-units` lost the claim that
one scale per operator reaches every bond). The review of a cluster that went
through a repair pass is a review of the first pass: the repairs were checked by
the main agent against the reviewer's report and by the test runs below, and
were not handed to a second reviewer. What the reviewers and the fix agents
found and did not fix is collected in "Left open, and new leads".

## The ten items at a glance

| # | Item | Recorded in | Numbers change on fix | Regressions |
|---|---|---|---|---|
| 1 | `init-kwargs` | 2026-08 record (a reviewer's sharpening near line 1723); 2026-09-24b New leads | yes, every caller that passed a setting | `test_audit_2026_09_25_construction.py` |
| 2 | `small-units` | 2026-09-24c New leads | yes, v2/v3 below a largest coefficient of about 1e-6 | `test_audit_2026_09_25_scale.py` |
| 3 | `clean-threshold` | 2026-09-24c New leads | yes, operators below 1e-8 | `test_audit_2026_09_25_scale.py` |
| 4 | `julia-warm-start` | 2026-09-24b finding 12 Status (open) | yes, julia_live setters | `test_audit_2026_09_25_session.py` |
| 5 | `get-gs-wf0` | 2026-09-24c New leads | yes, get_gs(wf0=x) on a solved chain | `test_audit_2026_09_25_construction.py` |
| 6 | `generalized-cache` | 2026-09-24c New leads | yes, readers after a correlator | `test_audit_2026_09_25_session.py` |
| 7 | `thermal-bypass` | 2026-09-24c New leads | yes, MBChain readers | `test_audit_2026_09_25_session.py` |
| 8 | `nh-injected` | 2026-09-24c New leads | yes, gs_energy(wf0=x, reconverge=False) | `test_audit_2026_09_25_session.py` |
| 9 | `maxde-per-site` | 2026-09-24c New leads | no | `test_audit_2026_09_25_session.py` |
| 10 | `rootn-ij` | 2026-09-24 finding 10 Status (open) | yes, lower-level ROOTN at sites other than (0,0) | `test_audit_2026_09_25_construction.py` |

The clusters and their files: `construction` (1, 5, 10: `manybodychain.py`'s
`__init__` and `get_gs`, `sites.py`'s new `check_settings`, `bosonchain.py`'s v2
refusal message, `rootndmrg.py`, `examples/utilities/multioperator_density`, and
`thermal.py`'s constructor through the addendum under item 1); `session` (4, 6,
7, 8, 9: `groundstate.py`, `thermal.py`'s `get_gs`, `nhdmrg.py`,
`nonhermitian/kpm.py`, `mpsjulialive/{groundstate,dynamics,excited,mps}.py`,
`manybodychain.py`'s `gs_energy`, setters and `gs_energy_fluctuation`
docstrings, `examples/groundstate/GS_enforce_maximum_fluctuation`); `scale` (3:
`multioperator.py`, `multioperatortk/canonical.py`, `meanfield.py`, a comment in
`mpsalgebra.py`); and `scale-cpp` (2: `mo_terms.h` and `chain_session.h` in both
`mpscpp2` and `mpscpp3`, with a rebuild of both extensions).

## Scope

Out of scope by construction: the vendored ITensor (`mpscpp2/ITensor/`,
`mpscpp3/ITensor/`), whose absolute thresholds item 2 works around from
dmrgpy's own headers without editing them; every open item of the earlier
records not among the ten (the `kpm_energy_truncate` window, the sliver below
`kpm_scale=1/2`, TDZ at `tdvp_gse_sweeps=0` under padding, `mode="ED"` TD's
swallowed keywords, the shared adjoint helper, the `vx_deterministic_start`
comment, the third pass's four deliberate "left as they were", and the rest of
the unreviewed leads); and `julia_live` except for item 4, which is Julia's.

Every Python process ran through the three-slot runner under "Shared helpers",
threads pinned and this checkout's `src` first on `PYTHONPATH`. The before runs
used a snapshot of `8dd2198`, made before any agent started with `git archive
8dd2198 -- src` (the vendored ITensor and TDVP folders excluded) and copies of
the two compiled extensions as they were, so every before output is the code
under audit and every after output the working tree. For item 2 the C++ stage
also ran a "hybrid" tree, the fixed Python with the old extensions, so that the
C++ change is measured on its own; `<scratch>` in the paths below is the
session's scratch folder, which does not outlive it, and is why every script and
output is carried inline.

## Items

### 1. Every chain constructor dropped its keywords, so `Spin_Chain(sites, maxm=4, nsweeps=10)` ran at maxm=30 and returned -5.1420906326 on a 12-site Heisenberg chain where the request gives -5.1323602278, and a misspelled keyword was accepted

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `construction`

**Status**: FIXED. Many_Body_Chain.__init__ now checks its keywords with sites.check_settings() before any session is built and assigns every one of them after initialize(), mode= included, so a constructor keyword is exactly an assignment made right after construction: a chain built with mode="ED" keeps its session, as sc.mode="ED" afterwards does, and sc.mode=None goes back to DMRG (-1.6282694864 on both paths of the 4-site "python" chain). A setting is what kpm_finite accepts for window_chain_kwargs (second-pass finding 6): a public attribute the chain already has that is not a method. Refused besides, each with a message naming its own entry point, are the chain's state (hamiltonian, conserved_sector, wf0, e0, computed_gs, ns, Id, the ground-state bookkeeping flags, the ED cache, and fermionic and use_ampo_hamiltonian, which the model class sets after the base constructor and would overwrite) and every attribute the model class built before calling the base constructor, so Fermionic_Chain(4, N=5) and Parafermionic_Chain(3, Sig=1) raise TypeError naming the key; mode= is checked for its value. Every forwarding constructor inherits this without an edit of its own. Thermal_Spin_Chain hands every setting to its MBChain, except that its own hardcoded mode="DMRG" overwrites MBChain.mode on the first get_gs(), as on the parent, so Thermal_Spin_Chain(..., mode="ED") still runs DMRG (left open, the edit is in the notes). The v2 boson refusal no longer names mode="ED" as a way out, since its check runs before a constructor keyword is applied. The first fix pass applied mode before initialize(), which left a mode="ED" chain with no session, so a later sc.mode=None and Thermal_Spin_Chain(..., mode="ED").get_gs() raised AttributeError where the parent ran; that ordering is reverted, and with it the claim that Bosonic_Chain(3, maxnb=[6]*3, itensor_version=2, mode="ED") builds (it raises the v2 ValueError, as on the parent). The one in-tree caller that passed a non-setting, examples/utilities/multioperator_density (spinful=False), was fixed and rerun (ED -3.876034843149142, DMRG -3.876034843149144). Pinned by tests/test_audit_2026_09_25_construction.py: test_a_constructor_setting_takes_effect and test_a_keyword_that_is_not_a_setting_raises_naming_every_one over nine constructors, test_private_names_methods_and_names_the_chain_lacks_are_not_settings, test_the_chain_state_is_not_a_setting, test_what_the_model_class_builds_is_not_a_setting, test_a_constructor_setting_goes_through_its_setter, test_a_constructor_setting_is_the_same_as_assigning_it, test_mode_at_construction_keeps_the_session over nine constructors, test_mode_at_construction_can_be_left_again, test_thermal_chain_with_mode_at_construction_still_solves (which passes on the parent and pins the absence of the first pass's regression), test_the_v2_boson_refusal_does_not_offer_mode_ed and test_an_unknown_mode_is_refused_at_construction. NUMBERS CHANGE for every caller that passed a setting to a chain constructor: on a 12-site open S=1/2 Heisenberg chain on itensor_version=3, Spin_Chain(["S=1/2"]*12, itensor_version=3, maxm=4, nsweeps=10).gs_energy() goes from -5.1420906326 (the maxm=30 answer) to -5.1323602278, the number the same request made by assignment returns. A mode="ED" keyword now takes effect for every read that does not pass mode= itself: on the 4-site open S=1/2 Heisenberg chain with 0.2*Sz_0, the maximum of get_dynamical_correlator(name="ZZ", es=linspace(0,3,7), delta=0.3) goes from 0.205144 ("python") and 0.205151 (v3), both of them DMRG, to 0.205190, the ED number. No in-tree number moves: the four-correlation example's Fermionic_Chain(n0, mode="ED") now keeps its mode, and every call there already passed mode="ED" (sum|ct_ed| 51.433292772677 and E0 -2.450108689391 on both trees).

**Status addendum** (the main agent, after the workflow): the `Thermal_Spin_Chain` edit this item left open was applied as specified: `Thermal_Spin_Chain.__init__` takes `mode="DMRG"`, validates it with `mode._check_mode`, keeps it as the wrapper's own `self.mode` and hands it to `MBChain` at construction too, so `Thermal_Spin_Chain(..., mode="ED")` solves by ED where it ran DMRG, on the parent as after the first pass. Pinned by `tests/test_audit_2026_09_25_construction.py::test_thermal_chain_mode_at_construction_is_the_mode_it_solves_with` (the mode on the wrapper and on `MBChain`, the default still `"DMRG"` (`None` since the 2026-09-25b fix pass, finding 3, where a `"DMRG"` pin stopped overriding a call's `mode="ED"`), and an unknown mode refused with `ValueError`); the construction, session and 2026-08 regression files pass with it (119 passed). NUMBERS CHANGE for every caller of `Thermal_Spin_Chain(..., mode="ED")`, which now gets the ED answer; no in-tree caller passes it.

**Recorded in**: 2026-08 record (a reviewer's sharpening near line 1723); 2026-09-24b New leads. Recorded as a reviewer's sharpening in docs/audit_2026_08_hole_hunt.md around line 1723 and as the first bullet of 'New leads, not reviewed' in docs/audit_2026_09_24b_hole_hunt.md; never a finding

**Where**: src/dmrgpy/manybodychain.py:289 on 8dd2198 (`self.initialize(**kwargs)`) and :291 (`initialize(self,**kwargs)`, which reads none of them); reached through every constructor that forwards its keywords there: spinchain.py:54, fermionchain.py:23,211,242,501,649, bosonchain.py:50,174, parafermionchain.py:19-23, mixedchain.py:106, spinfermionchain.py:10, thermal.py:20

Many_Body_Chain.__init__ took **kwargs and handed them to initialize(), whose own **kwargs has no consumer, so every keyword a caller gave a chain constructor, besides itensor_version and the model's own (maxnb, Z, T), was discarded and the chain ran at the defaults: Spin_Chain(['S=1/2']*4, maxm=50, bogus_key=1) built with maxm 30 and no error, and so did all nine forwarding constructors, Thermal_Spin_Chain included. The layer above believes the opposite, since the Spinful_Fermionic_Chain docstring (fermionchain.py:225-236) says the forwarding is there so that 'itensor_version (and anything else Many_Body_Chain accepts) can be set at construction', and itensor_version is exactly the one keyword bound before the drop. The reviewer measured the cost on every backend: on a 12-site Heisenberg chain with 0.2*Sz_0 asked for maxm=4, nsweeps=10, E(at construction) - E(assigned) = -9.635e-03 on v2 and v3 and -9.632e-03 on "python". No test passed a setting at construction and every example assigns sc.maxm afterwards, which is why nothing caught it; an AST scan of src/, tests/, examples/ and benchmarks/ found only two in-tree calls passing a keyword that is neither itensor_version nor the model's own, Fermionic_Chain(n0, mode="ED") in examples/staticcorrelators/four_correlation_tensor_sweep_VS_full and Fermionic_Chain(n, spinful=False) in examples/utilities/multioperator_density, both silently dropped. Infinite_Many_Body_Chain has no **kwargs and already raised TypeError on maxm=, so the infinite chains were never affected and are left as they are.

**Expected**: A keyword that names a chain setting, mode= included, takes effect exactly as if it had been assigned right after construction, through its setter where it has one (maxm=0 raises, as sc.maxm=0 does), and anything else, the chain's state and the operator lists the model class builds included, raises TypeError naming every offending key, sorted.

Repro (<scratch>/construction/01_init_kwargs.py (before: DMRGPY_SRC=<parent>/src run3.sh 01_init_kwargs.py | tee 01_init_kwargs.before.out; after: run3.sh 01_init_kwargs.py | tee 01_init_kwargs.after.out). The first pass's shorter script and its outputs are kept as 01_init_kwargs.firstpass.{py,before.out,after.out}.):

````python
# init-kwargs: a keyword handed to a chain constructor is dropped.
import functools
import numpy as np
import dmrgpy
from dmrgpy import (spinchain, fermionchain, bosonchain, parafermionchain,
                    mixedchain, thermal, infinitechain)
print = functools.partial(print, flush=True)
print("dmrgpy from", dmrgpy.__file__)


def build(label, f):
    try:
        c = f()
    except (TypeError, ValueError) as e:
        print("%-26s %s: %s" % (label, type(e).__name__, e))
        return None
    return c


# 1. the call of the record
c = build("Spin_Chain", lambda: spinchain.Spin_Chain(["S=1/2"]*4, maxm=50,
                                                    bogus_key=1))
if c is not None:
    print("Spin_Chain(maxm=50, bogus_key=1): built, maxm =", c.maxm,
          " hasattr(bogus_key) =", hasattr(c, "bogus_key"))

# 2. every constructor that forwards to Many_Body_Chain.__init__
S = dict(maxm=7, nsweeps=3, kpmmaxm=11)
ctors = [
    ("Spin_Chain", lambda: spinchain.Spin_Chain(["S=1/2"]*4, **S)),
    ("Fermionic_Chain", lambda: fermionchain.Fermionic_Chain(4, **S)),
    ("Spinful_Fermionic_Chain", lambda: fermionchain.Spinful_Fermionic_Chain(2, **S)),
    ("Majorana_Chain", lambda: fermionchain.Majorana_Chain(4, **S)),
    ("Bosonic_Chain", lambda: bosonchain.Bosonic_Chain(3, maxnb=[3]*3, **S)),
    ("SpinBoson_Chain", lambda: bosonchain.SpinBoson_Chain(["B", "S=1/2"], **S)),
    ("Parafermionic_Chain", lambda: parafermionchain.Parafermionic_Chain(4, **S)),
    ("Mixed_Spin_Fermion_Chain", lambda: mixedchain.Mixed_Spin_Fermion_Chain(
        ["S=1/2", "fermion"], **S)),
    ("Thermal_Spin_Chain", lambda: thermal.Thermal_Spin_Chain(["S=1/2"]*2,
                                                             T=0.5, **S).MBChain),
]
for label, f in ctors:
    c = build(label, f)
    if c is not None:
        print("%-26s maxm=%s nsweeps=%s kpmmaxm=%s   (asked 7, 3, 11)"
              % (label, c.maxm, c.nsweeps, c.kpmmaxm))

# 3. a setting that is not maxm, and a keyword that is not a setting at all
c = build("Fermionic_Chain", lambda: fermionchain.Fermionic_Chain(4, mode="ED"))
if c is not None: print("Fermionic_Chain(mode='ED').mode =", repr(c.mode))
c = build("Fermionic_Chain", lambda: fermionchain.Fermionic_Chain(4, spinful=False))
if c is not None: print("Fermionic_Chain(spinful=False): built")
c = build("Spin_Chain", lambda: spinchain.Spin_Chain(["S=1/2"]*4, maxm=50,
                                                    bogus_key=1, kpm_nscale=3))
if c is not None: print("Spin_Chain(maxm=50, bogus_key=1, kpm_nscale=3): built")
c = build("Bosonic_Chain", lambda: bosonchain.Bosonic_Chain(3, maxnb=[6]*3,
                                                        itensor_version=2, mode="ED"))
if c is not None: print("Bosonic_Chain(maxnb=[6]*3, itensor_version=2, mode='ED'): built, mode =", repr(c.mode))
c = build("Spin_Chain", lambda: spinchain.Spin_Chain(["S=1/2"]*4, maxm=0))
if c is not None: print("Spin_Chain(maxm=0): built, maxm =", c.maxm)

# 4. what it costs: a 12-site Heisenberg chain asked for maxm=4, nsweeps=10
def heisenberg(**kw):
    sc = spinchain.Spin_Chain(["S=1/2"]*12, itensor_version=3, **kw)
    h = 0
    for i in range(11):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc

a = heisenberg(maxm=4, nsweeps=10)          # asked at construction
b = heisenberg(); b.maxm, b.nsweeps = 4, 10  # set afterwards
d = heisenberg()                             # the defaults
ea, eb, ed = a.gs_energy(), b.gs_energy(), d.gs_energy()
print("12-site Heisenberg, v3: E(maxm=4 at construction) = %.10f  maxm=%d"
      % (ea, a.maxm))
print("12-site Heisenberg, v3: E(maxm=4 set afterwards)  = %.10f  maxm=%d"
      % (eb, b.maxm))
print("12-site Heisenberg, v3: E(defaults)               = %.10f  maxm=%d"
      % (ed, d.maxm))
print("exact (ED)                                        = %.10f"
      % d.gs_energy(mode="ED"))

# 5. mode="ED" at construction against mode="ED" assigned afterwards, and
# the way back to DMRG; Thermal_Spin_Chain writes its own mode onto MBChain
import io, contextlib
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()
def heis4(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h + 0.2*sc.Sz[0]
for tag in ("ctor", "assigned"):
    if tag == "ctor":
        sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python", mode="ED")
    else:
        sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
        sc.mode = "ED"
    m = sc.mode
    sc.set_hamiltonian(heis4(sc))
    e1 = sc.gs_energy()
    sc.mode = None
    try:
        e2 = "%.10f" % quiet(sc.gs_energy)
    except Exception as e:
        e2 = "%s: %s" % (type(e).__name__, e)
    print("mode='ED' %-8s  sc.mode = %-6r session built: %-5s"
          " gs_energy() = %.10f; then mode=None: gs_energy() = %s"
          % (tag, m, sc._session is not None, e1, e2))
for T in (0.0, 0.5):
    try:
        tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=T,
                                        itensor_version="python", mode="ED")
        h = 0
        for i in range(2):
            h = h + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
        tc.set_hamiltonian(h)
        wf = quiet(tc.get_gs)
        print("Thermal_Spin_Chain(T=%g, mode='ED').get_gs(): MBChain.mode = %r,"
              " <Sz0 Sz1> = %.6f" % (T, tc.MBChain.mode,
              np.real(tc.MBChain.vev(tc.Sz[0]*tc.Sz[1], wf=wf))))
    except Exception as e:
        print("Thermal_Spin_Chain(T=%g, mode='ED').get_gs(): %s: %s"
              % (T, type(e).__name__, e))

# 6. what the model class builds before the base constructor runs
c = build("Fermionic_Chain", lambda: fermionchain.Fermionic_Chain(
    4, itensor_version="python", N=5))
def show(x):
    return repr(x) if not isinstance(x, list) else \
        "a list of %d %s" % (len(x), type(x[0]).__name__)
if c is not None: print("Fermionic_Chain(N=5).N =", show(c.N))
try:
    c = parafermionchain.Parafermionic_Chain(3, itensor_version="python", Sig=1)
    print("Parafermionic_Chain(Sig=1).Sig =", show(c.Sig))
except Exception as e:
    print("%-26s %s: %s" % ("Parafermionic_Chain", type(e).__name__, e))

# 7. the infinite chain has no **kwargs, for comparison
build("Infinite_Spin_Chain", lambda: infinitechain.Infinite_Spin_Chain(
    ["S=1/2"]*2, maxm=50))
````

Observed, before:

````
[run3] slot 2 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
Spin_Chain(maxm=50, bogus_key=1): built, maxm = 30  hasattr(bogus_key) = False
Spin_Chain                 maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Fermionic_Chain            maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Spinful_Fermionic_Chain    maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Majorana_Chain             maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Bosonic_Chain              maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
SpinBoson_Chain            maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Parafermionic_Chain        maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Mixed_Spin_Fermion_Chain   maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Thermal_Spin_Chain         maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
Fermionic_Chain(mode='ED').mode = None
Fermionic_Chain(spinful=False): built
Spin_Chain(maxm=50, bogus_key=1, kpm_nscale=3): built
Bosonic_Chain              ValueError: itensor_version=2 only implements the fixed 4-level boson site (ITensor's BosonFourSite), but this chain asks for local boson dimension(s) [6]. Use itensor_version=3, itensor_version="python" or mode="ED", which all build a boson site of the requested dimension.
Spin_Chain(maxm=0): built, maxm = 30
12-site Heisenberg, v3: E(maxm=4 at construction) = -5.1420906326  maxm=30
12-site Heisenberg, v3: E(maxm=4 set afterwards)  = -5.1323602278  maxm=4
12-site Heisenberg, v3: E(defaults)               = -5.1420906326  maxm=30
exact (ED)                                        = -5.1420906328
mode='ED' ctor      sc.mode = None   session built: True  gs_energy() = -1.6282694864; then mode=None: gs_energy() = -1.6282694864
mode='ED' assigned  sc.mode = 'ED'   session built: True  gs_energy() = -1.6282694864; then mode=None: gs_energy() = -1.6282694864
Thermal_Spin_Chain(T=0, mode='ED').get_gs(): MBChain.mode = 'DMRG', <Sz0 Sz1> = -0.166667
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs(): MBChain.mode = 'DMRG', <Sz0 Sz1> = -0.125704
Fermionic_Chain(N=5).N = a list of 4 MultiOperator
Parafermionic_Chain(Sig=1).Sig = a list of 3 MultiOperator
Infinite_Spin_Chain        TypeError: Infinite_Many_Body_Chain.__init__() got an unexpected keyword argument 'maxm'
````

Observed, after:

````
[run3] slot 2 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
Spin_Chain                 TypeError: Spin_Chain() got unexpected keyword argument(s) bogus_key. A constructor keyword names a setting of the chain (maxm, nsweeps, noise, cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes effect as if assigned right after construction; a name the chain does not have would be stored where nothing reads it
Spin_Chain                 maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Fermionic_Chain            maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Spinful_Fermionic_Chain    maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Majorana_Chain             maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Bosonic_Chain              maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
SpinBoson_Chain            maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Parafermionic_Chain        maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Mixed_Spin_Fermion_Chain   maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Thermal_Spin_Chain         maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
Fermionic_Chain(mode='ED').mode = 'ED'
Fermionic_Chain            TypeError: Fermionic_Chain() got unexpected keyword argument(s) spinful. A constructor keyword names a setting of the chain (maxm, nsweeps, noise, cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes effect as if assigned right after construction; a name the chain does not have would be stored where nothing reads it
Spin_Chain                 TypeError: Spin_Chain() got unexpected keyword argument(s) bogus_key, kpm_nscale. A constructor keyword names a setting of the chain (maxm, nsweeps, noise, cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes effect as if assigned right after construction; a name the chain does not have would be stored where nothing reads it
Bosonic_Chain              ValueError: itensor_version=2 only implements the fixed 4-level boson site (ITensor's BosonFourSite), but this chain asks for local boson dimension(s) [6]. Use itensor_version=3 or itensor_version="python", which build a boson site of the requested dimension (and run by exact diagonalization under mode="ED" as well).
Spin_Chain                 ValueError: maxm must be >= 1, got 0
12-site Heisenberg, v3: E(maxm=4 at construction) = -5.1323602278  maxm=4
12-site Heisenberg, v3: E(maxm=4 set afterwards)  = -5.1323602278  maxm=4
12-site Heisenberg, v3: E(defaults)               = -5.1420906326  maxm=30
exact (ED)                                        = -5.1420906328
mode='ED' ctor      sc.mode = 'ED'   session built: True  gs_energy() = -1.6282694864; then mode=None: gs_energy() = -1.6282694864
mode='ED' assigned  sc.mode = 'ED'   session built: True  gs_energy() = -1.6282694864; then mode=None: gs_energy() = -1.6282694864
Thermal_Spin_Chain(T=0, mode='ED').get_gs(): MBChain.mode = 'DMRG', <Sz0 Sz1> = -0.166667
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs(): MBChain.mode = 'DMRG', <Sz0 Sz1> = -0.125704
Fermionic_Chain            TypeError: Fermionic_Chain() cannot set N at construction: that is the chain's state, not a setting (N: the Fermionic_Chain class builds it)
Parafermionic_Chain        TypeError: Parafermionic_Chain() cannot set Sig at construction: that is the chain's state, not a setting (Sig: the Parafermionic_Chain class builds it)
Infinite_Spin_Chain        TypeError: Infinite_Many_Body_Chain.__init__() got an unexpected keyword argument 'maxm'
````

**NUMBERS CHANGE**: YES. Spin_Chain(["S=1/2"]*12, itensor_version=3, maxm=4, nsweeps=10) with the open Heisenberg Hamiltonian: gs_energy() -5.1420906326 (ran at maxm=30) -> -5.1323602278 (runs at maxm=4, equal to assigning maxm/nsweeps afterwards). A constructor mode="ED" now takes effect: 4-site open Heisenberg + 0.2*Sz_0, max of the KPM ZZ correlator at delta=0.3 on es=linspace(0,3,7), 0.205144 ("python", DMRG) and 0.205151 (v3, DMRG) -> 0.205190 (ED). No in-tree caller moves (the four-correlation example passes mode="ED" per call: sum|ct_ed| 51.433292772677, E0 -2.450108689391 on both trees).

**Tests**: `tests/test_audit_2026_09_25_construction.py::test_a_constructor_setting_takes_effect`; `tests/test_audit_2026_09_25_construction.py::test_a_keyword_that_is_not_a_setting_raises_naming_every_one`; `tests/test_audit_2026_09_25_construction.py::test_private_names_methods_and_names_the_chain_lacks_are_not_settings`; `tests/test_audit_2026_09_25_construction.py::test_the_chain_state_is_not_a_setting`; `tests/test_audit_2026_09_25_construction.py::test_what_the_model_class_builds_is_not_a_setting`; `tests/test_audit_2026_09_25_construction.py::test_a_constructor_setting_goes_through_its_setter`; `tests/test_audit_2026_09_25_construction.py::test_a_constructor_setting_is_the_same_as_assigning_it`; `tests/test_audit_2026_09_25_construction.py::test_mode_at_construction_keeps_the_session`; `tests/test_audit_2026_09_25_construction.py::test_mode_at_construction_can_be_left_again`; `tests/test_audit_2026_09_25_construction.py::test_thermal_chain_with_mode_at_construction_still_solves`; `tests/test_audit_2026_09_25_construction.py::test_the_v2_boson_refusal_does_not_offer_mode_ed`; `tests/test_audit_2026_09_25_construction.py::test_an_unknown_mode_is_refused_at_construction`

**Reviewer (CONFIRMED, fix INCOMPLETE)**: Reproduced on the parent tree with the fix agent's own script (copied to the review folder): Spin_Chain(['S=1/2']*4, maxm=50, bogus_key=1) builds with maxm 30 and no error. All nine forwarding constructors (Spin, Fermionic, Spinful_Fermionic, Majorana, Bosonic, SpinBoson, Parafermionic, Mixed_Spin_Fermion, Thermal's MBChain) read maxm=30 nsweeps=15 kpmmaxm=50 when asked for 7/3/11. On a 12-site v3 Heisenberg chain, maxm=4/nsweeps=10 at construction returns -5.1420906326, the maxm=30 answer, where the same values assigned afterwards give -5.1323602278. I measured it on every backend too, with 0.2*Sz_0 added and seed 3: E(ctor)-E(assigned) = -9.635e-03 on v2 and v3 and -9.632e-03 on "python". Nothing documents this or excuses it, and fermionchain.py's Spinful docstring states the opposite. My own AST scan of src/, tests/, examples/ and benchmarks/ agrees with the fix agent's: the only in-tree non-own keyword is Fermionic_Chain(n0, mode="ED") in examples/staticcorrelators/four_correlation_tensor_sweep_VS_full. The multioperator_density call (spinful=False) no longer shows up because that example was edited. The rest is **kw forwarding: thermal.py and test helpers that pass itensor_version only. Nothing struck.

The settings half holds. On the working tree all nine constructors take maxm/nsweeps/kpmmaxm. Unknown, private, method and state names raise TypeError naming every key, sorted. maxm=0 raises ValueError. E(ctor)-E(assigned) is 3.553e-15 on v2 and v3 and 0.000e+00 on "python", and the 12-site v3 NUMBERS CHANGE line re-measures exactly (-5.1420906326 to -5.1323602278). Many_Body_Chain([1]*4, maxm=5) directly also works.

The mode= half is not "as if assigned right after construction", although the docstring, the check_settings error text and the expected-behaviour statement all promise that for every keyword. The narrowed claim to record: a constructor keyword naming a chain setting takes effect, on all nine forwarding constructors and every backend, as if assigned after construction, and anything else raises TypeError. mode= is the exception: it is applied before initialize(), so a chain built with mode="ED" has no DMRG session, where sc.mode="ED" after construction keeps one. Every ordinary reader I tried is identical on both paths and both backends (gs_energy, vev, get_excited, KPM correlator, deepcopy, bandwidth, site entropy, random_state). The difference shows on two routes. First, the mode=None round trip, which the fix agent recorded as a consequence: it now fails with an opaque AttributeError where assignment works. Second, an ordinary call on one of the nine constructors: Thermal_Spin_Chain(..., mode="ED").get_gs() ran on the parent (DMRG, mode silently ignored) and now raises AttributeError: 'NoneType' object has no attribute 'set_sweep_params'. The cause is that thermal.py:33 hardcodes the wrapper's own self.mode="DMRG" and get_gs() writes it onto MBChain at line 45. So the status line "Thermal_Spin_Chain needs no edit ... its settings now take effect on MBChain" is false for mode=. None of the new tests catch it, because SETTINGS has no mode and the Thermal row only checks MBChain attributes.

Suggested completion:
(a) Thermal_Spin_Chain.__init__ consumes mode= into its own self.mode, with default "DMRG", and does not forward it to MBChain. thermal.py is in the session cluster's diff, so coordinate with that owner.
(b) Either drop "mode" from the "as if assigned" wording and say that a chain built with mode="ED" is ED-only, or make the DMRG path raise a clear error when _session is None instead of the AttributeError.

Kept as leads, not verdict changes, as the fix agent listed them: Fermionic_Chain(4, N=5).N comes back as 5, overwriting the operator list, and Parafermionic_Chain(3, Sig=1) crashes inside its own __init__ with "TypeError: 'int' object is not subscriptable". Both names pass check_settings because the lists exist before the base constructor runs. cvm_maxm stays 30 under maxm=4, the same as assignment and documented.

On the tests: SETTINGS['verbose']=False equals the default and so cannot discriminate. That is harmless, since the other seven keys do. The NUMBERS CHANGE statement is complete as a general rule ("any caller that passed a setting"). It should add that a mode="ED" keyword now takes effect, which moves calls that do not pass mode= per call: my 4-site KPM ZZ peak went from 0.205144 (python) and 0.205151 (v3) on DMRG to 0.205190 on ED. It should also say that such a chain builds no session.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
Review folder <scratch>/review/construction/.

01_init_kwargs.before.out:
```
Spin_Chain(maxm=50, bogus_key=1): built, maxm = 30  hasattr(bogus_key) = False
Thermal_Spin_Chain         maxm=30 nsweeps=15 kpmmaxm=50   (asked 7, 3, 11)
12-site Heisenberg, v3: E(maxm=4 at construction) = -5.1420906326  maxm=30
12-site Heisenberg, v3: E(maxm=4 set afterwards)  = -5.1323602278  maxm=4
```
01_init_kwargs.after.out:
```
Spin_Chain                 TypeError: Spin_Chain() got unexpected keyword argument(s) bogus_key. ...
Thermal_Spin_Chain         maxm=7 nsweeps=3 kpmmaxm=11   (asked 7, 3, 11)
12-site Heisenberg, v3: E(maxm=4 at construction) = -5.1323602278  maxm=4
```
04_attack_init.before.out:
```
[2] 12-site E(ctor) - E(assigned)                          -> -9.635e-03   E(ctor)=-5.1573204922
[3] 12-site E(ctor) - E(assigned)                          -> -9.635e-03   E(ctor)=-5.1573204922
[python] 12-site E(ctor) - E(assigned)                     -> -9.632e-03   E(ctor)=-5.1573204922
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs()              -> MBChain.mode='DMRG', <Sz0 Sz1>=-0.125704
[python ctor] get_dynamical_correlator KPM                 -> 0.205144
```
04_attack_init.after.out:
```
[2] 12-site E(ctor) - E(assigned)                          -> 3.553e-15   E(ctor)=-5.1476856698
[3] 12-site E(ctor) - E(assigned)                          -> 3.553e-15   E(ctor)=-5.1476856698
[python] 12-site E(ctor) - E(assigned)                     -> 0.000e+00   E(ctor)=-5.1476882590
[python ctor] mode='ED' session=NoneType
[python ctor] mode=None; gs_energy()                       -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
[python assigned] mode='ED' session=Chain
[python assigned] mode=None; gs_energy()                   -> -1.6282694864
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs()              -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
Thermal_Spin_Chain(T=0, mode='ED').get_gs()                -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
Fermionic_Chain(4, N=5).N                                  -> 5
Parafermionic_Chain(3, Sig=1).Sig                          -> TypeError: 'int' object is not subscriptable
Spin_Chain(maxm=4): maxm=4 cvm_maxm=30
```
06_scan.out: only examples/staticcorrelators/four_correlation_tensor_sweep_VS_full/main.py:86 Fermionic_Chain ['mode'], outside ** forwarding. 07_tests.after.out: 43 passed. 07_tests.before.out (file copied out, DMRGPY_SRC=parent): 41 failed, 2 passed, the passes being test_get_gs_without_wf0_still_returns_the_stored_state[python] and [v3].
````

**Left open**: Thermal_Spin_Chain(..., mode="ED") still runs DMRG: thermal.py:33 sets the wrapper's own self.mode="DMRG" and get_gs() writes it onto MBChain (thermal.py:45), overwriting the forwarded keyword, as on the parent; thermal.py is the session cluster's, and the exact edit is in the notes. An unknown keyword to Thermal_Spin_Chain raises with the message naming Spin_Chain() rather than Thermal_Spin_Chain(). Infinite_Many_Body_Chain accepts no settings at construction (it raises TypeError on maxm=, which is honest, not silent); left as it is. cvm_maxm is set from the default maxm when the chain is built, so it stays 30 under maxm=4 unless it is given too, exactly as when maxm is assigned afterwards. (The first sentence is superseded by the Status addendum above.)

### 2. On v2 and v3 a Hamiltonian written in small energy units got a wrong ground state and wrong spectra from three absolute thresholds inside ITensor: toMPO's svdMPO dropped the exchange channels of every Heisenberg bond once their squared singular values summed below 1e-13 (vev read H back as 1/3 of itself and gs_energy()/s came out at the Neel -1.25 against -2.493577, from s=4e-7 on v3 and s=1e-7 on v2), davidson() replaced every Krylov direction whose residual was below 1e-10 by a random vector (2e-10 to 2.1e-9 off at s=1e-6 over the runs measured, and 1.8e-5 to 4.4e-5 at s=1e-8 on a chain whose MPO has nothing to truncate), and svdMPO's isZero skipped the KPM band-centre shift below 1e-14 (1.2 of the ED peak off from s=1.5e-14 once the first two were fixed)

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; cluster `scale-cpp`

**Status**: FIXED, on the claim as the review narrowed it: for a Hamiltonian written in small units, a term-built s*H, on v2 and v3 the ground state, the excited states, both band edges, `gs_energy_generalized` and every KPM spectrum now come out as the unit-scale calculation, measured to s=1e-12 here (the review measured gs_energy and get_excited to 1e-40) and KPM to s=1e-20. What changed: `mo_terms.h` on both backends gains `unit_scale_up()`, the power of two that brings a largest coefficient below 1 into [1,2) and 1.0 otherwise, and `to_mpo_unit()`, which hands ITensor's `toMPO` the AutoMPO multiplied by that factor and divides the MPO back, exactly; `build_mpo()` goes through it, and so do the Hamiltonians the four v3 and two v2 real-time routines build from an AutoMPO with an appended `-EGS*Id`, the KPM band-centre shift (`scaled_hamiltonian` on both, `scaled_hamiltonian_gs_anchored` on v3) and CVM's `z*Id` (both, reached by no public caller). Each `chain_session.h` keeps `hscale_up_`, read from the caller's terms in `set_hamiltonian(terms)` and reset to 1.0 by v3's `set_hamiltonian_mpo`, and every `dmrg()` on the Hamiltonian solves on `hscale_up_*H` through `solver_hamiltonian()`: `gs_energy` divides the energy back, `excited_states` multiplies the overlap weight by the same factor, v3's `gs_energy_generalized` solves on the scaled `Heff`, and `maximum_energy` solves on `-hscale_up_*H`. At a largest coefficient of 1 or more the term MPO and the solver run the unscaled code, byte for byte; the KPM shift and CVM's z, whose one coefficient is usually below 1, take the scaled build at ordinary scales too, which is the same operator exactly (the 6-site KPM row at J=1 is identical to ten printed digits on the parent, the first pass and this pass). Both extensions were rebuilt, and `make -n pybind` reports nothing to do in either. Not fixed, and a different defect: one scale per operator does not reach a bond whose strongest crossing term is far below the largest coefficient of the operator, since svdMPO truncates bond by bond against the absolute 1e-13, so an O(1) energy offset or one-site field next to exchange at 1e-7, and a weak link J'=1e-7 in a J=1 chain, read the exchange at 1/3 of itself on v3 and v2 before the fix and after it, at any units, while "python" drops them outright (0.0000 and 0.0892) through the relative return-sweep cutoff of `pyitensor/mpobuilder.py`; the review's gate on the largest multi-site coefficient was measured against and not adopted, since it leaves the weak link, and the per-bond design is in cpp_plan. Pinned by `tests/test_audit_2026_09_25_scale.py`: `test_ground_state_energy_is_scale_invariant_on_every_mps_backend` (v3 and v2 at s from 1 to 1e-8 held at 1e-11, "python" as the control at 1e-6), `test_small_units_ground_state_on_the_compiled_backends` (s=1e-10 and 1e-12), `test_the_mpo_of_a_small_unit_operator_keeps_every_channel` (the MPO alone, vev on a fixed state), `test_the_local_solver_is_scale_invariant` (the solver alone, on the one-channel Ising chain), `test_excited_states_in_small_units_match_ed`, `test_kpm_correlator_in_small_units_matches_ed`, `test_kpm_shift_survives_below_the_absolute_coefficient_floor` (s=1.5e-14 and 1e-20) and `test_generalized_ground_state_in_small_units`. Of these 39 tests, 29 fail on a tree with the fixed Python and the parent's extensions, and the 4 KPM-floor ones fail on the first-pass build. The open part is held by strict xfails that hold on all three trees, `test_a_bond_far_below_the_largest_coefficient_still_open` (v3, v2 and "python", each for the offset, the field and the weak link, on a fixed state) and `test_small_units_next_to_an_order_one_term_still_open` (gs_energy of s*H + 1 and s*H + Sz0 at s=1e-7 on v3 and v2), and python-1e-10 stays for the Lanczos of `pyitensor/dmrg.py`. NUMBERS CHANGE, measured with the fixed Python and the old extensions as the before, on the 6-site S=1/2 Heisenberg chain written as s*H (`maxm=30`, `nsweeps=10`, fresh chain per s): gs_energy()/s from the Neel -1.25 (v3 at s from 3e-7 to 1e-9, v2 from 1e-7 to 1e-9), -1.27 and -1.40 (v3 at 6e-7 and 5e-7 in one run, -1.12 at 6e-7 in the review's), -1.747 (v2 at 3e-7) and -1.22 to -0.71 (both at 1e-10 to 1e-12) to -2.4935771339 at every s, errors at or below 5.8e-15 over the three runs measured. At s=1e-6 the error goes from 2e-10 to 2.1e-9 over the runs measured to at most 3.6e-15. vev(s*H)/s on a fixed state goes from 0.6667 and 0.3333 of vev(H) to 1.0000. At s=1e-8 and 1e-10, `get_excited(n=4)/s` goes from -1.24999797 -1.24999436 -0.74999323 -0.74998836 (v3 at 1e-8) to -2.49357713 -2.00199536 x3, the KPM correlator from 2.8 (v3) and 1.1 (v2) of the ED peak off to 4.1e-4, and v3's `gs_energy_generalized/s` with A = 1 + 0.2*Sz0 from -1.3888862460 to -2.5748369148. The KPM correlator at s from 1.5e-14 to 1e-20 goes from 1.2 of the ED peak off on the first-pass build (a RuntimeError, nan or 0.78 to 2.0 on the older ones) to 4.1e-4 on both backends. The transverse-field Ising chain s*(sum ZZ + 0.7 sum X) goes from 1.8e-5 to 3.6e-5 off at s=1e-8 and 0.46 to 0.64 off at 1e-10 to 3.9e-13 on both. As a side effect of `build_mpo`, v3's NH-DMRG on s*(Heisenberg + 0.3j*Sz0) at s=1e-8 goes from -1.2 to -2.46104102, 3.1e-13 to 7.1e-12 off, and v3 TDVP real-time evolution at s=1e-8 from 4.0e-2 to 6.4e-4 off. More widely, every v2 or v3 calculation on a term-built Hamiltonian whose largest coefficient is below 1 now builds its MPO and runs its local solves on 2^k H, and every C++ KPM run whose band-centre shift is below 1 in magnitude builds that shift at 2^k. Where the MPS is not exact, the change is below the run-to-run noise: on 20 sites at `maxm=10` the J=0.5 Heisenberg energy is -4.341137444378 before and after on both backends, and the v3 J=1 XX mean over ten runs moves from -6.190466101716 to -6.190466092079 against a std of 7e-8. No number changes for the bond-local cases: the offset, field and weak-link rows of 20 are the same on the hybrid, the first-pass and the repaired build (0.3333 on v3, 0.6667 then 0.3333 on v2, 0.0000 and 0.0892 on "python").

**Status note on the gate** (the main agent): the reviewer of the locate stage recommended gating the scale at a power of two near 1e-2 (2^-7) rather than at a largest coefficient of 1, so that every J=0.5 model and the realified XX chain at J=1, whose largest AutoMPO coefficient is 0.5, stay on the unscaled path bit for bit; that review was not forwarded to the C++ stage, and the gate stayed at 1. It is kept there: 1 is the one threshold at which byte identity is exact by construction for every Hamiltonian at or above it, and below it the measured cost is what the Status gives, the J=0.5 Heisenberg energy on 20 sites at maxm=10 unchanged to every printed digit on both backends and the v3 J=1 XX mean over ten runs moving by 1e-8 against a run-to-run standard deviation of 7e-8, both within the noise of the random start. Results at those couplings are therefore not bit-identical across this fix, only statistically identical.

**Recorded in**: 2026-09-24c New leads. Vendored ITensor behaviour in both versions; dmrgpy has built every MPO through toMPO and solved through ITensor's dmrg() since the pybind11 sessions were written. It was masked below s=1e-8 until the clean-threshold fix of the same cluster, which is why the parent-tree rows at s <= 1e-8 read 0, and the KPM shift floor was masked behind the other two mechanisms until the first fix pass, whose reviewer measured it.

**Where**: src/dmrgpy/mpscpp2/mo_terms.h and src/dmrgpy/mpscpp3/mo_terms.h `build_mpo` (every term-built MPO: set_hamiltonian, vev, correlators, KPM vertices, build_operator) through `to_mpo_unit`; the evolution Hamiltonians built from an AutoMPO with an appended -EGS*Id (v3 quench, evolve_and_measure, quench_tdvp, quench_tdvp_gse; v2 quench, evolve_and_measure); the single-term shift AutoMPOs of the KPM rescaling (v3 `scaled_hamiltonian` and `scaled_hamiltonian_gs_anchored`, v2 `scaled_hamiltonian`) and CVM's z*Id (`cvm_dynamical_correlator`, both); every dmrg() on the Hamiltonian in both chain_session.h (gs_energy, excited_states, maximum_energy, and v3's gs_energy_generalized) through `solver_hamiltonian`. The thresholds themselves are vendored: mpscppN/ITensor/itensor/mps/autompo.cc's compressMPO calls `truncate(D,maxdim,mindim,cutoff)` with `cutoff = args.getReal("Cutoff",1E-13)` (v3 :1176, v2 :1157) and doRelCutoff=false (decomp.cc divides the cutoff by P(0) along with P), and skips a coefficient through the hardcoded `isZero(t.coef,eps)` with `eps = 1E-14` (v3 :1172 and :1287); iterativesolvers.h davidson(), the hardcoded `if(qnrm < 1E-10)` (v3:324, v2:295).

**Severity**: high: silently wrong ground states on the default backend (the Neel energy, or a non-Hermitian half-hopping MPO) for any term-built Hamiltonian whose largest coefficient was below about 6e-7 on v3 and 3e-7 on v2, a continuous accuracy loss from about 1e-6 down on both, carried into excited states, band edges, KPM spectra and the generalized solve, and a KPM spectrum off by its whole position below about 1.5e-14

Three numbers inside ITensor are absolute, calibrated for an operator of order one. svdMPO's truncate() is called with the Cutoff it reads from Args (`args.getReal("Cutoff",1E-13)`) and with doRelCutoff=false, and it divides that cutoff by the bond's largest squared singular value along with the weights, so the discarded weight is compared in absolute terms; the other two are hardcoded, svdMPO's `isZero(coef,1E-14)` on every coefficient it places and davidson()'s `qnrm < 1E-10` (the first pass's claim that neither threshold is reachable through Args was struck by its reviewer, and only these two are). On v3, which realifies the Heisenberg bond into ZZ at weight s^2 and S+S-, S-S+ at (s/2)^2, one hopping channel goes at s <= 6.32e-7 and both at s <= 4.47e-7, and on v2, whose bond carries XX, YY, ZZ at s^2 each, one goes below 3.16e-7 and two below 2.24e-7, the third kept by the minimum bond dimension, which is exactly what 03 part A reads with no eigensolver running (H at 0.6667 then 0.3333 of itself, XX+YY at 0.5). davidson() is the second, gradual mechanism: it is what the one-channel transverse-field Ising chain meets (03 part C, 1.8e-5 to 4.4e-5 off at s=1e-8, 0.24 to 0.64 at 1e-10), and the locate stage's standalone build of ITensor's own dmrg() had already shown that making that one test relative is exact to 2.7e-14 at every s from 1 to 1e-12, while making Approx0 and ErrGoal relative instead changes nothing. The third, isZero, bites a single-term AutoMPO first, and the one that matters is the shift*Id with which KPM moves the band centre of H to zero: once 0.62*s fell below 1e-14 the shift was skipped and the spectrum lost its position, 1.2 of the ED peak off from s=1.5e-14 against 4.1e-4 at 1.8e-14 on both backends (20 part B on the first-pass build). The fix does not touch the vendored ITensor. mo_terms.h gains `unit_scale_up(cmax)`, the power of two that brings a largest |coefficient| below 1 into [1,2) and 1.0 otherwise, and `to_mpo_unit(ampo,args)`, which hands toMPO a copy of the AutoMPO multiplied by it and divides the MPO back; `build_mpo`, the six evolution AutoMPOs, the three KPM shift AutoMPOs and the two CVM z*Id ones go through it. The scale is read from the AutoMPO ITensor actually sees, after v3's realification, so the v3 XX chain at J=1 (every coefficient 0.5 once written with S+/S-) takes the scaled MPO build too. AutoMPO keeps its terms ordered by operator string alone (LessNoCoef), so the scaled copy holds the same terms in the same order, and a real power of two keeps an exactly real coefficient exactly real, which svdMPO's is_real branch relies on. Each chain_session.h keeps `hscale_up_`, read from the caller's term list in set_hamiltonian(terms) and reset to 1.0 by v3's set_hamiltonian_mpo and by forget_everything_built_on_sites, and `solver_hamiltonian(H)` returns H itself at 1.0 and hscale_up_*H otherwise: gs_energy solves on it and divides the energy back before caching it, excited_states multiplies the overlap weight by hscale_up_ (the penalty is an energy), gs_energy_generalized solves on the scaled Heff and reads lambda back from innerC as before, and maximum_energy solves on -hscale_up_*H and divides back. The MPS dmrg() returns is normalized, so nothing downstream sees the factor, and the noise term, which denmatDecomp adds as noise*deltaRho, quadratic in H, before renormalizing by the trace, gets back its unit-scale strength with it. The scale-back is exact in both builds: ITensor's MPO operator*=(Real) multiplies the stored elements of the tensor at leftLim()+1 (v3 mpo.cc:104, v2 mpo.h:153, with itensor.cc:1015 and itensor_operators.cc:188 multiplying elements, USESCALE being defined in neither build), so multiplying by 2^-k is an exponent shift, and the SVD is homogeneous under a power of two away from LAPACK's underflow rescaling, so what changes is only which side of the absolute thresholds a number falls on. At a largest coefficient of 1 or more unit_scale_up returns 1.0, to_mpo_unit returns toMPO(ampo,args) itself and solver_hamiltonian a copy of H, so no floating-point operation changes in the term MPO or in the solver. The KPM shift and CVM's z are the exception: their single coefficient is usually below 1 (0.62 on the 6-site chain at J=1), so they take the scaled build at ordinary scales too, which gives the same operator exactly, with tensors that differ only in the unused not-yet-started channel of site 1, which carries 2^-k in place of 1; the 6-site KPM row of 10 is identical to all ten printed digits on the parent, the first-pass and the repaired build. Neither backend starts deterministically (v3 randomMPS, v2 MPS(sites) is a random product state, and the two parent runs of 10 differ in the last digits), so the ordinary-scale guarantee is shown by the gate plus agreement with ED at s=1 (10) and the golden files, not by bit comparison of two runs. Below 1 the change was measured where the MPS is not exact (14, 15): on 20 sites at maxm=10 the J=0.5 Heisenberg energy is the same to all 12 printed digits before and after on both backends, and the v3 XX chain's ten-run mean moved by 9.6e-9 against a run-to-run std of 7e-8, both 2.79e-4 above the exact free-fermion energy from truncation. What one scale per operator cannot reach, and the review found, is a hierarchy inside the Hamiltonian. truncate() runs bond by bond, each bond against its own largest weight, so a bond whose strongest crossing term is far below the largest coefficient of the operator still compares its channels against the absolute 1e-13 in units of that coefficient, whatever the units of H: an O(1) energy offset (sent as ('Id',site)) or one-site field next to exchange at 1e-7 reads the exchange at 1/3 of itself on v3 and v2 (20 part A), and so does a weak link J'=1e-7 in an ordinary J=1 chain (20 part C), on the hybrid, the first-pass and the repaired build alike. "python" fails the same probes harder (0.0000, 0.0892 and 0.0000), through the relative cutoff of pyitensor/mpobuilder.py's return sweep, which measures a weak bond against the Schmidt weight of the whole operator. That is a different defect from this item, since it happens at s=1 as well, and it is recorded open with strict xfails on all three backends, with the per-bond cure in cpp_plan. The review's suggestion of gating the scale on the largest multi-site coefficient was not adopted: it would fix the offset and field MPOs but not the weak link, whose largest coefficient already sits on a multi-site term, and moving the solver's scale onto it would give the field of s*H + Sz0 a coefficient of 2^24 in the solver, and the noise term, quadratic in it, a factor 2^48 (by reading, not measured). Because the parent tree reads 0 below s=1e-8 for the clean-threshold reason, the C++ half is isolated with a hybrid tree, the working tree's Python with the parent's compiled extensions (cmp-identical to the parent's), and the repair's own before is a first-pass tree, the working tree's Python with the first pass's extensions (cmp-identical to the ones this pass replaced).

**Expected**: E0(s*H) = s*E0(H), vev(s*X) = s*vev(X), the excited levels, the band edges and every KPM spectrum scale with s on v2 and v3 down to where double precision runs out, as they already do on ED and julia_live, while every calculation whose largest coefficient is 1 or more builds its term MPO and runs its local solver on exactly the code it ran before.

#### Locate stage

(a) is the cliff the lead measured, and it is located quantitatively. With no eigensolver at all, vev(s*X)/s on the scale-1 ground state (03, part A) is exact on "python" at every s, while on v3 the whole Hamiltonian reads 0.6667 of itself at s=6e-7 and 5e-7 and 0.3333 from 4e-7 down, the XX+YY part alone 0.5 from 6e-7 down and ZZ 1.0000 throughout, and on v2 H reads 0.6667 at 3e-7 and 0.3333 at 1e-7, XX+YY 0.5 from 3e-7 down, and a lone XX or YY 1.0000. That is exactly what svdMPO's truncation predicts: it normalizes the squared singular values by the largest AND divides the 1e-13 cutoff by it, so the discarded weight is compared in absolute terms, dropping from the smallest while the summed discarded weight stays at or below 1e-13, with MinDim=1 always keeping one channel. On v3 the Hamiltonian is realified, so a bond carries ZZ at weight s^2 and S+S-, S-S+ at (s/2)^2 each: one hopping channel goes at s <= sqrt(1e-13/0.25) = 6.32e-7, both at s <= sqrt(1e-13/0.5) = 4.47e-7. On v2 the bond carries XX, YY, ZZ at s^2 each: one goes at s < 3.16e-7, two at s < 2.24e-7, and the third is kept by MinDim. A single-channel operator (XX alone, a one-site term) takes truncate()'s origm==1 early return and is never cut, which is why vev(1e-7*Sz0) was always right on v2 and v3. With one of S+S- and S-S+ gone the MPO is not Hermitian, which is the v3 -1.28 to -1.33 at 5e-7 and 6e-7 whose state measures -1.30 to -1.73 with the true H (02). (b) is a second, gradual mechanism that the cliff hid: with the MPO built at scale 1 and multiplied by s exactly (set_hamiltonian(s*toMPO(H)), 03 part B) v3 is 1.8e-10 off at s=1e-6, 1.4e-8 at 3e-7, 5.5e-6 at 1e-8, 6.0e-4 at 1e-9 and 0.13 at 1e-10, and the transverse-field Ising chain, whose MPO has one channel per bond and so is exact at any s, degrades the same way on both v3 and v2 (part C, 2.1e-7 and 2.5e-7 at 1e-7, 4.4e-5 and 2.8e-5 at 1e-8). A standalone program linking the vendored libitensor.a (04) runs ITensor's own dmrg() on 6-site Heisenberg times s with a copy of iterativesolvers.h in which one threshold at a time is made relative to s: making the randomization test `qnrm < 1E-10` relative alone gives 2.4e-14 to 2.7e-14 at every s from 1 to 1e-12, while making the convergence tests (Approx0 = 1e-12 and ErrGoal) relative alone changes nothing (4.5e-6 at 1e-8, 0.35 at 1e-10), so it is the orthogonalization-failure branch, which replaces the Krylov direction by a random vector whenever the residual is below an absolute 1e-10, i.e. below a relative 1e-10/s; Approx0 bites only at s=1e-12. v2 has the same line at iterativesolvers.h:295 (by reading; its degradation is measured on the Ising chain). (b') "python"'s MPO is scale-free (part A), and running the same public gs_energy() with the local Lanczos handed the unit-normalized operator (monkeypatched in that process only, part D) takes it from 2.5e-9, 1.2e-6, 2.8e-5, 1.0e-2 and 2.6 off at s = 1e-8, 1e-9, 1e-10, 1e-11 and 1e-12 to 2.4e-14, 3.9e-14, 2.9e-12, 2.0e-12 and 4.4e-8, so it is `_lanczos_ground_state`'s absolute tests, both `beta < tol` and the value test `tol*max(1, |E|)`, which the one patch covers together; the remaining 4.4e-8 at 1e-12 was not investigated. julia_live (05, one run) returns -2.4935771339 at s=3e-7 and -2.4935771338 at 1e-9, and its vev(s*H)/s reads 1.0000 down to 1e-9, so ITensors.jl's OpSum-to-MPO is relative. Not measured: the Taylor MPO products of `evoloperator` and the MPO sums at Cutoff=cutoff_ (ITensor's svd() defaults to DoRelCutoff=true there, by reading, so they are expected to be scale-free).

Locate-stage status: LOCATED, not fixed. Three mechanisms, each separated by measurement. (a) v2 and v3 build every MPO through ITensor's toMPO, whose svdMPO truncates squared singular values on an absolute 1e-13 (and drops matrix elements below an absolute 1e-14), so the XX+YY channels of a Heisenberg bond go at s <= 6.32e-7 and 4.47e-7 on v3 (realified, one channel then both) and s < 3.16e-7 and 2.24e-7 on v2, predicted and measured to the digit by vev(s*X)/s on a fixed state with no eigensolver running; this is the lead's -1.25 and -1.747. (b) ITensor's davidson() replaces every Krylov direction whose residual is below an absolute 1e-10 by a random vector (v3 iterativesolvers.h:324, v2 :295), which costs accuracy continuously from about 1e-6 down even on an exactly scaled MPO (5.5e-6 off at 1e-8, 0.13 at 1e-10 on v3; the one-channel Ising chain degrades the same way on v2 and v3); a standalone build with that single test made relative is exact to 2.7e-14 at every s from 1 to 1e-12, and making Approx0 and ErrGoal relative instead changes nothing. (b') "python"'s MPO is scale-free, and its `_lanczos_ground_state` stops on an absolute breakdown and value test below |E|=1 (1.2e-6 off at 1e-9, 2.8e-5 at 1e-10), removed by handing it the unit-normalized operator (3.9e-14 at 1e-9). julia_live is unaffected (-2.4935771339 at 3e-7, -2.4935771338 at 1e-9). No code changed for this item; the plan is in cpp_plan, and `tests/test_audit_2026_09_25_scale.py::test_small_units_ground_state_still_open` holds strict xfails for v3 and v2 at s=1e-8, which the C++ stage must turn into plain tests, and for "python" at s=1e-10, which stays until `pyitensor/dmrg.py` is fixed. No number changes in this stage.

Repro (<scratch>/scale/02_small_units_onset.py, 03_small_units_mpo.py, cpp/04_davidson_exits.cc (+ cpp/patched/itensor/iterativesolvers.h), 05_julia_live.py):

````python
# ==== 02_small_units_onset.py ====
# small-units, where each backend starts failing: gs_energy()/s of s*H for
# the 6-site S=1/2 Heisenberg chain, a fresh chain per s (a v3 session warm
# starts from its previous state), against ED at s=1.  Next to it, the energy
# of the returned state measured with the UNSCALED H (vev(H) builds its own
# MPO at scale 1, where nothing is truncated), which separates "the state is
# wrong" from "the reported number is wrong".  Pinned sweep schedule.
import sys
import numpy as np
from dmrgpy import spinchain, cppext

L = 6
def heis(sc):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

ss = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 6e-7, 5e-7, 3e-7, 1e-7, 1e-8,
      1e-9, 1e-10, 1e-11, 1e-12]
backends = ["ED", "python"] + [v for v in (3, 2) if cppext.available(v)]
if len(sys.argv) > 1: backends = [a if a in ("ED", "python") else int(a) for a in sys.argv[1:]]
ref = None
for v in backends:
    for s in ss:
        np.random.seed(1)
        sc = spinchain.Spin_Chain(["S=1/2"]*L,
                                  itensor_version=("python" if v == "ED" else v))
        h = heis(sc)
        sc.set_hamiltonian(s*h)
        sc.maxm = 30; sc.nsweeps = 10
        mode = "ED" if v == "ED" else "DMRG"
        try:
            e = np.real(sc.gs_energy(mode=mode))/s
            eh = np.real(sc.vev(h, mode=mode))
            if ref is None: ref = e
            print("v=%-6s s=%.0e  gs_energy()/s=%.10f  <wf|H|wf>=%.10f  error=%.1e" % (
                  v, s, e, eh, abs(e - ref)), flush=True)
        except Exception as ex:
            print("v=%-6s s=%.0e  raised %s: %s" % (v, s, type(ex).__name__, str(ex)[:70]), flush=True)

# ==== 03_small_units_mpo.py ====
# small-units, the discriminant between (a) the MPO being built wrong and (b)
# the local eigensolver stopping on an absolute criterion.
#  A. No eigensolver at all: solve H at scale 1 once, then measure
#     vev(s*X)/s / vev(X) on that one state for X = H, the XX+YY part, the ZZ
#     part and one bond's XX+YY.  Every vev builds its MPO from the scaled
#     terms through the same builder set_hamiltonian uses.
#  B. v3 only: the Hamiltonian as an MPO built at scale 1 and multiplied by s
#     afterwards (set_hamiltonian(s*toMPO(H)), exact), so the solver sees a
#     correct small-unit operator; gs_energy()/s against exact.
#  C. v2 and v3: the transverse-field Ising chain, whose bond carries a single
#     operator string (ZZ), so its MPO has one channel per bond and nothing
#     for the compression to truncate; gs_energy()/s against ED.
#  D. "python": the same s*H with the local Lanczos solve rescaled to unit
#     norm (monkeypatched in this process only), against the stock solver.
import numpy as np
from dmrgpy import spinchain, cppext

L = 6
def parts(sc):
    xy = 0; zz = 0
    for i in range(L-1):
        xy = xy + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
        zz = zz + sc.Sz[i]*sc.Sz[i+1]
    return xy, zz

def chain(v, H=None):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 10
    return sc

cpp = [v for v in (3, 2) if cppext.available(v)]
ss = [1e-3, 1e-5, 1e-6, 6e-7, 5e-7, 4e-7, 3e-7, 1e-7, 1e-8, 1e-10, 1e-12]

print("A. vev(s*X)/s / vev(X) on the scale-1 ground state")
for v in cpp + ["python"]:
    sc = chain(v)
    xy, zz = parts(sc)
    b01 = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1]
    sc.set_hamiltonian(xy + zz)
    ops = (("H", xy + zz), ("XX+YY", xy), ("ZZ", zz), ("bond01 XX+YY", b01))
    if v == 2:
        ops = ops + (("XX", sum(sc.Sx[i]*sc.Sx[i+1] for i in range(L-1))),
                     ("YY", sum(sc.Sy[i]*sc.Sy[i+1] for i in range(L-1))))
    refs = dict((n, np.real(sc.vev(X))) for n, X in ops)
    print("  v=%s  reference vev: %s" % (v, "  ".join("%s=%.6f" % (n, refs[n]) for n, X in ops)))
    for n, X in ops:
        print("    %-13s %s" % (n, " ".join("%.0e:%.4f" % (s, np.real(sc.vev(s*X))/s/refs[n]) for s in ss)), flush=True)

if 3 in cpp:
    print("B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly")
    for s in [1.0, 1e-6, 6e-7, 5e-7, 3e-7, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]:
        sc = chain(3)
        xy, zz = parts(sc)
        H = xy + zz
        sc.set_hamiltonian(s*sc.toMPO(H))
        e = np.real(sc.gs_energy())/s
        print("  s=%.0e  gs_energy()/s=%.10f  <wf|H|wf>=%.10f  error=%.1e" % (
              s, e, np.real(sc.vev(H)), abs(e + 2.4935771338879)), flush=True)

print("C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond")
sc = chain("python")
tf = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1)) + 0.7*sum(sc.Sx)
sc.set_hamiltonian(tf)
eref = np.real(sc.gs_energy(mode="ED"))
print("  ED at s=1: %.10f" % eref)
for v in cpp:
    row = []
    for s in [1.0, 1e-6, 3e-7, 1e-7, 1e-8, 1e-10, 1e-11, 1e-12]:
        sc = chain(v)
        sc.set_hamiltonian(s*(sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1)) + 0.7*sum(sc.Sx)))
        e = np.real(sc.gs_energy())/s
        row.append("%.0e:%.1e" % (s, abs(e - eref)))
    print("  v=%s  |gs_energy()/s - ED|  %s" % (v, " ".join(row)), flush=True)

print("D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator")
from dmrgpy.pyitensor import dmrg as pydmrg
_stock = pydmrg._lanczos_ground_state
def _unit(matvec, v0, niter=30, tol=1e-12, residual_tol=None):
    c = np.linalg.norm(np.asarray(matvec(v0)))/np.linalg.norm(np.asarray(v0))
    if c == 0.: return _stock(matvec, v0, niter=niter, tol=tol, residual_tol=residual_tol)
    val, vec = _stock(lambda x: matvec(x)/c, v0, niter=niter, tol=tol,
                      residual_tol=residual_tol)
    return val*c, vec
for label, fn in (("stock", _stock), ("unit-normalized", _unit)):
    pydmrg._lanczos_ground_state = fn
    row = []
    for s in [1.0, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]:
        sc = chain("python")
        xy, zz = parts(sc)
        sc.set_hamiltonian(s*(xy + zz))
        e = np.real(sc.gs_energy())/s
        row.append("%.0e:%.1e" % (s, abs(e + 2.4935771338879)))
    print("  %-16s |gs_energy()/s - exact|  %s" % (label, " ".join(row)), flush=True)
pydmrg._lanczos_ground_state = _stock

# ==== cpp/04_davidson_exits.cc ====
// small-units (b): which of ITensor v3 davidson()'s absolute thresholds
// fires for a Hamiltonian in small units.  6-site S=1/2 Heisenberg built by
// AutoMPO at scale 1 (so the MPO is exact) and multiplied by s afterwards,
// then ITensor's own dmrg() at DebugLevel 3, whose davidson() prints why it
// leaves each local solve and when it replaces a Krylov direction by a
// random vector.  Standalone: links the vendored libitensor.a read-only.
#include "itensor/all.h"
#include <cstdlib>
using namespace itensor;
int main(int argc, char* argv[])
    {
    double s = std::atof(argv[1]);
    int N = 6;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo += "Sz",j,"Sz",j+1;
        }
    auto H = toMPO(ampo);
    H.ref(1) *= s;
    g_davidson_scale = s;
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 30;
    sweeps.cutoff() = 1E-12;
    sweeps.noise() = 0.;
    auto psi0 = randomMPS(sites,4);
    auto E = dmrg(psi0,H,sweeps,{"Quiet",true,"DebugLevel",3});
    printfln("RESULT s=%.0e E/s=%.12f error=%.1e",s,E/s,std::abs(E/s+2.4935771338879));
    return 0;
    }

# cpp/patched/itensor/iterativesolvers.h is a copy of mpscpp3/ITensor/itensor/iterativesolvers.h with this diff (diff patched original):
24,25d23
< inline double g_davidson_scale = 1.0;
< 
112,115d109
< #ifdef REL_CONVERGED
<     Real Approx0 = 1E-12*g_davidson_scale;
<     errgoal_ *= g_davidson_scale;
< #else
117d110
< #endif
331,333d323
< #ifdef REL_RANDOMIZE
<             if(qnrm < 1E-10*g_davidson_scale)
< #else
335d324
< #endif

# build and run (native, not Python, OMP_NUM_THREADS=1):
#   IT=<repo>/src/dmrgpy/mpscpp3/ITensor
#   for v in stock:"" randomize:"-DREL_RANDOMIZE" converged:"-DREL_CONVERGED" both:"-DREL_RANDOMIZE -DREL_CONVERGED"; do
#     /usr/bin/g++ -m64 -std=c++17 -fconcepts -O2 -DNDEBUG ${v#*:} -I./patched -I"$IT" 04_davidson_exits.cc -o 04_${v%%:*} -L"$IT/lib" -litensor -lpthread -lblas -llapack; done
#   for name in stock randomize converged both; do for s in 1 1e-6 1e-8 1e-10 1e-12; do OMP_NUM_THREADS=1 ./04_$name $s > run_${name}_$s.log 2>&1;
#     echo "$name s=$s randomizing=$(grep -c 'randomizing' run_${name}_$s.log) small_residual_exits=$(grep -c 'small residual' run_${name}_$s.log) maxiter_exits=$(grep -c 'ii == actual_maxiter' run_${name}_$s.log) $(grep RESULT run_${name}_$s.log)"; done; done | tee ../04_davidson_exits.out

# ==== 05_julia_live.py ====
# small-units on julia_live, one process: gs_energy()/s of s*H for the
# 6-site S=1/2 Heisenberg chain at s = 1, 3e-7 (where v2 and v3 fail) and
# 1e-9, fresh chain per s, and vev(s*H)/s on the scale-1 state (the MPO
# alone), against the exact -2.493577.
import numpy as np
from dmrgpy import spinchain

L = 6
def heis(sc):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

try:
    for s in (1.0, 3e-7, 1e-9):
        sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="julia_live")
        sc.maxm = 30; sc.nsweeps = 10
        sc.set_hamiltonian(s*heis(sc))
        print("julia_live s=%.0e  gs_energy()/s=%.10f" % (s, np.real(sc.gs_energy())/s), flush=True)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="julia_live")
    sc.maxm = 30; sc.nsweeps = 10
    h = heis(sc)
    sc.set_hamiltonian(h)
    ref = np.real(sc.vev(h))
    print("julia_live vev(s*H)/s / vev(H) on the scale-1 state: %s" % " ".join(
          "%.0e:%.4f" % (s, np.real(sc.vev(s*h))/s/ref) for s in (1e-6, 5e-7, 3e-7, 1e-7, 1e-9)))
except Exception as ex:
    print("julia_live raised %s: %s" % (type(ex).__name__, str(ex)[:200]))
````

Observed, before:

````
==== 02_small_units_onset.before.out (parent tree; rows at s <= 1e-8 are the clean-threshold drop) ====
[run3] slot 2 acquired
v=ED     s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=ED     s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.3e-15
v=ED     s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.0e-15
v=ED     s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=ED     s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=6.2e-15
v=ED     s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=ED     s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=ED     s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=1e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=ED     s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=1.2500000000  error=2.5e+00
v=ED     s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=1.2500000000  error=2.5e+00
v=ED     s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=1.2500000000  error=2.5e+00
v=ED     s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=1.2500000000  error=2.5e+00
v=ED     s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=1.2500000000  error=2.5e+00
v=python s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=python s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=python s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.6e-15
v=python s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=python s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=python s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=python s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.7e-13
v=python s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.2e-12
v=python s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.8e-12
v=python s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.8e-12
v=python s=1e-07  gs_energy()/s=-2.4935771338  <wf|H|wf>=-2.4935771338  error=6.2e-11
v=python s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1430080260  error=2.5e+00
v=python s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1430080260  error=2.5e+00
v=python s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1430080260  error=2.5e+00
v=python s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1430080260  error=2.5e+00
v=python s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1430080260  error=2.5e+00
v=3      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=3      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=3      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=3      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=6.0e-14
v=3      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-11
v=3      s=1e-06  gs_energy()/s=-2.4935771336  <wf|H|wf>=-2.4935771336  error=2.7e-10
v=3      s=6e-07  gs_energy()/s=-1.2841239778  <wf|H|wf>=-1.6186325401  error=1.2e+00
v=3      s=5e-07  gs_energy()/s=-1.1678737333  <wf|H|wf>=-1.3005744613  error=1.3e+00
v=3      s=3e-07  gs_energy()/s=-1.2499999976  <wf|H|wf>=-1.2499792096  error=1.2e+00
v=3      s=1e-07  gs_energy()/s=-1.2499999926  <wf|H|wf>=-1.2499937733  error=1.2e+00
v=3      s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1978425323  error=2.5e+00
v=3      s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=-0.0839250469  error=2.5e+00
v=3      s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=0.0855958671  error=2.5e+00
v=3      s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=-0.1444794419  error=2.5e+00
v=3      s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=0.1360720490  error=2.5e+00
v=2      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=2      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=2      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.0e-15
v=2      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.1e-13
v=2      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.5e-11
v=2      s=1e-06  gs_energy()/s=-2.4935771323  <wf|H|wf>=-2.4935771323  error=1.6e-09
v=2      s=6e-07  gs_energy()/s=-2.4935771298  <wf|H|wf>=-2.4935771298  error=4.1e-09
v=2      s=5e-07  gs_energy()/s=-2.4935771299  <wf|H|wf>=-2.4935771299  error=4.0e-09
v=2      s=3e-07  gs_energy()/s=-1.7469795528  <wf|H|wf>=-2.3972468362  error=7.5e-01
v=2      s=1e-07  gs_energy()/s=-1.2499999993  <wf|H|wf>=-1.2499884758  error=1.2e+00
v=2      s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=0.9759927757  error=2.5e+00
v=2      s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=0.6450172051  error=2.5e+00
v=2      s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=0.3412516392  error=2.5e+00
v=2      s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=0.3906018362  error=2.5e+00
v=2      s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=0.4976087393  error=2.5e+00

==== 03_small_units_mpo.before.out (parent tree; rows at s <= 1e-8 in A, C and D are the clean-threshold drop) ====
[run3] slot 1 acquired
A. vev(s*X)/s / vev(X) on the scale-1 ground state
  v=3  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.6667 5e-07:0.6667 4e-07:0.3333 3e-07:0.3333 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
  v=2  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974  XX=-0.831192  YY=-0.831192
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    YY            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
  v=python  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly
  s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.0e-14
  s=1e-06  gs_energy()/s=-2.4935771338  <wf|H|wf>=-2.4935771338  error=1.2e-10
  s=6e-07  gs_energy()/s=-2.4935771335  <wf|H|wf>=-2.4935771335  error=4.2e-10
  s=5e-07  gs_energy()/s=-2.4935771303  <wf|H|wf>=-2.4935771303  error=3.6e-09
  s=3e-07  gs_energy()/s=-2.4935771328  <wf|H|wf>=-2.4935771328  error=1.1e-09
  s=1e-07  gs_energy()/s=-2.4935770071  <wf|H|wf>=-2.4935770071  error=1.3e-07
  s=1e-08  gs_energy()/s=-2.4935734241  <wf|H|wf>=-2.4935734241  error=3.7e-06
  s=1e-09  gs_energy()/s=-2.4914855584  <wf|H|wf>=-2.4914855584  error=2.1e-03
  s=1e-10  gs_energy()/s=-2.2219245438  <wf|H|wf>=-2.2219245438  error=2.7e-01
  s=1e-11  gs_energy()/s=-2.2531646301  <wf|H|wf>=-2.2531646301  error=2.4e-01
  s=1e-12  gs_energy()/s=-1.4235726492  <wf|H|wf>=-1.4235726492  error=1.1e+00
C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond
  ED at s=1: -2.3275944192
  v=3  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:2.0e-09 3e-07:2.2e-08 1e-07:3.4e-08 1e-08:2.3e+00 1e-10:2.3e+00 1e-11:2.3e+00 1e-12:2.3e+00
  v=2  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:4.2e-09 3e-07:3.4e-08 1e-07:2.4e-07 1e-08:2.3e+00 1e-10:2.3e+00 1e-11:2.3e+00 1e-12:2.3e+00
D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator
  stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e+00 1e-09:2.5e+00 1e-10:2.5e+00 1e-11:2.5e+00 1e-12:2.5e+00
  unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.5e+00 1e-09:2.5e+00 1e-10:2.5e+00 1e-11:2.5e+00 1e-12:2.5e+00
````

Observed, the other runs of the locate stage:

````
==== 02_small_units_onset.after.out (working tree, clean-threshold fixed) ====
[run3] slot 2 acquired
v=ED     s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=ED     s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.3e-15
v=ED     s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.0e-15
v=ED     s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=ED     s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=6.2e-15
v=ED     s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=ED     s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=ED     s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=ED     s=1e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=ED     s=1e-08  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=ED     s=1e-09  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=ED     s=1e-10  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=ED     s=1e-11  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=ED     s=1e-12  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-15
v=python s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=python s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=python s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.6e-15
v=python s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=python s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=python s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=python s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.7e-13
v=python s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.2e-12
v=python s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.8e-12
v=python s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.8e-12
v=python s=1e-07  gs_energy()/s=-2.4935771338  <wf|H|wf>=-2.4935771338  error=6.2e-11
v=python s=1e-08  gs_energy()/s=-2.4935771314  <wf|H|wf>=-2.4935771314  error=2.5e-09
v=python s=1e-09  gs_energy()/s=-2.4935759285  <wf|H|wf>=-2.4935759285  error=1.2e-06
v=python s=1e-10  gs_energy()/s=-2.4935489641  <wf|H|wf>=-2.4935489655  error=2.8e-05
v=python s=1e-11  gs_energy()/s=-2.4833640664  <wf|H|wf>=-2.4833639894  error=1.0e-02
v=python s=1e-12  gs_energy()/s=0.1431039529  <wf|H|wf>=0.1430080260  error=2.6e+00
v=3      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=3      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=3      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.0e-15
v=3      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.1e-14
v=3      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.1e-12
v=3      s=1e-06  gs_energy()/s=-2.4935771334  <wf|H|wf>=-2.4935771334  error=5.4e-10
v=3      s=6e-07  gs_energy()/s=-1.3297714242  <wf|H|wf>=-1.7045581124  error=1.2e+00
v=3      s=5e-07  gs_energy()/s=-1.3151803810  <wf|H|wf>=-1.7258121032  error=1.2e+00
v=3      s=3e-07  gs_energy()/s=-1.2499999988  <wf|H|wf>=-1.2500120906  error=1.2e+00
v=3      s=1e-07  gs_energy()/s=-1.2499999956  <wf|H|wf>=-1.2499987381  error=1.2e+00
v=3      s=1e-08  gs_energy()/s=-1.2499953998  <wf|H|wf>=-1.2487404883  error=1.2e+00
v=3      s=1e-09  gs_energy()/s=-1.2498014082  <wf|H|wf>=-1.2526617712  error=1.2e+00
v=3      s=1e-10  gs_energy()/s=-1.1046520844  <wf|H|wf>=-0.8677643695  error=1.4e+00
v=3      s=1e-11  gs_energy()/s=-1.1742360199  <wf|H|wf>=-1.1597213669  error=1.3e+00
v=3      s=1e-12  gs_energy()/s=-0.8391956045  <wf|H|wf>=-0.4550315190  error=1.7e+00
v=2      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=2      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=2      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=2      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=2      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.2e-13
v=2      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.8e-12
v=2      s=1e-06  gs_energy()/s=-2.4935771329  <wf|H|wf>=-2.4935771329  error=1.0e-09
v=2      s=6e-07  gs_energy()/s=-2.4935771320  <wf|H|wf>=-2.4935771320  error=1.9e-09
v=2      s=5e-07  gs_energy()/s=-2.4935771312  <wf|H|wf>=-2.4935771312  error=2.7e-09
v=2      s=3e-07  gs_energy()/s=-1.7469795185  <wf|H|wf>=-2.3972444300  error=7.5e-01
v=2      s=1e-07  gs_energy()/s=-1.2499999792  <wf|H|wf>=-1.2501070218  error=1.2e+00
v=2      s=1e-08  gs_energy()/s=-1.2499967227  <wf|H|wf>=-1.2503754626  error=1.2e+00
v=2      s=1e-09  gs_energy()/s=-1.2499318808  <wf|H|wf>=-1.2607640889  error=1.2e+00
v=2      s=1e-10  gs_energy()/s=-1.2353273309  <wf|H|wf>=-1.3150341534  error=1.3e+00
v=2      s=1e-11  gs_energy()/s=-1.1878739023  <wf|H|wf>=-1.3241841309  error=1.3e+00
v=2      s=1e-12  gs_energy()/s=-0.7443806000  <wf|H|wf>=-0.6005520982  error=1.7e+00

==== 03_small_units_mpo.after.out (working tree) ====
[run3] slot 1 acquired
A. vev(s*X)/s / vev(X) on the scale-1 ground state
  v=3  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.6667 5e-07:0.6667 4e-07:0.3333 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333 1e-12:0.3333
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
  v=2  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974  XX=-0.831192  YY=-0.831192
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333 1e-12:0.3333
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    XX            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    YY            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
  v=python  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly
  s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.5e-14
  s=1e-06  gs_energy()/s=-2.4935771337  <wf|H|wf>=-2.4935771337  error=1.8e-10
  s=6e-07  gs_energy()/s=-2.4935771319  <wf|H|wf>=-2.4935771319  error=2.0e-09
  s=5e-07  gs_energy()/s=-2.4935771330  <wf|H|wf>=-2.4935771330  error=9.2e-10
  s=3e-07  gs_energy()/s=-2.4935771197  <wf|H|wf>=-2.4935771197  error=1.4e-08
  s=1e-07  gs_energy()/s=-2.4935771239  <wf|H|wf>=-2.4935771239  error=1.0e-08
  s=1e-08  gs_energy()/s=-2.4935716713  <wf|H|wf>=-2.4935716713  error=5.5e-06
  s=1e-09  gs_energy()/s=-2.4929802589  <wf|H|wf>=-2.4929802589  error=6.0e-04
  s=1e-10  gs_energy()/s=-2.3654432362  <wf|H|wf>=-2.3654432362  error=1.3e-01
  s=1e-11  gs_energy()/s=-2.2890382756  <wf|H|wf>=-2.2890382756  error=2.0e-01
  s=1e-12  gs_energy()/s=-1.8777708221  <wf|H|wf>=-1.8777708221  error=6.2e-01
C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond
  ED at s=1: -2.3275944192
  v=3  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:1.4e-09 3e-07:4.0e-08 1e-07:2.1e-07 1e-08:4.4e-05 1e-10:5.7e-01 1e-11:2.5e-01 1e-12:1.6e+00
  v=2  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:9.9e-10 3e-07:2.0e-09 1e-07:2.5e-07 1e-08:2.8e-05 1e-10:2.4e-01 1e-11:5.5e-02 1e-12:1.3e+00
D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator
  stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e-09 1e-09:1.2e-06 1e-10:2.8e-05 1e-11:1.0e-02 1e-12:2.6e+00
  unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.4e-14 1e-09:3.9e-14 1e-10:2.9e-12 1e-11:2.0e-12 1e-12:4.4e-08

==== 04_davidson_exits.out (standalone C++ against the vendored v3 libitensor.a; stock = the copied header with no macro, i.e. the shipped davidson) ====
stock s=1 randomizing=141 small_residual_exits=0 maxiter_exits=100 RESULT s=1e+00 E/s=-2.493577133888 error=2.8e-14
stock s=1e-6 randomizing=175 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-06 E/s=-2.493577133444 error=4.4e-10
stock s=1e-8 randomizing=186 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-08 E/s=-2.493569746516 error=7.4e-06
stock s=1e-10 randomizing=200 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-10 E/s=-1.969214891329 error=5.2e-01
stock s=1e-12 randomizing=100 small_residual_exits=100 maxiter_exits=0 RESULT s=1e-12 E/s=-1.272809114449 error=1.2e+00
randomize s=1 randomizing=137 small_residual_exits=0 maxiter_exits=100 RESULT s=1e+00 E/s=-2.493577133888 error=2.7e-14
randomize s=1e-6 randomizing=61 small_residual_exits=10 maxiter_exits=18 RESULT s=1e-06 E/s=-2.493577133888 error=2.6e-14
randomize s=1e-8 randomizing=61 small_residual_exits=9 maxiter_exits=10 RESULT s=1e-08 E/s=-2.493577133888 error=2.7e-14
randomize s=1e-10 randomizing=0 small_residual_exits=9 maxiter_exits=7 RESULT s=1e-10 E/s=-2.493577133888 error=2.4e-14
randomize s=1e-12 randomizing=0 small_residual_exits=21 maxiter_exits=0 RESULT s=1e-12 E/s=-2.493577133888 error=2.5e-14
converged s=1 randomizing=138 small_residual_exits=0 maxiter_exits=100 RESULT s=1e+00 E/s=-2.493577133888 error=2.8e-14
converged s=1e-6 randomizing=176 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-06 E/s=-2.493577133400 error=4.9e-10
converged s=1e-8 randomizing=178 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-08 E/s=-2.493572635959 error=4.5e-06
converged s=1e-10 randomizing=199 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-10 E/s=-2.145067468771 error=3.5e-01
converged s=1e-12 randomizing=200 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-12 E/s=-2.229333666609 error=2.6e-01
both s=1 randomizing=138 small_residual_exits=0 maxiter_exits=100 RESULT s=1e+00 E/s=-2.493577133888 error=2.4e-14
both s=1e-6 randomizing=143 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-06 E/s=-2.493577133888 error=2.7e-14
both s=1e-8 randomizing=138 small_residual_exits=0 maxiter_exits=100 RESULT s=1e-08 E/s=-2.493577133888 error=2.6e-14
both s=1e-10 randomizing=0 small_residual_exits=71 maxiter_exits=27 RESULT s=1e-10 E/s=-2.493577133888 error=2.3e-14
both s=1e-12 randomizing=0 small_residual_exits=76 maxiter_exits=23 RESULT s=1e-12 E/s=-2.493577133888 error=2.7e-14

==== 05_julia_live.out (working tree, one run; the file opens with 173 lines of juliapkg's dependency-resolution log, then) ====
julia_live s=1e+00  gs_energy()/s=-2.4935771339
julia_live s=3e-07  gs_energy()/s=-2.4935771339
julia_live s=1e-09  gs_energy()/s=-2.4935771338
julia_live vev(s*H)/s / vev(H) on the scale-1 state: 1e-06:1.0000 5e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-09:1.0000
````

**Reviewer of the locate stage (CONFIRMED)**: The location stands, and each of the three mechanisms is separated by a measurement I reran myself.

(a) The MPO. Part A of 03 (vev(s*X)/s on the scale-1 ground state, no solver involved) reproduces the fix agent's table exactly. On v3, H reads 0.6667 of itself at 6e-7 and 5e-7 and 0.3333 from 4e-7 down; on v2 it reads 0.6667 at 3e-7 and 0.3333 from 1e-7 down; ZZ, a lone XX and a lone YY stay at 1.0000; "python" is 1.0000 everywhere. I read the cited code. autompo.cc:1176 takes Cutoff=1E-13 from args, :1223 calls truncate(D,maxdim,mindim,cutoff) with the header defaults absoluteCutoff=false and doRelCutoff=false (decomp.h:422), and decomp.cc divides both P and the cutoff by P0, so the discarded sum of squared singular values is compared against an absolute 1e-13. In addition, autompo.cc:1171 sets a local `Real eps = 1E-14` and :1287 skips any term with |coef| below it. That last threshold I measured for the first time here: vev(eps*Sz0) is exactly 0 on v2 and v3 at eps <= 1e-15 and 1.000000 at 1e-14, while ED and "python" stay linear down to 1e-30.

(b) Davidson. I reran the prebuilt diagnostic binaries natively. The randomization test made relative, alone, gives 2.7e-14 at s=1e-8 and 2.8e-14 at 1e-10. The convergence tests made relative, alone, give 8.2e-6 and 6.3e-2, against stock 5.7e-6 and 1.2e-1. I read the lines: `if(qnrm < 1E-10)` at v3 iterativesolvers.h:324 and v2 :295.

(b') "python". Part D reproduces to every digit: stock 2.5e-9 / 1.2e-6 / 2.8e-5 / 1.0e-2 / 2.6 against unit-normalized 2.4e-14 / 3.9e-14 / 2.9e-12 / 2.0e-12 / 4.4e-8.

Narrowed sizes. Parts B and C are draws: my B gives 1.1e-5 at 1e-8 and 5.1e-1 at 1e-10, against the record's 5.5e-6 and 0.13. Quote them as ranges.

Not re-verified: julia_live, which I skipped deliberately. The fix agent's single run already triggered a juliapkg resolution nobody asked for.

No code changed for this item. My comments on cpp_plan:

(1) Put the gate lower, at a power of two near 1e-2 (for example 2^-7), not at cmax >= 1. The measured cost of not scaling is 4.4e-16 (v3) and 8.9e-16 (v2) at s=1e-2, and 4.0e-15 and 0.0 at 1e-3, so scaling between 1e-2 and 1 buys nothing. What the cmax >= 1 gate does is move every model whose largest AutoMPO coefficient lies in [2^-7,1) onto a path where Davidson's decisions change by design, so the bits change. That includes every J=0.5 model and the J=1 XX chain, whose realified terms reach AutoMPO as 0.5*S+S- and 0.5*S-S+ (cmax=0.5 after AutoMPO merges them). The plan itself offers this as an option; I recommend taking it. Above the gate, byte identity holds by construction, since to_mpo_unit returns toMPO(ampo,args) itself.

(2) Passing a small Cutoff through build_mpo is not an alternative. It would cover truncate(), but not the local eps=1E-14 isZero, which I measured biting at eps <= 1e-15, and not Davidson. The rescale is the right route.

(3) Single-term AutoMPOs (the shift and Id MPOs) are exempt from truncate() but not from isZero(1E-14). This matters only below 1e-14 (by reading).

(4) Real-time evolution is outside the plan, and I measured it to be unit-dependent through a separate absolute threshold (new issue 1). v3's applyExp (iterativesolvers.h) compares an error estimate that carries nrm, the state norm, against ErrGoal=1E-10, and has NormCutoff=1e-7 on the Lanczos beta, both absolute. The measured TD correlator of s*H, with es, delta and dt scaled, is 1.2e-9, 8.7e-8 and 7.0e-6 off on v3 at s=1e-2, 1e-4 and 1e-6, and on "python" 6.7e-9, 2.4e-6, 2.8e-5 and 8.8e-2 at 1e-2 to 1e-8. So if real-time is meant to be as covariant as E0, hscale_up_ should reach tdvp_step as well (evolve under up*H with dt/up), and plan (3) should cover pyitensor/tdvp.py::_lanczos_expm_multiply.

(5) The plan's accessors exist as it assumes: AutoMPO::terms()/sites() in both versions (v3 autompo.h:200-203, v2 :200/:203), and operator*=(MPO&,Real) scaling one tensor (v2 mpo.h:153, v3 mpo.h:133). The dmrg() call list is complete: v3 592, 710, 746, 11707; v2 176, 222, 1145. The Weight*up treatment of excited_states and the -dmrg(...,-up*H)/up form of maximum_energy read correctly.

(6) The claim that SVD is homogeneous under a power-of-two scaling is plausible but not measured, and it only matters below the gate.

Reviewer's evidence:

````
R=<scratch>/review/scale
R/03_small_units_mpo.after.out:
  v=3 H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.6667 5e-07:0.6667 4e-07:0.3333 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333 1e-12:0.3333
  v=2 H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.6667 1e-07:0.3333 ...
  v=2 XX            ... 3e-07:1.0000 1e-07:1.0000 ...
  B. s=1e-08  gs_energy()/s=-2.4935662807 ... error=1.1e-05 ; s=1e-10 ... error=5.1e-01
  D. stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e-09 1e-09:1.2e-06 1e-10:2.8e-05 1e-11:1.0e-02 1e-12:2.6e+00
     unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.4e-14 1e-09:3.9e-14 1e-10:2.9e-12 1e-11:2.0e-12 1e-12:4.4e-08
R/cpp/04_rerun.out (native, OMP_NUM_THREADS=1):
  stock s=1e-8 randomizing=183 ... RESULT s=1e-08 E/s=-2.493571437323 error=5.7e-06
  stock s=1e-10 randomizing=199 ... RESULT s=1e-10 E/s=-2.373812229138 error=1.2e-01
  randomize s=1e-8 ... RESULT s=1e-08 E/s=-2.493577133888 error=2.7e-14
  randomize s=1e-10 ... RESULT s=1e-10 E/s=-2.493577133888 error=2.8e-14
  converged s=1e-8 ... error=8.2e-06 ; converged s=1e-10 ... error=6.3e-02
Read: mpscpp3/ITensor/itensor/iterativesolvers.h:324 `if(qnrm < 1E-10)`; mpscpp2 :295 same; mpscpp3/ITensor/itensor/mps/autompo.cc:1171 `Real eps = 1E-14;`, :1176 `Real cutoff = args.getReal("Cutoff",1E-13);`, :1223 `truncate(D,maxdim,mindim,cutoff);`, :1287 `if(isZero(t.coef,eps)) continue;`; decomp.cc truncate: `if(absoluteCutoff || !doRelCutoff) { cutoff /= P0; }`; applyExp: `auto tol = args.getReal("ErrGoal",1E-10); ... auto beta_tol = args.getReal("NormCutoff",1e-7);`
R/r4_downstream_floors.after.out: v=2      1e-13:1.000000 1e-14:1.000000 1e-15:-0.000000 1e-16:-0.000000 1e-18:-0.000000 1e-20:-0.000000 1e-30:-0.000000
R/r6_td_small_units.after.out:
  v=python submode=TD  max|s*C_s(s*w) - C_1(w)|/max|C_1|  1e-02:6.66e-09  1e-04:2.42e-06  1e-06:2.76e-05  1e-08:8.76e-02
  v=3      submode=TD  max|s*C_s(s*w) - C_1(w)|/max|C_1|  1e-02:1.21e-09  1e-04:8.66e-08  1e-06:6.97e-06
````

The plan the locate stage handed to the C++ stage:

````
Two changes to dmrgpy's own C++ (mo_terms.h and chain_session.h of both backends), none to vendored ITensor, both gated on the largest |coefficient| of the terms so that every operator whose largest coefficient is at least 1 runs exactly today's code, and both using an exact power of two.

Shared helper, in both mo_terms.h (next to build_mpo: mpscpp3/mo_terms.h:306, mpscpp2/mo_terms.h:36):

    // Power of two bringing the largest |coef| into [1,2); 1.0 when it is
    // already >= 1 or the list is empty, which leaves the caller on today's path.
    inline double unit_scale_up(double cmax)
        {
        if (!(cmax > 0.0) || cmax >= 1.0) return 1.0;
        int e = 0; std::frexp(cmax,&e);      // cmax in [2^(e-1),2^e), e <= 0
        return std::ldexp(1.0,1-e);          // cmax*up in [1,2), exact
        }
    // toMPO with the operator at unit scale, so svdMPO's absolute cutoffs
    // (truncate() at Cutoff=1E-13 with doRelCutoff=false, autompo.cc ~1216,
    // and isZero(coef,1E-14), autompo.cc:1287) act relative to it.
    MPO inline
    to_mpo_unit(AutoMPO const& ampo, Args const& args)
        {
        double cmax = 0.0;
        for (auto const& t : ampo.terms()) cmax = std::max(cmax,std::abs(t.coef));
        double up = unit_scale_up(cmax);
        if (up == 1.0) return toMPO(ampo,args);   // byte-identical path
        auto scaled = AutoMPO(ampo.sites());
        for (auto t : ampo.terms()) { t.coef *= up; scaled.add(t); }
        auto W = toMPO(scaled,args);
        W *= 1.0/up;          // operator*=(MPO&,Real): one tensor, exact
        return W;
        }
(v2: the same with toMPO<ITensor>(...) and MPO; AutoMPO::terms() and sites() exist in both, autompo.h:201-204.)

(1) Mechanism (a), the MPO. build_mpo's `return toMPO(ampo,{"MaxDim",mpomaxm,"Exact",false});` (v3) and `return toMPO<ITensor>(ampo,{"Maxm",mpomaxm,"Exact",false});` (v2) become `to_mpo_unit(ampo,{...})`. That covers every term-built MPO: on v3 through mpo_from_terms (chain_session.h:4928, including sector_terms), so set_hamiltonian (536), vev/correlators/KPM vertices, build_operator (1223) and the rest; on v2 the 18 direct build_mpo calls (chain_session.h 153, 259, 260, 407, 427, 446, 533, 554, 660, 666, 667, 702, 756, 757, 777, 875, 876, 897). The evolution Hamiltonians built from an AutoMPO with an appended `-EGS,"Id",1` go to toMPO directly and need the same call, with no other change since the appended term is in ampo.terms() by then: v3 1311-1313 (quench), 1343-1344 (evolve_and_measure), 1398-1400 and 1482-1484 (Hshift); v2 663-665 and 700-701, where `MPO(ampo)` becomes `to_mpo_unit(ampo,{})`. The single-term AutoMPOs (Id, shifts, one-site operators: v3 2025-2144, 11750-12038, v2 602-604, 760-846, 1206-1310) have one channel per bond, which truncate() returns early on, and can stay.

(2) Mechanism (b), the local eigensolver. Add a member `double hscale_up_ = 1.0;` set in Chain::set_hamiltonian(terms) (v3 534-537, v2 151-158) to unit_scale_up(max |coef| of terms), and to 1.0 in v3's set_hamiltonian_mpo (an already-built MPO, whose scale the session cannot read; left as today and documented). At every dmrg() on H_, under `if (hscale_up_ != 1.0)`, solve on `hscale_up_*H_` and divide the energy back by hscale_up_: gs_energy (v3 592, v2 176); excited_states (v3 746, v2 222), where the overlap penalty `Weight` is in energy units and becomes `weight*hscale_up_` (the energies are then measured with innerC on the unscaled H_, unchanged); gs_energy_generalized (v3 710), solving on `hscale_up_*Heff` (lam comes from innerC and needs nothing); maximum_energy (v3 11707, v2 1145), `-dmrg(psi,(-hscale_up_)*H_,...)/hscale_up_`. minimum_energy goes through gs_energy. The MPS dmrg() returns is normalized, so nothing downstream sees the factor; the noise term, whose deltaRho scales as the square of the Hamiltonian, also gets back its scale-1 strength.

Why byte-identical at ordinary scales: when the largest coefficient is at least 1, unit_scale_up returns 1.0, to_mpo_unit returns toMPO(ampo,args) itself and every dmrg() branch is skipped, so no floating-point operation changes; Heisenberg at J=1, Hubbard at t=1, the Ising chain at J=1, every example and test at those scales, are untouched. Below 1, multiplying by a power of two is exact in IEEE arithmetic (no underflow at these scales), so every coefficient and every matrix svdMPO decomposes is exactly 2^k times today's and, LAPACK's SVD being homogeneous under such a scaling away from its underflow guard (norms below about 1e-146), so is every singular value; the only decisions that change are the comparisons against the absolute 1e-13, 1e-14 and 1e-10, which is the fix. If byte identity is also wanted for largest coefficients in [0.5,1) (J=0.5 models), gate at a lower power of two instead, at the cost of scale covariance between that gate and 1.

Evidence that this is sufficient: vev(s*X)/s on a fixed state isolates (a) and matches svdMPO's truncation to the digit; the standalone diagnostic 04 (ITensor's own dmrg() on an MPO multiplied by s, with only davidson's randomization test made relative) is 2.4e-14 to 2.7e-14 off exact at every s from 1 to 1e-12, and solving on H*up is by construction the unit-scale calculation. Tests for that stage: flip the strict xfails v3-1e-08 and v2-1e-08 in tests/test_audit_2026_09_25_scale.py::test_small_units_ground_state_still_open into plain tests, and add vev(s*H)/s = vev(H) on v2/v3 at s from 6e-7 to 1e-10, gs_energy()/s of the Ising chain to 1e-10 at s=1e-8, and excited_states at s=1e-8 against ED.

(3) Not C++, the \"python\" half, same shape: pyitensor/dmrg.py::_lanczos_ground_state's `if beta < tol: break` and `abs(cur_eval - prev_eval) < tol * max(1.0, abs(cur_eval))` are absolute once |eval| < 1. The mirror of (2) is to solve on the Hamiltonian MPO times unit_scale_up(cmax) in pyitensor/chain.py's ground-state, excited-state and band-edge entry points and divide the energy back, byte-identical at cmax >= 1; the measured stand-in (the unit-normalized Lanczos in 03 part D) takes s=1e-9 from 1.2e-6 to 3.9e-14 and s=1e-10 from 2.8e-5 to 2.9e-12. That then flips the python-1e-10 strict xfail.
````

#### Fix stage

Repro (<scratch>/scale-cpp/02_small_units_onset.py and 03_small_units_mpo.py (byte-identical to the locate stage's, cmp), 10_ordinary_scale.py, 11_small_units_consumers.py, 13_left_open.py, 14_truncated_below_one.py, 15_xx_spread.py; <scratch>/scale-cpp/repair/20_mixed_scale.py (and cmp-identical copies of 02, 03, 10, 11, 13 there); trees: parent <scratch>/parent/src, hybrid <scratch>/scale-cpp/hybrid/src, first-pass <scratch>/scale-cpp/repair/firstpass/src):

````python
# 02_small_units_onset.py and 03_small_units_mpo.py are byte-identical (cmp) to the locate stage's scripts of the same name, carried in full in that stage's record; 02 was run with the arguments `3 2`.

# ==== 10_ordinary_scale.py ====
# small-units, the ordinary-scale guard: one fixed set of calculations at
# energy scale 1 on v3 and v2, which the fix must leave as they were.  A
# 6-site S=1/2 Heisenberg chain at J=1 (largest AutoMPO coefficient 1 on both
# backends, so the scaled path is never taken): gs_energy, two vev, a KPM
# dynamical correlator and a short real-time evolution of <Sz1Sz2> from
# Sz0|GS> (an even observable, so spin-flip symmetry does not zero it), against
# ED.  Then two models that ARE at an ordinary scale but whose largest AutoMPO
# coefficient is below 1, so they take the scaled path after the fix: the
# Heisenberg chain at J=0.5, and the XX chain at J=1, which v3 realifies into
# S+S-/S-S+ at coefficient 0.5.  Printed as repr and as the error against ED,
# so a run-to-run comparison and a before/after comparison read off the same
# lines.  Pinned sweep schedule.
import numpy as np
from dmrgpy import spinchain, cppext, timedependent

L = 6
def heis(sc, J=1.0):
    h = 0
    for i in range(L-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    return h
def xx(sc):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    return h

def chain(v):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 10
    return sc

es = np.linspace(-0.5, 4.0, 7)
for v in [v for v in (3, 2) if cppext.available(v)]:
    sc = chain(v)
    sc.set_hamiltonian(heis(sc))
    e0 = np.real(sc.gs_energy())
    e0ed = np.real(sc.gs_energy(mode="ED"))
    a = np.real(sc.vev(sc.Sz[0]*sc.Sz[1]))
    aed = np.real(sc.vev(sc.Sz[0]*sc.Sz[1], mode="ED"))
    b = np.real(sc.vev(sc.Sx[2]*sc.Sx[3] + sc.Sy[2]*sc.Sy[3]))
    bed = np.real(sc.vev(sc.Sx[2]*sc.Sx[3] + sc.Sy[2]*sc.Sy[3], mode="ED"))
    print("v=%s J=1 Heisenberg  gs_energy=%r  err=%.1e" % (v, e0, abs(e0 - e0ed)))
    print("v=%s J=1 Heisenberg  vev(Sz0Sz1)=%r  err=%.1e  vev(XY23)=%r  err=%.1e" % (
          v, a, abs(a - aed), b, abs(b - bed)))
    x, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.2)
    x, yed = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.2,
                                         mode="ED")
    y, yed = np.asarray(y), np.asarray(yed)
    print("v=%s J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = %s  max|DMRG-ED|/peak=%.1e" % (
          v, " ".join("%.10f" % t for t in np.real(y)),
          np.max(np.abs(y - yed))/np.max(np.abs(yed))))
    O = sc.Sz[1]*sc.Sz[2]
    wf = sc.applyoperator(sc.Sz[0], sc.get_gs())
    t, m = timedependent.evolve_and_measure(sc, operator=O, nt=10, dt=0.05, wf=wf)
    t, med = timedependent.evolve_and_measure(sc, mode="ED", operator=O, nt=10,
                                              dt=0.05, wf=sc.applyoperator(
                                              sc.Sz[0], sc.get_gs(mode="ED"), mode="ED"))
    m, med = np.real(np.asarray(m)), np.real(np.asarray(med))
    print("v=%s J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: %s  max|DMRG-ED|=%.1e" % (
          v, " ".join("%.10f" % q for q in m), np.max(np.abs(m - med))))
    for label, build in (("J=0.5 Heisenberg", lambda c: heis(c, 0.5)),
                         ("J=1 XX chain", xx)):
        sc = chain(v)
        sc.set_hamiltonian(build(sc))
        e = np.real(sc.gs_energy())
        eed = np.real(sc.gs_energy(mode="ED"))
        print("v=%s %s  gs_energy=%r  ED=%r  err=%.1e" % (v, label, e, eed, abs(e - eed)),
              flush=True)

# ==== 11_small_units_consumers.py ====
# small-units, the consumers of the band edges and of the local solver on
# v3 and v2: s*H for the 6-site S=1/2 Heisenberg chain, fresh chain per s,
# each against mode="ED" on the same chain (the ED route is scale-free since
# the clean-threshold fix).  (i) The KPM dynamical correlator C[Sz0,Sz0] at
# es*s and delta*s, which reads both band edges (the upper one from the -H
# solve) and sums a shift into H as an MPO; the density scales as 1/s, so
# s*C is compared with the ED s*C, relative to the ED peak.  (ii) The first
# four levels from get_excited(), whose overlap penalty is an energy, /s
# against ED /s.  (iii) v3 only, gs_energy_generalized(A)/s with the metric
# A = 1 + 0.2*Sz0, against the dense generalized eigenproblem.  Pinned sweep
# schedule.
import numpy as np
import scipy.linalg as sla
from dmrgpy import spinchain, cppext

L = 6
def heis(sc):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

def chain(v, s):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 10
    sc.set_hamiltonian(s*heis(sc))
    return sc

def dense(M):
    return np.asarray(M.todense() if hasattr(M, "todense") else M)

es = np.linspace(-0.5, 4.0, 46)
for v in [v for v in (3, 2) if cppext.available(v)]:
    for s in (1.0, 1e-8, 1e-10):
        sc = chain(v, s)
        try:
            kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=s*es, delta=0.2*s)
            x, y = sc.get_dynamical_correlator(**kw)
            x, yed = sc.get_dynamical_correlator(mode="ED", **kw)
            y, yed = s*np.asarray(y), s*np.asarray(yed)
            kpm = "%.1e" % (np.max(np.abs(y - yed))/np.max(np.abs(yed)))
        except Exception as ex:
            kpm = "raised %s: %s" % (type(ex).__name__, str(ex)[:60])
        sc = chain(v, s)
        ex4 = np.real(np.asarray(sc.get_excited(n=4)))/s
        ed4 = np.real(np.asarray(sc.get_excited(n=4, mode="ED")))/s
        row = "v=%s s=%.0e  KPM max|DMRG-ED|/peak=%s  get_excited(n=4)/s=%s  max|DMRG-ED|=%.1e" % (
              v, s, kpm, " ".join("%.8f" % e for e in ex4), np.max(np.abs(ex4 - ed4)))
        if v == 3:
            sc = chain(v, s)
            A = 1.0 + 0.2*sc.Sz[0]
            lam = np.real(sc.gs_energy_generalized(A))/s
            Hm = dense(sc.get_ED_obj().get_operator(heis(sc)))
            Am = dense(sc.get_ED_obj().get_operator(A))
            lref = np.min(sla.eigh(Hm, Am, eigvals_only=True))
            row += "  gs_energy_generalized/s=%.10f  err=%.1e" % (lam, abs(lam - lref))
        print(row, flush=True)

# ==== 13_left_open.py ====
# small-units, what the C++ fix does not cover, measured rather than read:
# (i) NH-DMRG, which gs_energy() dispatches to for a non-Hermitian H and
# which solves with its own Arnoldi (a 1e-10*(1+|E|) residual test), on
# s*(6-site Heisenberg + 0.3j*Sz0), E0/s against ED; (ii) real-time
# evolution, <Sz1Sz2>(t) from Sz0|GS> under s*H over the same physical
# times as at s=1 (dt/s), against ED at the same s.  Fresh chain per s,
# pinned sweep schedule.
import numpy as np
from dmrgpy import spinchain, cppext, timedependent

L = 6
def heis(sc):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

def chain(v):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 10
    return sc

for v in [v for v in (3, 2) if cppext.available(v)]:
    for s in (1.0, 1e-8):
        sc = chain(v)
        sc.set_hamiltonian(s*(heis(sc) + 0.3j*sc.Sz[0]))
        try:
            e = complex(sc.gs_energy())/s
            eed = complex(sc.gs_energy(mode="ED"))/s
            nh = "E0/s=%.8f%+.8fj  |DMRG-ED|=%.1e" % (e.real, e.imag, abs(e - eed))
        except Exception as ex:
            nh = "raised %s: %s" % (type(ex).__name__, str(ex)[:60])
        sc = chain(v)
        sc.set_hamiltonian(s*heis(sc))
        O = sc.Sz[1]*sc.Sz[2]
        try:
            t, m = timedependent.evolve_and_measure(
                sc, operator=O, nt=10, dt=0.05/s, wf=sc.applyoperator(sc.Sz[0], sc.get_gs()))
            t, med = timedependent.evolve_and_measure(
                sc, mode="ED", operator=O, nt=10, dt=0.05/s,
                wf=sc.applyoperator(sc.Sz[0], sc.get_gs(mode="ED"), mode="ED"))
            ev = "evolution max|DMRG-ED|=%.1e" % np.max(np.abs(np.real(m) - np.real(med)))
        except Exception as ex:
            ev = "evolution raised %s: %s" % (type(ex).__name__, str(ex)[:60])
        print("v=%s s=%.0e  NH-DMRG %s  %s" % (v, s, nh, ev), flush=True)

# ==== 14_truncated_below_one.py ====
# small-units, what the fix moves at an ordinary scale where the MPS is NOT
# exact: a 20-site S=1/2 chain at maxm=10, nsweeps=10, where truncation sets
# the error and the local solver's path can matter.  Three fresh chains per
# row, to put a before/after difference next to the run-to-run spread.  The
# J=1 Heisenberg chain takes the unscaled path (largest coefficient 1); the
# J=0.5 Heisenberg chain and, on v3, the J=1 XX chain (realified to S+S-/S-S+
# at 0.5) take the scaled one after the fix.
import numpy as np
from dmrgpy import spinchain, cppext

L = 20
def model(sc, J, zz):
    h = 0
    for i in range(L-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + zz*sc.Sz[i]*sc.Sz[i+1])
    return h

for v in [v for v in (3, 2) if cppext.available(v)]:
    for label, J, zz in (("J=1 Heisenberg", 1.0, 1.0), ("J=0.5 Heisenberg", 0.5, 1.0),
                         ("J=1 XX chain", 1.0, 0.0)):
        es = []
        for run in range(3):
            np.random.seed(run)
            sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
            sc.maxm = 10; sc.nsweeps = 10
            sc.set_hamiltonian(model(sc, J, zz))
            es.append(np.real(sc.gs_energy()))
        es = np.array(es)
        print("v=%s %-16s gs_energy over 3 runs: %s  mean=%.12f  spread=%.1e" % (
              v, label, " ".join("%.12f" % e for e in es), es.mean(), es.max() - es.min()),
              flush=True)

# ==== 15_xx_spread.py ====
# small-units, follow-up of 14: the J=1 XX chain on v3 at 20 sites, maxm=10,
# nsweeps=10, is the one row of 14 whose MPO build takes the scaled path
# after the fix (realified to S+S-/S-S+ at coefficient 0.5; the solver scale
# stays 1, since the caller's largest coefficient is 1).  Ten fresh chains,
# to tell a shift of the mean from the run-to-run spread, next to the exact
# free-fermion energy (hopping 1/2, so single-particle energies cos(pi k/21)).
import numpy as np
from dmrgpy import spinchain

L = 20
exact = np.sum(np.sort(np.cos(np.pi*np.arange(1, L+1)/(L+1)))[:L//2])
es = []
for run in range(10):
    np.random.seed(run)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=3)
    sc.maxm = 10; sc.nsweeps = 10
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    sc.set_hamiltonian(h)
    es.append(np.real(sc.gs_energy()))
es = np.array(es)
print("exact %.12f" % exact)
print("v=3 J=1 XX chain, 10 runs: mean=%.12f  std=%.1e  min=%.12f  max=%.12f  mean-exact=%.2e" % (
      es.mean(), es.std(), es.min(), es.max(), es.mean() - exact))

# ==== repair/20_mixed_scale.py ====
# small-units, repair pass: the two gaps the review measured on the first
# fix, plus one neighbour that no global scale can reach.
#  A. The gate. svdMPO's absolute cutoff acts on the channels a term opens
#     across a bond, and a term confined to one site (an energy offset, sent
#     as ('Id',site), or a field) opens none, yet it counted towards the
#     largest coefficient the first fix scaled by. 6-site S=1/2 Heisenberg
#     chain, pinned sweeps. (i) No eigensolver: on the scale-1 ground state,
#     [vev(1 + s*H) - 1]/s and [vev(Sz0 + s*H) - vev(Sz0)]/s, each over
#     vev(H). (ii) gs_energy of s*H + 1 and of s*H + Sz0, fresh chain per
#     row, against mode="ED" on the same chain, both /s after removing the
#     offset (the ED column carries the roundoff of E ~ 1, 2.2e-16/s).
#  B. The single-term shift AutoMPO of the KPM rescaling, whose one
#     coefficient svdMPO's isZero skips below an absolute 1e-14: KPM
#     C[Sz0,Sz0] of s*H at es*s and delta*s against ED, s*C relative to the
#     ED peak.
#  C. Not the gate: a weak link inside the exchange itself. Bonds J=1 except
#     the middle one, J'*S2.S3; on the scale-1 uniform ground state,
#     [vev(Hs + J'*B23) - vev(Hs)]/J' / vev(B23), Hs being the chain without
#     that bond. Only J' crosses the middle bond, so svdMPO compares its
#     channels against the absolute 1e-13 whatever the rest of H is.
import numpy as np
from dmrgpy import spinchain, cppext

L = 6
def heis(sc, skip=None):
    h = 0
    for i in range(L-1):
        if i == skip: continue
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

def chain(v):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 10
    return sc

EXACT = -2.4935771338879
cpp = [v for v in (3, 2) if cppext.available(v)]

print("A(i). no eigensolver, ratio to vev(H) on the scale-1 ground state")
for v in cpp + ["python"]:
    sc = chain(v); h = heis(sc); sc.set_hamiltonian(h)
    ref = np.real(sc.vev(h)); z0 = np.real(sc.vev(sc.Sz[0]))
    row1 = " ".join("%.0e:%.4f" % (s, (np.real(sc.vev(s*h + 1.0)) - 1.0)/s/ref)
                    for s in (1e-6, 3e-7, 1e-7, 1e-8, 1e-10))
    row2 = " ".join("%.0e:%.4f" % (s, (np.real(sc.vev(sc.Sz[0] + s*h)) - z0)/s/ref)
                    for s in (1e-6, 3e-7, 1e-7, 1e-8, 1e-10))
    print("  v=%-6s offset 1 + s*H   %s" % (v, row1))
    print("  v=%-6s field Sz0 + s*H  %s" % (v, row2), flush=True)

print("A(ii). gs_energy, fresh chain per row, against ED on the same chain")
for v in cpp:
    for s in (1e-7, 1e-8, 1e-10):
        sc = chain(v); sc.set_hamiltonian(s*heis(sc) + 1.0)
        e = (np.real(sc.gs_energy()) - 1.0)/s
        eed = (np.real(sc.gs_energy(mode="ED")) - 1.0)/s
        print("  v=%s s=%.0e  s*H+1    (E-1)/s=%.10f  ED=%.10f  |DMRG-ED|=%.1e  |DMRG-exact|=%.1e" % (
              v, s, e, eed, abs(e - eed), abs(e - EXACT)), flush=True)
    for s in (1e-7, 1e-8, 1e-10):
        sc = chain(v); sc.set_hamiltonian(s*heis(sc) + sc.Sz[0])
        e = np.real(sc.gs_energy()); eed = np.real(sc.gs_energy(mode="ED"))
        print("  v=%s s=%.0e  s*H+Sz0  (E+0.5)/s=%.10f  ED=%.10f  |DMRG-ED|/s=%.1e" % (
              v, s, (e + 0.5)/s, (eed + 0.5)/s, abs(e - eed)/s), flush=True)

print("B. KPM C[Sz0,Sz0] of s*H, max|DMRG-ED|/ED peak")
es = np.linspace(-0.5, 4.0, 46)
for v in cpp:
    row = []
    for s in (1.0, 1e-12, 1.8e-14, 1.5e-14, 1e-14, 1e-16, 1e-20):
        sc = chain(v); sc.set_hamiltonian(s*heis(sc))
        kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=s*es, delta=0.2*s)
        try:
            x, y = sc.get_dynamical_correlator(**kw)
            x, yed = sc.get_dynamical_correlator(mode="ED", **kw)
            y, yed = s*np.asarray(y), s*np.asarray(yed)
            row.append("%.1e:%.1e" % (s, np.max(np.abs(y - yed))/np.max(np.abs(yed))))
        except Exception as ex:
            row.append("%.1e:raised %s" % (s, type(ex).__name__))
    print("  v=%s  %s" % (v, " ".join(row)), flush=True)

print("C. weak link J'*S2.S3 in a J=1 chain, [vev(Hs+J'*B23)-vev(Hs)]/J'/vev(B23)")
for v in cpp + ["python"]:
    sc = chain(v); h = heis(sc); sc.set_hamiltonian(h)
    hs = heis(sc, skip=2)
    b23 = sc.Sx[2]*sc.Sx[3] + sc.Sy[2]*sc.Sy[3] + sc.Sz[2]*sc.Sz[3]
    ref = np.real(sc.vev(b23)); e0 = np.real(sc.vev(hs))
    row = " ".join("%.0e:%.4f" % (jp, (np.real(sc.vev(hs + jp*b23)) - e0)/jp/ref)
                   for jp in (1e-5, 1e-6, 3e-7, 1e-7, 1e-8))
    print("  v=%-6s %s" % (v, row), flush=True)

# trees:
#   hybrid (working-tree Python, clean-threshold fixed, with the parent's compiled extensions, cmp-identical to parent/src/dmrgpy/mpscppN/_dmrgcpp*.so):
#     rsync -a --exclude='ITensor/' --exclude='__pycache__/' --exclude='*.so' --exclude='TDVP/' <repo>/src/dmrgpy <scratch>/scale-cpp/hybrid/src/
#     cp <scratch>/parent/src/dmrgpy/mpscpp{2,3}/_dmrgcpp*.so into the matching hybrid folders
#   first-pass (working-tree Python with the first pass's compiled extensions, cmp-identical to the .so files this pass replaced; every changed .py cmp-identical to hybrid's):
#     rsync -a --exclude='ITensor/' --exclude='__pycache__/' --exclude='TDVP/' --exclude='*.o' <repo>/src/dmrgpy <scratch>/scale-cpp/repair/firstpass/src/   (before the rebuild)
# run, from the script's folder:
#   DMRGPY_SRC=<scratch>/parent/src <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.before.out                     (parent)
#   DMRGPY_SRC=<scratch>/scale-cpp/hybrid/src <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.hybrid.out          (hybrid)
#   DMRGPY_SRC=<scratch>/scale-cpp/repair/firstpass/src <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.firstpass.out  (first pass)
#   <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.after.out                                                      (working tree, rebuilt)
# the test file against a tree: copied alone to a folder with no conftest (<scratch>/scale-cpp/repair/{firstpass,hybrid}_check/), so DMRGPY_SRC decides, and run with DMRGPY_SRC=<tree> run3.sh -m pytest test_audit_2026_09_25_scale.py -q -p no:cacheprovider
````

Observed, before:

````
==== 02_small_units_onset.before.out (parent tree 8dd2198, old extensions; the rows at s <= 1e-8 read 0 for the clean-threshold reason) ====
[run3] slot 1 acquired
v=3      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=3      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.3e-15
v=3      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.5e-14
v=3      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.1e-12
v=3      s=1e-06  gs_energy()/s=-2.4935771331  <wf|H|wf>=-2.4935771331  error=7.8e-10
v=3      s=6e-07  gs_energy()/s=-1.0577261243  <wf|H|wf>=-1.3313330636  error=1.4e+00
v=3      s=5e-07  gs_energy()/s=-1.2316327789  <wf|H|wf>=-1.4975884304  error=1.3e+00
v=3      s=3e-07  gs_energy()/s=-1.2499999985  <wf|H|wf>=-1.2500073270  error=1.2e+00
v=3      s=1e-07  gs_energy()/s=-1.2499999614  <wf|H|wf>=-1.2499114001  error=1.2e+00
v=3      s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=0.2407647172  error=2.5e+00
v=3      s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=-0.6676725099  error=2.5e+00
v=3      s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=-0.2127880030  error=2.5e+00
v=3      s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=0.2036014815  error=2.5e+00
v=3      s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=-0.1214190977  error=2.5e+00
v=2      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=2      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=2      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=2      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=2      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=7.2e-14
v=2      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.2e-11
v=2      s=1e-06  gs_energy()/s=-2.4935771330  <wf|H|wf>=-2.4935771330  error=8.6e-10
v=2      s=6e-07  gs_energy()/s=-2.4935771320  <wf|H|wf>=-2.4935771320  error=1.9e-09
v=2      s=5e-07  gs_energy()/s=-2.4935771316  <wf|H|wf>=-2.4935771316  error=2.3e-09
v=2      s=3e-07  gs_energy()/s=-1.7469795955  <wf|H|wf>=-2.3972679729  error=7.5e-01
v=2      s=1e-07  gs_energy()/s=-1.2499999962  <wf|H|wf>=-1.2499631156  error=1.2e+00
v=2      s=1e-08  gs_energy()/s=0.0000000000  <wf|H|wf>=0.4878061035  error=2.5e+00
v=2      s=1e-09  gs_energy()/s=0.0000000000  <wf|H|wf>=0.5432223460  error=2.5e+00
v=2      s=1e-10  gs_energy()/s=0.0000000000  <wf|H|wf>=0.9201522642  error=2.5e+00
v=2      s=1e-11  gs_energy()/s=0.0000000000  <wf|H|wf>=0.7829088225  error=2.5e+00
v=2      s=1e-12  gs_energy()/s=0.0000000000  <wf|H|wf>=0.5172500693  error=2.5e+00

==== 02_small_units_onset.hybrid.out (hybrid tree: working-tree Python with clean-threshold fixed, parent extensions; the unmasked before) ====
[run3] slot 1 acquired
v=3      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=3      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=6.8e-14
v=3      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.4e-11
v=3      s=1e-06  gs_energy()/s=-2.4935771337  <wf|H|wf>=-2.4935771337  error=2.0e-10
v=3      s=6e-07  gs_energy()/s=-1.2715204438  <wf|H|wf>=-1.7435466864  error=1.2e+00
v=3      s=5e-07  gs_energy()/s=-1.4010959621  <wf|H|wf>=-1.9592853070  error=1.1e+00
v=3      s=3e-07  gs_energy()/s=-1.2499999935  <wf|H|wf>=-1.2500472682  error=1.2e+00
v=3      s=1e-07  gs_energy()/s=-1.2499999971  <wf|H|wf>=-1.2499958456  error=1.2e+00
v=3      s=1e-08  gs_energy()/s=-1.2499940188  <wf|H|wf>=-1.2501478403  error=1.2e+00
v=3      s=1e-09  gs_energy()/s=-1.2496378620  <wf|H|wf>=-1.2529409033  error=1.2e+00
v=3      s=1e-10  gs_energy()/s=-1.0621729389  <wf|H|wf>=-0.7647090540  error=1.4e+00
v=3      s=1e-11  gs_energy()/s=-1.2245464562  <wf|H|wf>=-1.2380607066  error=1.3e+00
v=3      s=1e-12  gs_energy()/s=-0.8956681314  <wf|H|wf>=-0.8048914864  error=1.6e+00
v=2      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=2      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=2      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=2      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.2e-14
v=2      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.6e-11
v=2      s=1e-06  gs_energy()/s=-2.4935771328  <wf|H|wf>=-2.4935771328  error=1.1e-09
v=2      s=6e-07  gs_energy()/s=-2.4935771329  <wf|H|wf>=-2.4935771329  error=9.7e-10
v=2      s=5e-07  gs_energy()/s=-2.4935771322  <wf|H|wf>=-2.4935771322  error=1.7e-09
v=2      s=3e-07  gs_energy()/s=-1.7469795607  <wf|H|wf>=-2.3972684049  error=7.5e-01
v=2      s=1e-07  gs_energy()/s=-1.2499999086  <wf|H|wf>=-1.2500056100  error=1.2e+00
v=2      s=1e-08  gs_energy()/s=-1.2499964036  <wf|H|wf>=-1.2478933732  error=1.2e+00
v=2      s=1e-09  gs_energy()/s=-1.2499123274  <wf|H|wf>=-1.2521103245  error=1.2e+00
v=2      s=1e-10  gs_energy()/s=-1.1723913522  <wf|H|wf>=-1.1503696192  error=1.3e+00
v=2      s=1e-11  gs_energy()/s=-1.1324993350  <wf|H|wf>=-0.8349127760  error=1.4e+00
v=2      s=1e-12  gs_energy()/s=-0.7142427764  <wf|H|wf>=-0.6931622453  error=1.8e+00

==== 03_small_units_mpo.before.out (parent tree; rows at s <= 1e-8 in A, C and D are the clean-threshold drop, part B is unmasked) ====
[run3] slot 1 acquired
A. vev(s*X)/s / vev(X) on the scale-1 ground state
  v=3  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.6667 5e-07:0.6667 4e-07:0.3333 3e-07:0.3333 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
  v=2  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974  XX=-0.831192  YY=-0.831192
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    YY            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
  v=python  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:-0.0000 1e-10:-0.0000 1e-12:-0.0000
B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly
  s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.5e-14
  s=1e-06  gs_energy()/s=-2.4935771324  <wf|H|wf>=-2.4935771324  error=1.5e-09
  s=6e-07  gs_energy()/s=-2.4935771284  <wf|H|wf>=-2.4935771284  error=5.5e-09
  s=5e-07  gs_energy()/s=-2.4935771326  <wf|H|wf>=-2.4935771326  error=1.3e-09
  s=3e-07  gs_energy()/s=-2.4935771318  <wf|H|wf>=-2.4935771318  error=2.1e-09
  s=1e-07  gs_energy()/s=-2.4935771150  <wf|H|wf>=-2.4935771150  error=1.9e-08
  s=1e-08  gs_energy()/s=-2.4935685066  <wf|H|wf>=-2.4935685066  error=8.6e-06
  s=1e-09  gs_energy()/s=-2.4929698055  <wf|H|wf>=-2.4929698055  error=6.1e-04
  s=1e-10  gs_energy()/s=-2.0173907156  <wf|H|wf>=-2.0173907156  error=4.8e-01
  s=1e-11  gs_energy()/s=-2.3754020970  <wf|H|wf>=-2.3754020970  error=1.2e-01
  s=1e-12  gs_energy()/s=-1.6143189666  <wf|H|wf>=-1.6143189666  error=8.8e-01
C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond
  ED at s=1: -2.3275944192
  v=3  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:2.0e-09 3e-07:4.8e-09 1e-07:6.0e-08 1e-08:2.3e+00 1e-10:2.3e+00 1e-11:2.3e+00 1e-12:2.3e+00
  v=2  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:2.6e-09 3e-07:4.5e-08 1e-07:3.9e-07 1e-08:2.3e+00 1e-10:2.3e+00 1e-11:2.3e+00 1e-12:2.3e+00
D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator
  stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e+00 1e-09:2.5e+00 1e-10:2.5e+00 1e-11:2.5e+00 1e-12:2.5e+00
  unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.5e+00 1e-09:2.5e+00 1e-10:2.5e+00 1e-11:2.5e+00 1e-12:2.5e+00

==== 03_small_units_mpo.hybrid.out (hybrid tree) ====
[run3] slot 1 acquired
A. vev(s*X)/s / vev(X) on the scale-1 ground state
  v=3  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.6667 5e-07:0.6667 4e-07:0.3333 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333 1e-12:0.3333
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:0.5000 5e-07:0.5000 4e-07:0.5000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
  v=2  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974  XX=-0.831192  YY=-0.831192
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333 1e-12:0.3333
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:0.5000 1e-07:0.5000 1e-08:0.5000 1e-10:0.5000 1e-12:0.5000
    XX            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    YY            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
  v=python  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly
  s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.5e-14
  s=1e-06  gs_energy()/s=-2.4935771337  <wf|H|wf>=-2.4935771337  error=2.3e-10
  s=6e-07  gs_energy()/s=-2.4935771331  <wf|H|wf>=-2.4935771331  error=8.3e-10
  s=5e-07  gs_energy()/s=-2.4935771262  <wf|H|wf>=-2.4935771262  error=7.7e-09
  s=3e-07  gs_energy()/s=-2.4935771285  <wf|H|wf>=-2.4935771285  error=5.3e-09
  s=1e-07  gs_energy()/s=-2.4935769678  <wf|H|wf>=-2.4935769678  error=1.7e-07
  s=1e-08  gs_energy()/s=-2.4935715745  <wf|H|wf>=-2.4935715745  error=5.6e-06
  s=1e-09  gs_energy()/s=-2.4924989275  <wf|H|wf>=-2.4924989275  error=1.1e-03
  s=1e-10  gs_energy()/s=-2.1444799424  <wf|H|wf>=-2.1444799424  error=3.5e-01
  s=1e-11  gs_energy()/s=-2.3526837274  <wf|H|wf>=-2.3526837274  error=1.4e-01
  s=1e-12  gs_energy()/s=-1.8176580712  <wf|H|wf>=-1.8176580712  error=6.8e-01
C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond
  ED at s=1: -2.3275944192
  v=3  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:8.7e-10 3e-07:2.4e-08 1e-07:1.1e-07 1e-08:3.6e-05 1e-10:6.4e-01 1e-11:4.4e-01 1e-12:1.1e+00
  v=2  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:3.6e-09 3e-07:3.4e-09 1e-07:3.7e-07 1e-08:1.8e-05 1e-10:4.6e-01 1e-11:2.1e-01 1e-12:1.3e+00
D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator
  stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e-09 1e-09:1.2e-06 1e-10:2.8e-05 1e-11:1.0e-02 1e-12:2.6e+00
  unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.4e-14 1e-09:3.9e-14 1e-10:2.9e-12 1e-11:2.0e-12 1e-12:4.4e-08

==== 11_small_units_consumers.before.out (parent tree; s <= 1e-8 masked by the clean-threshold drop, ED reads 0 too) ====
[run3] slot 1 acquired
v=3 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=1.3e-15  gs_energy_generalized/s=-2.5748369148  err=8.9e-16
v=3 s=1e-08  KPM max|DMRG-ED|/peak=raised RuntimeError: Error condition in diagHermitian  get_excited(n=4)/s=0.00000000 0.00000000 0.00000000 0.00000000  max|DMRG-ED|=0.0e+00  gs_energy_generalized/s=0.0000000000  err=2.6e+00
v=3 s=1e-10  KPM max|DMRG-ED|/peak=raised RuntimeError: Error condition in diagHermitian  get_excited(n=4)/s=0.00000000 0.00000000 0.00000000 0.00000000  max|DMRG-ED|=0.0e+00  gs_energy_generalized/s=0.0000000000  err=2.6e+00
v=2 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=4.4e-16
v=2 s=1e-08  KPM max|DMRG-ED|/peak=nan  get_excited(n=4)/s=0.00000000 0.00000000 0.00000000 0.00000000  max|DMRG-ED|=0.0e+00
v=2 s=1e-10  KPM max|DMRG-ED|/peak=nan  get_excited(n=4)/s=0.00000000 0.00000000 0.00000000 0.00000000  max|DMRG-ED|=0.0e+00

==== 11_small_units_consumers.hybrid.out (hybrid tree) ====
[run3] slot 1 acquired
v=3 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=8.9e-16  gs_energy_generalized/s=-2.5748369148  err=0.0e+00
v=3 s=1e-08  KPM max|DMRG-ED|/peak=2.8e+00  get_excited(n=4)/s=-1.24999797 -1.24999436 -0.74999323 -0.74998836  max|DMRG-ED|=1.3e+00  gs_energy_generalized/s=-1.3888862460  err=1.2e+00
v=3 s=1e-10  KPM max|DMRG-ED|/peak=2.6e+00  get_excited(n=4)/s=-1.20582272 -0.94719711 -0.72347689 -0.70179304  max|DMRG-ED|=1.3e+00  gs_energy_generalized/s=-1.2476144765  err=1.3e+00
v=2 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=1.3e-15
v=2 s=1e-08  KPM max|DMRG-ED|/peak=1.1e+00  get_excited(n=4)/s=-1.24999963 -1.24999927 -0.74999962 -0.74999948  max|DMRG-ED|=1.3e+00
v=2 s=1e-10  KPM max|DMRG-ED|/peak=1.1e+00  get_excited(n=4)/s=-1.22515464 -1.10684999 -0.72157162 -0.68398136  max|DMRG-ED|=1.3e+00

==== 13_left_open.hybrid.out (hybrid tree) ====
[run3] slot 1 acquired
v=3 s=1e+00  NH-DMRG E0/s=-2.46104102+0.00000000j  |DMRG-ED|=2.4e-14  evolution max|DMRG-ED|=4.2e-13
v=3 s=1e-08  NH-DMRG E0/s=-1.20000000-0.00000000j  |DMRG-ED|=1.3e+00  evolution max|DMRG-ED|=4.0e-02
v=2 s=1e+00  NH-DMRG E0/s=-2.46104102+0.00000000j  |DMRG-ED|=2.4e-14  evolution max|DMRG-ED|=5.1e-06
v=2 s=1e-08  NH-DMRG E0/s=0.30000000-0.00000000j  |DMRG-ED|=2.8e+00  evolution max|DMRG-ED|=4.7e+49

==== 10_ordinary_scale.before.out (parent tree, run 1) ====
[run3] slot 1 acquired
v=3 J=1 Heisenberg  gs_energy=-2.493577133887928  err=1.8e-15
v=3 J=1 Heisenberg  vev(Sz0Sz1)=-0.22248720747516915  err=1.8e-14  vev(XY23)=-0.40504827583360703  err=9.0e-13
v=3 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=3 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229536936 -0.0229296109 -0.0228895422 -0.0228335911 -0.0227619017 -0.0226746584 -0.0225720849 -0.0224544431 -0.0223220321  max|DMRG-ED|=4.5e-13
v=3 J=0.5 Heisenberg  gs_energy=-1.2467885669439627  ED=np.float64(-1.2467885669439631)  err=4.4e-16
v=3 J=1 XX chain  gs_energy=-1.7469796037174667  ED=np.float64(-1.7469796037174663)  err=4.4e-16
v=2 J=1 Heisenberg  gs_energy=-2.493577133887927  err=8.9e-16
v=2 J=1 Heisenberg  vev(Sz0Sz1)=-0.22248720747588013  err=6.9e-13  vev(XY23)=-0.405048275832528  err=1.8e-13
v=2 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=2 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229550002 -0.0229320335 -0.0228928896 -0.0228376720 -0.0227665258 -0.0226796368 -0.0225772306 -0.0224595722 -0.0223269642  max|DMRG-ED|=5.1e-06
v=2 J=0.5 Heisenberg  gs_energy=-1.2467885669439636  ED=np.float64(-1.2467885669439631)  err=4.4e-16
v=2 J=1 XX chain  gs_energy=-1.746979603717468  ED=np.float64(-1.7469796037174663)  err=1.8e-15

==== 10_ordinary_scale.before_run2.out (parent tree, run 2 of the same script: the last digits move, so neither backend is deterministic run to run) ====
[run3] slot 1 acquired
v=3 J=1 Heisenberg  gs_energy=-2.493577133887927  err=8.9e-16
v=3 J=1 Heisenberg  vev(Sz0Sz1)=-0.2224872074746424  err=5.4e-13  vev(XY23)=-0.40504827583244746  err=2.6e-13
v=3 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=3 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229536936 -0.0229296109 -0.0228895422 -0.0228335911 -0.0227619017 -0.0226746584 -0.0225720849 -0.0224544431 -0.0223220321  max|DMRG-ED|=3.9e-13
v=3 J=0.5 Heisenberg  gs_energy=-1.246788566943962  ED=np.float64(-1.2467885669439631)  err=1.1e-15
v=3 J=1 XX chain  gs_energy=-1.7469796037174676  ED=np.float64(-1.7469796037174663)  err=1.3e-15
v=2 J=1 Heisenberg  gs_energy=-2.4935771338879293  err=3.1e-15
v=2 J=1 Heisenberg  vev(Sz0Sz1)=-0.22248720747571532  err=5.3e-13  vev(XY23)=-0.40504827583215053  err=5.6e-13
v=2 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=2 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229550002 -0.0229320335 -0.0228928896 -0.0228376720 -0.0227665258 -0.0226796368 -0.0225772306 -0.0224595722 -0.0223269642  max|DMRG-ED|=5.1e-06
v=2 J=0.5 Heisenberg  gs_energy=-1.2467885669439627  ED=np.float64(-1.2467885669439631)  err=4.4e-16
v=2 J=1 XX chain  gs_energy=-1.7469796037174667  ED=np.float64(-1.7469796037174663)  err=4.4e-16

==== 14_truncated_below_one.before.out (parent tree) ====
[run3] slot 1 acquired
v=3 J=1 Heisenberg   gs_energy over 3 runs: -8.682274888756 -8.682274888756 -8.682274888756  mean=-8.682274888756  spread=5.5e-14
v=3 J=0.5 Heisenberg gs_energy over 3 runs: -4.341137444378 -4.341137444378 -4.341137444378  mean=-4.341137444378  spread=7.1e-15
v=3 J=1 XX chain     gs_energy over 3 runs: -6.190466206576 -6.190466210518 -6.190466155824  mean=-6.190466190973  spread=5.5e-08
v=2 J=1 Heisenberg   gs_energy over 3 runs: -8.682274888756 -8.682274888756 -8.682274888756  mean=-8.682274888756  spread=1.4e-13
v=2 J=0.5 Heisenberg gs_energy over 3 runs: -4.341137444378 -4.341137444378 -4.341137444378  mean=-4.341137444378  spread=4.4e-15
v=2 J=1 XX chain     gs_energy over 3 runs: -6.190466036921 -6.190466174074 -6.190466096045  mean=-6.190466102346  spread=1.4e-07

==== 15_xx_spread.before.out (parent tree) ====
[run3] slot 1 acquired
exact -6.190744999827
v=3 J=1 XX chain, 10 runs: mean=-6.190466101716  std=6.8e-08  min=-6.190466222403  max=-6.190466030838  mean-exact=2.79e-04

==== repair/20_mixed_scale.before.out (parent tree; the s <= 1e-8 rows are the clean-threshold drop, and B raises or returns nan at s <= 1e-12 through the older mechanisms) ====
[run3] slot 1 acquired
A(i). no eigensolver, ratio to vev(H) on the scale-1 ground state
  v=3      offset 1 + s*H   1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000
  v=3      field Sz0 + s*H  1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000
  v=2      offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.0000 1e-10:0.0000
  v=2      field Sz0 + s*H  1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:-0.0000 1e-10:-0.0000
  v=python offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.0000 1e-08:-0.0000 1e-10:-0.0000
  v=python field Sz0 + s*H  1e-06:1.0000 3e-07:1.0000 1e-07:0.0892 1e-08:-0.0000 1e-10:-0.0000
A(ii). gs_energy, fresh chain per row, against ED on the same chain
  v=3 s=1e-07  s*H+1    (E-1)/s=-1.2500000057  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-08  s*H+1    (E-1)/s=-0.0000000666  ED=0.0000000000  |DMRG-ED|=6.7e-08  |DMRG-exact|=2.5e+00
  v=3 s=1e-10  s*H+1    (E-1)/s=-0.0000122125  ED=0.0000000000  |DMRG-ED|=1.2e-05  |DMRG-exact|=2.5e+00
  v=3 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.2499999913  ED=-2.0925528521  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-08  s*H+Sz0  (E+0.5)/s=0.0000000500  ED=0.0000000000  |DMRG-ED|/s=5.0e-08
  v=3 s=1e-10  s*H+Sz0  (E+0.5)/s=-0.0000011102  ED=0.0000000000  |DMRG-ED|/s=1.1e-06
  v=2 s=1e-07  s*H+1    (E-1)/s=-1.2499999869  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-08  s*H+1    (E-1)/s=0.0000000000  ED=0.0000000000  |DMRG-ED|=0.0e+00  |DMRG-exact|=2.5e+00
  v=2 s=1e-10  s*H+1    (E-1)/s=-0.0000033307  ED=0.0000000000  |DMRG-ED|=3.3e-06  |DMRG-exact|=2.5e+00
  v=2 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.0000000039  ED=-2.0925528521  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-08  s*H+Sz0  (E+0.5)/s=-0.0000000666  ED=0.0000000000  |DMRG-ED|/s=6.7e-08
  v=2 s=1e-10  s*H+Sz0  (E+0.5)/s=-0.0000022204  ED=0.0000000000  |DMRG-ED|/s=2.2e-06
B. KPM C[Sz0,Sz0] of s*H, max|DMRG-ED|/ED peak
  v=3  1.0e+00:4.1e-04 1.0e-12:raised RuntimeError 1.8e-14:raised RuntimeError 1.5e-14:raised RuntimeError 1.0e-14:raised RuntimeError 1.0e-16:raised RuntimeError 1.0e-20:raised RuntimeError
  v=2  1.0e+00:4.1e-04 1.0e-12:nan 1.8e-14:nan 1.5e-14:nan 1.0e-14:nan 1.0e-16:nan 1.0e-20:nan
C. weak link J'*S2.S3 in a J=1 chain, [vev(Hs+J'*B23)-vev(Hs)]/J'/vev(B23)
  v=3      1e-05:1.0000 1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:-0.0000
  v=2      1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:-0.0000
  v=python 1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:-0.0000 1e-08:-0.0000

==== repair/20_mixed_scale.hybrid.out (hybrid tree) ====
[run3] slot 1 acquired
A(i). no eigensolver, ratio to vev(H) on the scale-1 ground state
  v=3      offset 1 + s*H   1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=3      field Sz0 + s*H  1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      field Sz0 + s*H  1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=python offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.0000 1e-08:-0.0000 1e-10:0.0000
  v=python field Sz0 + s*H  1e-06:1.0000 3e-07:1.0000 1e-07:0.0892 1e-08:0.0892 1e-10:0.0892
A(ii). gs_energy, fresh chain per row, against ED on the same chain
  v=3 s=1e-07  s*H+1    (E-1)/s=-1.2499999857  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-08  s*H+1    (E-1)/s=-1.2499986823  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-10  s*H+1    (E-1)/s=-1.1975676006  ED=-2.4935808973  |DMRG-ED|=1.3e+00  |DMRG-exact|=1.3e+00
  v=3 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.2499999924  ED=-2.0925528521  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-08  s*H+Sz0  (E+0.5)/s=-1.2499999813  ED=-2.0925528088  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-10  s*H+Sz0  (E+0.5)/s=-1.2435541485  ED=-2.0925539079  |DMRG-ED|/s=8.5e-01
  v=2 s=1e-07  s*H+1    (E-1)/s=-1.2499999602  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-08  s*H+1    (E-1)/s=-1.2499997704  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-10  s*H+1    (E-1)/s=-1.1864076388  ED=-2.4935808973  |DMRG-ED|=1.3e+00  |DMRG-exact|=1.3e+00
  v=2 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.0000000050  ED=-2.0925528521  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-08  s*H+Sz0  (E+0.5)/s=-0.9999998607  ED=-2.0925528088  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-10  s*H+Sz0  (E+0.5)/s=-0.9938838641  ED=-2.0925539079  |DMRG-ED|/s=1.1e+00
B. KPM C[Sz0,Sz0] of s*H, max|DMRG-ED|/ED peak
  v=3  1.0e+00:4.1e-04 1.0e-12:raised RuntimeError 1.8e-14:2.0e+00 1.5e-14:raised RuntimeError 1.0e-14:raised RuntimeError 1.0e-16:raised RuntimeError 1.0e-20:raised RuntimeError
  v=2  1.0e+00:4.1e-04 1.0e-12:raised RuntimeError 1.8e-14:raised RuntimeError 1.5e-14:7.8e-01 1.0e-14:1.3e+00 1.0e-16:nan 1.0e-20:nan
C. weak link J'*S2.S3 in a J=1 chain, [vev(Hs+J'*B23)-vev(Hs)]/J'/vev(B23)
  v=3      1e-05:1.0000 1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333
  v=2      1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333
  v=python 1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:-0.0000 1e-08:-0.0000

==== repair/20_mixed_scale.firstpass.out (first-pass tree: the repair's own before) ====
[run3] slot 1 acquired
A(i). no eigensolver, ratio to vev(H) on the scale-1 ground state
  v=3      offset 1 + s*H   1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=3      field Sz0 + s*H  1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      field Sz0 + s*H  1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=python offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.0000 1e-08:-0.0000 1e-10:0.0000
  v=python field Sz0 + s*H  1e-06:1.0000 3e-07:1.0000 1e-07:0.0892 1e-08:0.0892 1e-10:0.0892
A(ii). gs_energy, fresh chain per row, against ED on the same chain
  v=3 s=1e-07  s*H+1    (E-1)/s=-1.2499999713  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-08  s*H+1    (E-1)/s=-1.2499996038  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-10  s*H+1    (E-1)/s=-1.2212630907  ED=-2.4935808973  |DMRG-ED|=1.3e+00  |DMRG-exact|=1.3e+00
  v=3 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.2499999880  ED=-2.0925528521  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-08  s*H+Sz0  (E+0.5)/s=-1.2499991486  ED=-2.0925528088  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-10  s*H+Sz0  (E+0.5)/s=-1.1075562689  ED=-2.0925539079  |DMRG-ED|/s=9.8e-01
  v=2 s=1e-07  s*H+1    (E-1)/s=-1.2499999991  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-08  s*H+1    (E-1)/s=-1.2499997148  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-10  s*H+1    (E-1)/s=-1.1903522612  ED=-2.4935808973  |DMRG-ED|=1.3e+00  |DMRG-exact|=1.3e+00
  v=2 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.0000000006  ED=-2.0925528521  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-08  s*H+Sz0  (E+0.5)/s=-0.9999999828  ED=-2.0925528088  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-10  s*H+Sz0  (E+0.5)/s=-0.9903067255  ED=-2.0925539079  |DMRG-ED|/s=1.1e+00
B. KPM C[Sz0,Sz0] of s*H, max|DMRG-ED|/ED peak
  v=3  1.0e+00:4.1e-04 1.0e-12:4.1e-04 1.8e-14:4.1e-04 1.5e-14:1.2e+00 1.0e-14:1.2e+00 1.0e-16:1.2e+00 1.0e-20:1.2e+00
  v=2  1.0e+00:4.1e-04 1.0e-12:4.1e-04 1.8e-14:4.1e-04 1.5e-14:1.2e+00 1.0e-14:1.2e+00 1.0e-16:1.2e+00 1.0e-20:1.2e+00
C. weak link J'*S2.S3 in a J=1 chain, [vev(Hs+J'*B23)-vev(Hs)]/J'/vev(B23)
  v=3      1e-05:1.0000 1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333
  v=2      1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333
  v=python 1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:-0.0000 1e-08:-0.0000

==== repair/22_tests.hybrid.out (the current test file, copied alone to repair/hybrid_check/, against the hybrid tree; supersedes the first pass's 12_tests_on_hybrid.out, which listed the first 25 of these) ====
E         comparison failed
E         comparison failed
E         comparison failed
E                 comparison failed
E                 comparison failed
E         comparison failed
E         comparison failed
E         comparison failed
E         comparison failed
E         comparison failed
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v3-1e-06]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v3-6e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v3-3e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v3-1e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v3-1e-08]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v2-1e-06]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v2-6e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v2-3e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v2-1e-07]
FAILED test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend[v2-1e-08]
FAILED test_audit_2026_09_25_scale.py::test_small_units_ground_state_on_the_compiled_backends[1e-10-v3]
FAILED test_audit_2026_09_25_scale.py::test_small_units_ground_state_on_the_compiled_backends[1e-10-v2]
FAILED test_audit_2026_09_25_scale.py::test_small_units_ground_state_on_the_compiled_backends[1e-12-v3]
FAILED test_audit_2026_09_25_scale.py::test_small_units_ground_state_on_the_compiled_backends[1e-12-v2]
FAILED test_audit_2026_09_25_scale.py::test_the_mpo_of_a_small_unit_operator_keeps_every_channel[v3]
FAILED test_audit_2026_09_25_scale.py::test_the_mpo_of_a_small_unit_operator_keeps_every_channel[v2]
FAILED test_audit_2026_09_25_scale.py::test_the_local_solver_is_scale_invariant[1e-08-v3]
FAILED test_audit_2026_09_25_scale.py::test_the_local_solver_is_scale_invariant[1e-08-v2]
FAILED test_audit_2026_09_25_scale.py::test_the_local_solver_is_scale_invariant[1e-10-v3]
FAILED test_audit_2026_09_25_scale.py::test_the_local_solver_is_scale_invariant[1e-10-v2]
FAILED test_audit_2026_09_25_scale.py::test_excited_states_in_small_units_match_ed[v3]
FAILED test_audit_2026_09_25_scale.py::test_excited_states_in_small_units_match_ed[v2]
FAILED test_audit_2026_09_25_scale.py::test_kpm_correlator_in_small_units_matches_ed[v3]
FAILED test_audit_2026_09_25_scale.py::test_kpm_correlator_in_small_units_matches_ed[v2]
FAILED test_audit_2026_09_25_scale.py::test_generalized_ground_state_in_small_units
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1.5e-14-v3]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1.5e-14-v2]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1e-20-v3]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1e-20-v2]
29 failed, 41 passed, 14 xfailed in 18.19s

==== repair/22_tests.firstpass.out (the same file against the first-pass tree; the saved tail opens with six lines of pytest's assertion introspection of the KPM-floor failure, then) ====
test_audit_2026_09_25_scale.py:358: AssertionError
=========================== short test summary info ============================
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1.5e-14-v3]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1.5e-14-v2]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1e-20-v3]
FAILED test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor[1e-20-v2]
4 failed, 66 passed, 14 xfailed in 14.91s
````

Observed, after:

````
==== repair/02_small_units_onset.after.out (working tree, both extensions rebuilt in the repair pass) ====
[run3] slot 1 acquired
v=3      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=3      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=3      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=3      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=3      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=3      s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=3      s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=3      s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-08  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-09  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=3      s=1e-10  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=3      s=1e-11  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=3      s=1e-12  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=2      s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=2      s=1e-01  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=1e-02  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=1e-03  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=2      s=1e-04  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=1e-05  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15
v=2      s=1e-06  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=4.4e-16
v=2      s=6e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.2e-15
v=2      s=5e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=2      s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.3e-15
v=2      s=1e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=3.1e-15
v=2      s=1e-08  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=1.8e-15
v=2      s=1e-09  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=0.0e+00
v=2      s=1e-10  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=2      s=1e-11  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=2      s=1e-12  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.7e-15

==== repair/03_small_units_mpo.after.out (working tree; part B, the Hamiltonian handed in as an already-built MPO, is unchanged by design) ====
[run3] slot 1 acquired
A. vev(s*X)/s / vev(X) on the scale-1 ground state
  v=3  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
  v=2  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974  XX=-0.831192  YY=-0.831192
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    YY            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
  v=python  reference vev: H=-2.493577  XX+YY=-1.662385  ZZ=-0.831192  bond01 XX+YY=-0.444974
    H             1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    XX+YY         1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    ZZ            1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
    bond01 XX+YY  1e-03:1.0000 1e-05:1.0000 1e-06:1.0000 6e-07:1.0000 5e-07:1.0000 4e-07:1.0000 3e-07:1.0000 1e-07:1.0000 1e-08:1.0000 1e-10:1.0000 1e-12:1.0000
B. v3, set_hamiltonian(s*toMPO(H)): the MPO built at scale 1, scaled exactly
  s=1e+00  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=2.5e-14
  s=1e-06  gs_energy()/s=-2.4935771334  <wf|H|wf>=-2.4935771334  error=4.4e-10
  s=6e-07  gs_energy()/s=-2.4935771316  <wf|H|wf>=-2.4935771316  error=2.3e-09
  s=5e-07  gs_energy()/s=-2.4935771307  <wf|H|wf>=-2.4935771307  error=3.2e-09
  s=3e-07  gs_energy()/s=-2.4935771213  <wf|H|wf>=-2.4935771213  error=1.3e-08
  s=1e-07  gs_energy()/s=-2.4935770053  <wf|H|wf>=-2.4935770053  error=1.3e-07
  s=1e-08  gs_energy()/s=-2.4935673645  <wf|H|wf>=-2.4935673645  error=9.8e-06
  s=1e-09  gs_energy()/s=-2.4923920370  <wf|H|wf>=-2.4923920370  error=1.2e-03
  s=1e-10  gs_energy()/s=-1.9771596164  <wf|H|wf>=-1.9771596164  error=5.2e-01
  s=1e-11  gs_energy()/s=-1.9573771481  <wf|H|wf>=-1.9573771481  error=5.4e-01
  s=1e-12  gs_energy()/s=-1.5058864492  <wf|H|wf>=-1.5058864492  error=9.9e-01
C. transverse-field Ising s*(sum ZZ + 0.7 sum X): one channel per bond
  ED at s=1: -2.3275944192
  v=3  |gs_energy()/s - ED|  1e+00:3.9e-13 1e-06:3.9e-13 3e-07:3.9e-13 1e-07:3.9e-13 1e-08:3.9e-13 1e-10:3.9e-13 1e-11:3.9e-13 1e-12:3.9e-13
  v=2  |gs_energy()/s - ED|  1e+00:4.0e-13 1e-06:3.9e-13 3e-07:4.0e-13 1e-07:3.9e-13 1e-08:3.9e-13 1e-10:3.9e-13 1e-11:3.9e-13 1e-12:3.9e-13
D. python, stock Lanczos against the same Lanczos on the unit-normalized local operator
  stock            |gs_energy()/s - exact|  1e+00:2.4e-14 1e-08:2.5e-09 1e-09:1.2e-06 1e-10:2.8e-05 1e-11:1.0e-02 1e-12:2.6e+00
  unit-normalized  |gs_energy()/s - exact|  1e+00:2.6e-14 1e-08:2.4e-14 1e-09:3.9e-14 1e-10:2.9e-12 1e-11:2.0e-12 1e-12:4.4e-08

==== repair/11_small_units_consumers.after.out (working tree) ====
[run3] slot 1 acquired
v=3 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=8.9e-16  gs_energy_generalized/s=-2.5748369148  err=4.4e-16
v=3 s=1e-08  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=3.6e-15  gs_energy_generalized/s=-2.5748369148  err=8.9e-16
v=3 s=1e-10  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=2.7e-15  gs_energy_generalized/s=-2.5748369148  err=8.9e-16
v=2 s=1e+00  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=1.8e-15
v=2 s=1e-08  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=1.3e-15
v=2 s=1e-10  KPM max|DMRG-ED|/peak=4.1e-04  get_excited(n=4)/s=-2.49357713 -2.00199536 -2.00199536 -2.00199536  max|DMRG-ED|=1.8e-15

==== repair/13_left_open.after.out (working tree) ====
[run3] slot 1 acquired
v=3 s=1e+00  NH-DMRG E0/s=-2.46104102-0.00000000j  |DMRG-ED|=2.4e-14  evolution max|DMRG-ED|=9.4e-13
v=3 s=1e-08  NH-DMRG E0/s=-2.46104102+0.00000000j  |DMRG-ED|=3.1e-13  evolution max|DMRG-ED|=6.4e-04
v=2 s=1e+00  NH-DMRG E0/s=-2.46104102-0.00000000j  |DMRG-ED|=2.3e-14  evolution max|DMRG-ED|=5.1e-06
v=2 s=1e-08  NH-DMRG E0/s=-0.61437572-0.14372027j  |DMRG-ED|=1.9e+00  evolution max|DMRG-ED|=9.1e+68

==== repair/10_ordinary_scale.after.out (working tree) ====
[run3] slot 1 acquired
v=3 J=1 Heisenberg  gs_energy=-2.493577133887928  err=1.8e-15
v=3 J=1 Heisenberg  vev(Sz0Sz1)=-0.2224872074760726  err=8.9e-13  vev(XY23)=-0.40504827583369485  err=9.8e-13
v=3 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=3 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229536936 -0.0229296109 -0.0228895422 -0.0228335911 -0.0227619017 -0.0226746584 -0.0225720849 -0.0224544431 -0.0223220321  max|DMRG-ED|=4.9e-13
v=3 J=0.5 Heisenberg  gs_energy=-1.246788566943963  ED=np.float64(-1.2467885669439631)  err=2.2e-16
v=3 J=1 XX chain  gs_energy=-1.7469796037174676  ED=np.float64(-1.7469796037174663)  err=1.3e-15
v=2 J=1 Heisenberg  gs_energy=-2.4935771338879236  err=2.7e-15
v=2 J=1 Heisenberg  vev(Sz0Sz1)=-0.22248720747571038  err=5.2e-13  vev(XY23)=-0.4050482758336055  err=8.9e-13
v=2 J=1 Heisenberg  KPM Re C[Sz0,Sz0](w) = 0.0000068136 0.0674737183 0.2308127466 0.0318732807 0.0003865168 0.0000096216 0.0000041675  max|DMRG-ED|/peak=4.2e-04
v=2 J=1 Heisenberg  <Sz1Sz2>(t) from Sz0|GS>, t=0..0.45: -0.0229617281 -0.0229550002 -0.0229320335 -0.0228928896 -0.0228376720 -0.0227665258 -0.0226796368 -0.0225772306 -0.0224595722 -0.0223269642  max|DMRG-ED|=5.1e-06
v=2 J=0.5 Heisenberg  gs_energy=-1.2467885669439622  ED=np.float64(-1.2467885669439631)  err=8.9e-16
v=2 J=1 XX chain  gs_energy=-1.7469796037174685  ED=np.float64(-1.7469796037174663)  err=2.2e-15

==== repair/20_mixed_scale.after.out (working tree: B fixed at every s; A and C unchanged, the open bond-local defect) ====
[run3] slot 1 acquired
A(i). no eigensolver, ratio to vev(H) on the scale-1 ground state
  v=3      offset 1 + s*H   1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=3      field Sz0 + s*H  1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=2      field Sz0 + s*H  1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333 1e-10:0.3333
  v=python offset 1 + s*H   1e-06:1.0000 3e-07:0.6667 1e-07:0.0000 1e-08:-0.0000 1e-10:0.0000
  v=python field Sz0 + s*H  1e-06:1.0000 3e-07:1.0000 1e-07:0.0892 1e-08:0.0892 1e-10:0.0892
A(ii). gs_energy, fresh chain per row, against ED on the same chain
  v=3 s=1e-07  s*H+1    (E-1)/s=-1.2499999991  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-08  s*H+1    (E-1)/s=-1.2499997704  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=3 s=1e-10  s*H+1    (E-1)/s=-1.1819500934  ED=-2.4935808973  |DMRG-ED|=1.3e+00  |DMRG-exact|=1.3e+00
  v=3 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.2499999646  ED=-2.0925528521  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-08  s*H+Sz0  (E+0.5)/s=-1.2499998703  ED=-2.0925528088  |DMRG-ED|/s=8.4e-01
  v=3 s=1e-10  s*H+Sz0  (E+0.5)/s=-1.2156076146  ED=-2.0925539079  |DMRG-ED|/s=8.8e-01
  v=2 s=1e-07  s*H+1    (E-1)/s=-1.2499999935  ED=-2.4935771281  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-08  s*H+1    (E-1)/s=-1.2499995039  ED=-2.4935771226  |DMRG-ED|=1.2e+00  |DMRG-exact|=1.2e+00
  v=2 s=1e-10  s*H+1    (E-1)/s=-0.7667100288  ED=-2.4935808973  |DMRG-ED|=1.7e+00  |DMRG-exact|=1.7e+00
  v=2 s=1e-07  s*H+Sz0  (E+0.5)/s=-1.0000000061  ED=-2.0925528521  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-08  s*H+Sz0  (E+0.5)/s=-0.9999996164  ED=-2.0925528088  |DMRG-ED|/s=1.1e+00
  v=2 s=1e-10  s*H+Sz0  (E+0.5)/s=-0.9981915294  ED=-2.0925539079  |DMRG-ED|/s=1.1e+00
B. KPM C[Sz0,Sz0] of s*H, max|DMRG-ED|/ED peak
  v=3  1.0e+00:4.1e-04 1.0e-12:4.1e-04 1.8e-14:4.1e-04 1.5e-14:4.1e-04 1.0e-14:4.1e-04 1.0e-16:4.1e-04 1.0e-20:4.1e-04
  v=2  1.0e+00:4.1e-04 1.0e-12:4.1e-04 1.8e-14:4.1e-04 1.5e-14:4.1e-04 1.0e-14:4.1e-04 1.0e-16:4.1e-04 1.0e-20:4.1e-04
C. weak link J'*S2.S3 in a J=1 chain, [vev(Hs+J'*B23)-vev(Hs)]/J'/vev(B23)
  v=3      1e-05:1.0000 1e-06:1.0000 3e-07:0.3333 1e-07:0.3333 1e-08:0.3333
  v=2      1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:0.3333 1e-08:0.3333
  v=python 1e-05:1.0000 1e-06:1.0000 3e-07:0.6667 1e-07:-0.0000 1e-08:-0.0000

==== 14_truncated_below_one.after.out (first-pass build; the repair changed no code path 14 or 15 reaches, since neither runs KPM or CVM) ====
[run3] slot 1 acquired
v=3 J=1 Heisenberg   gs_energy over 3 runs: -8.682274888756 -8.682274888756 -8.682274888756  mean=-8.682274888756  spread=4.1e-13
v=3 J=0.5 Heisenberg gs_energy over 3 runs: -4.341137444378 -4.341137444378 -4.341137444378  mean=-4.341137444378  spread=8.0e-14
v=3 J=1 XX chain     gs_energy over 3 runs: -6.190466028219 -6.190466092935 -6.190466057017  mean=-6.190466059390  spread=6.5e-08
v=2 J=1 Heisenberg   gs_energy over 3 runs: -8.682274888756 -8.682274888756 -8.682274888756  mean=-8.682274888756  spread=5.7e-14
v=2 J=0.5 Heisenberg gs_energy over 3 runs: -4.341137444378 -4.341137444378 -4.341137444378  mean=-4.341137444378  spread=3.4e-14
v=2 J=1 XX chain     gs_energy over 3 runs: -6.190466074618 -6.190466050839 -6.190466143214  mean=-6.190466089557  spread=9.2e-08

==== 15_xx_spread.after.out (first-pass build, as 14) ====
[run3] slot 1 acquired
exact -6.190744999827
v=3 J=1 XX chain, 10 runs: mean=-6.190466092079  std=7.3e-08  min=-6.190466222400  max=-6.190466028516  mean-exact=2.79e-04

==== repair/21_tests.after.out (tests/test_audit_2026_09_25_scale.py on the working tree) ====
[run3] slot 1 acquired
................................................................x......x [ 85%]
xxxxxxxxxxxx                                                             [100%]
70 passed, 14 xfailed in 16.92s
````

**NUMBERS CHANGE**: 6-site S=1/2 Heisenberg s*H, maxm=30, nsweeps=10, fresh chain per s, before = working-tree Python with the parent extensions: gs_energy()/s v3 s=6e-7 -1.2715204438 -> -2.4935771339, 5e-7 -1.4010959621 -> -2.4935771339, 3e-7..1e-9 about -1.25 -> -2.4935771339, 1e-10 -1.0621729389, 1e-11 -1.2245464562, 1e-12 -0.8956681314 -> -2.4935771339; v2 s=3e-7 -1.7469795607, 1e-7..1e-9 about -1.25, 1e-10 -1.1723913522, 1e-11 -1.1324993350, 1e-12 -0.7142427764 -> -2.4935771339 (errors <= 5.8e-15 after over three runs); s=1e-6 error v3 2.0e-10 (7.8e-10 on the parent, 2.1e-9 in the review's hybrid run) -> <= 3.6e-15, v2 1.1e-9 (8.6e-10, 1.6e-9) -> <= 4.4e-16. vev(s*H)/s / vev(H) on the scale-1 state: v3 0.6667 (6e-7, 5e-7) and 0.3333 (<= 4e-7), v2 0.6667 (3e-7) and 0.3333 (<= 1e-7) -> 1.0000; XX+YY 0.5000 -> 1.0000. get_excited(n=4)/s at s=1e-8: v3 -1.24999797 -1.24999436 -0.74999323 -0.74998836, v2 -1.24999963 -1.24999927 -0.74999962 -0.74999948 -> -2.49357713 -2.00199536 -2.00199536 -2.00199536. KPM C[Sz0,Sz0] (es*s, delta=0.2*s), max|DMRG-ED|/ED peak: s=1e-8 v3 2.8, v2 1.1 -> 4.1e-4; s=1.5e-14, 1e-14, 1e-16, 1e-20 1.2 on both (first-pass build; RuntimeError, nan or 0.78..2.0 on the hybrid) -> 4.1e-4. v3 gs_energy_generalized/s, A = 1 + 0.2*Sz0, s=1e-8: -1.3888862460 -> -2.5748369148. Ising s*(sum ZZ + 0.7 sum X), |gs_energy()/s - ED|: s=1e-8 v3 3.6e-5, v2 1.8e-5 -> 3.9e-13; s=1e-10 v3 0.64, v2 0.46 -> 3.9e-13. v3 NH-DMRG on s*(Heisenberg + 0.3j*Sz0), s=1e-8: -1.20000000 -> -2.46104102 (3.1e-13 to 7.1e-12 off). v3 TDVP <Sz1Sz2>(t) at s=1e-8, dt=0.05/s: max|DMRG-ED| 4.0e-2 -> 6.4e-4. Ordinary scale: 6-site J=1 Heisenberg, J=0.5 Heisenberg, J=1 XX on v2/v3 against ED unchanged within run-to-run spread, and the J=1 KPM row identical to ten printed digits (10, parent, first pass and repair); 20-site maxm=10: J=1 and J=0.5 Heisenberg identical to 12 digits, v3 J=1 XX ten-run mean -6.190466101716 -> -6.190466092079 against std 7e-8 (15). Unchanged, still wrong: [vev(s*H+1)-1]/s and [vev(s*H+Sz0)-vev(Sz0)]/s over vev(H) at s=1e-7 (0.3333 on v3 and v2, 0.0000 and 0.0892 on python), the weak link J'=1e-7 (0.3333 on v3 and v2, 0.0000 on python), gs_energy of s*H+1 ((E-1)/s about -1.25) and s*H+Sz0 ((E+0.5)/s -1.25 on v3, -1.00 on v2, against ED -2.0925528521) at s=1e-7.

**Tests**: `tests/test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_invariant_on_every_mps_backend`; `tests/test_audit_2026_09_25_scale.py::test_small_units_ground_state_on_the_compiled_backends`; `tests/test_audit_2026_09_25_scale.py::test_the_mpo_of_a_small_unit_operator_keeps_every_channel`; `tests/test_audit_2026_09_25_scale.py::test_the_local_solver_is_scale_invariant`; `tests/test_audit_2026_09_25_scale.py::test_excited_states_in_small_units_match_ed`; `tests/test_audit_2026_09_25_scale.py::test_kpm_correlator_in_small_units_matches_ed`; `tests/test_audit_2026_09_25_scale.py::test_kpm_shift_survives_below_the_absolute_coefficient_floor`; `tests/test_audit_2026_09_25_scale.py::test_generalized_ground_state_in_small_units`; `tests/test_audit_2026_09_25_scale.py::test_a_bond_far_below_the_largest_coefficient_still_open (strict xfail: v3, v2, python x offset, field, weak-link)`; `tests/test_audit_2026_09_25_scale.py::test_small_units_next_to_an_order_one_term_still_open (strict xfail: v3, v2 x offset, field)`; `tests/test_audit_2026_09_25_scale.py::test_small_units_ground_state_still_open (strict xfail, python-1e-10 only)`

**Reviewer (CONFIRMED, fix INCOMPLETE)**: I reproduced it on the hybrid tree (working-tree Python, parent extensions, cmp-identical to parent/src's .so, and every changed .py cmp-identical to the working tree at review time), which is the right before, since the parent tree masks every row at s <= 1e-8 as 0 for the clean-threshold reason. 02 with `3 2` gives the Neel -1.25 on v3 from s=3e-7 and on v2 from s=1e-7, -1.747 on v2 at 3e-7, and the gradual davidson loss at s=1e-6 (2.1e-9 on v3, 1.6e-9 on v2), so both mechanisms stand as stated. Two sub-claims are struck or corrected. First, the explanation says of the two thresholds that "neither is reachable through Args": that is false for the MPO half, since svdMPO reads its truncation cutoff as `args.getReal("Cutoff",1E-13)` (mpscpp3/ITensor/itensor/mps/autompo.cc:1176, mpscpp2 :1157); only the matrix-element `isZero(coef,eps)` with `eps=1E-14` (v3 :1172/:1287) and davidson's `qnrm < 1E-10` are hardcoded. This does not undermine the fix, since scaling covers the hardcoded isZero as well, but the sentence should not go into documentation.md as proposed in docs_needed item 2. Second, the s=1e-6 range quoted in the heading (2.0e-10 to 1.1e-9) is run-to-run: my hybrid run measured 2.1e-9 on v3, so the range is 2e-10 to 2e-9. It is not intended behaviour, not documented as a known issue, and the anchor is independent (ED and the scale-1 run of the same backend).

What holds: on the working tree (both .so newer than the headers, `make -n pybind` reports nothing to do in both), 02 is exact to at most 5.8e-15 at every s from 1 to 1e-12 on v3 and v2. So are spinless fermions and a 4-site Hubbard at s=1e-8 and 1e-10, v3 sector mode Sz=0 (ground state and three excited levels), one chain reused across s=1, 1e-8, 1, 1e-10, s down to 1e-40 for gs_energy and get_excited, and a complex-phase Hubbard (complex MPO through the scaled path on v2 and v3). By reading, v2 has one `H_ =` (set_hamiltonian) and v3 one (set_hamiltonian_mpo, reset to 1 there), so no rebuild path escapes `hscale_up_`. At a largest coefficient of 1 or more every changed call reduces to the old one (`solver_hamiltonian` returns H, `weight*1.0`, `(-1.0)*H_`, `/1.0`, `toMPO(ampo,Args::global())` being what `MPO(ampo)` and `toMPO(ampo)` already did), so I accept the byte-identity argument. The tests pin the property: 66 passed and 1 xfailed on the working tree, and the same file against the hybrid gives 25 failed, 41 passed, 1 xfailed, which is exactly the fix agent's count. Two things fail, which is why the fix is INCOMPLETE rather than HOLDS. (1) The gate reads the wrong quantity. `unit_scale_up(max|coef|)` takes the maximum over the whole term list, but svdMPO's absolute cutoff acts only on the bond-crossing channels, so an O(1) term that crosses no bond, whether an energy offset `+1.0` (sent as `('Id',site)` with coefficient 1) or a field on one site, keeps the whole Hamiltonian on the unscaled path, and a sub-cliff exchange next to it is truncated exactly as before. On the fixed tree s*Heis+1.0 gives (E-1)/s = -1.2499999913 (v3) and -1.2499999902 (v2) against ED -2.4935771281 at s=1e-7, and s*Heis+Sz0 gives -1.25 (v3) and -1.00 (v2) against ED -2.0925528521. The vev ratio with c=1 reads 0.3333 on both the hybrid and the fixed tree, while the c=0 column went from 0.3333 to 1.0000, so this is left unchanged by the fix, not a regression. It needs a six to seven decade hierarchy inside one Hamiltonian with the small side below about 6e-7 (v3) or 3e-7 (v2), and the most reachable construction is s*H plus an O(1) constant. So the status line "FIXED on v2 and v3 for every Hamiltonian set from terms, which is every set_hamiltonian(MultiOperator)" is struck. The narrowed claim: fixed for every term-built Hamiltonian whose largest coefficient sits on a multi-site term (every uniformly scaled model), not for one where a single-site or constant term of order 1 accompanies exchange below the cliff. (2) The single-term shift AutoMPOs the fix deliberately left on plain `toMPO` (v3 scaled_hamiltonian :11786, scaled_hamiltonian_gs_anchored :12072, cvm's `z*Id` :2039, and the rest listed under evidence) still meet `isZero(coef,1E-14)`. The KPM correlator is right at s=1.8e-14 and 1.2 of the ED peak off at 1.5e-14 on both backends, which brackets the band-centre shift 0.62*s crossing 1e-14. So "every KPM spectrum ... down to where double precision runs out" in the item's expected line is struck below s of about 1.5e-14 (left_open item 5 was "not measured", and it is measured now); it is practically irrelevant but cheap to close. Suggested changes, for the fix agent to measure rather than adopt on my reading. (a) For the MPO build, gate on the largest coefficient of the terms spanning more than one site, since those are the entries the bond SVD decomposes, falling back to the overall largest when there is none; an Id-only offset could alternatively be stripped and added back analytically. One consequence has to be measured, not read: for s*Heis+Sz0 the onsite term becomes about 2^22, and davidson's absolute ErrGoal then never fires, which the locate stage's 04 diagnostic suggests is harmless (stock s=1 already exits every local solve on maxiter, 100 of 100). (b) Route the single-term shift/Id AutoMPOs through `to_mpo_unit` too, which is exact for one term. (c) Add a strict xfail on s*Heis+Sz0 or s*Heis+1.0 at s=1e-7 on v2 and v3 if (a) is not done, the same pattern the file uses for python-1e-10. NUMBERS CHANGE is accurate for what it covers; add the non-change that mixed-scale Hamiltonians (an O(1) single-site or constant term next to exchange below the cliff) are unchanged before and after, still wrong, and add that limit to docs_needed item 1 next to the four surviving routes.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
Review folder <scratch>/review/scale-cpp/.

02.hybrid.out (hybrid tree, `02_small_units_onset.py 3 2`, cmp-identical to the locate stage's script):
```
v=3      s=1e-06  gs_energy()/s=-2.4935771318  <wf|H|wf>=-2.4935771318  error=2.1e-09
v=3      s=6e-07  gs_energy()/s=-1.1232830787  <wf|H|wf>=-1.5749902442  error=1.4e+00
v=3      s=3e-07  gs_energy()/s=-1.2499999995  <wf|H|wf>=-1.2499966105  error=1.2e+00
v=3      s=1e-10  gs_energy()/s=-1.2254632158  <wf|H|wf>=-1.1250622618  error=1.3e+00
v=2      s=1e-06  gs_energy()/s=-2.4935771323  <wf|H|wf>=-2.4935771323  error=1.6e-09
v=2      s=3e-07  gs_energy()/s=-1.7469795991  <wf|H|wf>=-2.3972559944  error=7.5e-01
v=2      s=1e-07  gs_energy()/s=-1.2499998420  <wf|H|wf>=-1.2500038983  error=1.2e+00
```
02.after.out (working tree): every row -2.4935771339, errors 0.0e+00 to 5.8e-15 on both backends, for example
```
v=3      s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
v=3      s=1e-12  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=5.8e-15
v=2      s=3e-07  gs_energy()/s=-2.4935771339  <wf|H|wf>=-2.4935771339  error=8.9e-16
```
tests.after.out: `66 passed, 1 xfailed in 22.96s`. tests.hybrid.out: the same 25 FAILED ids as the fix agent's 12_tests_on_hybrid.out, `25 failed, 41 passed, 1 xfailed in 22.69s`.

Gate hole, pA_gate.after.out (fixed tree):
```
v=3 (i) [vev(c+s*H)-c]/s/vev(H): s=1e-06,c=0:1.0000 s=1e-06,c=1:1.0000 s=3e-07,c=0:1.0000 s=3e-07,c=1:0.3333 s=1e-07,c=0:1.0000 s=1e-07,c=1:0.3333 s=1e-08,c=0:1.0000 s=1e-08,c=1:0.3333
v=3 (ii) s=1e-07 c=1  (gs_energy-c)/s=-1.2499999913  ED=-2.4935771281  err=1.2e+00
v=3 (iii) s=1e-07  H=s*Heis+Sz0  (E+0.5)/s=-1.2499999991  ED (E+0.5)/s=-2.0925528521  err/s=8.4e-01
v=2 (ii) s=1e-07 c=1  (gs_energy-c)/s=-1.2499999902  ED=-2.4935771281  err=1.2e+00
v=2 (iii) s=1e-07  H=s*Heis+Sz0  (E+0.5)/s=-1.0000000006  ED (E+0.5)/s=-2.0925528521  err/s=1.1e+00
```
pA_gate.hybrid.out, the same c=1 rows before the fix:
```
v=3 (i) [vev(c+s*H)-c]/s/vev(H): s=1e-06,c=0:1.0000 s=1e-06,c=1:1.0000 s=3e-07,c=0:0.3333 s=3e-07,c=1:0.3333 ...
v=3 (ii) s=1e-07 c=1  (gs_energy-c)/s=-1.2499999968  ED=-2.4935771281  err=1.2e+00
v=2 (iii) s=1e-07  H=s*Heis+Sz0  (E+0.5)/s=-1.0000000017  ED (E+0.5)/s=-2.0925528521  err/s=1.1e+00
```
p0_terms.out: a constant reaches the session as an Id term, `[((1e-07+0j), [('Sx', 1), ('Sx', 2)]), ((1+0j), [('Id', 2)])]`.

KPM shift onset, pC_kpm_extreme.after.out and pC2_kpm_threshold.after.out (the second's rows are s=1.8e-14 and 1.5e-14, printed by %.0e):
```
v=3 s=1e-12  ED peak*s=0.3134  DMRG peak*s=0.3133  max|DMRG-ED|/peak=4.1e-04
v=3 s=1e-14  ED peak*s=0.3134  DMRG peak*s=0.3905  max|DMRG-ED|/peak=1.2e+00
v=2 s=1e-20  ED peak*s=0.3134  DMRG peak*s=0.3905  max|DMRG-ED|/peak=1.2e+00
v=3 s=2e-14  ED peak*s=0.3134  DMRG peak*s=0.3133  max|DMRG-ED|/peak=4.1e-04
v=3 s=1e-14  ED peak*s=0.3134  DMRG peak*s=0.3905  max|DMRG-ED|/peak=1.2e+00
v=2 s=2e-14  ED peak*s=0.3134  DMRG peak*s=0.3133  max|DMRG-ED|/peak=4.1e-04
v=2 s=1e-14  ED peak*s=0.3134  DMRG peak*s=0.3905  max|DMRG-ED|/peak=1.2e+00
```
The single-term Id AutoMPOs still on plain toMPO in v3 are at chain_session.h lines 1578, 2039, 11786, 11868, 11885 and 12072 (`grep -n '"Id",1'`).

Neighbours that hold, pB_neighbours.after.out:
```
v=3 B1 spinless s=1e-10 gs/s  DMRG=-3.1169710301  ED=-3.1169710301  err=4.9e-15
v=3 B2 Hubbard s=1e-10 gs/s  DMRG=-6.8759428090  ED=-6.8759428090  err=2.7e-15
v=3 B3 sector Sz=0 s=1e-10 excited/s=[-2.49357713 -2.00199536 -1.42146154] ED=[-2.49357713 -2.00199536 -1.42146154] maxerr=4.5e-13
v=3 B4 reuse one chain: s=1e+00:-2.4935771339 s=1e-08:-2.4935771339 s=1e+00:-2.4935771339 s=1e-10:-2.4935771339
v=3 B5 s=1e-40 gs/s=-2.4935771339 excited(n=2)/s=[-2.49357713 -2.00199536] (exact -2.4935771339, -2.0019953605)
v=2 B2 Hubbard s=1e-10 gs/s  DMRG=-6.8759428090  ED=-6.8759428090  err=8.9e-16
v=2 B5 s=1e-40 gs/s=-2.4935771339 excited(n=2)/s=[-2.49357713 -2.00199536] (exact -2.4935771339, -2.0019953605)
```
pD_complex.after.out (phase-hopping Hubbard, complex MPO): gs err 3.6e-15 to 6.7e-15 on v3 and 4.4e-16 to 4.4e-15 on v2 at s=1, 1e-8, 1e-10. The third excited level is off by 7.7e-6 to 3.7e-4 at every s including s=1 (v3 s=1 maxerr=3.7e-04), which is the penalty solver's own convergence and not scale-dependent.

By reading: `grep -n "\bH_ ="` gives only mpscpp2/chain_session.h:153 (set_hamiltonian) and mpscpp3/chain_session.h:572 (set_hamiltonian_mpo). svdMPO reads `Real cutoff = args.getReal("Cutoff",1E-13);` at mpscpp3/ITensor/itensor/mps/autompo.cc:1176 and mpscpp2 :1157, with `Real eps = 1E-14;` hardcoded at v3 :1172.
````

**Left open**: (1) Bond-local truncation, a different defect found by this item's review and pinned by strict xfails: svdMPO truncates bond by bond against the bond's own largest weight and the absolute 1e-13, so a bond whose strongest crossing term is below about 4e-7 (v3) or 3e-7 (v2) of the largest coefficient of the operator loses channels at any units, s=1 included. s*H + 1 and s*H + Sz0 at s=1e-7 give (E-1)/s of about -1.25 against -2.4936 and (E+0.5)/s of -1.25 (v3) or -1.00 (v2) against -2.0926, and a weak link J'=1e-7 in a J=1 chain reads 1/3 of itself; the per-bond cure is in cpp_plan. (2) "python" fails the same probes through the return sweep of `pyitensor/mpobuilder.py::to_mpo`, whose relative cutoff (`_BUILD_CUTOFF` from pyitensor/chain.py) measures a bond against the Schmidt weight of the whole operator: 0.0000 for the offset and the weak link and 0.0892 for the field at s=1e-7, 0.6667 for the offset and the weak link at 3e-7. A new lead, located by reading, owned by no cluster, not fixed. (3) A Hamiltonian handed to v3 as an already-built MPO (`set_hamiltonian(s*toMPO(H))`, i.e. `set_hamiltonian_mpo`) carries no term list to read a scale from, so `hscale_up_` is 1.0 there by design and davidson()'s absolute 1e-10 still acts (03 part B in this pass: 9.8e-6 off at s=1e-8, 1.2e-3 at 1e-9, 0.52 at 1e-10). (4) Real-time evolution at small units: v3 TDVP is 6.4e-4 off at s=1e-8 against 9.4e-13 at s=1 (from 4.0e-2 before the MPO half, the rest not investigated, the TDVP Krylov exponentiation being the obvious suspect); the v2 MPO-Taylor stepper diverges at s=1e-8 before and after (4.7e+49 and 9.1e+68), not investigated. (5) NH-DMRG on v2 at s=1e-8 is still wrong (-0.61437572-0.14372027j against -2.46104102, from 0.3 before) while v3's is right; both stop their Arnoldi at the same 1e-10*(1+|E|), so the difference lies elsewhere in v2's NH path, not located. (6) "python": the absolute breakdown and value tests of `pyitensor/dmrg.py::_lanczos_ground_state` (1.2e-6 off at s=1e-9, 2.8e-5 at 1e-10), owned by no cluster, pinned by the python-1e-10 strict xfail. (7) Changed by reading and not measured: CVM's `z*Id` in `Chain::cvm_dynamical_correlator` (exact for one term; no Python caller reaches the C++ CVM, which cvm.py solves itself, by grep), and the small-unit behaviour of `quench_tebd`, whose shift term goes into `bond_hamiltonians` rather than toMPO. (8) Not measured: `metts_vev`/`metts_dynamical_correlator` (imaginary-time TDVP on H_), sector-mode (QN) solves at small units beyond the review's Sz=0 check, and the iDMRG/VUMPS solvers (their own dense rows and vx_lanczos with absolute 1e-12 tolerances, never through build_mpo).

### 3. `multioperator.clean_threshold` dropped every term whose |coefficient| was at or below an absolute 1e-8, at every consumption point and inside the canonical form, so vev(1e-9*Sz0)/1e-9 was 0.000000 against -0.189361 on ED, "python", v3 and v2, gs_energy()/s of a 6-site Heisenberg chain written at s=1e-8 was 0.000000 against -2.493577 on all four, the canonical form proved 1e-9*(Sz0+1j*Sx0) Hermitian and 1e-9*Sz0 zero, and meanfield's own absolute 1e-10 emptied the mean-field Hamiltonian below that scale

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; cluster `scale`

**Status**: FIXED. `multioperator.clean_threshold` is 1e-12 and relative: `_filter_small` drops a term at or below 1e-12 times the largest |coefficient| of its own list, so the `0*identity()` placeholder and every other exact zero still go while an operator in small units keeps its terms, and `canonical.canonical_dict` drops a summed coefficient at or below 1e-12 times the larger of the operator's largest coefficient and the sum of the magnitudes collected into that signature, the second being what the rounding dust of a long accumulation scales with (2000 copies per signature leave 2.11e-12 of the largest coefficient but 1.07e-15 of their own sum). `meanfield` compares against 1e-10 of the largest coefficient of the Hamiltonian it decouples, and `mpsalgebra.is_hermitian`'s comment now says the rescaling by 1/max|c| is needed only by the probe. Pinned by `tests/test_audit_2026_09_25_scale.py` (31 tests, 26 of which fail on 8dd2198): `test_an_operator_in_small_units_keeps_its_terms`, `test_exact_zeros_are_still_dropped`, `test_the_drop_is_relative_to_the_largest_coefficient`, `test_rounding_dust_of_a_cancellation_is_zero`, `test_ordinary_hamiltonians_are_proven_hermitian_at_every_scale`, `test_a_long_accumulation_on_one_signature_is_proven_hermitian`, `test_vev_is_linear_in_the_scale_of_the_operator` (ED, "python", v3, v2), `test_correlator_is_quadratic_in_the_scale_of_the_operators` (ED, "python"), `test_ground_state_energy_is_scale_covariant` and `test_meanfield_does_not_depend_on_units`; the three strict xfails in `test_small_units_ground_state_still_open` belong to the small-units item, the v3 and v2 ones to be flipped by the C++ stage and the "python" one to stay until `pyitensor/dmrg.py` is fixed. NUMBERS CHANGE: vev(eps*Sz0)/eps at eps = 1e-8, 1e-9 and 1e-12 from 0.000000 to -0.189361 on ED, "python", v3 and v2 (6-site S=1/2 Heisenberg + 0.3*Sz0, `maxm=30`, `nsweeps=10`), and the KPM correlator C[1e-9*Sz0, 1e-9*Sz3] from identically 0 to 1e-18*C[Sz0,Sz3] on ED and "python"; gs_energy()/s of s*(6-site Heisenberg) at s=1e-8 from 0.000000 to -2.493577 on ED and "python" and to -1.249992 on v3 and -1.250000 on v2, and at s=1e-9 to -2.493577 (ED), -2.4935759 ("python", 1.2e-6 off by its own solver) and about -1.2497 to -1.24999 (v3, v2), the v2 and v3 values still wrong through the small-units mechanisms rather than masked as zero; `spinchain_meanfield(p=0, mode="ED")` on s*(4-site Heisenberg + 0.3*Sz_tot) at s=1e-9 from E_MF/s = 0.000000 and <Sz_i> = -0.5 to -0.627273 and +0.5, the s=1 answer, and at s=1e-11 from ValueError to the same; any term between 1e-12 and 1e-8 of an operator's largest coefficient is now kept where it was dropped, and a term above 1e-8 in absolute size but at or below 1e-12 of the largest coefficient (possible only when that exceeds 1e4) is now dropped, each moving a result by at most that term. Behaviour change: (1e-9*(Sz0+1j*Sx0)).is_hermitian() and (1e-9*Sz0).is_zero() go from True to False, the Heisenberg chain's simplify() at s=1e-9 keeps 15 of 16 terms instead of 1, and an O(1) Hamiltonian with an anti-Hermitian part between about 1e-10 and 1e-8 of its largest coefficient is no longer called Hermitian by the chain, so gs_energy() sends it to NH-DMRG (Heisenberg + 1j*1e-9*Sz0: True to False).

**Recorded in**: 2026-09-24c New leads. The absolute 1e-8 is a literal in `MultiOperator.clean()` since 539821a (2019-04-22), named `clean_threshold` in 2bfc9cc (2026-07-17) and read by the canonical form since 593b394 (2026-09-19); meanfield's absolute 1e-10 since 43d1a35 (2026-09-06).

**Where**: src/dmrgpy/multioperator.py `clean_threshold`/`_filter_small` (consumed by `to_terms`, `MO2list`/`write`, `MO2matrix`, `clean`); src/dmrgpy/multioperatortk/canonical.py `canonical_dict` (behind `simplify`, `is_zero`, `is_hermitian`, `is_antihermitian`, `is_dagger_pair`); src/dmrgpy/meanfield.py `_tol` in `_mean_field_hamiltonian` and `spinchain_meanfield`

**Severity**: high: silently wrong numbers (exact zeros) on every backend, ED included, for any operator or Hamiltonian written below 1e-8, and a false Hermiticity proof below it

The drop is now relative, at 1e-12, the factor mpscpp3/mo_terms.h's combine_terms and pyitensor/sector.py already use for the same job. At a consumption point the reference is the largest |coefficient| of the list itself, so an exact zero always goes, an operator written in small units keeps every term, and a small term next to a large one is dropped only below what the large one's rounding resolves ((1e9*Sz0+Sz1) keeps both, (1e9*Sz0+1e-9*Sz1) drops the second). In the canonical form the test is relative to what went INTO each sum, not to what came out: a summed coefficient is dropped at or below 1e-12 times the larger of the operator's largest coefficient and the sum of the magnitudes collected into that signature. Taking the scale from the sums would keep everything in the one case this exists for, H - H^dagger of a Hermitian H, whose sums are exact zeros or rounding dust with nothing larger left to compare against. The per-signature sum is needed and was measured, not assumed: 2000 copies of each of two terms minus their adjoints leave 2.11e-12 of the largest coefficient, which a floor on that alone would keep (not proven Hermitian), but 1.07e-15 of their own sum. The global floor keeps the canonical form consistent with the consumption-point filter, so a term every backend drops cannot keep an operator from being proven zero. meanfield compares against 1e-10 of the largest coefficient of the Hamiltonian it decouples (p, dimensionless, stays against 1e-10 itself). mpsalgebra.is_hermitian's rescaling by 1/max|c| stays, since the probe's 1e-20 is absolute, and its comment no longer claims the rescaling is what makes the proof scale-free. One behaviour moves on O(1) operators: an anti-Hermitian part between about 1e-10 and 1e-8 of the largest coefficient used to be dropped inside the proof and never reached the probe; now the probe sees it and calls it non-Hermitian, so gs_energy() sends such a Hamiltonian to NH-DMRG (06: g=1e-9 True to False); between 1e-12 and about 1e-10 the proof says not proven and the probe, at its absolute 1e-20 on the squared witness norm, still says Hermitian (g=1e-11), and below 1e-12 the proof holds (g=1e-13). The small-unit rows of v2 and v3 now reach the backend and come back at the Neel energy, which is the small-units item, not this one.

**Expected**: An operator means the same thing in any unit of energy: vev(eps*A) = eps*vev(A), C[eps*A, eps*B] = eps^2*C[A,B], E0(s*H) = s*E0(H), and the Hermiticity proof and is_zero() give the same verdict at every scale, while exact zeros (the 0*identity() placeholder) are still dropped and the rounding dust of H - H^dagger for an ordinary Hermitian H is still provably zero.

Repro (<scratch>/scale/01_clean_threshold.py and <scratch>/scale/06_hermiticity_window.py):

````python
# ==== 01_clean_threshold.py ====
# clean-threshold: multioperator.clean_threshold = 1e-8 drops every term
# whose |coef| is at or below 1e-8 in absolute value, at every consumption
# point (to_terms, MO2matrix, write) and in the canonical form.  Measure (i)
# the drop itself, symbolically, (ii) what it does to vev and gs_energy
# through the public API on ED, "python", v3 and v2, and (iii) the traps a
# relative rule has to clear: roundoff dust left by an exact cancellation
# (H - H^dagger of an ordinary Heisenberg or Hubbard H, and a long
# accumulation of copies of one term), and a large coefficient next to a
# small one.
import numpy as np
from dmrgpy import spinchain, fermionchain, cppext, multioperator, meanfield
from dmrgpy.multioperatortk import canonical

L = 6
def heis(sc, n):
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

def raw_sums(MO):
    """{signature: (summed coefficient, sum of |c| that went in)}, with no
    filter at all, from the canonical signatures."""
    out = {}
    for term in MO.op:
        sig, c = canonical._canonical_signature(term)
        s, a = out.get(sig, (0.0, 0.0))
        out[sig] = (s + c, a + abs(c))
    return out

print("clean_threshold =", multioperator.clean_threshold)
sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")

print("(i) the drop, symbolically")
for eps in (1e-6, 1e-8, 1e-9, 1e-12):
    A = eps*sc.Sz[0]
    N = eps*(sc.Sz[0] + 1j*sc.Sx[0]) # not Hermitian at any eps
    print("  eps=%.0e  (eps*Sz0).to_terms()=%s  (eps*Sz0).is_zero()=%s  "
          "(eps*(Sz0+1j*Sx0)).is_hermitian()=%s" % (
          eps, A.to_terms(), A.is_zero(), N.is_hermitian()))
B = 1e9*sc.Sz[0] + sc.Sz[1]
print("  (1e9*Sz0+Sz1).to_terms() =", B.to_terms())
B = 1e9*sc.Sz[0] + 1e-9*sc.Sz[1]
print("  (1e9*Sz0+1e-9*Sz1).to_terms() =", B.to_terms())

print("(iii) roundoff dust after an exact cancellation")
D = 0.1*sc.Sz[0] + 0.2*sc.Sz[0] - 0.3*sc.Sz[0]
(d, a), = raw_sums(D).values()
print("  0.1*Sz0+0.2*Sz0-0.3*Sz0: raw sum %.3e (inputs %.1f), is_zero()=%s" % (
      abs(d), a, D.is_zero()))
fs = fermionchain.Spinful_Fermionic_Chain(4)
phi = 0.3
hub = 0
for i in range(3):
    t = np.exp(1j*phi*(i+1))
    hub = hub + t*fs.Cdagup[i]*fs.Cup[i+1] + np.conj(t)*fs.Cdagup[i+1]*fs.Cup[i]
    hub = hub + t*fs.Cdagdn[i]*fs.Cdn[i+1] + np.conj(t)*fs.Cdagdn[i+1]*fs.Cdn[i]
for i in range(4): hub = hub + 2.0*fs.Nup[i]*fs.Ndn[i]
h = heis(sc, L)
np.random.seed(7)
cs = np.random.random(2000)
rep = multioperator.msum([c*sc.Sz[0]*sc.Sz[1] for c in cs])
rep = rep + multioperator.msum([c*sc.Sx[1]*sc.Sx[2] for c in cs[::-1]])
for label, H in (("Heisenberg", h), ("Hubbard, complex hopping", hub),
                 ("2000 copies per signature", rep)):
    for s in (1.0, 1e-6, 1e-9, 1e-12):
        X = s*H
        sums = raw_sums(X - X.get_dagger())
        cmax = max(abs(t[0]) for t in X.op)
        dust = max(abs(v[0]) for v in sums.values())
        per = max(abs(v[0])/v[1] for v in sums.values() if v[1] > 0)
        print("  %-26s s=%.0e  max|dust|/cmax=%.2e  max|dust|/sum|c| of its signature=%.2e  "
              "(X-X^dag).is_zero()=%s  X.is_hermitian()=%s  len(X.simplify().op)=%d/%d" % (
              label, s, dust/cmax, per, (X - X.get_dagger()).is_zero(),
              X.is_hermitian(), len(X.simplify().op), len(X.op)))
    NH = 0.1*(fs.Cdagup[0]*fs.Cup[1] if label.startswith("Hub") else 1j*sc.Sz[0])
    for s in (1.0, 1e-9):
        print("  %-26s s=%.0e  (s*(H+0.1*non-Hermitian)).is_hermitian()=%s" % (
              label, s, (s*(H + NH)).is_hermitian()))

print("(ii) vev(eps*Sz0)/eps against vev(Sz0), 6-site Heisenberg + 0.3*Sz0")
backends = ["ED", "python"] + [v for v in (3, 2) if cppext.available(v)]
for v in backends:
    c = spinchain.Spin_Chain(["S=1/2"]*L,
                             itensor_version=("python" if v == "ED" else v))
    c.set_hamiltonian(heis(c, L) + 0.3*c.Sz[0])
    c.maxm = 30; c.nsweeps = 10
    mode = "ED" if v == "ED" else "DMRG"
    ref = np.real(c.vev(c.Sz[0], mode=mode))
    row = " ".join("%.0e:%.6f" % (eps, np.real(c.vev(eps*c.Sz[0], mode=mode))/eps)
                   for eps in (1e-6, 1e-8, 1e-9, 1e-12))
    print("  %-6s vev(Sz0)=%.6f  vev(eps*Sz0)/eps  %s" % (v, ref, row))

print("(ii) gs_energy()/s of s*(6-site Heisenberg), fresh chain per s, exact -2.493577")
for v in backends:
    row = []
    for s in (1.0, 1e-8, 1e-9):
        c = spinchain.Spin_Chain(["S=1/2"]*L,
                                 itensor_version=("python" if v == "ED" else v))
        c.set_hamiltonian(s*heis(c, L))
        c.maxm = 30; c.nsweeps = 10
        try:
            e = "%.6f" % (np.real(c.gs_energy(mode=("ED" if v == "ED" else "DMRG")))/s)
        except Exception as ex:
            e = "raised %s: %s" % (type(ex).__name__, str(ex)[:60])
        row.append("s=%.0e:%s" % (s, e))
    print("  %-6s %s" % (v, "  ".join(row)))

print("(ii) spinchain_meanfield(p=0, mode='ED') on s*(Heisenberg + 0.3*Sz_tot), 4 sites")
for s in (1.0, 1e-9, 1e-11):
    c = spinchain.Spin_Chain(["S=1/2"]*4)
    c.set_hamiltonian(s*(heis(c, 4) + 0.3*sum(c.Sz)))
    np.random.seed(3)
    try:
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            out = meanfield.spinchain_meanfield(c, p=0.0, m0=[[0., 0., 0.5]]*4,
                                                maxite=100, mode="ED")
        mz = np.array(out.get_magnetization(mode="ED")).T[:, 2].real
        print("  s=%.0e  E_MF/s=%.6f  <Sz_i>=%s" % (
              s, np.real(out.gs_energy(mode="ED"))/s, np.round(mz, 6)))
    except Exception as ex:
        print("  s=%.0e  raised %s: %s" % (s, type(ex).__name__, str(ex)[:80]))

# ==== 06_hermiticity_window.py ====
# clean-threshold, the one behaviour a relative drop changes on an O(1)
# operator: an anti-Hermitian part between 1e-12 and 1e-8 of the largest
# coefficient.  The absolute 1e-8 dropped it in the symbolic proof, which
# mpsalgebra.is_hermitian runs on op/max|c|, so the operator was proven
# Hermitian and never reached the numerical probe.  Measure the proof, the
# chain's verdict, and the canonical coefficient that survives.
import numpy as np
from dmrgpy import spinchain

L = 6
sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for g in (1e-6, 1e-9, 1e-11, 1e-13):
    H = h + 1j*g*sc.Sz[0]
    np.random.seed(1)
    print("g=%.0e  (H+1j*g*Sz0).is_hermitian()=%s  chain.is_hermitian=%s  "
          "terms in simplify(H-H^dag)=%d" % (
          g, H.is_hermitian(), sc.is_hermitian(H),
          len([t for t in (H - H.get_dagger()).simplify().op if t[0] != 0])))

# run: cd <scratch>/scale && DMRGPY_SRC=<scratch>/parent/src <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.before.out   (parent)
#      cd <scratch>/scale && <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.after.out                               (working tree)
````

Observed, before:

````
==== 01_clean_threshold.before.out (parent tree, 8dd2198) ====
[run3] slot 2 acquired
clean_threshold = 1e-08
(i) the drop, symbolically
  eps=1e-06  (eps*Sz0).to_terms()=[((1e-06+0j), [('Sz', 1)])]  (eps*Sz0).is_zero()=False  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  eps=1e-08  (eps*Sz0).to_terms()=[]  (eps*Sz0).is_zero()=True  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  eps=1e-09  (eps*Sz0).to_terms()=[]  (eps*Sz0).is_zero()=True  (eps*(Sz0+1j*Sx0)).is_hermitian()=True
  eps=1e-12  (eps*Sz0).to_terms()=[]  (eps*Sz0).is_zero()=True  (eps*(Sz0+1j*Sx0)).is_hermitian()=True
  (1e9*Sz0+Sz1).to_terms() = [((1000000000+0j), [('Sz', 1)]), ((1+0j), [('Sz', 2)])]
  (1e9*Sz0+1e-9*Sz1).to_terms() = [((1000000000+0j), [('Sz', 1)])]
(iii) roundoff dust after an exact cancellation
  0.1*Sz0+0.2*Sz0-0.3*Sz0: raw sum 5.551e-17 (inputs 0.6), is_zero()=True
  Heisenberg                 s=1e+00  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e-06  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e-09  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=1/16
  Heisenberg                 s=1e-12  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=1/16
  Heisenberg                 s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  Heisenberg                 s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=True
  Hubbard, complex hopping   s=1e+00  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e-06  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e-09  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=1/17
  Hubbard, complex hopping   s=1e-12  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=1/17
  Hubbard, complex hopping   s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  Hubbard, complex hopping   s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=True
  2000 copies per signature  s=1e+00  max|dust|/cmax=2.11e-12  max|dust|/sum|c| of its signature=1.07e-15  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-06  max|dust|/cmax=5.95e-13  max|dust|/sum|c| of its signature=3.01e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-09  max|dust|/cmax=1.94e-12  max|dust|/sum|c| of its signature=9.83e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-12  max|dust|/cmax=1.52e-12  max|dust|/sum|c| of its signature=7.69e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=1/4000
  2000 copies per signature  s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  2000 copies per signature  s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=True
(ii) vev(eps*Sz0)/eps against vev(Sz0), 6-site Heisenberg + 0.3*Sz0
  ED     vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
  python vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
  3      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
  2      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
(ii) gs_energy()/s of s*(6-site Heisenberg), fresh chain per s, exact -2.493577
  ED     s=1e+00:-2.493577  s=1e-08:0.000000  s=1e-09:0.000000
  python s=1e+00:-2.493577  s=1e-08:0.000000  s=1e-09:0.000000
  3      s=1e+00:-2.493577  s=1e-08:0.000000  s=1e-09:0.000000
  2      s=1e+00:-2.493577  s=1e-08:0.000000  s=1e-09:0.000000
(ii) spinchain_meanfield(p=0, mode='ED') on s*(Heisenberg + 0.3*Sz_tot), 4 sites
  s=1e+00  E_MF/s=-0.627273  <Sz_i>=[0.5 0.5 0.5 0.5]
  s=1e-09  E_MF/s=0.000000  <Sz_i>=[-0.5 -0.5 -0.5 -0.5]
  s=1e-11  raised ValueError: meanfield: the mean-field Hamiltonian is empty (every coefficient vanished)

==== 06_hermiticity_window.before.out (parent tree) ====
[run3] slot 1 acquired
g=1e-06  (H+1j*g*Sz0).is_hermitian()=False  chain.is_hermitian=False  terms in simplify(H-H^dag)=1
g=1e-09  (H+1j*g*Sz0).is_hermitian()=True  chain.is_hermitian=True  terms in simplify(H-H^dag)=0
g=1e-11  (H+1j*g*Sz0).is_hermitian()=True  chain.is_hermitian=True  terms in simplify(H-H^dag)=0
g=1e-13  (H+1j*g*Sz0).is_hermitian()=True  chain.is_hermitian=True  terms in simplify(H-H^dag)=0
````

Observed, after:

````
==== 01_clean_threshold.after.out (working tree) ====
[run3] slot 1 acquired
clean_threshold = 1e-12
(i) the drop, symbolically
  eps=1e-06  (eps*Sz0).to_terms()=[((1e-06+0j), [('Sz', 1)])]  (eps*Sz0).is_zero()=False  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  eps=1e-08  (eps*Sz0).to_terms()=[((1e-08+0j), [('Sz', 1)])]  (eps*Sz0).is_zero()=False  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  eps=1e-09  (eps*Sz0).to_terms()=[((1e-09+0j), [('Sz', 1)])]  (eps*Sz0).is_zero()=False  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  eps=1e-12  (eps*Sz0).to_terms()=[((1e-12+0j), [('Sz', 1)])]  (eps*Sz0).is_zero()=False  (eps*(Sz0+1j*Sx0)).is_hermitian()=False
  (1e9*Sz0+Sz1).to_terms() = [((1000000000+0j), [('Sz', 1)]), ((1+0j), [('Sz', 2)])]
  (1e9*Sz0+1e-9*Sz1).to_terms() = [((1000000000+0j), [('Sz', 1)])]
(iii) roundoff dust after an exact cancellation
  0.1*Sz0+0.2*Sz0-0.3*Sz0: raw sum 5.551e-17 (inputs 0.6), is_zero()=True
  Heisenberg                 s=1e+00  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e-06  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e-09  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e-12  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=15/16
  Heisenberg                 s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  Heisenberg                 s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  Hubbard, complex hopping   s=1e+00  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e-06  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e-09  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e-12  max|dust|/cmax=0.00e+00  max|dust|/sum|c| of its signature=0.00e+00  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=16/17
  Hubbard, complex hopping   s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  Hubbard, complex hopping   s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  2000 copies per signature  s=1e+00  max|dust|/cmax=2.11e-12  max|dust|/sum|c| of its signature=1.07e-15  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-06  max|dust|/cmax=5.95e-13  max|dust|/sum|c| of its signature=3.01e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-09  max|dust|/cmax=1.94e-12  max|dust|/sum|c| of its signature=9.83e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e-12  max|dust|/cmax=1.52e-12  max|dust|/sum|c| of its signature=7.69e-16  (X-X^dag).is_zero()=True  X.is_hermitian()=True  len(X.simplify().op)=2/4000
  2000 copies per signature  s=1e+00  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
  2000 copies per signature  s=1e-09  (s*(H+0.1*non-Hermitian)).is_hermitian()=False
(ii) vev(eps*Sz0)/eps against vev(Sz0), 6-site Heisenberg + 0.3*Sz0
  ED     vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:-0.189361 1e-09:-0.189361 1e-12:-0.189361
  python vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:-0.189361 1e-09:-0.189361 1e-12:-0.189361
  3      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:-0.189361 1e-09:-0.189361 1e-12:-0.189361
  2      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:-0.189361 1e-09:-0.189361 1e-12:-0.189361
(ii) gs_energy()/s of s*(6-site Heisenberg), fresh chain per s, exact -2.493577
  ED     s=1e+00:-2.493577  s=1e-08:-2.493577  s=1e-09:-2.493577
  python s=1e+00:-2.493577  s=1e-08:-2.493577  s=1e-09:-2.493577
  3      s=1e+00:-2.493577  s=1e-08:-1.249992  s=1e-09:-1.249659
  2      s=1e+00:-2.493577  s=1e-08:-1.250000  s=1e-09:-1.249986
(ii) spinchain_meanfield(p=0, mode='ED') on s*(Heisenberg + 0.3*Sz_tot), 4 sites
  s=1e+00  E_MF/s=-0.627273  <Sz_i>=[0.5 0.5 0.5 0.5]
  s=1e-09  E_MF/s=-0.627273  <Sz_i>=[0.5 0.5 0.5 0.5]
  s=1e-11  E_MF/s=-0.627273  <Sz_i>=[0.5 0.5 0.5 0.5]

==== 06_hermiticity_window.after.out (working tree) ====
[run3] slot 1 acquired
g=1e-06  (H+1j*g*Sz0).is_hermitian()=False  chain.is_hermitian=False  terms in simplify(H-H^dag)=1
g=1e-09  (H+1j*g*Sz0).is_hermitian()=False  chain.is_hermitian=False  terms in simplify(H-H^dag)=1
g=1e-11  (H+1j*g*Sz0).is_hermitian()=False  chain.is_hermitian=True  terms in simplify(H-H^dag)=1
g=1e-13  (H+1j*g*Sz0).is_hermitian()=True  chain.is_hermitian=True  terms in simplify(H-H^dag)=0
````

**NUMBERS CHANGE**: vev(eps*Sz0)/eps at eps=1e-8,1e-9,1e-12: 0.000000 -> -0.189361 on ED, python, v3, v2 (6-site S=1/2 Heisenberg + 0.3*Sz0, maxm=30, nsweeps=10). gs_energy()/s of s*(6-site Heisenberg), s=1e-8: 0.000000 -> -2.493577 (ED, python), -1.249992 (v3), -1.250000 (v2); s=1e-9: 0.000000 -> -2.493577 (ED), -2.4935759 (python), -1.249659 (v3), -1.249986 (v2). meanfield p=0 mode=ED, s*(4-site Heisenberg + 0.3*Sz_tot): s=1e-9 E_MF/s 0.000000 -> -0.627273 and <Sz_i> -0.5 -> +0.5; s=1e-11 ValueError -> -0.627273. KPM C[1e-9*Sz0,1e-9*Sz3]: identically 0 -> 1e-18*C[Sz0,Sz3] on ED and python. Chain Hermiticity verdict of Heisenberg + 1j*1e-9*Sz0: True -> False (NH-DMRG dispatch).

**Tests**: `tests/test_audit_2026_09_25_scale.py::test_an_operator_in_small_units_keeps_its_terms`; `tests/test_audit_2026_09_25_scale.py::test_exact_zeros_are_still_dropped`; `tests/test_audit_2026_09_25_scale.py::test_the_drop_is_relative_to_the_largest_coefficient`; `tests/test_audit_2026_09_25_scale.py::test_rounding_dust_of_a_cancellation_is_zero`; `tests/test_audit_2026_09_25_scale.py::test_ordinary_hamiltonians_are_proven_hermitian_at_every_scale`; `tests/test_audit_2026_09_25_scale.py::test_a_long_accumulation_on_one_signature_is_proven_hermitian`; `tests/test_audit_2026_09_25_scale.py::test_vev_is_linear_in_the_scale_of_the_operator`; `tests/test_audit_2026_09_25_scale.py::test_correlator_is_quadratic_in_the_scale_of_the_operators`; `tests/test_audit_2026_09_25_scale.py::test_ground_state_energy_is_scale_covariant`; `tests/test_audit_2026_09_25_scale.py::test_meanfield_does_not_depend_on_units`

**Reviewer (CONFIRMED, fix HOLDS)**: I reproduced it with my own run of the fix agent's 01 on the parent tree (8dd2198). vev(eps*Sz0)/eps is 0.000000 at eps=1e-8, 1e-9 and 1e-12 on ED, "python", v3 and v2; gs_energy()/s of the 6-site Heisenberg chain written at s=1e-8 and 1e-9 is 0.000000 on all four; spinchain_meanfield at s=1e-9 gives E_MF/s=0 with <Sz_i>=-0.5, and at s=1e-11 it raises ValueError; and the canonical form calls 1e-9*(Sz0+1j*Sx0) Hermitian and 1e-9*Sz0 zero. None of this is intended or documented anywhere beyond the third pass's own 'New leads' bullet, and it is not a probe artifact: ED is exact and has no solver in the way. Nothing is struck from the defect itself. One size statement is narrowed. The v3 and v2 values after the fix at small s are single draws from ITensor's randomMPS, which np.random.seed does not reach. My run gave -1.249995 and -1.249908 on v3 and -1.250000 and -1.249763 on v2 at s=1e-8 and 1e-9, against the record's -1.249992/-1.249659 and -1.250000/-1.249986. So the record should say 'about -1.2497 to -1.2500' and not quote the digits.

The fix holds for the property it claims, and nothing that should stay put moved. The cluster file gives 31 passed and 3 xfailed on the fixed tree; my copy run against the parent gives 26 failed, 5 passed and 3 xfailed. A regression slice (test_multioperator_canonical, test_audit_2026_09_24c_hermiticity, test_meanfield, test_jordan_wigner) gives 56 passed.

Ordinary case. I built a hybrid tree (the parent plus only this cluster's four files) to separate this cluster's changes from other clusters' uncommitted edits. On an 8-site Heisenberg chain + 0.3*Sz0 at J=1, ED and "python" agree on E0, <Sz3> and an 8-point KPM C(w) to every printed digit (15 for E0) across parent, hybrid and fixed, with and without a 6e-17 dust term. v2 and v3 cannot be compared below their own run-to-run noise: the parent against itself, on identical input (its dust row drops the term), differs at about 1e-15 in E0, 1e-11 in <Sz3> and 1e-6 in C(w), and the cross-tree differences show no direction beyond that.

The one change on O(1) Hamiltonians is a correctness gain, and the record should say so in numbers. For Heisenberg + 0.3*Sz0 + 1j*1e-9*Sz0, NH-DMRG now returns Im E0 = -1.894e-10 = g*<Sz0> (the exact first order) on ED, "python", v3 and v2, where the parent returned exactly 0; Re E0 is unchanged to 12 digits. The window between 1e-12 and 1e-10 (g=1e-11) still drops the decay rate on both trees, through the probe's absolute 1e-20. That is disclosed and older than this fix.

The NUMBERS CHANGE statement is incomplete in three places:
- The KPM and EX correlators C[1e-9*Sz0, 1e-9*Sz3] now equal 1e-18*C to within 3.7e-15 (KPM) and 4.7e-16 (EX) on v3 and v2 as well, not only on ED and "python".
- "python" submode="TD" at eps <= 1e-8 goes from a ZeroDivisionError on the parent to a silently wrong spectrum (8.76e-2 relative at 1e-8, 2.42e-1 at 1e-9); see new issue 1. v3 TD comes out 1.0e-4 off.
- On v2 and v3, vev(eps*Sz0) is still exactly 0 at eps <= 1e-15. That is ITensor's isZero(1E-14), which belongs to small-units.

Corners the relative floor introduced (new issues 2 to 4). They need a coefficient ratio of at least 1e12 (1e10 in meanfield), so they do not overturn HOLDS, but each one is a regression against the parent:
- the scale is taken before summation, so cancelled pairs and identity offsets count;
- meanfield counts the constant offset in its scale;
- non-finite coefficients empty the whole operator.

A simpler design removes all three. Raw terms at a consumption point are never rounding dust, because dust appears only when canonical_dict sums a signature, and criterion (b) already handles that case. So _filter_small could drop exact zeros only (and raise on non-finite coefficients), and canonical_dict could keep criterion (b) alone. This is an option, not a requirement. It would make (1e9*Sz0+1e-9*Sz1).simplify() keep both terms, so test_the_drop_is_relative_to_the_largest_coefficient would change.

Test pinning. That test passes on the parent for the wrong reason (1e-9 is also below the absolute 1e-8); (1e9*Sz0 + 1e-5*Sz1) would discriminate between the trees. The remaining 26 parent failures do pin the property. The mpsalgebra.py:330 comment is now accurate, and there is no double rescale.

Reviewer's evidence:

````
Review folder R=<scratch>/review/scale

R/01_clean_threshold.before.out (parent):
  ED     vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
  3      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:0.000000 1e-09:0.000000 1e-12:0.000000
  ED     s=1e+00:-2.493577  s=1e-08:0.000000  s=1e-09:0.000000
  s=1e-11  raised ValueError: meanfield: the mean-field Hamiltonian is empty (every coefficient vanished)
R/01_clean_threshold.after.out (fixed):
  2      vev(Sz0)=-0.189361  vev(eps*Sz0)/eps  1e-06:-0.189361 1e-08:-0.189361 1e-09:-0.189361 1e-12:-0.189361
  python s=1e+00:-2.493577  s=1e-08:-2.493577  s=1e-09:-2.493577
  3      s=1e+00:-2.493577  s=1e-08:-1.249995  s=1e-09:-1.249908
  2      s=1e+00:-2.493577  s=1e-08:-1.250000  s=1e-09:-1.249763
  s=1e-11  E_MF/s=-0.627273  <Sz_i>=[0.5 0.5 0.5 0.5]
R/test_fixed.out: 31 passed, 3 xfailed in 8.86s ; R/test_parent.out: 26 failed, 5 passed, 3 xfailed in 7.40s ; R/test_regression_surface.out: 56 passed in 6.71s
R/r1b_nh_window_field.{before,after}.out (Heisenberg+0.3*Sz0+1j*g*Sz0):
  parent g=1e-09 v=3      chain.is_hermitian=True  gs_energy=-2.523188643451+0.000e+00j
  fixed  g=1e-09 v=3      chain.is_hermitian=False  gs_energy=-2.523188643451-1.894e-10j
  fixed  g=1e-09 v=ED     chain.is_hermitian=False  gs_energy=-2.523188643451-1.894e-10j
  fixed  g=1e-11 v=ED     chain.is_hermitian=True  gs_energy=-2.523188643451+0.000e+00j
R/r3_ordinary_identity.{parent,hybrid,fixed}.out, identical on all three trees:
  dust=0.0e+00 v=ED     E0=-3.406163131268207  <Sz3>=0.116290261219943  C(w)=0.000179650538 0.042467933006 ...
  dust=0.0e+00 v=python E0=-3.406163131268208  <Sz3>=0.116290261220425  C(w)=0.000179472871 0.042472819171 ...
  v3 noise on the parent against itself (same input): E0=-3.406163131268207 <Sz3>=0.116290261216837 C=0.000180274271 vs E0=-3.406163131268213 <Sz3>=0.116290261222451 C=0.000179579511
R/r4_downstream_floors.after.out:
  v=3      submode=KPM  max|C1|=0.2557  rel.dev=3.74e-15
  v=2      submode=EX   max|C1|=0.1467  rel.dev=4.73e-16
  v=python submode=TD   max|C1|=0.1464  rel.dev=2.42e-01
  v=3      submode=TD   max|C1|=0.1464  rel.dev=1.04e-04
  v=3      1e-13:1.000000 1e-14:1.000000 1e-15:-0.000000 1e-16:-0.000000 ...
R/r5b_td_krylov_parent.before.out: stock ... 1e-06:2.58e-05  1e-08:raised ZeroDivisionError  1e-09:raised ZeroDivisionError
````

### 4. On `itensor_version="julia_live"`, `set_initial_wf` and `set_initial_wf_guess` never reach the solver, which starts from its own random state (|<target|gs>|^2 = 0.0008 for a guess that is an exact eigenstate, 0.0538 for a state set as it is, on a 3-site chain), `set_gs` raises `AttributeError` on every Julia MPS, and `gs_energy(wf0=x)` on a solved chain raises `TypeError`

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `session`

**Status**: FIXED. `mark_injected()` now marks a julia_live state the way it marks one on the session backends, and both julia_live branches of `groundstate.gs_energy()` go through a new `groundstate._gs_energy_julia()`, which reads the mark with `gs_energy_single()`'s precedence: an explicit `wf0=` is swept from, or taken unswept with `reconverge=False`; a state set with `set_gs()`/`set_initial_wf()` is taken unswept through `_take_injected_state()`, which on a backend with no session computes e0 = <x|H|x> with the Julia MPS's own algebra and records no solver key, as julia_live's own solves do not; and a `set_initial_wf_guess()` state is swept from. `mpsjulialive/groundstate.get_gs_dmrg()` always solves now (its `computed_gs` short circuit returned the stored MPS alone, which its only caller unpacked as a pair), and hands the solver a Julia-side copy of the start with the stale prime level on its Link indices cleared, which operator algebra leaves there and which made the first warm start raise a DimensionMismatch inside eigsolve; the Julia MPS carries `mode="DMRG"`, which `set_gs()` dispatches on. The julia_live KPM measures a set state from its own energy with the window on the band edges, taking the lower edge from a solve of its own when the state was supplied (`mpsjulialive/dynamics._min_energy()`) and restoring the supplied mark that `restart()` clears, as the session backends do since the 2026-09-24c audit's finding 1: anchored on the state's energy instead, the curve is 1.28 of the peak off `mode="ED"`'s, against 8.2e-4. Pinned by `tests/test_audit_2026_09_25_session.py::test_warm_start_setters_reach_the_solver[julia_live]` (with `python` and `v3` rows for the same contract) and `::test_julia_kpm_measures_a_set_state_from_its_own_energy`. NUMBERS CHANGE for every julia_live caller of `set_initial_wf`/`set_initial_wf_guess`: on the 3-site Heisenberg chain the guess of an exact doublet member goes from |<target|gs>|^2 = 0.0008 (0.0008 to 0.68 over runs) to 1.0000, `set_initial_wf` of the other member from 0.0538 to 1.0000, and on a fresh chain `set_initial_wf(x)` then `gs_energy()` from the solved -1.000000 to <x|H|x> (-0.102571 for that run's x). `set_gs()` and `gs_energy(wf0=x)` on a solved chain raised before, so they have no old number; a julia_live KPM after a plain solve is unchanged, since there e0 is the band edge.

**Recorded in**: 2026-09-24b finding 12 Status (open). 867e2b4 built the injection mark for the session backends and left julia_live out; recorded as open in the Status of docs/audit_2026_09_24b_hole_hunt.md finding 12. The wf.mode read in set_gs and the get_gs_dmrg short circuit are older.

**Where**: src/dmrgpy/groundstate.py:443-455 on 8dd2198 (gs_energy's two julia_live branches call get_gs_dmrg without reading pending_injection), groundstate.py:122 (mark_injected returns unmarked when there is no session), groundstate.py:663 (set_gs reads wf.mode, which mpsjulialive/mps.py's MPS never had), src/dmrgpy/mpsjulialive/groundstate.py:10-11 (get_gs_dmrg's computed_gs short circuit returns the stored MPS alone, which its only caller unpacks as a pair), src/dmrgpy/mpsjulialive/dynamics.py:155-157 (the KPM window takes emin = e0, which for a set state is its own energy, not the band edge)

The injection machinery marks a state set by hand and lets gs_energy_single() take it unswept, with e0 = <x|H|x>, or sweep from it. julia_live has no session, so mark_injected() stored the state and returned unmarked, and gs_energy()'s julia_live branches called get_gs_dmrg() with no start, which solved from a random state and replaced the one that was set: a guess that is an exact eigenstate came back with an overlap anywhere from 0.0008 to 0.68 over runs, and a state set as it is was re-solved. set_gs() never got that far, since the Julia MPS has no mode attribute and set_gs() dispatches on it. The explicit gs_energy(wf0=x) on a solved chain reached get_gs_dmrg() with computed_gs set, which returned the stored MPS alone, and the caller's unpacking raised. Taking the mark into the Julia solve also exposed that an MPS built by operator algebra carries a stale prime level on its Link indices, which dmrg() cannot take, and that the julia_live KPM anchored its window on the state's energy, which is the band edge only for a solved state.

**Expected**: The contract gs_energy_single() has on the other backends: set_gs()/set_initial_wf() taken unswept with e0 = <x|H|x>, set_initial_wf_guess() swept from, gs_energy(wf0=x) swept from on a solved chain, the caller's state not moved, and vev and the correlators measuring the set state, KPM from its own energy on the band-edge window, as mode="ED" does.

Repro (<scratch>/session/01_julia_warm_start.py (outputs 01_julia_warm_start.{before,after}.out and .filtered.out) and 07_julia_kpm_window.py (07_julia_kpm_window.{before,after}.out)):

````python
# 01_julia_warm_start.py
# julia-warm-start: do set_gs, set_initial_wf and set_initial_wf_guess reach
# the julia_live solver, under the contract gs_energy_single has on the other
# backends? 3-site Heisenberg chain, whose ground doublet has two exact
# members up/dn (the 2Sz=+-1 components of the solved state). Julia's own
# RNG is not seeded by numpy, so the solved state varies from run to run.
import io, contextlib
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()

def attempt(label, f):
    try:
        return f()
    except Exception as err:
        print("%s raised %s: %s" % (label, type(err).__name__, err))

def heis(n=3):
    sc = spinchain.Spin_Chain([2]*n, itensor_version="julia_live")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 10, 10
    return sc

def energy(sc, wf):
    return float(np.real(wf.aMb(sc.hamiltonian, wf)/wf.dot(wf)))

def fid(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))

np.random.seed(1)
sc = heis()
e_solved = quiet(sc.gs_energy)
s = sc.get_gs().copy()
Szt = sc.Sz[0] + sc.Sz[1] + sc.Sz[2]
up = (Szt + 0.5)*s ; up = up*(1/np.sqrt(up.dot(up).real))
dn = (0.5 - 1*Szt)*s ; dn = dn*(1/np.sqrt(dn.dot(dn).real))
sz = float(np.real(sc.vev(Szt)))
target = dn if sz > 0 else up   # the member the solved state is far from
other = up if target is dn else dn
x = s + 0.4*(sc.Sx[0]*s) ; x = x*(1/np.sqrt(x.dot(x).real)) # not an eigenstate
print("solved e0 = %.6f, <Sz_tot> of the solved state = %+.4f" % (e_solved, sz))

# set_initial_wf_guess: the next solve should sweep from the guess, an exact
# eigenstate, so the result stays on it
sc.set_initial_wf_guess(target)
e1 = quiet(sc.gs_energy)
print("set_initial_wf_guess(target): gs_energy = %.6f, |<target|gs>|^2 = %.4f"
      % (e1, fid(sc.get_gs(), target)))

# set_initial_wf: taken as it is, no sweep
sc.set_initial_wf(other)
e2 = quiet(sc.gs_energy)
print("set_initial_wf(other):         gs_energy = %.6f, |<other|gs>|^2 = %.4f"
      % (e2, fid(sc.get_gs(), other)))

# set_gs of a state that is not an eigenstate: gs_energy should be <x|H|x>,
# and vev and a correlator should read x
def set_gs_x():
    sc.set_gs(x)
    e3 = quiet(sc.gs_energy)
    print("set_gs(x):                     gs_energy = %.6f, <x|H|x> = %.6f, |<x|gs>|^2 = %.4f"
          % (e3, energy(sc, x), fid(sc.get_gs(), x)))
    print("  vev(Sz0 Sz1) = %+.6f, <x|Sz0 Sz1|x> = %+.6f"
          % (np.real(sc.vev(sc.Sz[0]*sc.Sz[1])),
             np.real(x.dot(sc.Sz[0]*(sc.Sz[1]*x)))))
    es = np.linspace(-6.0, 6.0, 601)
    es, d = quiet(lambda: sc.get_dynamical_correlator(
        name=(sc.Sz[0], sc.Sz[1]), submode="KPM", es=es, delta=0.2))
    print("  KPM (Sz0,Sz1) after set_gs: int C = %+.6f, e0 = %.6f, |<x|gs>|^2 = %.4f"
          % (float(np.real(np.trapezoid(d, es))), sc.e0, fid(sc.get_gs(), x)))
attempt("set_gs(x)", set_gs_x)

# the caller's state is not moved by a sweep from it
keep = x.copy() ; e_x = energy(sc, x)
sc.set_initial_wf_guess(x)
e4 = quiet(sc.gs_energy)
print("set_initial_wf_guess(x): gs_energy = %.6f; the caller's x after: <H> %.12f -> %.12f, "
      "|<keep|x>|^2 = %.12f" % (e4, e_x, energy(sc, x), fid(x, keep)))

# the explicit warm start on a solved chain
def explicit():
    e5 = quiet(lambda: sc.gs_energy(wf0=target))
    print("gs_energy(wf0=target) on a solved chain: %.6f, |<target|gs>|^2 = %.4f"
          % (e5, fid(sc.get_gs(), target)))
attempt("gs_energy(wf0=target) on a solved chain", explicit)

# the same on a chain that was never solved: set_initial_wf(x) then gs_energy()
sc2 = heis()
x2 = sc2.random_state() ; x2 = x2*(1/np.sqrt(x2.dot(x2).real))
sc2.set_initial_wf(x2)
e6 = quiet(sc2.gs_energy)
print("fresh chain, set_initial_wf(random x): gs_energy = %.6f, <x|H|x> = %.6f, |<x|gs>|^2 = %.4f"
      % (e6, energy(sc2, x2), fid(sc2.get_gs(), x2)))

# 07_julia_kpm_window.py
# julia-warm-start, the KPM half: after set_gs of the first excited
# eigenstate |1> of the 3-site chain in Bz=0.3, does the julia_live KPM
# agree with mode="ED" KPM, which measures from E_1 on the band-edge window?
# Compared with the window anchored on the state's own energy instead (what
# emin = e0 would give), by switching state_supplied off.
import io, contextlib
import numpy as np
import dmrgpy
from dmrgpy import spinchain, groundstate
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()

def heis(sc, B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h + B*(sc.Sz[0] + sc.Sz[1] + sc.Sz[2])

ES = np.linspace(-0.8, 2.4, 641)
KW = dict(submode="KPM", delta=0.05, es=ES)
ed = spinchain.Spin_Chain(["S=1/2"]*3)
ed.set_hamiltonian(heis(ed, 0.3))
ed.get_gs(mode="ED")
ee, ww = ed.get_excited_states(n=2, mode="ED")
ed.set_gs(ww[1])
_, yed = quiet(lambda: ed.get_dynamical_correlator(mode="ED",
        name=(ed.Sx[0] + 1j*ed.Sy[0], ed.Sx[0] - 1j*ed.Sy[0]), **KW))
peak = np.max(np.abs(yed))

sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="julia_live")
sc.set_hamiltonian(heis(sc, 0.3))
sc.maxm, sc.nsweeps = 20, 12
quiet(sc.gs_energy)
Sp = [sc.Sx[i] + 1j*sc.Sy[i] for i in range(3)]
one = (Sp[0] + Sp[1] + Sp[2])*sc.get_gs()
one = one*(1/np.sqrt(one.dot(one).real))
pair = (Sp[0], sc.Sx[0] - 1j*sc.Sy[0])
try:
    sc.set_gs(one)
    _, y = quiet(lambda: sc.get_dynamical_correlator(name=pair, **KW))
    print("julia_live after set_gs(|1>): e0 = %.6f, max|y - y_ED| = %.3e of the %.3f peak"
          % (sc.gs_energy(), np.max(np.abs(y - yed))/peak, peak))
    real = groundstate.state_supplied
    groundstate.state_supplied = lambda self: False # window on e0 = E_1
    try:
        _, y0 = quiet(lambda: sc.get_dynamical_correlator(name=pair, **KW))
    finally:
        groundstate.state_supplied = real
    print("  with the window anchored on e0 instead: max|y - y_ED| = %.3e of the peak"
          % (np.max(np.abs(y0 - yed))/peak))
except Exception as err:
    print("julia_live after set_gs(|1>) raised %s: %s" % (type(err).__name__, err))

# Invocation, from the cluster folder (before: parent tree; after: working tree). The
# grep drops juliapkg's resolver lines, and the .filtered.out copies of 01 also drop the
# Julia Pkg lines ([uuid] package lists, Updating/Resolving/Info/Project/Manifest):
# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 01_julia_warm_start.py 2>&1 | grep -v "^\[juliapkg\]\|^           |" | tee 01_julia_warm_start.before.out
# cd $WF/session && ../run3.sh -u 01_julia_warm_start.py 2>&1 | grep -v "^\[juliapkg\]\|^           |" | tee 01_julia_warm_start.after.out
# (07_julia_kpm_window.py the same way, into 07_julia_kpm_window.{before,after}.out)
````

Observed, before:

````
01_julia_warm_start.before.filtered.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
solved e0 = -1.000000, <Sz_tot> of the solved state = +0.4640
set_initial_wf_guess(target): gs_energy = -1.000000, |<target|gs>|^2 = 0.0008
set_initial_wf(other):         gs_energy = -1.000000, |<other|gs>|^2 = 0.0538
set_gs(x) raised AttributeError: 'MPS' object has no attribute 'mode'
set_initial_wf_guess(x): gs_energy = -1.000000; the caller's x after: <H> -0.971648960582 -> -0.971648960582, |<keep|x>|^2 = 1.000000000000
gs_energy(wf0=target) on a solved chain raised TypeError: cannot unpack non-iterable MPS object
fresh chain, set_initial_wf(random x): gs_energy = -1.000000, <x|H|x> = -0.381718, |<x|gs>|^2 = 0.4194

07_julia_kpm_window.before.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
julia_live after set_gs(|1>) raised AttributeError: 'MPS' object has no attribute 'mode'

(Two earlier parent runs of the first version of 01, before the set_gs lines were wrapped, gave |<target|gs>|^2 = 0.6736 and 0.6793 for the guess and 0.1560 and 0.0544 for set_initial_wf; Julia's RNG is not seeded.)
````

Observed, after:

````
01_julia_warm_start.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
solved e0 = -1.000000, <Sz_tot> of the solved state = +0.4684
set_initial_wf_guess(target): gs_energy = -1.000000, |<target|gs>|^2 = 1.0000
set_initial_wf(other):         gs_energy = -1.000000, |<other|gs>|^2 = 1.0000
set_gs(x):                     gs_energy = -0.976470, <x|H|x> = -0.976470, |<x|gs>|^2 = 1.0000
  vev(Sz0 Sz1) = -0.141181, <x|Sz0 Sz1|x> = -0.141181
  KPM (Sz0,Sz1) after set_gs: int C = -0.140912, e0 = -0.976470, |<x|gs>|^2 = 1.0000
set_initial_wf_guess(x): gs_energy = -1.000000; the caller's x after: <H> -0.976469806592 -> -0.976469806592, |<keep|x>|^2 = 1.000000000000
gs_energy(wf0=target) on a solved chain: -1.000000, |<target|gs>|^2 = 1.0000
fresh chain, set_initial_wf(random x): gs_energy = -0.102571, <x|H|x> = -0.102571, |<x|gs>|^2 = 1.0000

07_julia_kpm_window.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
julia_live after set_gs(|1>): e0 = -0.850000, max|y - y_ED| = 8.187e-04 of the 6.130 peak
  with the window anchored on e0 instead: max|y - y_ED| = 1.278e+00 of the peak
````

**NUMBERS CHANGE**: julia_live only. 3-site Heisenberg chain: set_initial_wf_guess(exact doublet member) then gs_energy(): |<target|gs>|^2 0.0008 (0.0008 to 0.68 over runs) -> 1.0000; set_initial_wf(other member): 0.0538 -> 1.0000; fresh chain set_initial_wf(random x) then gs_energy(): the solved -1.000000 -> <x|H|x> (-0.102571 in the after run). set_gs() and gs_energy(wf0=x) on a solved chain raised before (AttributeError, TypeError). KPM after a plain solve unchanged.

**Tests**: `tests/test_audit_2026_09_25_session.py::test_warm_start_setters_reach_the_solver[julia_live]`; `tests/test_audit_2026_09_25_session.py::test_warm_start_setters_reach_the_solver[python]`; `tests/test_audit_2026_09_25_session.py::test_warm_start_setters_reach_the_solver[v3]`; `tests/test_audit_2026_09_25_session.py::test_julia_kpm_measures_a_set_state_from_its_own_energy[julia_live]`

**Reviewer (CONFIRMED, fix HOLDS)**: The defect is real as stated. I took the parent tree by reading and did not rerun it, since the JIT is costly and importing juliacall from the parent tree re-resolves the shared Julia project. In the parent, mark_injected() returns False whenever there is no session (`if getattr(self,"_session",None) is None: ... return False`). Both julia_live branches of gs_energy() call `get_gs_dmrg(self,**kwargs)` with no start. get_gs_dmrg short-circuits with `if self.computed_gs: return self.wf0`, which its only caller unpacks as a pair, and it otherwise starts from `self.random_state()`. mpsjulialive/mps.py's MPS has no `mode` attribute, which set_gs() dispatches on. So every setter was dropped, set_gs() raised AttributeError and gs_energy(wf0=x) on a solved chain raised TypeError, as the fix agent's before output shows. The sizes attached (|<target|gs>|^2 = 0.0008, 0.0538) are single draws from an unseeded Julia RNG; the fix agent says as much (0.0008 to 0.68 over runs). The property does not depend on the draw, since no start state was passed at all.

I reran 01 and 07 on the fixed tree in one process, and they reproduce the fix agent's after output. Both julia_live tests pass (2 passed in 109.31s). I also checked edges the fix agent did not measure, and all behave: get_gs(best=True,n=2), which now ends in a marked set_initial_wf; on a non-Hermitian julia chain, set_gs(x) (taken at <x|H|x>) and set_initial_wf_guess(x) (swept to the NH energy, x unmoved); set_initial_wf(None) after set_initial_wf(x) (dropped, solve runs); and get_excited_states after set_gs, which matches ED.

Three corrections to the record. (1) docs_needed lists "gs_energy(reconverge=...) without wf0= on julia_live raises TypeError" as a behaviour change, and it is not one. On a fresh chain the parent raised TypeError too, since get_gs_dmrg's signature has no reconverge. On a solved chain both trees return the stored e0 silently through Many_Body_Chain.gs_energy's gs_is_current short circuit (measured -1.61602540). Strike it from the list, or word it as "on a chain with no current state". (2) mpsjulialive/dynamics._min_energy re-solves the chain at full maxm/nsweeps on every KPM call of a supplied state, uncached. This is cost, not correctness. (3) The excited-state run printed an ITensorMPS deprecation warning (`inner(x::MPS, A::MPO, y::MPS)` with mismatched site indices), a warn_once from excited_states_dmrg's own Julia code. I did not spend a parent JIT run on attributing it. The energies match ED.

The left-open julia solver key (a maxm ramp returns the first energy) is gs_is_current's documented choice for julia_live and does not count against this fix.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
Fixed tree, 09_julia_edges.after.out (runs 01 and 07 verbatim, then edges):
solved e0 = -1.000000, <Sz_tot> of the solved state = -0.4993
set_initial_wf_guess(target): gs_energy = -1.000000, |<target|gs>|^2 = 1.0000
set_initial_wf(other):         gs_energy = -1.000000, |<other|gs>|^2 = 1.0000
set_gs(x):                     gs_energy = -0.974704, <x|H|x> = -0.974704, |<x|gs>|^2 = 1.0000
  vev(Sz0 Sz1) = -0.151774, <x|Sz0 Sz1|x> = -0.151774
  KPM (Sz0,Sz1) after set_gs: int C = -0.151490, e0 = -0.974704, |<x|gs>|^2 = 1.0000
set_initial_wf_guess(x): gs_energy = -1.000000; the caller's x after: <H> -0.974704406321 -> -0.974704406321, |<keep|x>|^2 = 1.000000000000
gs_energy(wf0=target) on a solved chain: -1.000000, |<target|gs>|^2 = 1.0000
fresh chain, set_initial_wf(random x): gs_energy = -0.282855, <x|H|x> = -0.282855, |<x|gs>|^2 = 1.0000
julia_live after set_gs(|1>): e0 = -0.850000, max|y - y_ED| = 8.187e-04 of the 6.130 peak
  with the window anchored on e0 instead: max|y - y_ED| = 1.278e+00 of the peak
(a) get_gs(best=True,n=2): gs_energy = -1.61602540, <wf|H|wf> = -1.61602540, |<wf|gs>|^2 = 1.000000 (exact 4-site -1.6160254)
(b) NH set_gs(x): gs_energy = (-1.526676+0j), <x|H|x> = (-1.526676+0j), NH e0 = (-1.596396+0j), |<x|gs>|^2 = 1.000000
(b) NH set_initial_wf_guess(x): gs_energy = (-1.596396-0j) (NH e0 (-1.596396+0j)), caller's x kept: 1.000000000000
(c) solved chain, gs_energy(reconverge=True) returned -1.61602540 (e0 -1.61602540)
(c) fresh chain gs_energy(reconverge=True) raised TypeError: gs_energy(reconverge=...) on julia_live needs a state to take or sweep from, given as wf0=
(d) set_initial_wf(x); set_initial_wf(None); gs_energy() = -1.61602540, <x|H|x> = 0.11623052
(e) excited states after set_gs(gs): [-1.616025 -0.957107 -0.957107], ED [-1.616025 -0.957107 -0.957107]
pytest -k julia_live: "2 passed, 20 deselected in 109.31s (0:01:49)"
Parent, by reading (parent/src/dmrgpy): groundstate.py:122 `if getattr(self,"_session",None) is None: self._gs_injected = None; return False`; groundstate.py:443-455 `e0,wf0 = get_gs_dmrg(self,**kwargs)` in both julia branches; mpsjulialive/groundstate.py:9-15 `if self.computed_gs: return self.wf0 ... if wf0 is None: psi0 = self.random_state()`; mpsjulialive/mps.py has no `mode`.
Files: <scratch>/review/session/09_julia_edges.py, 09_julia_edges.after.out, pytest_julia.out
````

**Left open**: julia_live still records no solver key, so a maxm ramp on one chain returns the first energy every time, measured by 08_julia_solver_key.py on an 8-site Heisenberg chain: python -3.186822, -3.371801, -3.374933 at maxm = 2, 4, 16 against julia_live -3.194321 three times. That is gs_is_current()'s documented choice for julia_live and was left as it is; it contradicts documentation.md's 'a convergence ramp over sc.maxm ... returned the first energy every time. The key now lives in groundstate.solver_key()' for this backend. `mpsjulialive/dynamics._min_energy()` re-solves the chain at its full maxm/nsweeps on every KPM call of a supplied state, uncached, which is cost rather than correctness (the reviewer's note). The non-Hermitian julia_live start (get_gs_nhdmrg) now receives the guess as well; the reviewer measured set_gs(x) taken at <x|H|x> and set_initial_wf_guess(x) swept to the NH energy with x unmoved there.

### 5. `get_gs(wf0=x)` on a chain whose ground state is current returned the stored state without reading x: with `reconverge=False` it returned the solved state (e0 = -2.5231886435, |<w|x>|^2 = 0.0107) instead of x itself (<x|H|x> = 0.0388284989), the answer `gs_energy` gives for the same keywords

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `construction`

**Status**: FIXED. The short circuit in Many_Body_Chain.get_gs now reads `groundstate.gs_is_current(self) and kwargs.get("wf0") is None`, the condition gs_energy() uses, so the two entry points cannot disagree on wf0= again; a repeated get_gs() with no start state still returns the stored object with no session call. The state reaches the session, so every later reader on the chain (get_gs(), gs_energy(), vev(), the correlators) measures x, which the reviewer confirmed on v2, v3 and "python". Pinned by tests/test_audit_2026_09_25_construction.py: test_get_gs_takes_wf0_as_it_is_on_a_current_chain and test_get_gs_warm_starts_from_wf0_on_a_current_chain on "python" and v3, and test_get_gs_without_wf0_still_returns_the_stored_state, which passes on both trees and pins the kept short circuit. NUMBERS CHANGE for get_gs(wf0=x, reconverge=False) on a chain with a current ground state, and for every reader after it: on the 6-site open S=1/2 Heisenberg chain with 0.3*Sz_0, itensor_version="python", maxm=20, nsweeps=10, np.random.seed(1), the returned state goes from the solved one (|<w|x>|^2 = 0.010694, e0 = -2.5231886435) to x itself (|<w|x>|^2 = 1, e0 = <x|H|x> = 0.0388284989). get_gs(wf0=x) with reconverge unspecified now sweeps from x instead of returning the stored object, and lands on the same ground state (e0 -2.5231886435 and |<w|gs>|^2 = 1 to the ten printed digits). get_gs(best=True, wf0=x) on a current chain now raises TypeError (best_gs() takes no wf0=), where it returned the stored MPS; on a fresh chain it raised the same on both trees. The v3 lines draw x from ITensor's own generator, which numpy's seed does not reach, so their x and the [3, maxm=3] energy in its ninth digit differ run to run; the "python" lines are the reproducible ones. No in-tree number moves.

**Recorded in**: 2026-09-24c New leads. Recorded by reading only in the 'New leads, not reviewed' of docs/audit_2026_09_24c_hole_hunt.md (third bullet); the shape 867e2b4 fixed for gs_energy(wf0=x) only (docs/audit_2026_09_24b_hole_hunt.md finding 12's Status)

**Where**: src/dmrgpy/manybodychain.py:1267-1269 on 8dd2198 (`if groundstate.gs_is_current(self): return self.wf0`, ahead of every keyword)

get_gs() opened its DMRG branch with `if groundstate.gs_is_current(self): return self.wf0`, before any keyword was read, so on a solved chain wf0= was never looked at: get_gs(wf0=x, reconverge=False), which takes x as it is, returned the solved state, and get_gs(wf0=x), which sweeps from x, returned the stored object without sweeping. Every later reader then measured the solved state rather than x (the reviewer's v2 run: vev(Sz0) = -0.1893606778 against <x|Sz0|x> = -0.1777307241, gs_energy() = E0). 867e2b4 fixed exactly this shape in gs_energy() and left get_gs() on the old condition, so the two public entry points to the same solver disagreed on the same keywords; on a chain with nothing stored both read x, which is the path every test and example takes (the two in-tree callers, examples/spin_models/anisotropic_heisenberg_S1 and examples/magnetization/magnetization_distribution_maxent, call set_hamiltonian() right before get_gs(wf0=wf), so their state is never current there). The reviewer reproduced it on v2, v3 and "python". Of the other keywords get_gs() forwards to gs_energy(), maxde= changes the answer too, and gs_energy() itself returns its stored energy ahead of it (the [3, maxm=3] lines, -3.3468165390 unrefined against an exact -3.4061631313, identical through get_gs and gs_energy); since the rule is to treat the keywords the way gs_energy() does, get_gs() now bypasses on wf0= only and the maxde= half is left to gs_energy()'s owner. With reconverge unspecified the fix makes x a warm start, meaning that the call sweeps from x and replaces the stored state with the result, which lands on the same ground state.

**Expected**: get_gs(**kwargs) returns the stored state under exactly the condition gs_energy(**kwargs) returns the stored energy, a current state and no wf0=, so get_gs(wf0=x, reconverge=False) returns x with e0 = <x|H|x> and get_gs(wf0=x) sweeps from x, on a solved chain as on a fresh one.

Repro (<scratch>/construction/02_get_gs_wf0.py (before: DMRGPY_SRC=<parent>/src run3.sh 02_get_gs_wf0.py | tee 02_get_gs_wf0.before.out; after: run3.sh 02_get_gs_wf0.py | tee 02_get_gs_wf0.after.out, rerun on the repaired tree)):

````python
# get-gs-wf0: get_gs(wf0=x) on a chain whose ground state is current
# returns the stored state without reading x.
import functools
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print = functools.partial(print, flush=True)
print("dmrgpy from", dmrgpy.__file__)


def heisenberg(v, n=6):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.3*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    return sc, h


def ov(sc, a, b):
    return abs(sc.overlap(a, b))**2/abs(sc.overlap(a, a)*sc.overlap(b, b))


def energy(sc, h, w):
    return np.real(sc.vev(h, wf=w))/np.real(sc.overlap(w, w))


for v in (3, "python"):
    np.random.seed(1)
    # the reference: gs_energy(wf0=x, reconverge=False) on a solved chain
    # takes x as it is (867e2b4)
    sc, h = heisenberg(v)
    sc.gs_energy()
    x = sc.random_state()
    print("[%s] gs_energy(wf0=x, reconverge=False) = %.10f   <x|H|x> = %.10f"
          % (v, sc.gs_energy(wf0=x, reconverge=False), energy(sc, h, x)))
    # the call under test, on a chain whose state is current
    sc, h = heisenberg(v)
    e0 = sc.gs_energy()
    gs = sc.get_gs()
    x = sc.random_state()
    print("[%s] E0 = %.10f   <x|H|x> = %.10f   |<x|gs>|^2 = %.3e"
          % (v, e0, energy(sc, h, x), ov(sc, x, gs)))
    w = sc.get_gs(wf0=x, reconverge=False)
    print("[%s] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = %.6f"
          "  |<w|gs>|^2 = %.6f  <w|H|w> = %.10f  sc.e0 = %.10f"
          % (v, ov(sc, w, x), ov(sc, w, gs), energy(sc, h, w), sc.e0))
    # with reconverge unspecified x is a warm start: a sweep from it
    sc, h = heisenberg(v)
    sc.gs_energy()
    gs = sc.get_gs()
    x = sc.random_state()
    w = sc.get_gs(wf0=x)
    print("[%s] get_gs(wf0=x): returned the stored object: %s"
          "  |<w|gs>|^2 = %.10f  sc.e0 = %.10f"
          % (v, w is gs, ov(sc, w, gs), sc.e0))
    # and on a chain whose state is not current, for comparison
    sc, h = heisenberg(v)
    x = sc.random_state()
    w = sc.get_gs(wf0=x, reconverge=False)
    print("[%s] same call, nothing stored:   |<w|x>|^2 = %.6f"
          "  <x|H|x> = %.10f  sc.e0 = %.10f"
          % (v, ov(sc, w, x), energy(sc, h, x), sc.e0))

# the keywords get_gs forwards besides wf0, on a solved chain at maxm=3:
# maxde= asks for a refinement
sc, h = heisenberg(3, n=8)
sc.maxm = 3
e = sc.gs_energy()
print("[3, maxm=3] E0 = %.10f   exact %.10f" % (e, sc.gs_energy(mode="ED")))
w = sc.get_gs(maxde=1e-4)
print("[3, maxm=3] get_gs(maxde=1e-4): sc.e0 = %.10f" % sc.e0)
print("[3, maxm=3] gs_energy(maxde=1e-4)    = %.10f" % sc.gs_energy(maxde=1e-4))
````

Observed, before:

````
[run3] slot 2 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
[3] gs_energy(wf0=x, reconverge=False) = -0.0789370491   <x|H|x> = -0.0789370491
[3] E0 = -2.5231886435   <x|H|x> = -0.0920979907   |<x|gs>|^2 = 6.527e-04
[3] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 0.000653  |<w|gs>|^2 = 1.000000  <w|H|w> = -2.5231886435  sc.e0 = -2.5231886435
[3] get_gs(wf0=x): returned the stored object: True  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
[3] same call, nothing stored:   |<w|x>|^2 = 1.000000  <x|H|x> = -0.1797510522  sc.e0 = -0.1797510522
[python] gs_energy(wf0=x, reconverge=False) = -0.0534127198   <x|H|x> = -0.0534127198
[python] E0 = -2.5231886435   <x|H|x> = 0.0388284989   |<x|gs>|^2 = 1.069e-02
[python] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 0.010694  |<w|gs>|^2 = 1.000000  <w|H|w> = -2.5231886435  sc.e0 = -2.5231886435
[python] get_gs(wf0=x): returned the stored object: True  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
[python] same call, nothing stored:   |<w|x>|^2 = 1.000000  <x|H|x> = 0.1522104199  sc.e0 = 0.1522104199
[3, maxm=3] E0 = -3.3468165390   exact -3.4061631313
[3, maxm=3] get_gs(maxde=1e-4): sc.e0 = -3.3468165390
[3, maxm=3] gs_energy(maxde=1e-4)    = -3.3468165390
````

Observed, after:

````
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
[3] gs_energy(wf0=x, reconverge=False) = -0.1601793432   <x|H|x> = -0.1601793432
[3] E0 = -2.5231886435   <x|H|x> = -0.1141279816   |<x|gs>|^2 = 3.049e-02
[3] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 1.000000  |<w|gs>|^2 = 0.030494  <w|H|w> = -0.1141279816  sc.e0 = -0.1141279816
[3] get_gs(wf0=x): returned the stored object: False  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
[3] same call, nothing stored:   |<w|x>|^2 = 1.000000  <x|H|x> = -0.2432283276  sc.e0 = -0.2432283276
[python] gs_energy(wf0=x, reconverge=False) = -0.0534127198   <x|H|x> = -0.0534127198
[python] E0 = -2.5231886435   <x|H|x> = 0.0388284989   |<x|gs>|^2 = 1.069e-02
[python] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 1.000000  |<w|gs>|^2 = 0.010694  <w|H|w> = 0.0388284989  sc.e0 = 0.0388284989
[python] get_gs(wf0=x): returned the stored object: False  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
[python] same call, nothing stored:   |<w|x>|^2 = 1.000000  <x|H|x> = 0.1522104199  sc.e0 = 0.1522104199
[3, maxm=3] E0 = -3.3468165360   exact -3.4061631313
[3, maxm=3] get_gs(maxde=1e-4): sc.e0 = -3.3468165360
[3, maxm=3] gs_energy(maxde=1e-4)    = -3.3468165360
````

**NUMBERS CHANGE**: YES. get_gs(wf0=x, reconverge=False) on a solved chain, 6-site open Heisenberg + 0.3*Sz_0, itensor_version="python", maxm=20, nsweeps=10, np.random.seed(1): returned state and e0 go from the solved state (e0 -2.5231886435, |<w|x>|^2 0.010694) to x (e0 = <x|H|x> = 0.0388284989, |<w|x>|^2 1), and every later reader (gs_energy(), vev, correlators) follows. get_gs(wf0=x) now sweeps from x; same energy to 10 digits. get_gs(best=True, wf0=x) on a current chain now raises TypeError (returned the stored MPS). No in-tree caller moves.

**Tests**: `tests/test_audit_2026_09_25_construction.py::test_get_gs_takes_wf0_as_it_is_on_a_current_chain`; `tests/test_audit_2026_09_25_construction.py::test_get_gs_warm_starts_from_wf0_on_a_current_chain`; `tests/test_audit_2026_09_25_construction.py::test_get_gs_without_wf0_still_returns_the_stored_state`

**Reviewer (CONFIRMED, fix HOLDS)**: Reproduced on the parent tree on every session backend, v2 included, which the fix agent did not run. On a solved chain, get_gs(wf0=x, reconverge=False) returns the solved state: |<w|x>|^2 = 8.959e-05 on v2, 0.0334 on v3 and 0.0123 on python, with e0 = -2.5231886435. Every later reader then measures the solved state rather than x: vev(Sz0) = -0.1893606778 against <x|Sz0|x> = -0.1777 (v2) and 0.0539 (python), gs_energy() = E0, and the KPM ZZ peak 0.1626. get_gs(wf0=x) with no reconverge argument returns the stored object without a sweep. gs_energy(wf0=x) reads x on the same chain, so the two entry points disagree, which is the shape 867e2b4 fixed in gs_energy only. The maxde= sub-claim, recorded as a lead for the session cluster, also holds: gs_energy(maxde=1e-4) at maxm=3 on 8 sites gives -3.4061631313 on a fresh chain and -3.3468165405 on a current one. Nothing struck.

The condition is now exactly gs_energy()'s: gs_is_current(self) and kwargs.get('wf0') is None. On the working tree, on v2, v3 and python, get_gs(wf0=x, reconverge=False) on a current chain returns x (|<w|x>|^2 = 1.00000000) with e0 = <x|H|x>. Every downstream reader follows: get_gs() is w, gs_energy() = <x|H|x>, vev(Sz0) = <x|Sz0|x> to all printed digits, and the KPM peak moves, so the state reached the session. get_gs(wf0=x) sweeps from x and lands on E0. The kept short circuit is pinned by a test that passes on both trees. get_gs(wf0=stored, reconverge=False) now returns a copy rather than the same object, with e0-E0 ~1e-15, which does not matter.

Sharpenings for the record:
(1) NUMBERS CHANGE should say that every subsequent reader on the chain now measures x (vev, gs_energy, correlators), not only the returned object.
(2) get_gs(best=True, wf0=x) on a current chain returned the stored MPS on the parent and now raises TypeError: best_gs() got an unexpected keyword argument 'wf0'. A fresh chain raised the same on both trees. This belongs in 'these now raise', not only in left_open.

The maxde= bypass is correctly left to gs_energy's owner. If that owner widens the condition, get_gs must widen identically; the fix agent's suggestion of one shared helper is the right shape.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
02_get_gs_wf0.before.out (python):
```
[python] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 0.010694  |<w|gs>|^2 = 1.000000  <w|H|w> = -2.5231886435  sc.e0 = -2.5231886435
[python] get_gs(wf0=x): returned the stored object: True  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
```
02_get_gs_wf0.after.out:
```
[python] get_gs(wf0=x, reconverge=False): |<w|x>|^2 = 1.000000  |<w|gs>|^2 = 0.010694  <w|H|w> = 0.0388284989  sc.e0 = 0.0388284989
[python] get_gs(wf0=x): returned the stored object: False  |<w|gs>|^2 = 1.0000000000  sc.e0 = -2.5231886435
```
05_attack_gs_rootn.before.out:
```
[2] E0=-2.5231886435  get_gs(wf0=x,rc=False): |<w|x>|^2=0.00008959  e0=-2.5231886435  <x|H|x>=0.4503481356
[2]   then get_gs() is w: True  |<get_gs()|x>|^2=0.00008959  gs_energy()=-2.5231886435  vev(Sz0)=-0.1893606778  <x|Sz0|x>=-0.1777307241
[2] current chain: get_gs(best=True, wf0=x)                  -> MPS
[3, maxm=3, fresh] gs_energy(maxde=1e-4) = -3.4061631313   maxm after = 3
[3, maxm=3, current] gs_energy(maxde=1e-4) = -3.3468165392   maxm after = 3
```
05_attack_gs_rootn.after.out:
```
[2] E0=-2.5231886435  get_gs(wf0=x,rc=False): |<w|x>|^2=1.00000000  e0=0.5071157023  <x|H|x>=0.5071157023
[2]   then get_gs() is w: True  |<get_gs()|x>|^2=1.00000000  gs_energy()=0.5071157023  vev(Sz0)=0.0467503506  <x|Sz0|x>=0.0467503506
[3]   then get_gs() is w: True  |<get_gs()|x>|^2=1.00000000  gs_energy()=0.1501866848  vev(Sz0)=-0.0336798185  <x|Sz0|x>=-0.0336798185
[python]   then get_gs() is w: True  |<get_gs()|x>|^2=1.00000000  gs_energy()=-0.0534127198  vev(Sz0)=0.0538021571  <x|Sz0|x>=0.0538021571
[2] current chain: get_gs(best=True, wf0=x)                  -> TypeError: best_gs() got an unexpected keyword argument 'wf0'
[3, maxm=3, current] gs_energy(maxde=1e-4) = -3.3468165405   maxm after = 3
```
The parent-tree test run fails test_get_gs_takes_wf0_as_it_is_on_a_current_chain[python,v3] and test_get_gs_warm_starts_from_wf0_on_a_current_chain[python,v3]; all pass on the working tree.
````

**Left open**: gs_energy(maxde=...) on a current chain returns the stored energy without refining (the [3, maxm=3] lines, and the reviewer's -3.3468165405 on a current chain against -3.4061631313 on a fresh one), the same shape as wf0=, in gs_energy's own short circuit, which belongs to the session cluster; get_gs follows gs_energy's condition and does not bypass on maxde= by itself, so if that condition widens, get_gs's must widen with it. get_gs(best=True, **kwargs) calls groundstate.best_gs(self, n=n, **kwargs), which takes no other keyword, so any keyword next to best=True raises TypeError (pre-existing on a fresh chain).

### 6. After `gs_energy_generalized`, the next dynamical correlator re-solves a plain ground state over the generalized one and measures that, on a Hermitian chain whose Hamiltonian-send cache is empty (|<wg|wf0>|^2 = 0.60, <Sz0> -0.48 to 0 and a (Sz0,Sz1) sum rule of -0.2225 against the generalized state's -0.1988, 6-site chain) and on a non-Hermitian chain whether or not it was solved first (|<wg|wf0>|^2 = 0.6877, e0 from lambda = -2.17581-0.213083j to -1.596396, <Sz0> -0.4529 to 0, 4-site chain), on `"python"` and v3; and after a plain NH-DMRG solve the first correlator solves NH-DMRG a second time

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `session`

**Status**: FIXED. On the Hermitian route `gs_energy_generalized()` sends the Hamiltonian through `send_hamiltonian()`, so the send cache records it and an identical re-send is skipped; nothing is lost by the skip, since both sessions drop their own energy cache after the generalized solve anyway and the band edges they keep are still the Hamiltonian's. On the non-Hermitian route, which the first pass had missed and the reviewer found, `nhdmrg.gs_energy_generalized_nhdmrg()` and `nhdmrg.gs_energy_nhdmrg()` now end with `_record_hamiltonian_sent()`, which puts H on the session through the same send cache on the session backends (julia_live has neither), so `ground_state_on_session()` finds the solve's state current; no NH reader needs the session to hold more, since NH-KPM hands the session its own terms and states and CVM_explicit reads `self.wf0`. A plain NH solve for an `H=` other than the chain's own Hamiltonian is stored as before but not recorded as sent. The CAVEAT in `gs_energy_generalized()`'s docstring now says that a correlator reads the generalized state on every chain, and that the remedy is `restart()` and then `gs_energy()`, since a bare `gs_energy()` returns lambda (on 8dd2198 too, measured by 10_generalized_edges.py). Pinned by `tests/test_audit_2026_09_25_session.py::test_correlator_after_gs_energy_generalized_reads_the_generalized_state` (`fresh` and `solved` on `python` and `v3`, the `fresh` rows failing against 8dd2198), `::test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state` (`fresh` and `solved` on `python` and `v3`, all four failing against 8dd2198) and `::test_first_nh_correlator_does_not_resolve_a_solved_chain` (`python`, `v3`, failing against 8dd2198). NUMBERS CHANGE for every reader after a correlator that followed `gs_energy_generalized()`: on a Hermitian chain never solved before (6-site Heisenberg chain, A = 1 + 0.8 Sz0), `e0` goes from -2.493577 to -3.597994 (lambda), `vev(Sz0)` from 0.0000 to -0.4801, and the integral of the (Sz0,Sz1) KPM spectrum from -0.2225 to -0.1988; on a non-Hermitian chain with no correlator run before the generalized solve, solved first or not (4-site Heisenberg chain plus 0.3j Sz0 + 0.2 Sx1, same A), `e0` after an NH-KPM goes from -1.596396 to lambda = -2.17581-0.213083j, |<wg|wf0>|^2 from 0.6877 to 1.0000 and `vev(Sz0)` from 0 to -0.4529, on `python` and v3. A Hermitian chain solved first and a non-Hermitian one on which a correlator had run are unchanged. Behaviour changes without number changes: the first correlator after a plain NH solve makes no second NH-DMRG solve (1 to 0), its NH-KPM spectrum unchanged to six digits at es = 0, 2, 4; `gs_energy_generalized()` on an MPO Hamiltonian on v3 runs (-3.59799449, the MultiOperator route's value) where it raised `AttributeError`.

**Recorded in**: 2026-09-24c New leads. The bare re-send in gs_energy_generalized is older than the send cache, and the NH entry points never sent H at all; the re-solve both cause is ground_state_on_session's, at the top of every correlator since 867e2b4.

**Where**: src/dmrgpy/groundstate.py:610 on 8dd2198 (a bare self._session.set_hamiltonian that leaves _session_ham_cache as it was); src/dmrgpy/nhdmrg.py gs_energy_nhdmrg and gs_energy_generalized_nhdmrg on 8dd2198 (neither puts H on the session, since NH-DMRG hands the session its terms per call); both read by groundstate.ground_state_on_session (hamiltonian_on_session False, so computed_gs is reset and a plain solve runs) at dynamics.py:210 and kpmdmrg.py:173

gs_energy_generalized() put H on the session itself, bypassing send_hamiltonian(), so the chain's send cache still said whatever it said before. On a Hermitian chain that was solved first the cache already matched, and the correlator read the generalized state, which is what the function's own CAVEAT promises; on a chain never solved the cache was empty, ground_state_on_session() found the stored state current but its Hamiltonian not on the session, reset computed_gs and solved a plain ground state from scratch, which then replaced wf0 and e0 for every later reader. On a non-Hermitian Hamiltonian gs_energy_generalized() dispatches to nhdmrg.gs_energy_generalized_nhdmrg() before that line, and neither it nor a plain NH solve sends H, so solving first did not help: every chain on which no correlator had run before the generalized solve lost the generalized state to a plain NH-DMRG re-solve, and after a plain NH solve the first correlator solved NH-DMRG a second time, which costs a solve and replaces the pair with an equivalent one.

**Expected**: Every reader after gs_energy_generalized() sees the generalized state and lambda, as the CAVEAT says, whether or not the chain had been solved before, on a Hermitian and a non-Hermitian Hamiltonian alike, and the first correlator after a plain NH solve reads that solve's pair.

Repro (<scratch>/session/02_generalized_cache.py (Hermitian route), 09_nh_generalized_cache.py (non-Hermitian route) and 10_generalized_edges.py (the CAVEAT's remedy and an MPO Hamiltonian), each with .before.out/.after.out):

````python
# 02_generalized_cache.py
# generalized-cache: after gs_energy_generalized(A) on a chain whose
# Hamiltonian-send cache is empty, does the next correlator read the
# generalized state, as the CAVEAT says, or re-solve a plain ground state
# over it? 6-site Heisenberg chain, metric A = 1 + 0.8*Sz0 (eigenvalues 0.6
# and 1.4, positive definite), whose generalized state has <Sz0> != 0.
import io, contextlib
import numpy as np
import dmrgpy
from dmrgpy import spinchain, cppext
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()

def heis(version, n=6):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    return sc

def fid(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))

ES = np.linspace(-6.0, 6.0, 601)
for version in ("python", 3):
    if version == 3 and not cppext.available(3): continue
    for solved_first in (False, True):
        np.random.seed(3)
        sc = heis(version)
        if solved_first: quiet(sc.gs_energy)
        A = 1 + 0.8*sc.Sz[0]
        lam = quiet(lambda: sc.gs_energy_generalized(A))
        wg = sc.wf0.copy()
        sz_g = float(np.real(wg.dot(sc.Sz[0]*wg)/wg.dot(wg)))
        zz_g = float(np.real(wg.dot(sc.Sz[0]*(sc.Sz[1]*wg))/wg.dot(wg)))
        es, d = quiet(lambda: sc.get_dynamical_correlator(
            name=(sc.Sz[0], sc.Sz[1]), submode="KPM", es=ES, delta=0.2))
        sumrule = float(np.real(np.trapezoid(d, es)))
        print("%-6s solved first=%-5s lambda=%.6f  after the correlator: "
              "e0=%.6f |<wg|wf0>|^2=%.4f vev(Sz0)=%+.4f (generalized %+.4f)  "
              "int C(Sz0,Sz1)=%+.4f (generalized <Sz0 Sz1>=%+.4f)"
              % (version, solved_first, lam, sc.e0, fid(sc.wf0, wg),
                 np.real(sc.vev(sc.Sz[0])), sz_g, sumrule, zz_g))

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh 02_generalized_cache.py 2>&1 | tee 02_generalized_cache.before.out
# cd $WF/session && ../run3.sh -u 02_generalized_cache.py 2>&1 | tee 02_generalized_cache.after.out

# 09_nh_generalized_cache.py
# generalized-cache, the non-Hermitian route (the reviewer's 11_nh_generalized.py,
# with a count of NH-DMRG solves added). On a non-Hermitian chain
# gs_energy_generalized() goes to nhdmrg.gs_energy_generalized_nhdmrg(), and
# neither it nor a plain NH solve puts H on the session, so the next correlator's
# ground_state_on_session() may find H missing and re-solve plain NH-DMRG over the
# generalized state. Three histories on the same chain:
#   fresh:        gs_energy_generalized(A), NH-KPM
#   solved:       gs_energy(), gs_energy_generalized(A), NH-KPM
#   solved+corr:  gs_energy(), NH-KPM, gs_energy_generalized(A), NH-KPM
# and a plain history, gs_energy() then NH-KPM, counting the NH-DMRG solves the
# correlator makes. 4-site Heisenberg chain plus 0.3j*Sz0 + 0.2*Sx1, A = 1 + 0.8*Sz0.
import io, contextlib, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain, cppext, groundstate, nhdmrg
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()

def fid(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))

def nh_chain(version, n=4):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3j*sc.Sz[0] + 0.2*sc.Sx[1])
    sc.maxm, sc.nsweeps = 20, 10
    return sc

solves = [0]
_nhdmrg = nhdmrg.nhdmrg
def counted(*args, **kwargs):
    solves[0] += 1
    return _nhdmrg(*args, **kwargs)
nhdmrg.nhdmrg = counted

KPM = dict(submode="KPM", es=np.linspace(0.0, 4.0, 5), delta=0.3, E_max=10.0, n=50)
for version in ("python", 3):
    if version == 3 and not cppext.available(3): continue
    for history in ("fresh", "solved", "solved+corr"):
        np.random.seed(9)
        sc = nh_chain(version)
        corr = lambda: quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **KPM))
        if history != "fresh": quiet(sc.gs_energy)
        if history == "solved+corr": corr()
        lam = quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
        wg = sc.wf0.copy()
        on = groundstate.hamiltonian_on_session(sc)
        corr()
        print("%-6s %-12s lambda=%s H on session after the generalized solve: %-5s "
              "after NH-KPM: e0=%s |<wg|wf0>|^2=%.4f vev(Sz0)=%s (generalized %s)"
              % (version, history, np.round(lam, 6), on, np.round(sc.e0, 6), fid(sc.wf0, wg),
                 np.round(sc.vev(sc.Sz[0]), 4), np.round(wg.dot(sc.Sz[0]*wg)/wg.dot(wg), 4)))
    np.random.seed(9)
    sc = nh_chain(version)
    solves[0] = 0
    e_nh = quiet(sc.gs_energy)
    after_solve = solves[0]
    es, d = quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **KPM))
    print("%-6s plain        gs_energy()=%s: NH-DMRG solves by gs_energy() %d, "
          "by the first NH-KPM after it %d; that NH-KPM at es=0,2,4: %s"
          % (version, np.round(e_nh, 6), after_solve, solves[0] - after_solve,
             np.round(d[::2], 6)))

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 09_nh_generalized_cache.py 2>&1 | tee 09_nh_generalized_cache.before.out
# cd $WF/session && ../run3.sh -u 09_nh_generalized_cache.py 2>&1 | tee 09_nh_generalized_cache.after.out

# 10_generalized_edges.py
# generalized-cache, two edges of the fix. (a) The CAVEAT's remedy: after
# gs_energy_generalized(A), does gs_energy() return lambda (the stored state is
# current) or the plain ground-state energy, and does restart() then
# gs_energy() give the plain one? Hermitian 6-site Heisenberg chain and the
# non-Hermitian 4-site chain, A = 1 + 0.8*Sz0. (b) An already-built MPO
# Hamiltonian (toMPO) on v3: gs_energy_generalized() against the MultiOperator
# route on the same chain.
import io, contextlib, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain, cppext
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()

def heis(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h

for version in ("python", 3):
    if version == 3 and not cppext.available(3): continue
    for nh in (False, True):
        np.random.seed(11)
        n = 4 if nh else 6
        sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
        h = heis(sc) + (0.3j*sc.Sz[0] + 0.2*sc.Sx[1] if nh else 0)
        sc.set_hamiltonian(h)
        sc.maxm, sc.nsweeps = 20, 10
        lam = quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
        e1 = quiet(sc.gs_energy)
        sc.restart()
        e2 = quiet(sc.gs_energy)
        print("%-6s %-13s lambda=%s  gs_energy() right after=%s  restart(), gs_energy()=%s"
              % (version, "non-Hermitian" if nh else "Hermitian", np.round(lam, 6),
                 np.round(e1, 6), np.round(e2, 6)))

if cppext.available(3):
    np.random.seed(12)
    sc = spinchain.Spin_Chain([2]*6, itensor_version=3)
    sc.maxm, sc.nsweeps = 20, 10
    h = heis(sc)
    sc.set_hamiltonian(h)
    lam_mo = quiet(lambda: sc.gs_energy_generalized(1 + 0.8*sc.Sz[0]))
    sc2 = spinchain.Spin_Chain([2]*6, itensor_version=3)
    sc2.maxm, sc2.nsweeps = 20, 10
    try:
        sc2.set_hamiltonian(sc2.toMPO(heis(sc2)))
        lam_mpo = quiet(lambda: sc2.gs_energy_generalized(1 + 0.8*sc2.Sz[0]))
        print("v3 MPO Hamiltonian: gs_energy_generalized = %.8f, MultiOperator route %.8f"
              % (lam_mpo, lam_mo))
    except Exception as err:
        print("v3 MPO Hamiltonian: gs_energy_generalized raised %s: %s (MultiOperator route %.8f)"
              % (type(err).__name__, str(err)[:200], lam_mo))

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 10_generalized_edges.py 2>&1 | tee 10_generalized_edges.before.out
# cd $WF/session && ../run3.sh -u 10_generalized_edges.py 2>&1 | tee 10_generalized_edges.after.out
````

Observed, before:

````
02_generalized_cache.before.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
python solved first=False lambda=-3.597994  after the correlator: e0=-2.493577 |<wg|wf0>|^2=0.6011 vev(Sz0)=+0.0000 (generalized -0.4801)  int C(Sz0,Sz1)=-0.2225 (generalized <Sz0 Sz1>=-0.1988)
python solved first=True  lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=False lambda=-3.597994  after the correlator: e0=-2.493577 |<wg|wf0>|^2=0.6011 vev(Sz0)=-0.0000 (generalized -0.4801)  int C(Sz0,Sz1)=-0.2225 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=True  lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)

09_nh_generalized_cache.before.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
python fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=(-0-0j) (generalized (-0.4529+0j))
python solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=-0j (generalized (-0.4529-0j))
python solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529+0j) (generalized (-0.4529+0j))
python plain        gs_energy()=(-1.596396-0j): NH-DMRG solves by gs_energy() 1, by the first NH-KPM after it 1; that NH-KPM at es=0,2,4: [-8.3803  -7.878347e+00j -1.372571+4.803210e-01j -0.02633 +1.240000e-04j]
3      fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396+0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=0j (generalized (-0.4529-0j))
3      solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=(-0-0j) (generalized (-0.4529+0j))
3      solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529+0j))
3      plain        gs_energy()=(-1.596396+0j): NH-DMRG solves by gs_energy() 1, by the first NH-KPM after it 1; that NH-KPM at es=0,2,4: [-8.3803  -7.878347e+00j -1.372571+4.803210e-01j -0.02633 +1.240000e-04j]

10_generalized_edges.before.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
python Hermitian     lambda=-3.597994  gs_energy() right after=-3.597994  restart(), gs_energy()=-2.493577
python non-Hermitian lambda=(-2.17581-0.213083j)  gs_energy() right after=(-2.17581-0.213083j)  restart(), gs_energy()=(-1.596396-0j)
3      Hermitian     lambda=-3.597994  gs_energy() right after=-3.597994  restart(), gs_energy()=-2.493577
3      non-Hermitian lambda=(-2.17581-0.213083j)  gs_energy() right after=(-2.17581-0.213083j)  restart(), gs_energy()=(-1.596396+0j)
v3 MPO Hamiltonian: gs_energy_generalized raised AttributeError: 'StaticOperator' object has no attribute 'to_terms' (MultiOperator route -3.59799449)
````

Observed, after:

````
02_generalized_cache.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
python solved first=False lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
python solved first=True  lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=False lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=True  lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)

09_nh_generalized_cache.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
python fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529+0j))
python solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529-0j))
python solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529-0j))
python plain        gs_energy()=(-1.596396-0j): NH-DMRG solves by gs_energy() 1, by the first NH-KPM after it 0; that NH-KPM at es=0,2,4: [-8.3803  -7.878347e+00j -1.372571+4.803210e-01j -0.02633 +1.240000e-04j]
3      fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529+0j))
3      solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529+0j))
3      solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529-0j))
3      plain        gs_energy()=(-1.596396+0j): NH-DMRG solves by gs_energy() 1, by the first NH-KPM after it 0; that NH-KPM at es=0,2,4: [-8.3803  -7.878347e+00j -1.372571+4.803210e-01j -0.02633 +1.240000e-04j]

10_generalized_edges.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
python Hermitian     lambda=-3.597994  gs_energy() right after=-3.597994  restart(), gs_energy()=-2.493577
python non-Hermitian lambda=(-2.17581-0.213083j)  gs_energy() right after=(-2.17581-0.213083j)  restart(), gs_energy()=(-1.596396-0j)
3      Hermitian     lambda=-3.597994  gs_energy() right after=-3.597994  restart(), gs_energy()=-2.493577
3      non-Hermitian lambda=(-2.17581-0.213083j)  gs_energy() right after=(-2.17581-0.213083j)  restart(), gs_energy()=(-1.596396+0j)
v3 MPO Hamiltonian: gs_energy_generalized = -3.59799449, MultiOperator route -3.59799449
````

**NUMBERS CHANGE**: Hermitian, 6-site Heisenberg chain, A = 1 + 0.8*Sz0, python and v3, never solved before gs_energy_generalized(): after a KPM correlator e0 -2.493577 -> -3.597994, vev(Sz0) 0.0000 -> -0.4801, integral of the (Sz0,Sz1) spectrum -0.2225 -> -0.1988; solved-first chains unchanged. Non-Hermitian, 4-site Heisenberg chain plus 0.3j*Sz0 + 0.2*Sx1, same A, python and v3, fresh or solved first: after an NH-KPM e0 -1.596396 -> lambda (-2.17581-0.213083j), |<wg|wf0>|^2 0.6877 -> 1.0000, vev(Sz0) 0 -> -0.4529; a chain with a correlator run before the generalized solve unchanged. No number change for a plain NH solve (NH-KPM at es=0,2,4 identical to six digits), only its second NH-DMRG solve by the first correlator goes (1 -> 0). gs_energy_generalized on a v3 MPO Hamiltonian: AttributeError -> -3.59799449 (MultiOperator route -3.59799449).

**Tests**: `tests/test_audit_2026_09_25_session.py::test_correlator_after_gs_energy_generalized_reads_the_generalized_state[fresh-python]`; `tests/test_audit_2026_09_25_session.py::test_correlator_after_gs_energy_generalized_reads_the_generalized_state[fresh-v3]`; `tests/test_audit_2026_09_25_session.py::test_correlator_after_gs_energy_generalized_reads_the_generalized_state[solved-python]`; `tests/test_audit_2026_09_25_session.py::test_correlator_after_gs_energy_generalized_reads_the_generalized_state[solved-v3]`; `tests/test_audit_2026_09_25_session.py::test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state[fresh-python]`; `tests/test_audit_2026_09_25_session.py::test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state[fresh-v3]`; `tests/test_audit_2026_09_25_session.py::test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state[solved-python]`; `tests/test_audit_2026_09_25_session.py::test_nh_correlator_after_gs_energy_generalized_reads_the_generalized_state[solved-v3]`; `tests/test_audit_2026_09_25_session.py::test_first_nh_correlator_does_not_resolve_a_solved_chain[python]`; `tests/test_audit_2026_09_25_session.py::test_first_nh_correlator_does_not_resolve_a_solved_chain[v3]`

**Reviewer (CONFIRMED, fix INCOMPLETE)**: Reproduced to every printed digit on the parent tree, on python and v3. On a chain whose send cache is empty, the correlator after gs_energy_generalized() re-solves the plain ground state: e0 -3.597994 -> -2.493577, |<wg|wf0>|^2 = 0.6011, vev(Sz0) -0.4801 -> 0, sum rule -0.2225 against the generalized -0.1988. A chain solved first reads the generalized state. Nothing in CLAUDE.md, the known_issue files or ROADMAP marks this as intended, and the function's own CAVEAT promises the opposite.

Narrowed: the fix holds on a Hermitian Hamiltonian, on python and v3, fresh and solved-first alike. My rerun matches the fix agent's after output to every digit, and the solved-first rows are unchanged between the trees. The comment's claim that nothing is lost by skipping an identical re-send is right for v3: Chain::gs_energy_generalized sets `have_wf0_energy_ = false; // stale` (mpscpp3/chain_session.h, 46 lines into the function), as pyitensor sets `_wf0_energy = None`.

The same public function on a non-Hermitian Hamiltonian still has the defect, unchanged by the fix. gs_energy_generalized() dispatches to gs_energy_generalized_nhdmrg() before the new send_hamiltonian() line. NH-DMRG never fills the send cache, so the next correlator's ground_state_on_session() finds H missing and re-solves plain NH-DMRG over the generalized state. This happens on the fresh chain and on the solved-first chain alike (a plain NH solve does not send H either), on both trees and on python and v3: |<wg|wf0>|^2 = 0.6877, e0 goes from lambda (-2.17581-0.213083j) to -1.596396, and <Sz0> from -0.4529 to 0. Only a chain that had run a correlator before the generalized solve keeps the generalized state. So the expected contract ("every reader after gs_energy_generalized() sees the generalized state and lambda ... whether or not the chain had been solved before") fails on that route, and solving first does not help there.

Two record consequences. The docs_needed CAVEAT sentence ("this holds whether or not the chain was solved before") needs "on a Hermitian Hamiltonian", or it is false as written. And the nh-injected left_open ("the first correlator on a solved NH chain re-solves NH-DMRG once; not wrong, a second solve") is wrong after a generalized NH solve, where that re-solve discards the state.

The suggested completion is send_hamiltonian(self) at the end of gs_energy_generalized_nhdmrg(), and in gs_energy_nhdmrg() too on the session backends, which would also remove the extra NH re-solve left open under nh-injected. send_hamiltonian() of a non-Hermitian H is already exercised on python and v3 by the new NH take path (_take_injected_state), which passes its tests. Alternatively, record the NH route as open.

The new tests pin the Hermitian route only: the fresh rows fail on 8dd2198 by the fix agent's run. I did not verify the note that an MPO (StaticOperator) Hamiltonian now goes through send_hamiltonian().

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
02 rerun, parent:
python solved first=False lambda=-3.597994  after the correlator: e0=-2.493577 |<wg|wf0>|^2=0.6011 vev(Sz0)=+0.0000 (generalized -0.4801)  int C(Sz0,Sz1)=-0.2225 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=False lambda=-3.597994  after the correlator: e0=-2.493577 |<wg|wf0>|^2=0.6011 vev(Sz0)=-0.0000 (generalized -0.4801)  int C(Sz0,Sz1)=-0.2225 (generalized <Sz0 Sz1>=-0.1988)
02 rerun, fixed:
python solved first=False lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
3      solved first=False lambda=-3.597994  after the correlator: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 (generalized -0.4801)  int C(Sz0,Sz1)=-0.1988 (generalized <Sz0 Sz1>=-0.1988)
(solved-first rows identical on both trees: e0=-3.597994 |<wg|wf0>|^2=1.0000 vev(Sz0)=-0.4801 int C=-0.1988)
11_nh_generalized.py, fixed tree (parent identical line for line):
python fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=(-0-0j) (generalized (-0.4529+0j))
python solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=-0j (generalized (-0.4529-0j))
python solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529+0j) (generalized (-0.4529+0j))
3      fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=0j (generalized (-0.4529-0j))
3      solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=-0j (generalized (-0.4529-0j))
3      solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529-0j) (generalized (-0.4529+0j))
Files: <scratch>/review/session/02.before.out, 02.after.out, 11_nh_generalized.py, 11.before.out, 11.after.out
````

**Left open**: The CAVEAT itself stands: a correlator after gs_energy_generalized() measures the generalized state with lambda as its origin, which is now at least the same answer on every chain. The MPO route was run on v3 only; on the other backends send_hamiltonian() raises NotImplementedError unless the session has set_hamiltonian_mpo, by reading, not run.

### 7. `Thermal_Spin_Chain.get_gs()` at T>1e-5 assigns `MBChain.wf0` and `MBChain.hamiltonian` behind the setters, so on a 3-site chain at T=1 `MBChain.gs_energy()` returns the singlet Hamiltonian's -2.25 against <wf|H|wf> = -0.4164, a KPM correlator on `"python"`, v2 and v3 re-solves the plain ground state over the annealed state (sum rule -0.1664, and <Sz0 Sz1> -0.1667 afterwards, against -0.0694), and on `mode="ED"` `MBChain.vev()` and the KPM sum rule read the singlet state (0.0 against -0.0694)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `session`

**Status**: FIXED. `Thermal_Spin_Chain.get_gs()` at T>1e-5 hands the physical Hamiltonian and the annealed state to `MBChain` through `set_hamiltonian()` and `set_gs()`, so the state is marked as set by hand, the send cache and the ED object are the physical Hamiltonian's, and the next read takes the state unswept with its own energy; the T<=1e-5 branch no longer writes a normalized copy back over `MBChain`'s own solved state. On `mode="ED"` `MBChain.gs_energy()` is the physical Hamiltonian's lowest eigenvalue, not <wf|H|wf>, which is the 2026-09-24c record's open `gs_energy(mode="ED")` after `set_gs`, left as it is. Pinned by `tests/test_audit_2026_09_25_session.py::test_every_reader_of_the_thermal_chain_measures_the_annealed_state` on `python`, `v2`, `v3` and `ED` (all four fail against 8dd2198; the ED row pins `vev()` and the KPM sum rule, not `gs_energy()`) and `::test_zero_temperature_thermal_chain_is_the_solved_ground_state`. NUMBERS CHANGE on `MBChain` after `tc.get_gs()` at T>1e-5, never on the state `tc.get_gs()` returns: on the 3-site Heisenberg chain at T=1, on `python`, v2 and v3, `MBChain.gs_energy()` goes from -2.250000 to -0.416388, the integral of the (Sz0,Sz1) KPM spectrum from -0.166370 to -0.069291 and `MBChain.vev(Sz0 Sz1)` after that correlator from -0.166667 to -0.069398; on `mode="ED"`, `MBChain.vev(Sz0 Sz1)` from 0.000000 to -0.069398, the KPM sum rule from 0.000000 to -0.069222 and `MBChain.gs_energy()` from -2.250000 to -1.000000, the lowest eigenvalue; at T=1.0001e-5 on `python`, `MBChain.gs_energy()` from -2.250000 to -0.999977, the sum rule from -0.166370 to -0.166366 and `vev` after it from -0.166667 to -0.166663. The exact Boltzmann value at T=1 is -0.071372, the difference being anneal()'s first-order steps.

**Recorded in**: 2026-09-24c New leads. The direct assignments predate 867e2b4; the correlator re-solve is ground_state_on_session's since 867e2b4.

**Where**: src/dmrgpy/thermal.py:59-62 on 8dd2198 (MBChain.wf0 = wf0 and MBChain.hamiltonian = self.hamiltonian after the singlet solve), read by groundstate.ground_state_on_session, Many_Body_Chain.gs_energy (stored e0 of the singlet solve, still current by its solver key) and the ED object built for the singlet Hamiltonian

The annealed purified state is built on MBChain after a solve of the singlet Hamiltonian, and then written onto MBChain by plain assignment together with the physical Hamiltonian. MBChain's bookkeeping still describes the singlet solve: computed_gs and the solver key say the stored state is current, e0 is the singlet energy, the send cache holds the singlet terms, and on mode="ED" the ED object is the singlet Hamiltonian's with the singlet ground state in it. So gs_energy() returns -2.25, vev() and the KPM correlator on ED read the singlet state, and on the DMRG backends the first correlator finds the physical Hamiltonian not on the session, re-solves a plain ground state and replaces the annealed state for every later reader. The same holds just above the threshold, at T=1.0001e-5. The state tc.get_gs() returns is right, which is what the existing example and tests read.

**Expected**: MBChain holds the physical Hamiltonian and the annealed state as a state set by hand: vev() and every correlator on MBChain measure the annealed state on every mode, and gs_energy() is <wf|H|wf> on the DMRG backends; on mode="ED" gs_energy() stays the physical Hamiltonian's lowest eigenvalue, which is how the 2026-09-24c record leaves gs_energy(mode="ED") after set_gs.

Repro (<scratch>/session/03_thermal_bypass.py (revised in this pass to add v2, the ED gs_energy and sum rule, and T=1.0001e-5; both trees rerun)):

````python
# 03_thermal_bypass.py
# thermal-bypass: Thermal_Spin_Chain.get_gs() at T>1e-5 puts the annealed
# purified state on MBChain by assigning wf0 and hamiltonian directly. Does a
# reader on MBChain see that state? 3-site S=1/2 Heisenberg chain at T=1,
# purified on 6 sites, on python, v2, v3 and mode="ED", and on python just
# above the T=1e-5 threshold. Reference: <wf|Sz0 Sz1|wf> on the annealed state
# (by hand) and the exact Boltzmann average from ED.
import io, contextlib
import numpy as np
import dmrgpy
from dmrgpy import spinchain, thermal, cppext
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()

def fid(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))

n, T = 3, 1.0
ref = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h + ref.Sx[i]*ref.Sx[i+1] + ref.Sy[i]*ref.Sy[i+1] + ref.Sz[i]*ref.Sz[i+1]
ref.set_hamiltonian(h)
ed = ref.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense())
ZZ = np.array(ed.MO2matrix(ref.Sz[0]*ref.Sz[1]).todense())
w, U = np.linalg.eigh(H)
p = np.exp(-(w-w[0])/T) ; p = p/p.sum()
exact = float(np.real(np.sum(p*np.diag(U.conj().T @ ZZ @ U))))
print("exact thermal <Sz0 Sz1> at T=%.1f: %+.6f" % (T, exact))

ES = np.linspace(-6.0, 6.0, 601)
for version, T in (("python", 1.0), (2, 1.0), (3, 1.0), ("ED", 1.0),
                   ("python", 1.0001e-5)):
    if version in (2, 3) and not cppext.available(version): continue
    np.random.seed(4)
    kw = dict(itensor_version=version) if version != "ED" else {}
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=T, **kw)
    ht = 0
    for i in range(n-1):
        ht = ht + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    tc.set_hamiltonian(ht)
    if version == "ED": tc.mode = "ED"
    mb = tc.MBChain
    if version != "ED": mb.maxm, mb.nsweeps = 30, 10
    wf = quiet(tc.get_gs)
    op = tc.Sz[0]*tc.Sz[1]
    by_hand = float(np.real(wf.dot(op*wf)/wf.dot(wf)))
    v0 = float(np.real(quiet(lambda: mb.vev(op))))
    line = "%-6s T=%-9.7g <wf|Sz0 Sz1|wf>=%+.6f  MBChain.vev=%+.6f" % (version, T, by_hand, v0)
    e_mb = float(np.real(quiet(mb.gs_energy)))
    e_wf = float(np.real(wf.dot(ht*wf)/wf.dot(wf)))
    line += "  MBChain.gs_energy()=%+.6f (<wf|H|wf>=%+.6f)" % (e_mb, e_wf)
    kwm = dict(mode="ED") if version == "ED" else {}
    es, d = quiet(lambda: mb.get_dynamical_correlator(
        name=(tc.Sz[0], tc.Sz[1]), submode="KPM", es=ES, delta=0.2, **kwm))
    v1 = float(np.real(quiet(lambda: mb.vev(op))))
    line += ("  after a KPM correlator: int C(Sz0,Sz1)=%+.6f, vev=%+.6f"
             % (float(np.real(np.trapezoid(d, es))), v1))
    if version != "ED": line += ", |<wf|MBChain.wf0>|^2=%.4f" % fid(mb.wf0, wf)
    print(line)

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 03_thermal_bypass.py 2>&1 | tee 03_thermal_bypass.before.out
# cd $WF/session && ../run3.sh -u 03_thermal_bypass.py 2>&1 | tee 03_thermal_bypass.after.out
````

Observed, before:

````
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
exact thermal <Sz0 Sz1> at T=1.0: -0.071372
python T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-2.250000 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.166370, vev=-0.166667, |<wf|MBChain.wf0>|^2=0.1032
2      T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-2.250000 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.166370, vev=-0.166667, |<wf|MBChain.wf0>|^2=0.4044
3      T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-2.250000 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.166370, vev=-0.166667, |<wf|MBChain.wf0>|^2=0.0363
ED     T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.000000  MBChain.gs_energy()=-2.250000 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.000000, vev=-0.000000
python T=1.0001e-05 <wf|Sz0 Sz1|wf>=-0.166663  MBChain.vev=-0.166663  MBChain.gs_energy()=-2.250000 (<wf|H|wf>=-0.999977)  after a KPM correlator: int C(Sz0,Sz1)=-0.166370, vev=-0.166667, |<wf|MBChain.wf0>|^2=0.0516
````

Observed, after:

````
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
exact thermal <Sz0 Sz1> at T=1.0: -0.071372
python T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-0.416388 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.069291, vev=-0.069398, |<wf|MBChain.wf0>|^2=1.0000
2      T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-0.416388 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.069291, vev=-0.069398, |<wf|MBChain.wf0>|^2=1.0000
3      T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-0.416388 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.069291, vev=-0.069398, |<wf|MBChain.wf0>|^2=1.0000
ED     T=1         <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398  MBChain.gs_energy()=-1.000000 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.069222, vev=-0.069398
python T=1.0001e-05 <wf|Sz0 Sz1|wf>=-0.166663  MBChain.vev=-0.166663  MBChain.gs_energy()=-0.999977 (<wf|H|wf>=-0.999977)  after a KPM correlator: int C(Sz0,Sz1)=-0.166366, vev=-0.166663, |<wf|MBChain.wf0>|^2=1.0000
````

**NUMBERS CHANGE**: 3-site Heisenberg chain, on MBChain after tc.get_gs(), tc.get_gs() itself unchanged. At T=1 on python, v2 and v3: gs_energy() -2.250000 -> -0.416388; integral of the (Sz0,Sz1) KPM spectrum -0.166370 -> -0.069291; vev(Sz0 Sz1) after it -0.166667 -> -0.069398. At T=1 on mode="ED": vev(Sz0 Sz1) 0.000000 -> -0.069398; KPM sum rule 0.000000 -> -0.069222; gs_energy() -2.250000 -> -1.000000 (the lowest eigenvalue, not <wf|H|wf> = -0.416388). At T=1.0001e-5 on python: gs_energy() -2.250000 -> -0.999977; sum rule -0.166370 -> -0.166366; vev after it -0.166667 -> -0.166663.

**Tests**: `tests/test_audit_2026_09_25_session.py::test_every_reader_of_the_thermal_chain_measures_the_annealed_state[python]`; `tests/test_audit_2026_09_25_session.py::test_every_reader_of_the_thermal_chain_measures_the_annealed_state[v3]`; `tests/test_audit_2026_09_25_session.py::test_every_reader_of_the_thermal_chain_measures_the_annealed_state[v2]`; `tests/test_audit_2026_09_25_session.py::test_every_reader_of_the_thermal_chain_measures_the_annealed_state[ED]`; `tests/test_audit_2026_09_25_session.py::test_zero_temperature_thermal_chain_is_the_solved_ground_state[python]`; `tests/test_audit_2026_09_25_session.py::test_zero_temperature_thermal_chain_is_the_solved_ground_state[v3]`

**Reviewer (CONFIRMED, fix HOLDS)**: Reproduced to every printed digit on python, v3 and ED, and also on v2, which the fix agent did not run. On the parent, v2 gives gs_energy -2.25 against <wf|H|wf> -0.416388 and a KPM sum rule of -0.166370. On ED the KPM sum rule on MBChain was -0.000000, the singlet state, as well as vev 0.0. The state tc.get_gs() returns was right on both trees.

Holds on python, v2, v3 and ED for vev() and the KPM correlator, and holds just above the threshold (T=1.0001e-5). The T=1e-5 branch is identical on both trees. examples/finite_temperature/thermal_purification_VS_exact, which asserts, prints byte-identical output on both trees and passes.

The expected statement ("gs_energy() is <wf|H|wf> ... on every mode") holds on the DMRG backends only. On mode="ED", MBChain.gs_energy() returns the physical Hamiltonian's lowest eigenvalue, -1.000000, not <wf|H|wf> = -0.416388. That is the 2026-09-24c record's left-open "gs_energy(mode=\"ED\") after set_gs still reports the lowest eigenvalue", so it is not a new defect, but it is a number that moved.

The NUMBERS CHANGE statement is incomplete. It should also list: on mode="ED", MBChain.gs_energy() -2.250000 -> -1.000000 (the lowest eigenvalue, not <wf|H|wf>) and the MBChain KPM sum rule -0.000000 -> -0.069222; v2 moves exactly as python and v3 do (gs_energy -2.250000 -> -0.416388, sum rule -0.166370 -> -0.069291, vev after it -0.166667 -> -0.069398); and at T=1.0001e-5 on python, gs_energy -2.250000 -> -0.999977.

The new test's ED row checks vev only, so the ED gs_energy is not pinned, which is consistent with the left-open choice.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
10_neighbours.py section A, parent:
2      T=1 <wf|Sz0Sz1|wf>=-0.069398 vev=-0.069398 <wf|H|wf>=-0.416388 gs_energy()=-2.250000  KPM int C=-0.166370, vev after=-0.166667
ED     T=1 <wf|Sz0Sz1|wf>=-0.069398 vev=-0.000000 <wf|H|wf>=-0.416388 gs_energy()=-2.250000  KPM int C=-0.000000, vev after=-0.000000
python T=1e-05 <wf|Sz0Sz1|wf>=-0.166667 vev=-0.166667 <wf|H|wf>=-1.000000 gs_energy()=-1.000000  KPM int C=-0.166370, vev after=-0.166667
python T=1.0001e-05 <wf|Sz0Sz1|wf>=-0.166663 vev=-0.166663 <wf|H|wf>=-0.999977 gs_energy()=-2.250000  KPM int C=-0.166370, vev after=-0.166667
fixed:
2      T=1 <wf|Sz0Sz1|wf>=-0.069398 vev=-0.069398 <wf|H|wf>=-0.416388 gs_energy()=-0.416388  KPM int C=-0.069291, vev after=-0.069398
ED     T=1 <wf|Sz0Sz1|wf>=-0.069398 vev=-0.069398 <wf|H|wf>=-0.416388 gs_energy()=-1.000000  KPM int C=-0.069222, vev after=-0.069398
python T=1e-05 <wf|Sz0Sz1|wf>=-0.166667 vev=-0.166667 <wf|H|wf>=-1.000000 gs_energy()=-1.000000  KPM int C=-0.166370, vev after=-0.166667
python T=1.0001e-05 <wf|Sz0Sz1|wf>=-0.166663 vev=-0.166663 <wf|H|wf>=-0.999977 gs_energy()=-0.999977  KPM int C=-0.166366, vev after=-0.166663
03 rerun fixed: "python ... MBChain.gs_energy()=-0.416388 (<wf|H|wf>=-0.416388)  after a KPM correlator: int C(Sz0,Sz1)=-0.069291, vev=-0.069398, |<wf|MBChain.wf0>|^2=1.0000" and "ED     <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.069398"; parent: "python ... MBChain.gs_energy()=-2.250000 ... int C(Sz0,Sz1)=-0.166370, vev=-0.166667, |<wf|MBChain.wf0>|^2=0.1032", "ED     <wf|Sz0 Sz1|wf>=-0.069398  MBChain.vev=-0.000000".
thermal_purification_VS_exact, both trees: "Purification thermal energy (T=1.00) = -0.6073708211696766 ... TEST PASSED" (all four T identical).
Files: <scratch>/review/session/03.before.out, 03.after.out, 10.before.out, 10.after.out, 12_thermal_example.before.out, 12_thermal_example.after.out
````

**Left open**: On mode="ED" MBChain.gs_energy() is the lowest eigenvalue (-1.000000), not <wf|H|wf> (-0.416388), the 2026-09-24c record's open item, not changed here. A dynamical correlator of the purified state with H acting on the physical sites only is not the thermal correlator (that needs H minus the ancilla copy), so only the static readers and the sum rule are meaningful here; not a defect of this fix.

### 8. After `set_gs(x)` on a non-Hermitian chain the NH-KPM pairs x with the left eigenvector of the last NH-DMRG solve and renormalizes by <psil|x> = 0.93, returning a spectrum for a pair that is not biorthogonal (or raises `AttributeError` when there was no solve), on `"python"` and v3, and `gs_energy(wf0=x)`, with or without `reconverge=False`, drops x and returns the NH-DMRG energy -1.596396 (|<x|wf0>|^2 = 0.01 to 0.11), on `"python"`, v2 and v3

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `session`

**Status**: FIXED, with the contract stated: a right state with no left state of its own gives no biorthogonal KPM. `gs_energy_nhdmrg()` and `gs_energy_generalized_nhdmrg()` record the right state their `nh_left_wf` pairs with, the object itself, as `mark_injected()` does, and `nonhermitian/kpm.py` raises `RuntimeError` naming what is missing whenever the chain's right state is not that one, which `set_gs()`, `set_initial_wf()` and `gs_energy(wf0=x, reconverge=False)` all produce; a left state for an arbitrary x is not computed, since a state that is not an eigenstate has no biorthogonal partner. `gs_energy(wf0=x)` on a non-Hermitian chain raises `TypeError`, NH-DMRG taking no start state, and `gs_energy(wf0=x, reconverge=False)` takes x as it is with e0 = <x|H|x>, as `set_gs()` does. `set_initial_wf_guess()` on a non-Hermitian chain still re-solves from a random start, as the existing comment in `groundstate.gs_energy` says. Pinned by `tests/test_audit_2026_09_25_session.py::test_nh_kpm_refuses_a_right_state_with_no_left_partner` and `::test_nh_kpm_runs_on_the_pair_of_a_solve` on `python` and `v3`, and `::test_nh_gs_energy_wf0_is_taken_or_refused` on `python`, `v2` and `v3`. NUMBERS CHANGE: `gs_energy(wf0=x, reconverge=False)` on the 4-site chain goes from the NH-DMRG -1.596396 to <x|H|x> on `python`, v2 and v3 (0.060813-0.000217j for the `python` run's seeded x; v2 and v3 draw x on the C++ side, 0.420727-0.035754j and -0.00544+0.024734j in the after run); the NH-KPM after `set_gs()` returns no spectrum where it returned one on `python` and v3 (d[0] = -4.7804-9.2243j at es=0), and `gs_energy(wf0=x)` raises where it returned -1.596396 on all three.

**Recorded in**: 2026-09-24c New leads. The NH-KPM pairing and gs_energy_nhdmrg's permissive keywords predate 867e2b4; 867e2b4 made set_gs() take x on a non-Hermitian chain, which is what leaves nh_left_wf paired with a different right state.

**Where**: src/dmrgpy/nonhermitian/kpm.py:38-46 on 8dd2198 (psil = self.nh_left_wf, whatever solve it came from, and norm = psil.dot(psir)), src/dmrgpy/nhdmrg.py gs_energy_nhdmrg (accepts and ignores unknown keywords, wf0 included), groundstate.gs_energy's non-Hermitian branch (the pending take only without kwargs)

The non-Hermitian KPM needs a biorthogonal pair, and the chain keeps the left state of its last NH-DMRG solve in nh_left_wf with nothing saying which right state it belongs to. After set_gs(x) the right state is x, the left state is still the solve's, and the correlator divides by <psil|x> so that the pair looks biorthogonal, which hides that it is not; with no previous solve it raised AttributeError instead. A right state that is not an eigenstate has no biorthogonal partner at all, so there is nothing consistent to pair it with. Separately, gs_energy(wf0=x) reached gs_energy_nhdmrg(), which accepts and ignores unknown keywords, so x was dropped without a word and NH-DMRG re-solved from random, reconverge=False included; this half is also on v2, whose NH-KPM is NotImplementedError on both trees.

**Expected**: A right state with no left state of its own gives no biorthogonal KPM, so the call raises naming what is missing; gs_energy(wf0=x) either honours x under a stated contract or raises TypeError.

Repro (<scratch>/session/04_nh_injected.py (revised in this pass to add v2; both trees rerun)):

````python
# 04_nh_injected.py
# nh-injected: after set_gs(x) on a non-Hermitian chain, what left state does
# the non-Hermitian KPM pair x with, and does gs_energy(wf0=x) read x?
# 4-site Heisenberg chain plus a non-Hermitian onsite term 0.3j*Sz0 + 0.2*Sx1,
# on python, v2 and v3 (the NH-KPM is NotImplementedError on v2 on both trees).
import io, contextlib
import numpy as np
import dmrgpy
from dmrgpy import spinchain, cppext
print("dmrgpy from", dmrgpy.__file__)

def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()):
        return f()

def nh_chain(version, n=4):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.3j*sc.Sz[0] + 0.2*sc.Sx[1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    return sc

def fid(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))

ES = np.linspace(0.0, 4.0, 5)
KPM = dict(name=None, submode="KPM", es=ES, delta=0.3, E_max=10.0, n=50)
for version in ("python", 2, 3):
    if version in (2, 3) and not cppext.available(version): continue
    np.random.seed(5)
    # (a) set_gs(x) after an NH-DMRG solve: which left state does KPM use?
    sc = nh_chain(version)
    e_nh = quiet(sc.gs_energy)
    psir, psil = sc.wf0.copy(), sc.nh_left_wf
    x = psir + 0.5*(sc.Sx[2]*psir) ; x = x*(1/np.sqrt(x.dot(x).real))
    sc.set_gs(x)
    KPM["name"] = (sc.Sz[0], sc.Sz[0])
    try:
        es, d = quiet(lambda: sc.get_dynamical_correlator(**KPM))
        print("%-6s set_gs(x) after NH-DMRG: NH-KPM ran, pairing x with the "
              "left state of the solve: |<x|psir>|^2=%.4f, <psil|x>=%s, "
              "d[0]=%s" % (version, fid(x, psir), np.round(psil.dot(x), 4),
                           np.round(d[0], 4)))
    except Exception as err:
        print("%-6s set_gs(x) after NH-DMRG: NH-KPM raised %s: %s"
              % (version, type(err).__name__, err))
    # (b) set_gs(x) on a chain that never ran NH-DMRG (a state of its own:
    # an MPS carries its chain's site indices)
    sc2 = nh_chain(version)
    x2 = sc2.random_state() ; x2 = x2*(1/np.sqrt(x2.dot(x2).real))
    sc2.set_gs(x2)
    KPM["name"] = (sc2.Sz[0], sc2.Sz[0])
    try:
        quiet(lambda: sc2.get_dynamical_correlator(**KPM))
        print("%-6s set_gs(x) on a fresh chain: NH-KPM ran" % version)
    except Exception as err:
        print("%-6s set_gs(x) on a fresh chain: NH-KPM raised %s: %s"
              % (version, type(err).__name__, err))
    # (c) gs_energy(wf0=x) on a non-Hermitian chain
    sc3 = nh_chain(version)
    x = sc3.random_state() ; x = x*(1/np.sqrt(x.dot(x).real))
    try:
        e = quiet(lambda: sc3.gs_energy(wf0=x))
        ex = x.dot(sc3.hamiltonian*x)/x.dot(x)
        print("%-6s gs_energy(wf0=x): returned %s, <x|H|x>=%s, NH-DMRG e0=%s, "
              "|<x|wf0>|^2=%.4f" % (version, np.round(e, 6), np.round(ex, 6),
                                    np.round(e_nh, 6), fid(sc3.wf0, x)))
    except Exception as err:
        print("%-6s gs_energy(wf0=x) raised %s: %s"
              % (version, type(err).__name__, err))
    try:
        e = quiet(lambda: sc3.gs_energy(wf0=x, reconverge=False))
        ex = x.dot(sc3.hamiltonian*x)/x.dot(x)
        print("%-6s gs_energy(wf0=x, reconverge=False): returned %s, "
              "<x|H|x>=%s, |<x|wf0>|^2=%.4f" % (version, np.round(e, 6),
                                                np.round(ex, 6), fid(sc3.wf0, x)))
    except Exception as err:
        print("%-6s gs_energy(wf0=x, reconverge=False) raised %s: %s"
              % (version, type(err).__name__, err))

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 04_nh_injected.py 2>&1 | tee 04_nh_injected.before.out
# cd $WF/session && ../run3.sh -u 04_nh_injected.py 2>&1 | tee 04_nh_injected.after.out
````

Observed, before:

````
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
python set_gs(x) after NH-DMRG: NH-KPM ran, pairing x with the left state of the solve: |<x|psir>|^2=0.9443, <psil|x>=(0.9309+0j), d[0]=(-4.7804-9.2243j)
python set_gs(x) on a fresh chain: NH-KPM raised AttributeError: 'Spin_Chain' object has no attribute 'nh_left_wf'
python gs_energy(wf0=x): returned (-1.596396-0j), <x|H|x>=(0.060813-0.000217j), NH-DMRG e0=(-1.596396-0j), |<x|wf0>|^2=0.0230
python gs_energy(wf0=x, reconverge=False): returned (-1.596396-0j), <x|H|x>=(0.060813-0.000217j), |<x|wf0>|^2=0.0230
2      set_gs(x) after NH-DMRG: NH-KPM raised NotImplementedError: NH-KPM dynamical correlator is only implemented for itensor_version 3 or "python" so far, got 2
2      set_gs(x) on a fresh chain: NH-KPM raised NotImplementedError: NH-KPM dynamical correlator is only implemented for itensor_version 3 or "python" so far, got 2
2      gs_energy(wf0=x): returned (-1.596396+0j), <x|H|x>=(0.403387+0.093055j), NH-DMRG e0=(-1.596396-0j), |<x|wf0>|^2=0.0117
2      gs_energy(wf0=x, reconverge=False): returned (-1.596396-0j), <x|H|x>=(0.403387+0.093055j), |<x|wf0>|^2=0.0117
3      set_gs(x) after NH-DMRG: NH-KPM ran, pairing x with the left state of the solve: |<x|psir>|^2=0.9443, <psil|x>=(0.9309-0j), d[0]=(-4.7804-9.2243j)
3      set_gs(x) on a fresh chain: NH-KPM raised AttributeError: 'Spin_Chain' object has no attribute 'nh_left_wf'
3      gs_energy(wf0=x): returned (-1.596396+0j), <x|H|x>=(0.109134-0.046879j), NH-DMRG e0=(-1.596396-0j), |<x|wf0>|^2=0.0485
3      gs_energy(wf0=x, reconverge=False): returned (-1.596396+0j), <x|H|x>=(0.109134-0.046879j), |<x|wf0>|^2=0.0485

(The first-pass parent runs of the python/v3 version of this script gave |<x|wf0>|^2 = 0.1071 and 0.0647 on v3, whose random states are drawn by the C++ side and are not seeded; the same holds for v2.)
````

Observed, after:

````
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
python set_gs(x) after NH-DMRG: NH-KPM raised RuntimeError: the non-Hermitian KPM correlator needs the left eigenvector that pairs with the chain's right state, and this chain has none: its state was set with set_gs(), set_initial_wf() or gs_energy(wf0=...), which give a right state only. Solve the biorthogonal pair with NH-DMRG instead, restart() and then gs_energy()
python set_gs(x) on a fresh chain: NH-KPM raised RuntimeError: the non-Hermitian KPM correlator needs the left eigenvector that pairs with the chain's right state, and this chain has none: its state was set with set_gs(), set_initial_wf() or gs_energy(wf0=...), which give a right state only. Solve the biorthogonal pair with NH-DMRG instead, restart() and then gs_energy()
python gs_energy(wf0=x) raised TypeError: gs_energy(wf0=...) on a non-Hermitian Hamiltonian: NH-DMRG takes no start state, so a state to sweep from cannot be honoured; pass reconverge=False (and nothing else) to take it as it is, as set_gs() does
python gs_energy(wf0=x, reconverge=False): returned (0.060813-0.000217j), <x|H|x>=(0.060813-0.000217j), |<x|wf0>|^2=1.0000
2      set_gs(x) after NH-DMRG: NH-KPM raised NotImplementedError: NH-KPM dynamical correlator is only implemented for itensor_version 3 or "python" so far, got 2
2      set_gs(x) on a fresh chain: NH-KPM raised NotImplementedError: NH-KPM dynamical correlator is only implemented for itensor_version 3 or "python" so far, got 2
2      gs_energy(wf0=x) raised TypeError: gs_energy(wf0=...) on a non-Hermitian Hamiltonian: NH-DMRG takes no start state, so a state to sweep from cannot be honoured; pass reconverge=False (and nothing else) to take it as it is, as set_gs() does
2      gs_energy(wf0=x, reconverge=False): returned (0.420727-0.035754j), <x|H|x>=(0.420727-0.035754j), |<x|wf0>|^2=1.0000
3      set_gs(x) after NH-DMRG: NH-KPM raised RuntimeError: the non-Hermitian KPM correlator needs the left eigenvector that pairs with the chain's right state, and this chain has none: its state was set with set_gs(), set_initial_wf() or gs_energy(wf0=...), which give a right state only. Solve the biorthogonal pair with NH-DMRG instead, restart() and then gs_energy()
3      set_gs(x) on a fresh chain: NH-KPM raised RuntimeError: the non-Hermitian KPM correlator needs the left eigenvector that pairs with the chain's right state, and this chain has none: its state was set with set_gs(), set_initial_wf() or gs_energy(wf0=...), which give a right state only. Solve the biorthogonal pair with NH-DMRG instead, restart() and then gs_energy()
3      gs_energy(wf0=x) raised TypeError: gs_energy(wf0=...) on a non-Hermitian Hamiltonian: NH-DMRG takes no start state, so a state to sweep from cannot be honoured; pass reconverge=False (and nothing else) to take it as it is, as set_gs() does
3      gs_energy(wf0=x, reconverge=False): returned (-0.00544+0.024734j), <x|H|x>=(-0.00544+0.024734j), |<x|wf0>|^2=1.0000
````

**NUMBERS CHANGE**: 4-site Heisenberg chain plus 0.3j*Sz0 + 0.2*Sx1: gs_energy(wf0=x, reconverge=False) -1.596396 (NH-DMRG) -> <x|H|x> on python (0.060813-0.000217j, seeded), v2 and v3 (unseeded C++ draws: 0.420727-0.035754j and -0.00544+0.024734j in the after run); the NH-KPM after set_gs(x) now raises RuntimeError where it returned a spectrum on python and v3 (d[0] = -4.7804-9.2243j at es=0); gs_energy(wf0=x) now raises TypeError where it returned -1.596396 on python, v2 and v3.

**Tests**: `tests/test_audit_2026_09_25_session.py::test_nh_kpm_refuses_a_right_state_with_no_left_partner[python]`; `tests/test_audit_2026_09_25_session.py::test_nh_kpm_refuses_a_right_state_with_no_left_partner[v3]`; `tests/test_audit_2026_09_25_session.py::test_nh_kpm_runs_on_the_pair_of_a_solve[python]`; `tests/test_audit_2026_09_25_session.py::test_nh_kpm_runs_on_the_pair_of_a_solve[v3]`; `tests/test_audit_2026_09_25_session.py::test_nh_gs_energy_wf0_is_taken_or_refused[python]`; `tests/test_audit_2026_09_25_session.py::test_nh_gs_energy_wf0_is_taken_or_refused[v3]`; `tests/test_audit_2026_09_25_session.py::test_nh_gs_energy_wf0_is_taken_or_refused[v2]`

**Reviewer (CONFIRMED, fix HOLDS)**: Reproduced on python and v3, to every printed digit where the state is seeded (python). v3 draws x on the C++ side, so its <x|H|x> differs between runs, as the fix agent notes. After set_gs(x) the NH-KPM ran with the solve's left state (<psil|x> = 0.9309). On a fresh chain it raised AttributeError. gs_energy(wf0=x), with or without reconverge=False, returned the NH-DMRG -1.596396. The gs_energy(wf0) half is also real on v2 (returned -1.596396, |<x|wf0>|^2 = 0.011409). The KPM half cannot occur on v2, whose NH-KPM raises NotImplementedError on both trees.

Holds on python and v3, and on v2 for the gs_energy(wf0=x) half: TypeError without reconverge=False, the take with it. On v2 the pairing check is unreachable, since NH-KPM is NotImplementedError there before and after. The identity-based pairing survives a clone by reading: Many_Body_Chain.__deepcopy__ passes both wf0 and _nh_left_for through deepcopy(v,memo). A clone also re-solves anyway, since its send cache is None. The check also refuses set_gs() of the solved NH state itself, a copy. That is conservative and consistent with the stated contract.

The NUMBERS CHANGE statement names python and v3 only; v2 changes too. There, gs_energy(wf0=x, reconverge=False) goes from the NH-DMRG -1.596396 to <x|H|x> ((0.38016+0.06014j) in my run), and gs_energy(wf0=x) now raises TypeError.

The fix agent's left_open line ("the first correlator on a solved non-Hermitian chain re-solves NH-DMRG once ... not wrong, a second solve") is not harmless after gs_energy_generalized() on a non-Hermitian chain, where that re-solve discards the generalized state. See generalized-cache.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
04 rerun, fixed:
python set_gs(x) after NH-DMRG: NH-KPM raised RuntimeError: the non-Hermitian KPM correlator needs the left eigenvector that pairs with the chain's right state, and this chain has none: ...
python gs_energy(wf0=x) raised TypeError: gs_energy(wf0=...) on a non-Hermitian Hamiltonian: NH-DMRG takes no start state, ...
python gs_energy(wf0=x, reconverge=False): returned (0.060813-0.000217j), <x|H|x>=(0.060813-0.000217j), |<x|wf0>|^2=1.0000
3      gs_energy(wf0=x, reconverge=False): returned (0.151733-0.006075j), <x|H|x>=(0.151733-0.006075j), |<x|wf0>|^2=1.0000
parent:
python set_gs(x) after NH-DMRG: NH-KPM ran, pairing x with the left state of the solve: |<x|psir>|^2=0.9443, <psil|x>=(0.9309+0j), d[0]=(-4.7804-9.2243j)
python set_gs(x) on a fresh chain: NH-KPM raised AttributeError: 'Spin_Chain' object has no attribute 'nh_left_wf'
python gs_energy(wf0=x, reconverge=False): returned (-1.596396-0j), <x|H|x>=(0.060813-0.000217j), |<x|wf0>|^2=0.0230
v2 (10_neighbours.py section B), parent:
v2 NH-KPM after set_gs(x) raised NotImplementedError: NH-KPM dynamical correlator is only implemented for itensor_version 3 or "python" so far, got 2
v2 gs_energy(wf0=x) returned (-1.596396+0j)
v2 gs_energy(wf0=x, reconverge=False) = (-1.596396-0j), <x|H|x> = (0.39877+0.090491j), |<x|wf0>|^2 = 0.011409
v2, fixed:
v2 gs_energy(wf0=x) raised TypeError: gs_energy(wf0=...) on a non-Hermitian Hamiltonian: NH-DMRG takes no start state, ...
v2 gs_energy(wf0=x, reconverge=False) = (0.38016+0.06014j), <x|H|x> = (0.38016+0.06014j), |<x|wf0>|^2 = 1.000000
Files: <scratch>/review/session/04.before.out, 04.after.out, 10.before.out, 10.after.out
````

**Left open**: set_initial_wf_guess(x) on a non-Hermitian chain still drops x and re-solves from random, the same shape as the dropped wf0, left as the existing comment states; making it raise would be the consistent choice. The first pass recorded here that the first correlator on a solved non-Hermitian chain re-solves NH-DMRG once; that was not harmless after gs_energy_generalized(), where the re-solve discarded the generalized state, and it is fixed under generalized-cache (the solve count goes 1 to 0). The pairing check refuses set_gs() of the solved right state itself, a copy, which is conservative and consistent with the contract.

### 9. No docstring says that `maxde` is a fluctuation per site while `gs_energy_fluctuation()` returns the total, and the loop prints its per-site reading unlabelled as `Energy fluctuation =` (0.026246 for a state whose `gs_energy_fluctuation()` is 0.2625, 10-site S=1 chain at maxm=10 on v3); on 8dd2198 the example and the user guide already read it per site

`bug` &middot; severity **LOW** &middot; NARROWED &middot; cluster `session`

**Status**: FIXED, narrowed on reproduction: the lead was measured on 867e2b4, and on 8dd2198 the user guide's fluctuation paragraph already says per site and `examples/groundstate/GS_enforce_maximum_fluctuation` already plots `gs_energy_fluctuation()/n` against `maxde`, so those halves did not reproduce; what did is that `Many_Body_Chain.gs_energy()` and `groundstate.gs_energy_single()` said nothing about `maxde`, and the loop printed its per-site reading with the same bare label as the total. `maxde` stays per site, the decision taken. Both docstrings now say that `maxde` is a tolerance on ||(H-<H>)|psi>||/ns, so that it is `gs_energy_fluctuation()/ns` that is compared with it, `gs_energy_fluctuation()`'s docstring says it returns the total, and the loop prints `Energy fluctuation per site =`. The example pins its schedule at the defaults it used to inherit (`maxm=10`, `nsweeps=15`, `noise=1e-7`, `cutoff=1e-12`), labels the request and the axis as per site and saves `GS_enforce_maximum_fluctuation.png`. Pinned by `tests/test_audit_2026_09_25_session.py::test_maxde_is_documented_per_site` and `::test_maxde_compares_the_fluctuation_per_site` on `python` and `v3` (the loop's first reading is `gs_energy_fluctuation()/ns` of the state it starts from, and it stops once that is below `maxde`; against 8dd2198 the second fails only on the label). No number changes, in the library or in the example, whose three points are -12.872525, -12.894560 and -12.894560 with 2.62e-02, 2.03e-04 and 1.76e-06 per site on both trees.

**Recorded in**: 2026-09-24c New leads. The per-site division is from 4731b5a; the silence is the docstrings'. The example and user-guide halves of the lead were fixed in 8dd2198.

**Where**: src/dmrgpy/groundstate.py:396-407 on 8dd2198 (de/self.ns compared with maxde, printed as 'Energy fluctuation ='), Many_Body_Chain.gs_energy's one-line docstring, groundstate.gs_energy_single's docstring, Many_Body_Chain.gs_energy_fluctuation's docstring

The loop divides the fluctuation by ns before comparing it with maxde, which is the intensive choice and the one the main agent kept, but a reader of gs_energy() had no way to know, and the loop's own print used the same bare label as the total. The lead was measured on 867e2b4; on 8dd2198 the user guide's fluctuation paragraph (user_guide.md:709-712 and its .tex twin at 819-821) already says per site, and examples/groundstate/GS_enforce_maximum_fluctuation already plots gs_energy_fluctuation()/n against maxde, so neither of those halves reproduces. The example did not pin its sweep schedule or save its figure, which the new-example rules ask for.

**Expected**: gs_energy()'s docstring says maxde is a tolerance on ||(H-<H>)|psi>||/ns, gs_energy_fluctuation()'s says it returns the total, and the example plots the per-site fluctuation against maxde with its schedule pinned.

Repro (<scratch>/session/05_maxde_per_site.py (example runs: 06_example.{before,after}.out, figure example_after/GS_enforce_maximum_fluctuation.png)):

````python
# 05_maxde_per_site.py
# maxde-per-site: the maxde loop compares a per-site fluctuation against
# maxde while gs_energy_fluctuation() returns the total. Does the code say
# which one maxde is, and does the example plot like with like? 10-site S=1
# Heisenberg chain at maxm=10 on v3, the example's own setting.
import io, contextlib, inspect, re
import numpy as np
import dmrgpy
from dmrgpy import spinchain, groundstate
from dmrgpy.manybodychain import Many_Body_Chain
print("dmrgpy from", dmrgpy.__file__)

n = 10
sc = spinchain.Spin_Chain(["S=1"]*n, itensor_version=3)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 10, 10
np.random.seed(6)
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    sc.gs_energy()
    first = sc.gs_energy_fluctuation()
    sc.set_hamiltonian(h)
    sc.gs_energy(maxde=1e-3)
reads = [l for l in buf.getvalue().splitlines() if "Energy fluctuation" in l]
print("gs_energy_fluctuation() of the first solve: %.3e total, %.3e per site"
      % (first, first/n))
print("the maxde loop's first reading:", reads[0] if reads else None)
for name, f in (("Many_Body_Chain.gs_energy", Many_Body_Chain.gs_energy),
                ("groundstate.gs_energy_single", groundstate.gs_energy_single)):
    doc = inspect.getdoc(f) or ""
    print("%s docstring mentions maxde: %s, says per site: %s"
          % (name, "maxde" in doc, bool(re.search(r"per.site", doc))))

# cd $WF/session && DMRGPY_SRC=$WF/parent/src ../run3.sh -u 05_maxde_per_site.py 2>&1 | tee 05_maxde_per_site.before.out
# cd $WF/session && ../run3.sh -u 05_maxde_per_site.py 2>&1 | tee 05_maxde_per_site.after.out
# The example, from a scratch copy of examples/groundstate/GS_enforce_maximum_fluctuation/main.py:
# cd $WF/session/example_before && DMRGPY_SRC=$WF/parent/src ../../run3.sh -u main.py 2>&1 | tee ../06_example.before.out
# cd $WF/session/example_after && ../../run3.sh -u main.py 2>&1 | tee ../06_example.after.out
````

Observed, before:

````
05_maxde_per_site.before.out:
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
gs_energy_fluctuation() of the first solve: 2.625e-01 total, 2.625e-02 per site
the maxde loop's first reading: Energy fluctuation =  0.026245793110247745 10
Many_Body_Chain.gs_energy docstring mentions maxde: False, says per site: False
groundstate.gs_energy_single docstring mentions maxde: False, says per site: False

06_example.before.out (the example as it stands on 8dd2198, already plotting per site):
[run3] slot 1 acquired
maxde 0.1 Energy -12.872524926071373 fluctuation per site 0.026245796394550318
Energy fluctuation =  0.026245799876326865 10
Energy fluctuation =  0.0022828272966165116 20
maxde 0.001 Energy -12.894559608430365 fluctuation per site 0.00020322613328381612
Energy fluctuation =  0.026245799382343744 10
Energy fluctuation =  0.002282827268299778 20
Energy fluctuation =  0.00020322612837852416 40
Energy fluctuation =  2.5217686847781315e-06 80
Energy fluctuation =  1.764827710902212e-06 160
maxde 1e-06 Energy -12.894560132185685 fluctuation per site 1.7648286797271286e-06
````

Observed, after:

````
05_maxde_per_site.after.out:
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
gs_energy_fluctuation() of the first solve: 2.625e-01 total, 2.625e-02 per site
the maxde loop's first reading: Energy fluctuation per site =  0.026245794783396947 10
Many_Body_Chain.gs_energy docstring mentions maxde: True, says per site: True
groundstate.gs_energy_single docstring mentions maxde: True, says per site: True

06_example.after.out (the revised example, schedule pinned):
[run3] slot 1 acquired
maxde 0.1 Energy -12.87252492605654 fluctuation per site 0.02624579970103009
Energy fluctuation per site =  0.026245799871723922 10
Energy fluctuation per site =  0.002282827296630291 20
maxde 0.001 Energy -12.894559608430347 fluctuation per site 0.00020322613328218882
Energy fluctuation per site =  0.026245794814816276 10
Energy fluctuation per site =  0.002282827088299012 20
Energy fluctuation per site =  0.00020322610118516448 40
Energy fluctuation per site =  2.5217825300480836e-06 80
Energy fluctuation per site =  1.764833021024837e-06 160
maxde 1e-06 Energy -12.894560132185717 fluctuation per site 1.7648261662615573e-06
````

**NUMBERS CHANGE**: No number changes. The loop's print label changes from 'Energy fluctuation =' to 'Energy fluctuation per site ='.

**Tests**: `tests/test_audit_2026_09_25_session.py::test_maxde_is_documented_per_site`; `tests/test_audit_2026_09_25_session.py::test_maxde_compares_the_fluctuation_per_site[python]`; `tests/test_audit_2026_09_25_session.py::test_maxde_compares_the_fluctuation_per_site[v3]`

**Reviewer (NARROWED, fix HOLDS)**: The narrowed claim, which is what goes into the record: at 8dd2198 neither Many_Body_Chain.gs_energy() nor groundstate.gs_energy_single() said anything about maxde, and the maxde loop printed its per-site reading with the bare label "Energy fluctuation =", the same wording as the total that gs_energy_fluctuation() returns. The value printed is 0.026245..., where gs_energy_fluctuation() is 2.625e-01, on a 10-site S=1 chain at maxm=10 on v3.

Struck: the lead's example half and user-guide half. `git show 8dd2198` shows the example already plotting gs_energy_fluctuation()/n against maxde, and docs/user_guide.md at 8dd2198 line 710 already reads "reads the same quantity **per site**". The fix agent narrowed this themselves and I confirm it.

Both docstrings now say that maxde is per site, the loop prints "Energy fluctuation per site =", and gs_energy_fluctuation()'s docstring says it returns the total. The example's pinned schedule (maxm=10, nsweeps=15, noise=1e-7, cutoff=1e-12) equals the library defaults (manybodychain.py:106, 140, 175). Its three points are unchanged on both trees to about 1e-12, energies -12.8725249, -12.8945596, -12.8945601 and per-site fluctuations 2.62e-2, 2.03e-4, 1.76e-6, and the saved figure plots per site against per site. No library number moves. test_maxde_is_documented_per_site greps the docstrings and is brittle, but it pins the property. test_maxde_compares_the_fluctuation_per_site fails on the parent only on the label.

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
05 rerun, parent: "the maxde loop's first reading: Energy fluctuation =  0.026245794235743664 10", "Many_Body_Chain.gs_energy docstring mentions maxde: False, says per site: False"; fixed: "the maxde loop's first reading: Energy fluctuation per site =  0.026245791226805298 10", "Many_Body_Chain.gs_energy docstring mentions maxde: True, says per site: True", "groundstate.gs_energy_single docstring mentions maxde: True, says per site: True"; both: "gs_energy_fluctuation() of the first solve: 2.625e-01 total, 2.625e-02 per site".
Example, 8dd2198 version on parent: "maxde 0.1 Energy -12.872524926065362 fluctuation per site 0.026245799541198756" ... "maxde 1e-06 Energy -12.894560132185717 fluctuation per site 1.764814195436194e-06"; working-tree version on fixed tree: "maxde 0.1 Energy -12.872524925878833 fluctuation per site 0.026245796342797573" ... "maxde 1e-06 Energy -12.894560132185699 fluctuation per site 1.7648143170821063e-06".
Files: <scratch>/review/session/05.before.out, 05.after.out, 06.before.out, 06.after.out, ex_after/GS_enforce_maximum_fluctuation.png
````

**Left open**: The example's last point still sits just above its request (1.76e-06 per site against 1e-06), since the loop stops after five doublings; its comment says so.

### 10. On the lower-level `get_dynamical_correlator_MB` route `submode="ROOTN"` returned C[Sz_0,Sz_0] for `name="ZZ"` at any `i=`/`j=`, 4.065e-02 off the (1,1) curve on a 0.103 peak and 1.213e-01 off the (0,2) one on 0.080, and accepted a misspelled keyword

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; cluster `construction`

**Status**: FIXED. rootndmrg.dynamical_correlator now takes i=0, j=0 and passes them to str2MO, the signature cvm.dynamical_correlator has, and has no **kwargs left, so a keyword it does not read raises TypeError naming it. Every in-tree ROOTN call (the submode-parametrized tests, the two ROOTN examples, the Kondo spectrum route) passes only es, delta, name, N and nkry, the test files that reach it all pass, and the reviewer found the Kondo DMRG route at submode="ROOTN" unchanged to ten digits. Pinned by tests/test_audit_2026_09_25_construction.py::test_rootn_lower_level_route_honours_the_sites at (1,1) and (0,2), the shape of tests/test_audit_2026_09_24_realtime.py::test_lower_level_route_honours_the_sites, and test_rootn_rejects_an_unknown_keyword. NUMBERS CHANGE for get_dynamical_correlator_MB(submode="ROOTN") with a string name and sites other than (0,0): on the 4-site open S=1/2 Heisenberg chain, itensor_version="python", maxm=20, nsweeps=8, es=linspace(0,3,5), delta=0.4, N=4, nkry=12, the (1,1) curve goes from the (0,0) one, 4.065e-02 from ED's (1,1) on a 0.1031 peak, to 4.2e-16 from it (5.1e-16 and 5.8e-16 in other runs, the DMRG start being random; the reviewer measured 1.296e-13 on v2 and 4.236e-13 on v3). The public route is unchanged.

**Recorded in**: 2026-09-24 finding 10 Status (open). docs/audit_2026_09_24_hole_hunt.md finding 10's Status: 'the ROOTN sibling (rootndmrg.py:54) still drops i/j on this route and stays open'

**Where**: src/dmrgpy/rootndmrg.py:7 on 8dd2198 (the consumerless **kwargs of dynamical_correlator's signature) and rootndmrg.py:54 (`operatornames.str2MO(self,name)` with no sites), reached from dynamics.py:254-255

rootndmrg.dynamical_correlator took **kwargs with no consumer and called str2MO(self, name) without sites, so on the lower-level route, where name= reaches the submode still as a string, i= and j= fell into that **kwargs and str2MO took its own i=j=0: get_dynamical_correlator_MB(submode="ROOTN", name="ZZ", i=1, j=1) returned the [Sz_0, Sz_0] curve bit for bit (max|MB(i,j) - pair(0,0)| = 0.000e+00 at both site pairs measured), identically on v2, v3 and "python" in the reviewer's runs. The public get_dynamical_correlator resolves name= together with i=/j= before dispatch, which is why the public route, and every test on it, was never affected. The 2026-09-24 audit's finding 10 fixed TD and TDZ on this route and left ROOTN open; the same **kwargs also took a misspelled nkry (nkyr=3) and returned the default-parameter spectrum, the documentation.md 4.10 shape of a **kwargs with no consumer.

**Expected**: The sites asked for: get_dynamical_correlator_MB(name="ZZ", i=, j=, submode="ROOTN") equals the explicit [Sz_i, Sz_j] pair on the public route, and a keyword ROOTN does not read raises TypeError, as on CVM.

Repro (<scratch>/construction/03_rootn_ij.py (before: DMRGPY_SRC=<parent>/src run3.sh 03_rootn_ij.py | tee 03_rootn_ij.before.out; after: run3.sh 03_rootn_ij.py | tee 03_rootn_ij.after.out, rerun on the repaired tree)):

````python
# rootn-ij: on the lower-level route ROOTN drops i=/j= of a string name
import functools
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print = functools.partial(print, flush=True)
print("dmrgpy from", dmrgpy.__file__)

n = 4
sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 20, 8
K = dict(es=np.linspace(0.0, 3.0, 5), delta=0.4, N=4, nkry=12)

for (i, j) in [(1, 1), (0, 2)]:
    y = np.asarray(sc.get_dynamical_correlator_MB(submode="ROOTN", name="ZZ",
                                                  i=i, j=j, **K)[1])
    y_ij = np.asarray(sc.get_dynamical_correlator(submode="ROOTN",
            name=[sc.Sz[i], sc.Sz[j]], **K)[1])
    y_00 = np.asarray(sc.get_dynamical_correlator(submode="ROOTN",
            name=[sc.Sz[0], sc.Sz[0]], **K)[1])
    print("i,j=%d,%d  max|MB(i,j) - pair(i,j)| = %.3e   max|MB(i,j) - pair(0,0)|"
          " = %.3e   max|pair(i,j)| = %.3e"
          % (i, j, np.max(np.abs(y-y_ij)), np.max(np.abs(y-y_00)),
             np.max(np.abs(y_ij))))

# the exact reference for (1,1), mode="ED" on the public route
y_ed = np.asarray(sc.get_dynamical_correlator(mode="ED", submode="ROOTN",
        name="ZZ", i=1, j=1, **K)[1])
y = np.asarray(sc.get_dynamical_correlator_MB(submode="ROOTN", name="ZZ",
                                              i=1, j=1, **K)[1])
print("i,j=1,1  max|MB(1,1) - ED(1,1)| = %.3e" % np.max(np.abs(y-y_ed)))

# a misspelled keyword on the same route
try:
    sc.get_dynamical_correlator_MB(submode="ROOTN", name="ZZ", nkyr=3, **K)
    print("nkyr=3 (for nkry): accepted")
except TypeError as e:
    print("nkyr=3 (for nkry): TypeError:", e)
````

Observed, before:

````
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
i,j=1,1  max|MB(i,j) - pair(i,j)| = 4.065e-02   max|MB(i,j) - pair(0,0)| = 0.000e+00   max|pair(i,j)| = 1.031e-01
i,j=0,2  max|MB(i,j) - pair(i,j)| = 1.213e-01   max|MB(i,j) - pair(0,0)| = 0.000e+00   max|pair(i,j)| = 8.039e-02
i,j=1,1  max|MB(1,1) - ED(1,1)| = 4.065e-02
nkyr=3 (for nkry): accepted
````

Observed, after:

````
[run3] slot 1 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
i,j=1,1  max|MB(i,j) - pair(i,j)| = 0.000e+00   max|MB(i,j) - pair(0,0)| = 4.065e-02   max|pair(i,j)| = 1.031e-01
i,j=0,2  max|MB(i,j) - pair(i,j)| = 0.000e+00   max|MB(i,j) - pair(0,0)| = 1.213e-01   max|pair(i,j)| = 8.039e-02
i,j=1,1  max|MB(1,1) - ED(1,1)| = 4.178e-16
nkyr=3 (for nkry): TypeError: dynamical_correlator() got an unexpected keyword argument 'nkyr'
````

**NUMBERS CHANGE**: YES. get_dynamical_correlator_MB(submode="ROOTN", name="ZZ", i=1, j=1) on the 4-site open Heisenberg chain, itensor_version="python", maxm=20, nsweeps=8, es=linspace(0,3,5), delta=0.4, N=4, nkry=12: max|MB - ED(1,1)| 4.065e-02 (it was C[Sz_0,Sz_0]) -> 4.2e-16, on a 0.1031 peak; (0,2) moves by 1.213e-01 on 0.08039. Public route unchanged.

**Tests**: `tests/test_audit_2026_09_25_construction.py::test_rootn_lower_level_route_honours_the_sites`; `tests/test_audit_2026_09_25_construction.py::test_rootn_rejects_an_unknown_keyword`

**Reviewer (CONFIRMED, fix HOLDS)**: Reproduced on the parent tree on "python", and I also ran v2 and v3, which the fix agent did not. get_dynamical_correlator_MB(submode='ROOTN', name='ZZ', i=1, j=1) equals the (0,0) curve bit for bit (max|MB-(0,0)| = 0.000e+00) and sits 4.065e-02 from ED's (1,1) on a 0.1031 peak, identically on all three backends. The (0,2) line is 1.213e-01 off on a 0.0804 peak. A misspelled nkyr=3 was accepted. This is the sibling that docs/audit_2026_09_24_hole_hunt.md finding 10 left open by name, and nothing documents it as intended. Nothing struck.

The ROOTN signature now matches cvm.dynamical_correlator (i=0, j=0, no **kwargs). On the working tree, MB(1,1) equals the public route to 0.000e+00 on v2, v3 and python. MB minus ED is 1.296e-13 on v2, 4.236e-13 on v3 and 5.818e-16 on python, and MB minus (0,0) is 4.065e-02. nkyr= raises TypeError naming it. The public route is unchanged across the trees (public[:3] = 0.02948021, 0.10314648, 0.09112285 on both). The Kondo DMRG route with submode='ROOTN' gives the same dIdV to ten digits on both trees, so nothing that route forwards was a keyword ROOTN now refuses. The tests fail on the parent (3 of 3) and pass on the working tree.

One residual, pre-existing and uniform across submodes, so not this fix's to carry: on the lower-level MB route a pair name= next to i=/j= is still dropped silently, on ROOTN and on CVM alike (0.000e+00 on both trees), where the public route raises TypeError (second-pass finding 16).

*This review is of the first fix pass; the cluster then went through a repair pass, whose result is the Status above.*

Reviewer's evidence:

````
03_rootn_ij.before.out:
```
i,j=1,1  max|MB(i,j) - pair(i,j)| = 4.065e-02   max|MB(i,j) - pair(0,0)| = 0.000e+00   max|pair(i,j)| = 1.031e-01
i,j=0,2  max|MB(i,j) - pair(i,j)| = 1.213e-01   max|MB(i,j) - pair(0,0)| = 0.000e+00   max|pair(i,j)| = 8.039e-02
i,j=1,1  max|MB(1,1) - ED(1,1)| = 4.065e-02
nkyr=3 (for nkry): accepted
```
03_rootn_ij.after.out:
```
i,j=1,1  max|MB(i,j) - pair(i,j)| = 0.000e+00   max|MB(i,j) - pair(0,0)| = 4.065e-02   max|pair(i,j)| = 1.031e-01
i,j=0,2  max|MB(i,j) - pair(i,j)| = 0.000e+00   max|MB(i,j) - pair(0,0)| = 1.213e-01   max|pair(i,j)| = 8.039e-02
i,j=1,1  max|MB(1,1) - ED(1,1)| = 5.143e-16
nkyr=3 (for nkry): TypeError: dynamical_correlator() got an unexpected keyword argument 'nkyr'
```
05_attack_gs_rootn.before.out:
```
[2] ROOTN (1,1): max|MB-public|=4.065e-02  max|MB-ED|=4.065e-02  max|MB-(0,0)|=0.000e+00  public[:3]=[0.02948021+0.j 0.10314648+0.j 0.09112285+0.j]
[3] ROOTN (1,1): max|MB-public|=4.065e-02  max|MB-ED|=4.065e-02  max|MB-(0,0)|=0.000e+00  public[:3]=[0.02948021+0.j 0.10314648+0.j 0.09112285+0.j]
Kondo DMRG order=2 submode=ROOTN dIdV                        -> [11.66089377  8.76691383  7.24856614  7.24856614  8.76691383 11.66089377]
```
05_attack_gs_rootn.after.out:
```
[2] ROOTN (1,1): max|MB-public|=0.000e+00  max|MB-ED|=1.296e-13  max|MB-(0,0)|=4.065e-02  public[:3]=[0.02948021+0.j 0.10314648+0.j 0.09112285+0.j]
[3] ROOTN (1,1): max|MB-public|=0.000e+00  max|MB-ED|=4.236e-13  max|MB-(0,0)|=4.065e-02  public[:3]=[0.02948021+0.j 0.10314648+0.j 0.09112285+0.j]
[python] ROOTN (1,1): max|MB-public|=0.000e+00  max|MB-ED|=5.818e-16  max|MB-(0,0)|=4.065e-02  public[:3]=[0.02948021-0.j 0.10314648-0.j 0.09112285-0.j]
[python] MB ROOTN name=pair i=1 j=1: max|y - pair|           -> 0.000e+00
[python] MB CVM name=pair i=1 j=1: max|y - pair|             -> 0.000e+00
Kondo DMRG order=2 submode=ROOTN dIdV                        -> [11.66089377  8.76691383  7.24856614  7.24856614  8.76691383 11.66089377]
```
````

**Left open**: On the lower-level get_dynamical_correlator_MB route a pair name= next to i=/j= is still dropped silently, on ROOTN and on CVM alike (max|y - pair| = 0.000e+00 on both trees in the reviewer's runs), where the public route raises TypeError for the same call (second-pass finding 16); pre-existing and uniform across submodes, so not this item's.

## Left open, and new leads

Nothing below was fixed in this pass. The first list is what the main agent
takes from the entries and the reviews, in order of how much it matters; the
reviewers' own reports follow verbatim, since several carry repros.

- `gs_energy(maxde=...)` on a chain whose ground state is current returns the
  stored energy without refining (-3.3468165405 against -3.4061631313 on a fresh
  chain, found by the construction reviewer): the same short circuit as `wf0=`,
  in `gs_energy`'s own condition, and `get_gs` must widen with it. (fixed in
  the 2026-09-25b fix pass, see tests/test_audit_2026_09_25b_construction.py)
- `submode="TD"` and every TDVP real-time route depend on units on `"python"` and
  v3 through the Krylov exponentiator's absolute error goal, in the operator norm
  as well as in the Hamiltonian's scale (the scale reviewer: 8.76e-02 of the
  peak at an operator scale of 1e-8 on `"python"`, 1.04e-04 at 1e-9 on v3). On
  `"python"` at operator scales at or below 1e-8 the parent raised
  `ZeroDivisionError`, since the operator had no terms; after item 3 the same
  call returns a spectrum that far off, so item 3 turned a loud failure into a
  quiet one there. The suggested cure is to evolve the normalized state in
  `quench_tdvp` and multiply the correlator back, and on v3 the same in
  `chain_session.h`, which needs a rebuild.
- The bond-local truncation of item 2: a bond whose strongest crossing term is
  far below the largest coefficient of the operator loses channels at any units
  on v2, v3 and `"python"`, pinned by strict xfails; the per-bond design is in
  the locate stage's plan.
- A Hamiltonian set as an already-built MPO on v3, real-time evolution at small
  units (v3 TDVP 6.4e-4 off at 1e-8; the v2 MPO-Taylor stepper diverging),
  NH-DMRG on v2 at small units, and the absolute tests of `"python"`'s own
  Lanczos, all recorded under item 2.
- Item 3's relative rule drops a term more than twelve decades below the largest
  coefficient of its operator, so `1e6*Id + 1e-7*Sz0` loses its field and
  `1e6*Sz1 + 1e-7*Sz0 - 1e6*Sz1` is proven zero, where the absolute rule kept
  both (the scale reviewer). This is kept by design: a factor near roundoff
  would let rounding dust in coefficients a caller computed make an ordinary
  Hamiltonian fail its Hermiticity proof and go to NH-DMRG.
- `julia_live` records no solver key, so a `maxm` ramp on one chain returns the
  first energy every time (-3.194321 at maxm 2, 4 and 16 on an 8-site chain)
  (fixed in the 2026-09-25b fix pass, see tests/test_audit_2026_09_25b_session.py),
  and its KPM re-solves the chain for the lower band edge on every call of a
  supplied state.
- On the lower-level `get_dynamical_correlator_MB` route a pair `name=` next to
  `i=`/`j=` is still dropped silently, on ROOTN and CVM alike, where the public
  route raises (second pass finding 16).
- `set_initial_wf_guess(x)` on a non-Hermitian chain still drops x and re-solves
  from random, as the existing comment says; raising would be the consistent
  choice.
- `get_gs(best=True, **kwargs)` forwards to a `best_gs` that takes no keyword, and
  an unknown keyword to `Thermal_Spin_Chain` raises naming `Spin_Chain()`. (both
  fixed in the 2026-09-25b fix pass, see
  tests/test_audit_2026_09_25b_construction.py)

### What the construction reviewer reported beyond its items

````
1. Introduced by the init-kwargs fix, severity LOW to MEDIUM. Thermal_Spin_Chain(sites, T=..., mode="ED").get_gs() now raises AttributeError: 'NoneType' object has no attribute 'set_sweep_params'. On the parent the same call ran DMRG and silently ignored mode. The cause: mode="ED" now reaches MBChain, which builds no session, but thermal.py:33 hardcodes the wrapper's own self.mode="DMRG" and thermal.py:45 writes that onto MBChain in get_gs(). Repro, from 04_attack_init.py in the review folder:
```
tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=0.5, itensor_version="python", mode="ED")
# Heisenberg on tc.Sx/Sy/Sz, tc.set_hamiltonian(h), then
tc.get_gs()
```
Parent (04_attack_init.before.out):
```
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs()              -> MBChain.mode='DMRG', <Sz0 Sz1>=-0.125704
```
Fixed tree (04_attack_init.after.out):
```
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs()              -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
Thermal_Spin_Chain(T=0, mode='ED').get_gs()                -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
```
Suggested fix: Thermal_Spin_Chain.__init__ pops mode= into its own self.mode, keeping "DMRG" as the default, and does not forward it to MBChain. The session cluster's own test already uses tc.mode="ED" as the knob. thermal.py is in the session cluster's diff, so coordinate with that owner.

2. Also introduced by the fix: mode="ED" at construction is not the same as assigning it afterwards. Spin_Chain(sites, mode="ED") has _session None, while sc.mode="ED" after construction keeps a session. Setting sc.mode=None later and asking for DMRG then fails with the same opaque AttributeError. The fix agent recorded the round trip, but the docstring and the check_settings error text still claim equivalence for every keyword.
```
[python ctor] mode=None; gs_energy()                       -> AttributeError: 'NoneType' object has no attribute 'set_sweep_params'
[python assigned] mode=None; gs_energy()                   -> -1.6282694864
```
Either document that mode="ED" at construction makes an ED-only chain, or have the DMRG path raise a clear error when _session is None.

3. Pre-existing lead, not introduced here and uniform across submodes: on the lower-level get_dynamical_correlator_MB route, a pair name= next to i=/j= is still dropped silently. ROOTN and CVM give 0.000e+00 difference on both trees, while the public route raises TypeError for the same call (second-pass finding 16).

4. Confirmed leads the fix agent already listed. Fermionic_Chain(4, N=5).N comes back as 5, and Parafermionic_Chain(3, Sig=1) crashes inside its own __init__ with TypeError: 'int' object is not subscriptable: check_settings admits operator lists that exist before the base constructor runs. gs_energy(maxde=1e-4) on a current chain returns the stored energy without refining (-3.3468165405 against -3.4061631313 on a fresh chain), and that one belongs to the session cluster.
````

### The construction agent's notes

````
What the repair changed, against the reviewer's verdicts. init-kwargs (INCOMPLETE): mode= is no longer applied before initialize(); every setting, mode included, is assigned after it, so a constructor keyword is exactly an assignment afterwards. That removes both regressions the reviewer found (a mode="ED" chain with no session, so sc.mode=None raised AttributeError, and Thermal_Spin_Chain(..., mode="ED").get_gs() raising where the parent ran), and it gives back the one thing mode-first bought: Bosonic_Chain(3, maxnb=[6]*3, itensor_version=2, mode="ED") raises the v2 ValueError, as on the parent, and that message (bosonchain.py, owned here) no longer offers mode="ED" as a way out. On the reviewer's leads, the operator lists a model class builds before calling the base constructor are now refused too: Many_Body_Chain.__init__ snapshots set(vars(self)) at its top, which is exactly what the model class set before handing over, and check_settings refuses those names; this is an attribute test, not a guess on the value type, and it closes Fermionic_Chain(4, N=5).N == 5 (the first pass overwrote the list, the parent dropped the keyword) and Parafermionic_Chain(3, Sig=1) failing inside its own constructor. get-gs-wf0 and rootn-ij (HOLDS): no code change; the record takes the reviewer's sharpenings (every later reader measures x, v2 confirmed, get_gs(best=True, wf0=x) now raising on a current chain, the pair-name lead on the MB route).

Thermal_Spin_Chain, the exact edit for the session cluster (thermal.py is theirs): change the signature to `def __init__(self,sites,T=0.1,mode="DMRG",**kwargs):`, so mode is consumed by the wrapper and not forwarded to Spin_Chain, validate it with `from .mode import _check_mode; _check_mode(mode,"Thermal_Spin_Chain mode (mode=)")`, and replace the hardcoded `self.mode = "DMRG"` (thermal.py:33) with `self.mode = mode`. get_gs() already writes self.mode onto MBChain, so Thermal_Spin_Chain(..., mode="ED") then solves by ED; the session cluster's own tests drive it through tc.mode="ED", which stays compatible. Until then, Thermal_Spin_Chain(..., mode="ED") runs DMRG, as on the parent (the Thermal lines of 01 and 05, <Sz0 Sz1> -0.125704 at T=0.5 and -0.166667 at T=0 on both trees).

Supplementary evidence, the reviewer's own probe rerun on both trees by me: <scratch>/construction/05_review_attack_init.py (a copy of review/construction/04_attack_init.py), outputs 05_review_attack_init.before.out and .after.out. On the repaired tree every reader is identical on the ctor and assigned paths on "python" and v3, the round trip works on both, and the 12-site E(ctor)-E(assigned) is -2.220e-14 (v2), -8.882e-15 (v3), 0.000e+00 ("python"), against -9.635e-03, -9.635e-03, -9.632e-03 on the parent. Key lines, after:
```
[python ctor] mode='ED' session=Chain
[python ctor] get_dynamical_correlator KPM                 -> 0.205190
[python ctor] mode=None; gs_energy()                       -> -1.6282694864
[3 ctor] mode='ED' session=Chain
[3 ctor] get_dynamical_correlator KPM                      -> 0.205190
[3 ctor] mode=None; gs_energy()                            -> -1.6282694864
Thermal_Spin_Chain(T=0.5, mode='ED').get_gs()              -> MBChain.mode='DMRG', <Sz0 Sz1>=-0.125704
Thermal_Spin_Chain(T=0, mode='ED').get_gs()                -> MBChain.mode='DMRG', <Sz0 Sz1>=-0.166667
Fermionic_Chain(4, N=5).N                                  -> TypeError: Fermionic_Chain() cannot set N at construction: that is the chain's state, not a setting (N: the Fermionic_Chain class builds it)
Parafermionic_Chain(3, Sig=1).Sig                          -> TypeError: Parafermionic_Chain() cannot set Sig at construction: that is the chain's state, not a setting (Sig: the Parafermionic_Chain class builds it)
[2] 12-site E(ctor) - E(assigned)                          -> -2.220e-14   E(ctor)=-5.1476856698
[3] 12-site E(ctor) - E(assigned)                          -> -8.882e-15   E(ctor)=-5.1476856698
[python] 12-site E(ctor) - E(assigned)                     -> 0.000e+00   E(ctor)=-5.1476882590
```
before (parent):
```
[python ctor] mode=None session=Chain
[python ctor] get_dynamical_correlator KPM                 -> 0.205144
[3 ctor] mode=None session=Chain
[3 ctor] get_dynamical_correlator KPM                      -> 0.205151
[2] 12-site E(ctor) - E(assigned)                          -> -9.635e-03   E(ctor)=-5.1573204922
[3] 12-site E(ctor) - E(assigned)                          -> -9.635e-03   E(ctor)=-5.1573204922
[python] 12-site E(ctor) - E(assigned)                     -> -9.632e-03   E(ctor)=-5.1573204922
```

In-tree caller check (04_intree_callers.py in the cluster folder), rerun on the repaired tree; the before output is unchanged:
Before:
```
[run3] slot 1 acquired
dmrgpy from <scratch>/parent/src/dmrgpy/__init__.py
four_correlation_tensor example: fced.mode = None  session built: True
  sum|ct_ed| = 51.433292772677   E0(ED) = -2.450108689391
multioperator_density example: Fermionic_Chain(6, spinful=False) built
```
After:
```
[run3] slot 2 acquired
dmrgpy from <repo>/src/dmrgpy/__init__.py
four_correlation_tensor example: fced.mode = 'ED'  session built: True
  sum|ct_ed| = 51.433292772677   E0(ED) = -2.450108689391
multioperator_density example: TypeError: Fermionic_Chain() got unexpected keyword argument(s) spinful. A constructor keyword names a setting of the chain (maxm, nsweeps, noise, cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes effect as if assigned right after construction; a name the chain does not have would be stored where nothing reads it
```
(The first pass's after output, kept as 04_intree_callers.firstpass.after.out, had "session built: False".) The caller scan was an AST walk over src/, tests/, examples/ and benchmarks/ (scan_ctor_kwargs.py); the reviewer's own scan agrees.

Leads for other clusters, not fixed here:
- gs_energy(maxde=...) on a chain whose state is current returns the stored energy without refining (the [3, maxm=3] lines of 02), gs_energy's own short circuit running ahead of maxde= as it used to run ahead of wf0=. If the session cluster widens gs_energy's bypass, Many_Body_Chain.get_gs's condition must widen identically; the robust shape is one groundstate helper of (self, kwargs) that both read.
- get_gs(best=True, **kwargs) calls groundstate.best_gs(self, n=n, **kwargs), and best_gs(sc, n=1) takes nothing else, so any keyword next to best=True raises TypeError.
- On the lower-level get_dynamical_correlator_MB route a pair name= next to i=/j= is dropped silently on every submode checked (ROOTN, CVM), where the public route raises (second-pass finding 16).
````

### What the session reviewer reported beyond its items

````
1. This is the non-Hermitian twin of generalized-cache, not fixed by this cluster, and pre-existing on 8dd2198. gs_energy_generalized() on a non-Hermitian Hamiltonian dispatches to nhdmrg.gs_energy_generalized_nhdmrg() before the new send_hamiltonian() line, and neither NH solve puts H on the session. So the next correlator's groundstate.ground_state_on_session() finds H missing, resets computed_gs and re-solves plain NH-DMRG over the generalized state. Solving first does not help. Only a chain on which an earlier correlator had already sent H keeps the generalized state.

Repro: <scratch>/review/session/11_nh_generalized.py. It builds a 4-site Heisenberg chain plus 0.3j*Sz0 + 0.2*Sx1, with metric A = 1 + 0.8*Sz0, and runs three histories (fresh, solved, solved plus a correlator), each followed by an NH-KPM on (Sz0,Sz0) with E_max=10, n=50. The parent and fixed trees give the same output line for line. Fixed tree:
```
python fresh        lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=(-0-0j) (generalized (-0.4529+0j))
python solved       lambda=(-2.17581-0.213083j) H on session after the generalized solve: False after NH-KPM: e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 vev(Sz0)=-0j (generalized (-0.4529-0j))
python solved+corr  lambda=(-2.17581-0.213083j) H on session after the generalized solve: True  after NH-KPM: e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 vev(Sz0)=(-0.4529+0j) (generalized (-0.4529+0j))
3      fresh        ... e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 ...
3      solved       ... e0=(-1.596396-0j) |<wg|wf0>|^2=0.6877 ...
3      solved+corr  ... e0=(-2.17581-0.213083j) |<wg|wf0>|^2=1.0000 ...
```
The suggested fix is send_hamiltonian(self) at the end of gs_energy_generalized_nhdmrg(), and of gs_energy_nhdmrg(), on the session backends. It would also remove the extra NH re-solve that nh-injected left open.

2. These are record inaccuracies, not code defects, each detailed under its item. On ED, MBChain.gs_energy() after the thermal fix is the lowest eigenvalue (-1.000000) and not <wf|H|wf> (-0.416388), which the fix agent's expected statement and NUMBERS CHANGE do not say. The julia_live gs_energy(reconverge=...) TypeError is listed as a behaviour change, and it is not one.
````

### The session agent's notes

````
What changed in this repair pass. The reviewer found generalized-cache INCOMPLETE: the non-Hermitian route of gs_energy_generalized() still lost the generalized state to a plain NH-DMRG re-solve on every chain on which no correlator had run, solved first or not. I reproduced it on the parent (09_nh_generalized_cache.py) and fixed it in src/dmrgpy/nhdmrg.py: a new _record_hamiltonian_sent() at the end of gs_energy_generalized_nhdmrg() and gs_energy_nhdmrg() puts H on the session through groundstate.send_hamiltonian() when there is a session, and gs_energy_nhdmrg() skips it for an H= other than the chain's own. It also removes the second NH-DMRG solve the first correlator made after a plain NH solve (counted 1 to 0; spectrum identical to six digits), which the first pass had recorded as harmless under nh-injected. I also rewrote the CAVEAT in gs_energy_generalized()'s docstring (groundstate.py): its remedy said gs_energy() recomputes a genuine ground state, which it does not on either tree (10_generalized_edges.py: gs_energy() right after returns lambda), so it now says restart() then gs_energy(). The record corrections the reviewer asked for are made from reruns on both trees: 03 now carries v2, the ED gs_energy and sum rule, and T=1.0001e-5; 04 now carries v2; the julia reconverge line is struck from docs_needed. Tests added: the NH generalized rows (4), the NH solve count (2), the v2 rows of the thermal and NH wf0 tests, and the ED row of the thermal test now pins the KPM sum rule.
The MPO edge: on 8dd2198 gs_energy_generalized() on a StaticOperator Hamiltonian raised AttributeError ('StaticOperator' object has no attribute 'to_terms'); on v3 it now runs and gives -3.59799449, equal to the MultiOperator route (10_generalized_edges.py). I did not run it on "python" or v2; by reading, send_hamiltonian() raises NotImplementedError there unless the session has set_hamiltonian_mpo.
Things I could not verify or left alone. The reviewer saw an ITensorMPS deprecation warning (inner(x::MPS, A::MPO, y::MPS) with mismatched site indices) from excited_states_dmrg's Julia code in a julia_live excited-state run after set_gs; it was not attributed to a tree and I did not spend a parent JIT run on it; the energies matched ED. Julia's RNG and the C++ RNG of v2/v3 are not seeded by numpy, so the before and after runs of 01 and 04 use different x on those backends, and each line compares a result with its own run's x. The outputs of 01, 02, 05, 06 and 07 are the first pass's, since nothing they cover changed; 03, 04, 09 and 10 were run in this pass on both trees. The first-pass repro files carry no invocation comment lines, which the record shows appended to the script text as before; 04's file now has them.
Scratch folder: <scratch>/session/ holds scripts 01 to 10 with their .before/.after outputs, the two example runs and example_after/GS_enforce_maximum_fluctuation.png. No git state was changed, nothing under mpscpp2/mpscpp3 was touched, and in this pass the only edits were to nhdmrg.py, groundstate.py (gs_energy_generalized's docstring) and the cluster's test file.
````

### What the scale reviewer reported beyond its items

````
1. submode="TD" (and so every TDVP real-time route) depends on units on "python" and v3, through the Krylov exponentiator's absolute error goal, in the operators' norm as well as in the Hamiltonian's scale. On "python", pyitensor/tdvp.py::_lanczos_expm_multiply stops on `beta0*beta*abs(exp_col0[-1]) < errgoal` with errgoal=1e-10 absolute. beta0 is the norm of the local tensor, and pyitensor/chain.py::quench_tdvp evolves psi1 = B|gs> without normalizing it (it renormalizes to norm0 after each step), while beta carries the units of H. On v3, ITensor's applyExp has the same shape: ErrGoal 1E-10 on an estimate carrying nrm, plus NormCutoff 1e-7 on beta (by reading).

Measured, first with small operators. max|C[eps*Sz0,eps*Sz3]/eps^2 - C[Sz0,Sz3]|/max|C| on "python" is 6.66e-09 (eps=1e-2), 2.42e-06 (1e-4), 2.58e-05 (1e-6), 8.76e-02 (1e-8) and 2.42e-01 (1e-9); eps=1 run twice gives 0.0e+00. With the Krylov problem unit-normalized (monkeypatched in that process) it is 9.88e-13 to 1.66e-12 at every eps. On v3 it is 1.04e-04 at 1e-9; on v2 (MPO-Taylor) 3.83e-13. Then with a Hamiltonian in small units (s*H, with es, delta and dt scaled): "python" gives 6.66e-09, 2.42e-06, 2.76e-05 and 8.76e-02 at s = 1e-2 to 1e-8, and unit-normalized 8.04e-10 at 1e-4 and 3.58e-06 at 1e-8, the remainder consistent with the ground-state solver's own error. v3 gives 1.21e-09, 8.66e-08 and 6.97e-06 at 1e-2 to 1e-6. The problem is older than this fix: the parent gives the same 6.66e-09 / 2.42e-06 / 2.58e-05 at eps=1e-2 to 1e-6, and raised ZeroDivisionError at eps <= 1e-8. So the clean-threshold fix turns that loud failure into a silently wrong spectrum. Suggested fix: in quench_tdvp, evolve the normalized state and multiply the correlator by norm0. For the Hamiltonian's scale, the errgoal test needs to be dimensionless (gated if ordinary TD is to stay byte-identical). Repro: R/r4_downstream_floors.py, R/r5_td_krylov.py (+ .after.out), R/r5b_td_krylov_parent.py (+ .before.out), R/r6_td_small_units.py, R/r7_td_small_units_krylov.py, with R=<scratch>/review/scale.

2. The relative floor takes its scale from the raw, pre-summation coefficients, identity offsets included, both in _filter_small and in canonical_dict's floor (a). As a result, an operator exactly equal to 1e-7*Sz0 is proven zero and measures 0. Fixed tree: `1e6*Sz1 + 1e-7*Sz0 - 1e6*Sz1  to_terms=[((1000000+0j), [('Sz', 2)]), ((-1000000+0j), [('Sz', 2)])]  is_zero=True  simplify=[0.0]` and `vev(1e6*Sz1 + 1e-7*Sz0 - 1e6*Sz1)/(1e-07*<Sz0>) = -0.000000`. Parent: `is_zero=False  simplify=[(1e-07+0j)]` and `= 1.000008`. Likewise `1e6*Id + 1e-7*Sz0` now has `to_terms=[((1000000+0j), [('Id', 1)])]`, where the parent kept the field. This needs ratios of at least 1e12, so it is low severity. Suggested fix: drop exact zeros only at consumption points (raw terms are never rounding dust) and keep criterion (b) alone in canonical_dict, or at least exclude Id terms from cmax. Repro: R/r2_floor_edges.py (+ .before/.after.out).

3. meanfield counts the constant offset in its scale: `tol = _tol*_scale(b,J,const)` and `_scale(b,hmf,const,...)`. On a 4-site Heisenberg chain + 0.3*Sz_tot + c*Id, the fixed tree at c=1e11 gives `E_MF-c=0.000000  <Sz_i>=[-0.5 -0.5 -0.5 -0.5]`, where the parent gives `E_MF-c=-0.627274  <Sz_i>=[0.5 0.5 0.5 0.5]`; c=1e8 is fine on both. Suggested fix: leave const out of _scale, since an offset is not a coupling. Repro: R/r2_floor_edges.py part E.

4. Non-finite coefficients now empty the operator. Fixed tree: `nan*Sz0 + Sz1    to_terms=[]`, `Sz1 + nan*Sz0    to_terms=[((1+0j), [('Sz', 2)])]` (order-dependent), `inf*Sz0 + Sz1    to_terms=[]`. The parent kept Sz1 in both nan cases and kept the inf term, so it failed loudly downstream. Suggested fix: raise on a non-finite coefficient, or at least guard cmax with math.isfinite. Repro: R/r2_floor_edges.py part D.

5. Belongs to small-units, measured here for the first time: v2 and v3 give vev(eps*Sz0) exactly 0 at eps <= 1e-15 (`1e-14:1.000000 1e-15:-0.000000`), through the local `eps = 1E-14` behind isZero at autompo.cc:1287. ED and "python" are linear down to 1e-30. The cpp_plan's to_mpo_unit covers it.

6. Minor, recorded as leads. The TD correlator's default dt=0.1 is absolute, with nt = damping_periods/(delta*dt), so a Hamiltonian in small units at a scaled delta takes 1/s times the steps (my first r6 run did not finish in 600 s until I passed dt=0.1/s). effectivehamiltonian.fit_matrix and get_effective_hamiltonian_couplings drop couplings below an absolute tol=1e-4 (by reading, not measured); that is a documented keyword, but its default is in energy units.
````

### The scale agent's notes

````
Scratch: <scratch>/scale/ holds 01, 02, 03, 05, 06 (.py with .before.out/.after.out; 05 only after, one julia run as the brief allowed), cpp/04_davidson_exits.cc with its four binaries and per-run logs, cpp/patched/itensor/iterativesolvers.h (a scratch copy, the vendored header is untouched), cpp/04_patch.diff and 04_davidson_exits.out, and parent_check/ (the new test file run against the parent tree). The 04 binaries were compiled with /usr/bin/g++ against mpscpp3/ITensor's headers and lib/libitensor.a read-only (nothing under mpscpp2/mpscpp3 was edited or rebuilt, no make) and ran natively with OMP_NUM_THREADS=1 outside run3.sh, since they are not Python; each run takes well under a second. The julia_live run triggered juliapkg's automatic dependency resolution on import (registry update and Pkg.add/precompile of the declared ITensors, ITensorMPS, ITensorNHDMRG, PythonCall, OpenSSL_jll into ~/.julia/environments/pyjuliapkg, reporting no packages added or removed at the end), which is juliacall's standard behaviour when it sees two juliapkg.json files (site-packages dmrgpy and this checkout) and was not asked for; worth knowing for the julia cluster. v3 values at 5e-7 and 6e-7 differ between runs (-1.284 against -1.330 at 6e-7) because ITensor's randomMPS is not seeded by np.random.seed; every conclusion rests on the deterministic vev ratios and on repeated trends, not on those digits. I measured the v2 Davidson line only by reading plus the one-channel Ising degradation, not with a v2 diagnostic build. The meanfield helper `_scale` is new and module-private; `_mean_field_hamiltonian` keeps its signature, since tests/test_meanfield.py calls it positionally.
````

### What the scale-cpp reviewer reported beyond its items

````
1. Gate on the wrong quantity (left unchanged by the fix, not a regression). `unit_scale_up(max_abs_coef(terms))` and `to_mpo_unit` scale by the largest coefficient of the whole list, but svdMPO truncates only the bond-crossing channels, so an O(1) constant (sent as `('Id',site)`) or one-site field keeps a sub-cliff exchange on the unscaled path. Repro: review/scale-cpp/pA_gate.py, a 6-site S=1/2 Heisenberg chain, fresh chain per row, maxm=30, nsweeps=10. Fixed tree: s*Heis+1.0 at s=1e-7 gives (E-1)/s = -1.2499999913 on v3 and -1.2499999902 on v2 against ED -2.4935771281; s*Heis+Sz0 gives (E+0.5)/s = -1.2499999991 on v3 and -1.0000000006 on v2 against ED -2.0925528521; [vev(1+s*H)-1]/s/vev(H) reads 0.3333 at s <= 1e-7 on the hybrid and on the fixed tree alike. It needs a six to seven decade hierarchy with the small side below about 6e-7 (v3) or 3e-7 (v2). 2. KPM below s of about 1.5e-14 on v2 and v3: `scaled_hamiltonian`'s single-term `shift*Id` AutoMPO (v3 chain_session.h:11786, similarly :12072, cvm's `z*Id` :2039, and :1578, :11868, :11885) goes through plain `toMPO`, whose `isZero(coef,1E-14)` drops the band-centre shift 0.62*s. Repro: review/scale-cpp/pC_kpm_extreme.py and pC2_kpm_threshold.py. max|DMRG-ED|/peak is 4.1e-04 at s=1e-12 and at s=1.8e-14, and 1.2e+00 at 1.5e-14, 1e-14, 1e-16 and 1e-20, on both backends. This is the fix agent's left_open item 5, now measured.
````

### The scale-cpp agent's notes

````
Scratch: <scratch>/scale-cpp/ holds the first pass (02, 03, 10, 11, 13, 14, 15 with .before/.hybrid/.after outputs, hybrid/src, hybrid_check/, parent_check/). The repair pass is in scale-cpp/repair/: 20_mixed_scale.py with .before.out (parent), .hybrid.out, .firstpass.out and .after.out; cmp-identical copies of 02, 03, 10, 11, 13 with their .after.out on the final build; 21_tests.after.out, 22_tests.firstpass.out and 22_tests.hybrid.out (the test file against the first-pass and hybrid trees, from firstpass_check/ and hybrid_check/); 23, 24, 25 regression outputs; make_v2.log and make_v3.log; and firstpass/src, the working tree's Python with the first pass's two _dmrgcpp .so files, snapshotted before the rebuild (cmp-identical to the ones replaced, and every changed .py cmp-identical to hybrid's). Four trees ran 20: the parent masks s <= 1e-8 for the clean-threshold reason, the hybrid is the unmasked before of the whole item, the first-pass tree is the repair's own before, and the working tree is the after. The review's two gaps were handled differently, following its own options: gap (2), the KPM shift below the isZero floor, is fixed; gap (1), the gate, is narrowed with strict xfails, because probe 20 part C showed the mechanism is per bond (a weak link at s=1 fails the same way) and part A showed "python" fails it too, so the multi-site gate the review proposed would have been a partial cure of a different defect. One sub-claim of my first pass was struck by the review and is corrected in the code comment and in the prose: svdMPO's truncation cutoff is read from Args. The 1.5e-14 KPM boundary is the band-centre shift 0.62*s crossing isZero's 1e-14, as the review bracketed it. The comment line of 20 says [vev(1 + s*H) - 1]/s while the code computes vev(s*h + 1.0), the same quantity; I changed the code before the first run and left the comment. Nothing else was running dmrgpy during the rebuild (the other Python processes on the host were pyqula tests and old notebook kernels). No git state was changed; the .so files are gitignored, so files_changed is the four headers and the test file.
````

## Statements this pass overturns

- `CLAUDE.md`'s paragraph for the 2026-09-24 third pass: "the v2/v3 failure on a
  Hamiltonian in small energy units, which is not the MPO builder". It is the MPO
  construction after all, ITensor's own `toMPO` (svdMPO's absolute truncation, a
  different builder from the `"python"` one that pass's finding 13 fixed), plus
  Davidson's absolute 1e-10 randomization threshold below it (item 2).
- The same paragraph's "`multioperator.clean_threshold` dropping every term below
  an absolute 1e-8", now fixed (item 3).
- `gs_energy_generalized()`'s CAVEAT and the user guide: "recompute a genuine
  ground state (`gs_energy()`/`nhdmrg()`) first". A bare `gs_energy()` afterwards
  returns lambda, on `8dd2198` too; the remedy is `restart()` and then
  `gs_energy()` (item 6).
- `documentation.md` 4.10's "a convergence ramp over `sc.maxm` ... The key now
  lives in `groundstate.solver_key()`" holds on the session backends only, not on
  `julia_live` (item 4).
- The 2026-09-24c lead on `maxde`: its example and user-guide halves were already
  fixed on `8dd2198` (item 9). Its lead on `gs_energy_generalized`: the
  non-Hermitian route loses the generalized state on every chain with no earlier
  correlator, solved first or not, not only on a chain whose cache is empty
  (item 6).
- The 2026-09-24b finding 12 Status, "`julia_live` has the same shape ... left as
  it is", and the 2026-09-24 finding 10 Status, "the ROOTN sibling ... stays
  open" (items 4 and 10).
- The locate stage of item 2 said svdMPO's cutoff is not reachable through Args;
  its reviewer struck that (the cutoff is read from Args, and only
  `isZero(coef,1E-14)` and Davidson's `qnrm < 1E-10` are hardcoded).

## Shared helpers

The three-slot runner every Python process went through (`<scratch>/run3.sh`):

````bash
#!/bin/bash
# Runs one python invocation under a workflow-wide cap of three concurrent runs.
# Usage, from the folder holding the script:
#   <scratch>/run3.sh NN_slug.py [args] 2>&1 | tee NN_slug.out
#   <scratch>/run3.sh -m pytest tests/test_x.py -q 2>&1 | tail -30
# Waits (checking every 5 s) for one of three lock slots, then runs python3 with
# threads pinned and this checkout's src on PYTHONPATH. The slot is released when
# the process exits, including when it is killed.
LOCKDIR=<scratch>/locks
while true; do
  for i in 1 2 3; do
    exec 9>"$LOCKDIR/slot$i.lock"
    if flock -n 9; then
      echo "[run3] slot $i acquired" >&2
      MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
        MPLBACKEND=Agg PYTHONPATH=${DMRGPY_SRC:-<repo>/src} \
        python3 "$@"
      rc=$?
      exec 9>&-
      exit $rc
    fi
    exec 9>&-
  done
  sleep 5
done
````

The pre-fix snapshot, made before any agent started:

````bash
cd <repo> && git archive 8dd2198 -- src ':!src/dmrgpy/mpscpp2/ITensor' \
    ':!src/dmrgpy/mpscpp3/ITensor' ':!src/dmrgpy/mpscpp3/TDVP' | tar -x -C <scratch>/parent
for v in 2 3; do cp src/dmrgpy/mpscpp$v/_dmrgcpp*.so <scratch>/parent/src/dmrgpy/mpscpp$v/; done
# a before run:
DMRGPY_SRC=<scratch>/parent/src <scratch>/run3.sh NN_slug.py 2>&1 | tee NN_slug.before.out
````
