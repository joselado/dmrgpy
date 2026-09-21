---
name: hole-hunt
description: Run a multi-lens hole hunt (audit) of the dmrgpy Python layer, or add a lens, a finding or a fix cluster to an existing one. Use this whenever the user asks to audit, hole-hunt, hunt for bugs, sweep for holes, cross-check the backends against each other, or look for silently-wrong numbers, and also when they ask to record or fix a finding in docs/audit_*_hole_hunt.md. This is the repository's established audit process, run in 2026-08 and 2026-09, and it has a fixed record shape and a fixed evidence standard that a free-form bug hunt will not reproduce.
---

# Hole hunt

A hole hunt is a parallel, multi-lens search for behaviour in the dmrgpy Python
layer that is *silently* wrong: a number that is plausible and incorrect, a
dispatch that answers a question nobody asked, a kwarg with no consumer. The two
previous hunts are `docs/audit_2026_08_hole_hunt.md` (five lenses, 21 findings)
and `docs/audit_2026_09_hole_hunt.md` (eight lenses, 36 findings). Read the
scope section and the lens table of the most recent one before starting: a
finding already recorded there is not a new finding.

What makes this process worth following rather than improvising: every claim in
the record was *executed*, and every claim was then handed to a second agent
whose only brief was to refute it. That is what makes the record trustworthy
enough that a later fix does not have to re-derive the evidence. Predicted
output, remembered output and output from a stale `.so` all break that, so they
are the failure mode to guard against throughout.

## 1. Fix the frame

Record, and keep true for the whole hunt:

- The commit (`git rev-parse --short HEAD`) and that the tree is clean.
- Whether both compiled extensions are current. If a lens will touch C++, rebuild
  first, because a fix that lands mid-hunt and replaces `_dmrgcpp*.so` invalidates
  every other lens's measurements. The 2026-09 hunt took its one C++ fix by hand,
  separately, for exactly this reason.
- The invocation every repro uses:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=<this worktree>/src python3 <script>
```

Both halves matter. Unpinned threads make a timing claim meaningless, and a bare
import resolves to whichever checkout `site-packages` is symlinked into, so an
ad-hoc probe can silently test code that is not the code under audit.

## 2. Choose the lenses

Four to eight, each a one-line brief naming one class of problem, and as
file-disjoint as you can make them so the fix clusters afterwards can run in
parallel. Previous sets, to vary rather than repeat:

- 2026-08: silent backend and method dispatch fallbacks, dropped keyword
  arguments, feature-by-feature combinations, garbage-in-garbage-out cases that
  should raise, and cross-backend numerical disagreement on chains small enough
  for ED to be exact.
- 2026-09: `python-backend-parity`, `wavefunction-consumers`, `dispatch-matrix`,
  `pyitensor-performance`, `cpp-v3-completeness`, `ed-and-operators`,
  `recent-commits`, `docs-examples-drift`.

Out of scope by construction, and stated in the record so the exclusion is on
the page rather than in someone's head: vendored ITensor (`mpscpp2/ITensor/`,
`mpscpp3/ITensor/`); the legacy bugs `CLAUDE.md` says are deliberately
reproduced (`evoloperator`'s z^3/6 term on `H2`, the `"moise"` key, the
unreachable `"tevol_fit_td"` branch); the open `docs/known_issue_*.md` items;
anything already in either audit record; and gaps `ROADMAP.md` already marks as
absent. Decide explicitly whether `itensor_version="julia_live"` is in scope,
since its juliacall JIT cost dominates any lens that touches it, and say so.

## 3. Hunt, then refute

Spawn one `hole-hunter` agent per lens, all in one message so they run
concurrently. Each returns candidate findings, each carrying a repro script that
was actually run and its verbatim output.

Then spawn one `finding-reviewer` agent per candidate, briefed to refute it. A
reviewer that reproduces the repro and finds the behaviour intended, already
documented, or an artifact of the probe returns `REFUTED`, and that candidate
never enters the record. A reviewer that narrows the claim returns the narrowed
version, and the narrowed version is what gets written down. The 2026-09 record
carries several findings whose sub-claims the reviewer struck; keeping the strike
visible is part of the point.

## 4. Write the record

`docs/audit_<YYYY_MM>_hole_hunt.md`, in the shape both existing records share:

````markdown
# Audit, <YYYY-MM>: <n>-lens hole hunt

<one paragraph: date, commit, tree state, that every repro was executed and
every finding handed to an independent reviewer briefed to refute it, and that
REFUTED findings are not reproduced here>

<one paragraph: this file is the evidence, not a task list; fixed entries keep
their repro and gain a **Status** line rather than being deleted>

## The <n> lenses

| Lens | Brief |
|---|---|

## Scope

<what is excluded by construction, and the pinned-threads invocation above>

## Findings

### 1. <the defect stated as a claim, with its measured size, in one sentence>

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `<lens-name>`

**Where**: `<file:line list, every site the defect reaches>`

<prose: what the code does, what the layer above believes it does, why every
existing test passes through it>

**Expected**: <what a correct implementation returns>

Repro:

```bash
<the exact pinned-threads invocation that was run>
```

Observed:

```
<verbatim output, not retyped>
```

**Reviewer (CONFIRMED)**: <the attempt to refute it, and what survived>

**Suggested fix**: <one paragraph>
````

A `**Status**` line goes directly under the classification line once the finding
has been acted on, and `**Reviewer on severity**` is the variant used where the
reviewer accepted the defect and disputed how bad it is.

Two things about the finding heading, because they are what makes the record
readable a year later: state the defect as a claim rather than a topic
("`vev(op, npow=n)` silently ignores `npow` on every ED route" beats "npow
handling"), and put the measured size in it where there is one.

## 5. Fix in file-disjoint clusters

Group the findings into clusters that do not touch the same files, and take one
cluster at a time or in parallel agents. Each cluster gets one regression file,
`tests/test_audit_<YYYY_MM>_<cluster>.py`, and each fix gets a `**Status**` line
appended to its finding: `FIXED`, `PARTIAL` (say which half landed), or the
reasoning if the behaviour turns out to be intended after all. Name the tests
that pin it, or say explicitly that none does.

Pin the property, not a golden number, wherever the property is what was wrong:
the 2026-09 TDVP fix is pinned by a test that asserts the *order* in `dt` rather
than a value. Where a fix removes a bug from an existing computation, keeping
the pre-fix reference construction verbatim inside the test is the cheapest way
to prove the new path agrees with the old one where it should.

## 6. Say when numbers change

A fix that makes a previously-returned number different is a different kind of
event from a fix that makes a crash stop, because saved results elsewhere are
now not comparable. Every such fix gets `NUMBERS CHANGE` in its `**Status**`
line, naming the old value, the new one and the exact chain it was measured on,
and the consolidated list goes into `CLAUDE.md`'s paragraph for that audit.
Other projects save results produced by this library, so this is not
bookkeeping.

## 7. Close the loop

Update `CLAUDE.md` with a paragraph for the hunt: how many lenses, how many
findings, where the regressions live, which fixes changed numbers, and which
items are open rather than fixed. An open item stays in the record with what is
known about it, the way the `kpm_energy_truncate` window problem and the
`submode="TD"`/`"TDZ"` convention items did, rather than being dropped because
it did not get fixed.
