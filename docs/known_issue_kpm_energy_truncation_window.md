# Known issue: `kpm_energy_truncate` with too small a `kpm_scale`

**Status**: partially guarded, not fixed. Affects `itensor_version=3` and
`itensor_version="python"`, only with the opt-in, off-by-default
`kpm_energy_truncate=True`. `itensor_version=2` has no energy truncation
at all (it now raises rather than ignoring the flag), and the ED path is
unaffected.

## What happens

KPM energy truncation (Holzner et al., PRB 83, 195115 (2011), Sec. III-B)
projects the high-energy components out of each Chebyshev vector, so the
recursion only has to resolve the part of the spectrum a correlator
actually needs. Its premise is that the correlator's spectral weight lies
inside the ground-state-anchored window `[E0, E0+Ws]` that `kpm_scale`
sets. `kpm_scale` is precisely the knob that shrinks that window, and the
whole point of the method is that shrinking it buys resolution.

Push it past the point where the weight fits, and the discarded weight is
not removed: it piles up at the window edges. The returned spectrum then
peaks at roughly 0.96x the top of the retained window with several times
the correct amplitude, and looks entirely plausible. On a 6-site
Heisenberg chain, `<Sz_0 ; Sz_0>` at `delta=0.1`:

| `kpm_scale` | peak | exact peak | `max|y|` |
|---|---|---|---|
| 0.7 (default) | 0.489 | 0.489 | 1.11 |
| 0.5 | 0.489 | 0.489 | 1.02 |
| 0.4 | 0.000 | 0.489 | 1.73 |
| 0.3 | 1.113 | 0.489 | 3.96 |

Nothing warns. The moments stay bounded by `mu0` (at `kpm_scale=0.4` they
become exactly periodic, `0.25, 0.141, 0.25, 0.141, ...`, i.e. all weight
sitting at the rescaled band edges `x=+-1`), so neither
`pyitensor/chain.py`'s `_check_kpm_moment` nor `mpscpp3`'s
`check_kpm_moment` can fire -- their job is to catch moments that grow
past the bound, and these do not.

It is a cliff, not a gradual degradation. Between `kpm_scale=0.60` and
`0.55` on the chain the audit measured, the fraction of spectral weight
inside the window goes from 92.4% to 90.3% -- two percentage points -- and
the answer goes from correct to garbage, at a constant number of moments.
Under-convergence is not involved: the untruncated path at the same
`kpm_scale` gives the right peak, and the numbers above are unchanged at
the recommended `kpm_truncate_dK=30` / `nsweeps=10`.

The algorithm itself is sound where its premise holds. On a 7-site chain
with a decoupled `20*Sz[6]` term (bandwidth inflated, correlator width
unchanged) truncation reproduces the exact peak down to `kpm_scale=0.15`.

## What is guarded now

Both backends compute Holzner's Eq. (40)/(41) diagnostics
(`avg_truncated_weight`, `state_change_norm`) inside the truncation and
used to discard them at the call site. They are now checked: if a single
truncation removes more than half of the Chebyshev vector it was given,
the calculation raises instead of returning a spectrum.

```
RuntimeError: kpm_energy_truncate: the energy truncation removed 100.0% of
a Chebyshev vector, i.e. essentially all of it -- the retained energy
window is far too small for this correlator ...
```

That catches the unambiguous end of the failure: healthy runs on the same
chain sit at a relative state change of 0.01-0.17, while the broken ones
annihilate the vector outright (0.67, 0.76, 1.00 in the table above).
`kpm_scale=0.4` and `0.3` now raise on both backends; `0.7` and `0.5` are
untouched.

## What is still silent

The band between "clearly fine" and "vector annihilated". A run that
discards, say, 10-15% of the weight can still shift the peak while the
relative state change stays modest, and no threshold on this diagnostic
alone separates that from a legitimately aggressive truncation. Fixing it
properly means comparing the retained window against where the
correlator's weight actually is -- a spectral estimate the method does not
currently make.

Practical guidance until then: **treat `kpm_energy_truncate=True` as
requiring a check against the untruncated result at the same
`kpm_scale`**, rather than as a free accuracy knob. If the peak moves,
the window is too small. Note also that truncation is measurably *slower*
than not truncating on both backends here (see the ROADMAP entry), so it
is a paper-fidelity/regime feature rather than a default worth reaching
for.

## Where the code is

- `src/dmrgpy/pyitensor/chain.py::_maybe_energy_truncate` (the guard) and
  `src/dmrgpy/pyitensor/kpm_energy_truncation.py::energy_truncate` (the
  diagnostics)
- `src/dmrgpy/mpscpp3/chain_session.h::check_kpm_truncation`, called from
  `kpm_moments_truncated_full` / `kpm_moments_truncated_accelerated`
- dispatch: `src/dmrgpy/kpmdmrg.py::dynamical_correlator_moments`

Found by the 2026-08 cross-backend audit; the full reproduction is in
`docs/audit_2026_08_hole_hunt.md`.
