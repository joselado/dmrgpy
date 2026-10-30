# Repro scripts, night of 2026-10-30

The measurements behind "Night of 2026-10-30" in `docs/clik_mps_lindblad_plan.md`.
Run from this folder with `MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 <script>`;
the 5-site runs take 40 to 70 minutes on one core. `cases.py` holds the chains
(driven XXZ, Bose-Hubbard with dissipative hopping) and is imported by the rest.

- `count.py xxz|bh [errgoal]`: small-chain spectra against ED, with the Arnoldi iteration counts.
- `bh5.py 64 20 0.5 0.2 [errgoal [substeps]]`: the paper's 5-site check against ED, with snapshot timings.
- `bh5_head.py 64 20 0.5 0.2`: the same, printing how many modes the sector filter keeps.
- `spurious.py L N`: the charge variance of each Ritz vector, against ED (5 4 is the 5-site check).
- `eg_check.py L N goals` and `eg_sub.py L N goal nsub`: log-norm, trace and charge variance of single steps, at several error goals or with sub-steps.
- `stiff.py kappa dt on|off`: the 3-site case with stronger dissipation, with and without the sub-step gate.
- `defl_sector.py none|sector L N`: the second branch with and without per-step projection of the sector identity.
- `os_leak.py dt n` and `os_leak2.py dt n cutoff rr`: the one-site TDVP leak, with a rank-revealing split and a fine reference.
- `os_dbg.py nh|herm`: one-site and two-site TDVP from a full-rank in-sector state.
- `hybrid.py xxz|bh refresh`: two-site TDVP until the bonds stop growing, then one-site.
- `km_check.py`: batched against looped matrix elements on 5 sites.
- `gs_micro.py`: the Gram-Schmidt bookkeeping, columns against rows.
- `sector_id.py`: the first version of the sector identity MPS (now `vectorize.sector_identity_mps`).
- `guide_snip.py`: the user guide's open-system example against ED.
