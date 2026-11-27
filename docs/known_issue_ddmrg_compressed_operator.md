# Known issue: the variational correction vector on long chains

**Status**: GUARDED on 2026-11-27, not fixed. `cvm_solver="variational"`
raises where it used to return wrong numbers
(`tests/test_ddmrg_operator_guard.py`). The fix below is open.

## What goes wrong

The variational correction vector (dynamical DMRG, Jeckelmann) minimizes
W(x) = <x|M|x> - 2 Re<b|x> with M = (H-E0-omega)^2 + eta^2 and
b = -eta B|GS>. Its minimum gives the spectral weight, which is not
negative for the pair (B^dagger, B) this route serves. Both v3
(`Chain::ddmrg_correction_vector` in `mpscpp3/chain_session.h`) and
`"python"` (`pyitensor/chain.py`) build M as an MPO: the sum H - c, then
the product of that with itself, then the sum with eta^2, all compressed
at the DMRG cutoff.

From about 16 sites this compressed M is not positive any more, so W has
no minimum and the local conjugate gradient meets negative curvature. On
two spins on two Wilson chains (`ysr.spin_chain_on_metal`, Js = 0.3,
Lambda = 2, maxm 64), the error of <GS|M|GS> in units of eta^2 = 4e-4,
at omega = 0 and cutoff 1e-10 on v3, measured by a second model (Fable):

| sites | 8 | 12 | 16 | 20 | 28 | 36 | 44 |
|---|---|---|---|---|---|---|---|
| error / eta^2 | 0 | 0 | -11 | -75 | -271 | -357 | -381 |

On 44 sites <GS|M|GS> is -0.152 against +4.0e-4, and `"python"` gives
-0.131. The 44-site references of `docs/new_methods/chains_on_metal/`
came out negative at every bias for this reason. A smaller cutoff is not
a fix: on 44 sites the error is -92, -35, -11.0 and -10.7 eta^2 at
cutoff 1e-11, 1e-12, 1e-14 and 1e-16, and it stays there at 1e-20.
Whether the product or the sums set this floor was not isolated.

## The guard

Each backend compares <GS|M|GS> with ||H|GS>||^2 - 2c<H> + c^2 + eta^2,
c = E0 + omega, from the uncompressed H (an exact apply on `"python"`,
ITensor's exact `innerC(A,x,B,y)` on v3), and raises above a tenth of
eta^2 plus 1e-10 of the terms that cancel. `cvm.py` raises on a weight
below -1e-6 of ||B|GS>||^2/(pi eta) on every backend, `julia_live`
included, whose `ddmrg.jl` has no operator check.

## The proposed fix, not done

Build H^2 as the uncompressed product, each site tensor multiplied over
the middle physical leg and the two links of each bond merged into one
(a combiner on v3). It does not depend on omega, so it is built once per
call. On `"python"` this gives <GS|M|GS> right to 1e-13 on 16 and 44
sites. Then keep the shift and eta out of the MPO, with c = E0 + omega:
M_eff v = H2_eff v - 2c H_eff v + (c^2 + eta^2) v, from the environments
of H2 and of H. This second step was not run.

It is necessary and not shown to be sufficient. The uncompressed M costs
8 times more per sweep on `"python"` (163 s against 20 s on 16 sites at
cvm_maxm 16). With it, 7 of 30 local solves on 16 sites stop at the cap of
200 CG iterations in one sweep. On 20 sites with M good to 0.02 eta^2,
cvm_maxm 32 stays 30% (omega 0) and 7 times (omega 0.3) below TASK after 4
sweeps, for a reason that is open. Not measured: more than one sweep with
the uncompressed M, cutoff 0 on v3, `ddmrg.jl`, and ED above 8 sites.
