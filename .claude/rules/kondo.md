---
paths:
  - "src/dmrgpy/kondospectrumtk/**"
  - "tests/test_kondo_spectrum*.py"
  - "tests/test_kondo_screening.py"
  - "src/dmrgpy/wilsonchain.py"
  - "src/dmrgpy/ysr.py"
  - "tests/test_ysr*.py"
  - "src/dmrgpy/ysrtk/**"
  - "src/dmrgpy/kondoscreening.py"
  - "src/dmrgpy/kondocloud.py"
  - "tests/test_kondo_cloud.py"
  - "tests/test_atom_iets.py"
  - "examples/kondo/**"
---

# STM/Kondo tunneling spectra (`kondospectrumtk/`)

These are the notes CLAUDE.md used to carry as its "STM/Kondo tunneling
spectra" section; they load when a file under `kondospectrumtk/`, a Kondo
test or a Kondo example is read or edited.

`get_kondo_spectrum` (`Spin_Chain` and the spinful fermion chains, see the
last paragraph) is the third-order perturbation theory of
Ternes, NJP 17 063016 (2015), arXiv:1505.04430; `mode="ED"` (full
diagonalization, any T) is the reference implementation and `mode="DMRG"`
(T=0 only, except `submode="TASK"`, see below) reproduces its two T=0-reducible pieces via the dynamical
correlator and the third one via a two-time TDVP construction. Read
`conductance.py`'s module docstring before touching a prefactor: every
normalization there was pinned against the paper's absolutely scaled Fig. 7
(arXiv v1 numbering; the docstrings' "Fig. 3"/"Fig. F" are that figure and
Fig. 5), and `tests/test_kondo_spectrum_paper_fig7.py` holds
pixel-digitized values from it (zero-bias peaks, +-4 mV tails, the 10 T step
asymmetry) that agree to ~0.005, except the in-field Kondo term, where the
code deliberately departs from the paper (see below). Two things there are
easy to get wrong and were, until 2026-09-12: the potential-interference
term's direct and exchange diagrams enter with OPPOSITE signs (a symmetric
electron trace leaves the reversed order's hole-like sign flip
uncompensated, unlike the Levi-Civita Kondo term), so it vanishes
identically at B=0; the summed form tilted the zero-field peak and made the
10 T asymmetry 5x the figure's. And F(eps,T) is the paper's symmetric closed
form eq. 22 thermally broadened with Theta', not the electron-like defining
integral eq. 20 applied to both diagrams, which differs at O(|eV|/omega0)
(3.5% at |eV|=omega0/5) and carries a sharp band-edge singularity at
eps=omega0. On the ED side the third-order sums run over thermally occupied
initial states only (O(n_occ dim^2), not dim^3) and F is tabulated once per
FBuilder (0.4 s; it used to be a per-point quadrature that took 28 s and 4
GB at 64 states).

**The Kondo term's exchange diagram departs from the paper (2026-09-26).**
The paper's eq. 25, and the formula under its Fig. 6, put the
reversed-order diagram's log at F(eV-eps_mi) = F(eV+eps_im). A brute-force
second-order T-matrix of the paper's own Hamiltonian, built as explicit
fermion x impurity operators with every intermediate state's denominator
read off H0 (`tests/test_kondo_spectrum_tmatrix.py`), puts it at
F(eV-(e_f-e_m)): the sample vertex acts first and creates the outgoing
electron, so the hole goes on shell at an outgoing energy of e_i-e_m, and the
outgoing electron has eV-eps_if. The two differ only when f != i, so the U
term (I_fi forces f=i), every zero-bias value and a free S=1/2 at B=0 are
unaffected, which is why none of the paper's figure values caught it; its
own 10 T curves were drawn from eq. 25, and
`test_kondo_spectrum_paper_fig7.py` now checks them against a paper-form
reconstruction (`_paper_form`) and pins the code's own values separately.
Do NOT "restore" eq. 25 to get back on the figure. The two-time DMRG/ED
route cannot express the new argument along G's t2 axis, so it builds a
second, sheared function Gx(s,tau) = G(-s,tau+s) =
<psi_l(-s-tau)|Sk|psi_j(-s)> from checkpoints of the same trajectories, and `kondo_term_from_two_time` takes (t, G, Gx) triples with
the one-sided kernel `K_F` (K_W is gone; the s->t direction reuses
conj(K_F), since its cell-averaged quadratures were most of the cost).

Since the 2026-09-24 audit: `mode="ED"` raises `TypeError` on any keyword
it does not read (a misspelled `Jrho_s` used to remove the Kondo peak
silently); `mode="DMRG"`'s T=0 is the single converged state unless `n_gs=`
asks for the equal-weight manifold average, for which `set_gs` alone
suffices, since it hands each member to the session unswept on the next
read (it used to push with `session.set_wavefunction`, which re-swept the
member, 2026-09-24 second-pass findings 13 and 14); `n_gs>1` takes only the
submodes that read the chain's state; and the potential term needs `es`
past the top of the S_k spectrum, not just past max|eV| as the second-order
term does, and uses trapezoid weights, so a non-uniform `es` is fine (the
Kondo potential term on a non-uniform `es` and the Kondo terms on a
non-increasing `es` both moved numbers when they were fixed, 2026-09-24
finding 13 and second-pass finding 17). `submode="CVM"`'s CG tolerance used
to be absolute, so the Kondo route at its default broadening returned a flat
curve of about 1e-8 where ED has a peak; it is relative to ||b|| since the
2026-09-25b pass.

**`submode="TASK"` (2026-10-22)** is a third T=0 DMRG route for all three
terms: `kondospectrumtk/taskkondo.py` builds a `KondoSpectrum`-shaped
object from the tangent-space Krylov spectrum (`pyitensor/task.py`,
arXiv:2510.07279) and feeds it to `conductance.py` unchanged, so every
prefactor and the exchange-diagram argument are inherited from the ED route
rather than reimplemented. Exact at saturation (`cutoff=0`) once the Krylov
space is exhausted (2e-11 on a 10-site S=1 impurity chain at `nkry=400`);
`itensor_version="python"`, 3 and `"julia_live"` (`taskdmrg.BACKENDS`;
v3 exports its site tensors to the same Python code,
`taskdmrg.session_spectrum`), the default on all three, and the correlator
route (order=2 only) on v2; the near-zero Ritz energies are snapped
to 0 or the zero-bias point is rounding-dependent on a degenerate ground
state. It is in `_N_GS_SUBMODES`, checked on a Kramers doublet. A Ritz
energy below the snapping window means the state handed over is not the
ground state (DMRG stopped at a mixture of a near-degenerate doublet), and
`TaskKondoSpectrum` warns (2026-10 audit, finding 3).

**The TASK spectrum is batched and cached (2026-11-10).** Four changes,
none of which moves a number (1.7e-10 on a 24-site chain, 117 s to 28 s):
`TangentSpace.apply_many` applies an operator to a stack of tangent
vectors, and `block_krylov` and `projected_spectrum` use it; the operator
MPOs go through `task.compress_mpo` before the P S P matrices are built
(AutoMPO gives a one-site operator bond dimension 2); `task.real_starts`
removes a global phase from each start when the state and H are real, so
Sy no longer makes the basis complex; and `taskdmrg.session_spectrum`
keeps the last `SPECTRUM_CACHE_SIZE` spectra on the chain
(`_task_spectrum_cache`) under a digest of the exported arrays. The key is
the content on purpose: do NOT replace it by a flag of the chain, since
`n_gs` sets each member with `set_gs` and a flag would return the first
member's spectrum for all of them. The cache serves `n_gs=1` only: with
`n_gs>1` every call gets new members from `get_excited_states`, so every
call builds again (measured, 2 builds per `n_gs=2` call). The MPOs of the
starts and of H are left uncompressed, which is not tested either way. At
`maxm=100` the factor is 3.9 (172 s to 45 s, with `BATCH_BYTES` at 256 MB;
1 GB was slower and doubled the peak memory). Two fresh v3 solves of that
chain at `maxm=40` differ by 3e-3 of the peak with the same code, so
compare old and new code on one exported state, never across runs. On
d22 (40 sites): 691 s to 217 s at bond dimension 100 and 3701 s to 1350 s
at 200, on one core, with less memory (3.1 to 1.9 GB, 10.7 to 5.8 GB); and
1477 s on the RTX A2000 at 100, so the card loses by 7 in float64. The gain is on NumPy
only: on JAX old and new take the same time, and the new code's JAX host
memory is higher (4.6 GB against 2.2 GB on 24 sites), which is open.
`tests/test_task_batched_spectrum.py`.

**Unconverged TASK poles are measured and warned about (2026-11-10).**
`TaskKondoSpectrum.residual` is the Lanczos residual of every Ritz state
(`block_krylov(residual=True)`), `convergence(emax)` the count, weight and
largest step of the poles that have one, and `get_kondo_spectrum` warns
above `taskkondo.UNCONVERGED_STEP`. At T=0 a pole with a residual is a
sharp step at an energy that moves with `nkry` and with rounding of the
solve: two solves of one 24-site chain gave dI/dV different by 3e-3 of the
peak (1.3e-2 from `nkry=200`), and JAX against NumPy by 2e-5 on one state,
old and new code alike, with the converged poles equal to 4e-17 eV. This
is the Krylov error, not a bug, and a larger `nkry` shrinks it without
closing it on a long chain. Do NOT smooth the steps by the residual
(tested: the difference between runs went down by 4 and the rms error
against `nkry=300` doubled). Seeding the solve was not tried; it would hide
the difference between runs, not the error. What was measured to remove it
is the thermal step of the leads with the sample in its ground state
(2.9e-6 between two solves at 0.3 K, full third order), valid while kT is
far below the lowest pole. That is `T>0` on the TASK route since
2026-11-10 (the author's pick: "T with a guard"): `TaskKondoSpectrum(T=)`
keeps `p` on the ground state and raises `ValueError` when the Boltzmann
factor of the lowest pole above the ground state exceeds
`THERMAL_OCCUPATION` = 1e-3. Do NOT turn the raise into a warning or
populate the Ritz states thermally: the tangent space holds the
excitations of the ground state, not those of an excited state, so its
sums from an excited row are incomplete. The guard (`_check_cold`) takes
the lowest pole minus its residual and wants it converged (a Ritz energy
is an upper bound: 0.47 meV at `nkry=3` for 0.29 meV), counts the poles
within 1e-3 kT of zero as the manifold of `n_gs` and raises when there
are more than `n_gs`-1 (a truncated partner sits at 1.5e-9 eV, not at 0),
and raises on a reached pole below the state. Keep all three. At T>0 the warning counts a pole only when its third nearest pole
is further than 5.4 kT (`THERMAL_WIDTH`), a rule fitted on one 24-site
chain, and weighs it by min(1, d/(5.4 kT)), with d = min(r, r^2/g) the
distance a Ritz energy with residual r can still move, g being the
distance to the nearest other Ritz energy (2026-11-12, the author's pick
after the 2026-11 hunt: a gapped 20-site chain at 0.1 K warned of 4.1e-3
for a true change of 6.5e-8; now 3.0e-5). Measured on two chains, the
size stays between the true change and the old size, and a uniform
24-site chain still warns where the true change is near 1e-3. T=0 is
unchanged; there the size is that of one step and the whole curve on
pole-resolved biases can move more (0.34 against 6.4e-3 on the gapped
chain), which is left as it is. The other DMRG submodes still raise at T>0. In `_overlap_columns`,
`row_overlaps` and `project_rows` a real store is never multiplied by a
complex array directly, since NumPy then copies the store as complex
(4.7 GB against 1.9 GB on 40 sites at bond dimension 100).
`tests/test_kondo_task_convergence.py`.

**The 2026-11 hunt (2026-11-12)**, `docs/audit_2026_11_hole_hunt.md`,
left these rules. `conductance.second_order_dIdV` sums over the occupied
initial states with the relative `P_CUT` of the third-order terms
(`p_cut=0` is the full sum) and builds its Theta arrays per block of
`CHUNK_ENTRIES`; it was O(n_eV dim^2) in time and memory (6.5 GB at 904
TASK states and 241 biases). `third_order_kondo_dIdV` builds its
coefficient tensor per block of occupied states under `COEFF_BYTES`, with
BOTH tunnelling directions inside each block and the exchange diagram
accumulated per group of biases under `EXCHANGE_BYTES`: one block is the
old code path, and building the blocks once per direction cost 36 to 46%,
so do NOT split them by direction. `FBuilder.__call__` evaluates F0 only
outside its table; `Theta_prime` is a series below |x|=0.05, where the
closed form lost every digit; the table is accurate for omega0 >~ kT only
(2e-4 of F off for omega0 << kT, outside the paper's regimes, left as it
is, since widening it moves numbers); `stepfunctions.cached_fbuilder`
keeps the last four tables across calls, shared, so never mutate one. The
ED eigenbasis spin operator is sparse times dense. The leading-log flow
runs in real block arithmetic when the block form of the bare vertex is
real to `REAL_TOL` (the complex code stays for a field along y), and
`_hole_sum` reads r and omega0+|x| from tables built once per flow, up to
`HOLE_TABLE_BYTES` (256 MB), so the 9-spin call with 241 biases is the
same as before; still one core. A site with no spin gives the potential
scattering 2*pi*4*U^2 on the LL route, as on order 2 and 3: the stop needs
a strict crossing and a zero bare vertex gives a coupling ratio of 0. The
two-time G comes from three bra trajectories of S_l|GS> run once, so a t2
row costs only the three trajectories of Gx; below saturation this moves
the term by the TDVP truncation error (up to 1.8e-3 of the ED peak at
`maxm=4`), with the error against ED the same size as before. The grid of
the two-time term is validated before any TDVP step, and the correlator
route warns when `es` cuts the elastic line or stops below max|eV|. In
`get_kondo_spectrum` the mode is validated first, and `n_gs` before
anything is built. Not applied: computing Sz once for the second-order
and potential terms of the correlator route (4% of the call).

**Spinful fermion chains (2026-10-22).** `get_kondo_spectrum` lives on
`kondospectrumtk/kondomethod.py::KondoSpectrumMixin`, inherited by
`Spin_Chain`, `Spinful_Fermionic_Chain` and `Spinful_Fermionic_Chain_Native`;
nothing in the route reads the site statistics, only `Sx/Sy/Sz[site]`, the
Hamiltonian and the session. `edkondo` takes its spin matrices from
`edobj.get_operator()`, so a conserved sector restricts them (Nf=n makes a
6-site Hubbard chain 924 states instead of 4096) and an Sz sector raises.
TASK agrees with ED to 1e-11 on both classes at saturation; at maxm=12 on 6
sites the interleaved class missed by 3.5e-2 and the native by 8e-3 (at
equal maxm, which is not the same truncation on twice the bonds, so the
two-site S_site there is a likely cause, not an isolated one). `tests/test_kondo_spectrum_spinful.py`,
`examples/kondo/kondo_third_order_hubbard_VS_heisenberg`. The potential
term's `es` coverage check compares against the measured
`sum_k <S_k S_k>` on a fermion orbital (below 3/4 by the empty and doubly
occupied weight), not S(S+1). On a chain too short for the session the
correlator route's Kondo term is the exact ED one, and `order=3` with a
correlator submode raises on v2 and `julia_live`, which have no TDVP for
the two-time term (2026-10 audit, findings 2 and 4).

**Kondo screening on a Wilson chain (2026-10-30)** is a different method
that shares only the folder of its examples (`examples/kondo/kondo_screening_*`):
`wilsonchain.py` (Wilson chain of a band, `anderson_chain()` for one or two
levels) and `kondoscreening.py` (the NRG-MPS diagnostics of
arXiv:2501.10746, `get_kondo_screening` on both spinful chains through
`KondoScreeningMixin`). Nothing in `kondospectrumtk/` is involved. Every
diagnostic is a vev (M_imp = 4<Sz_imp^2>, M_n = 4<Sz_imp Sz_n>, Schmidt
coefficients from the single-occupancy block of the impurity density
matrix), so ED cross-checks it; what limits it is `maxm`, see
`docs/development.md` 4.8d. `tests/test_kondo_screening.py`.

**Open Wilson chains (2026-10-30)**, arXiv:1611.05291, are in `wilsonchain.py`
(`open_wilson_chain`, `renormalized_wilson_chain`, `boson_bath_chain`,
`spin_boson_chain`, `oscillator_chain`, `anderson_chain(kind=)`): C2 must give
Re G_bath(0) of the continuum at every length for a spectrum that vanishes
at omega=0, which `tests/test_open_wilson_chain.py` pins to 1e-10; for a
finite density of states at the Fermi level the shifts are regularized
principal values and C2 is only approximate. See `docs/development.md` 4.8e.

**The real-space Kondo cloud (2026-10-30)**, `kondocloud.py`, is Holzner
et al., arXiv:0906.2933: `anderson_lead()` (dot at the end of a uniform
lead, L sites in total, native class by default) and `get_kondo_cloud`
on `KondoScreeningMixin` (Sigma(x), xi_a, T_K of eq. 4). Sigma(L-1) = 0
is exact only for a singlet, so L even at N = L; `half_filling=True`
(the Nf sector) is needed away from ed = -U/2 and did not converge on
`"python"` at L = 40. `tests/test_kondo_cloud.py`, `docs/development.md`
4.8d.

**The thermal Krylov route (2026-11-12)**, `get_kondo_spectrum(mode="ED",
submode="KRYLOV")` in `kondospectrumtk/krylovkondo.py`, asked for by the
TaSe2 session and moved from its scripts (`edk_thermal.py`,
`si13_cmf.py`): the initial states are the low states of the sparse ED
Hamiltonian, down to the relative Boltzmann weight `p_cut` (1e-8), and
each one gets ONE block Krylov space from S_x|i>, S_y|i>, S_z|i> of the
tip site, orthonormalized as a whole. Do NOT split it into three Lanczos
spaces: once Sz is not conserved they overlap and double count, 2.7
relative error on 5 spins in a 120 degree field (pinned by
`test_three_separate_spaces_double_count_with_a_field`). The space is kept
orthogonal to |i>, so [i] + Ritz states is orthonormal, and a
`StateSpectrum` with p one-hot on i feeds the three sums of
`conductance.py` UNCHANGED; do not reimplement them (the TaSe2 version had
U=0 only). Ritz energies within the degeneracy tolerance of E_i are
snapped to it, as on TASK. Exact once each space is exhausted: 1e-15 to
1e-12 against ED on every model of `tests/test_kondo_krylov.py`.
`su2=True` takes one member per multiplet with weight 2S+1 in the sector
of the smallest |Sz_tot|, checks [H, S_tot]=0 on a random state and
raises otherwise. mode="ED" accepts no other submode (it stays an unknown
keyword, TypeError), and order="LL" with it raises. Full ED refuses
above 2000 states (`algebra.maxsize`), so above 10 spins 1/2 this is the
only ED route. README: `docs/new_methods/kondo_krylov_thermal/`.

**Leading-log resummation (2026-11-06)**, `kondospectrumtk/leadinglog.py`,
is a new method documented in `docs/new_methods/kondo_leading_log/README.md`
(derivation, references, measured checks): one-loop poor man's scaling with
frequency-dependent couplings (Rosch, Paaske, Kroha, Woelfle,
cond-mat/0408506 eq. 6; cond-mat/0301106) in the impurity eigenbasis, every
vertex four impurity matrices (mu=0 potential, 1..3 exchange). The sample
vertex G(0) and the on-shell tip vertex A_on (element (x,y) at e_x-e_y) flow
together, as in the papers, and A(eV) is a quadrature over them (the first
commit kept the running A(eV) on the right-hand side, a different closure at
O(g^2) without a log; do not go back to it); the smooth
kernel 1/(D+r_x) - 1/(D+omega0+|x|) integrates exactly to F0, so its first
order IS `conductance.py`'s third order (tested as g^2 and odd-part g^3
scaling), and `SIGN=+1` maps g to `Jrho_s`. G(0) is not Hermitian; do not
symmetrize it. Fixed Gamma0, no level renormalization; NaN plus a
warning when the sample coupling or the per-bias tip coupling
(`info["coupling"]`) reaches order one, a warning above 0.5.
`tests/test_kondo_leading_log.py`. Public as `get_kondo_spectrum(order="LL")`,
`mode="ED"` only, at any T; `mode="DMRG"` raises because the TASK space does
not close under the products of spin operators the flow needs.

**The leading-log flow at T>0 (2026-11-12)** replaces the cutoff kernel by
its thermal average with the weight of `FBuilder`,
`thermalkernel.ThermalKernel`, passed to `_joint_flow` as `kernel` (tip)
and `onshell_kernel` (G and A_on); with neither the flow is the T=0 code
line for line, hole tables included. Its first order is the thermal third
order (odd part 8.011 when g is halved, against 1.95 for the shortcut
Gamma0 -> sqrt(Gamma0^2+(c kT)^2): do NOT replace the kernel by a broader
Gamma0). The amplitude is averaged and then squared, the on-shell flow is
thermal by default (`onshell="thermal"`; "T0" differs at second order), and
at T>0 a transition is open for the coupling check when p_i Theta is above
`open_cut` = 1e-2 of the largest p. The kernel rows of the exchange sum are
filled on one core, like `_hole_sum`, which makes the kernel 2.6 times the
T=0 one on 6 spins; what T really adds is occupied states, linear in their
number (`p_cut`), and memory (`memory_model`, `mem_gb` chunks the biases).
It came from a script on top of dmrgpy in a 1T-TaSe2 project (`llT.py`) and
agrees with it to 1e-15. README section 8,
`tests/test_kondo_leading_log_thermal.py`.

**The leading-log flow is one solve on the occupied columns (2026-11-11).**
`leadinglog._joint_flow` integrates G, A_on and the tip vertex at every
bias together, with no dense output, and `leading_log_dIdV` asks only for
the columns of the occupied states (`tip_vertex(cols=occ)`): d^2 per bias
and occupied state instead of d^3. A vertex is a 2d x 2d matrix
(`_block`), the `_C` contraction is its matrix product, and the hole
channel is the product with the electron index transposed (`_swap`), so
the flow of G and A_on is four 2d x 2d products per step. The tip
potential scattering does not flow (direct and exchange cancel), so the
flow runs at U=0 and U is added at the end. Below D_stop the quadrature
is the closed-form integral of the kernel. Numbers are those of the
two-pass code to 1e-11 (pinned as `_TWO_PASS` in the tests); 5 spins 1/2
went from 110 s to 1.6 s, and 9 spins with 241 biases take 392 s on d22
(658 s on one core, 1.7 GB). Do NOT go back to `onshell_flow` plus one
flow per bias: its dense output is 350 MB per step at 512 states;
`onshell_flow` stays as the reference the tests read. The card was
measured and loses in float64 (0.235 s against 0.121 s on six cores for
the four products), and complex64 is 3e-4 off. `_hole_sum` (numba, one
core) is what is left; numba threads gave 0.38 s to 0.04 s and were not
added, since they do not follow the BLAS thread settings. The author
chose to leave it on one core (2026-11-11); do not thread it unasked.

**The leading-log flow on a device (2026-11-12)**, `leadinglog_device.py`,
is the same flow with the state on a GPU, taken by `_joint_flow` when
`pyitensor.backend.set_backend("jax")` is on (`leadinglog._on_device`,
which reads the backend module only if it is already imported, since
importing the pyitensor package costs 0.8 s). It is a port, not a new
method: the right-hand side in jax.numpy with the layout of the state of
`_joint_flow` (real block form or complex (4,d,d)), and SciPy's DOP853 as
`_Dop853` with the state on the device (13 stages, at most 20 copies
against about 24 in SciPy; the scalars, the brentq of the stop and the
thermal tables stay on the host). So the numbers are the NumPy ones to
roundoff: below 1.1e-15 on JAX's CPU device and 8e-15 on the RTX A2000
against one core. The step sizes agree only to about 1e-7 of themselves,
since the DOP853 error estimate is a difference of nearly equal stages and
its last digits depend on the order of the sum; do not test the device
integrator for identical steps. Keep the two right-hand sides in step: a
change to the equations of `_joint_flow` or `_tip_bracket` goes into
`_DeviceFlow` too, and `tests/test_kondo_leading_log_device.py` compares
them on the real and complex routes, at T>0 with both on-shell kernels, at
a stop and with chunks. The tip part has two forms (`TIP_FORM`): "columns"
is `_tip_bracket` line for line (the electron index summed first, an
operand of 4 d^2 numbers per occupied column), and "block" applies the
kernel to the vertex and does one product per bias for all the occupied
columns, twice the arithmetic but the exchange kernel built once per bias
and nothing that grows with their number. "auto" takes "block" from
`BLOCK_FROM` occupied columns. With `mem_gb` the biases are chunked by
`leadinglog_device.memory_model` on the device. The exchange sum runs on
the card there, which does not touch the author's one-core rule for
`_hole_sum` on the host.

**The tip coupling check reads open transitions, with a floor
(2026-11-11).** `info["coupling"]` is |Jrho_s| times the largest
flowed/bare ratio over the transitions that are open at the bias
(`Th > 0`, the ones in dI/dV), with the bare vertex bounded below by
`BARE_FLOOR` = 0.1 of the largest one. The author asked for it ("make it
work and give finite results", "with some regularization if needed")
after the old check, over every final state down to 1e-8, made 211 of
241 biases NaN at 9 spins and Jrho_s=-0.05. Both parts are needed and
both are pinned: a closed transition with a large bare vertex set the
check on the S=1 test model, and open transitions with a bare vertex
near zero still gave a NaN on 4 spins 1. Do NOT lower the floor to a cut
like 1e-3: a forbidden transition gets a vertex of order g ln times the
large matrix elements (4 to 10% of the largest bare vertex on the dimer
at Jrho_s=-0.05), which is under control and far above 1e-3. The values
of dI/dV do not change, only which biases are NaN or warned. The NaN
mask around D_stop still takes every energy difference of the spectrum,
so a stopped flow on a dense spectrum is all NaN; left as it is. Section
7 of the method's README has the measurements.

**A Kondo spin on a superconductor (2026-11-06)**, `ysr.py`, is the
composite-operator dI/dV with a BCS Wilson chain (uniform onsite pairing
on every site), plan in `docs/superconducting_substrate_dIdV_plan.md`,
notes in `docs/development.md` 4.8f and 4.8g. Three things are
load-bearing. The impurity is native spinful sites with a charge penalty
(2S of them for spin S, bound by a ferromagnetic Hund coupling, and
`fc.impurity` is a list), because only that class has Sz sectors, and
E_ysr = E_free - E_screened is a few sector solves (Sz=0 and Sz=1 for
S=1/2: even total parity, an odd bath with the bound quasiparticle, and
odd total parity; every 2Sz up to 2S+1 with Dz or B);
at Delta=0 an Sz-only solve keeps the N it started with, so
`subgap_energies` raises there. On the default mesh (`mesh="mean"`) pass
J*`lambda_factor(Lambda)` to stand for the continuum J; on `mesh="zp"`
(Zitko-Pruschke, arXiv:0810.4765) pass J bare, and `fc.dos_factor` is the
factor that applies, which `ysrtk/zaverage.py` reads for its continuum
normalization (it was A_Lambda too large on zp before it did). On the mean
mesh a coupling v on the chain acts as v/A_Lambda (checked on
the classical YSR energy to 5e-4 of Delta, and on the two-loop slope of
the measured T_K); the first version divided, and was 25 times off.
And use Lambda=4: at Lambda=2 on 38 sites maxm=128 was not converged at the
scale of Delta. Averaged over four twists z, Delta_c/T_K = 4.41 at the
paper's J=0.2 (rho J=0.1, NRG 3.7) with T_K from the local susceptibility,
and 3.94 at every rho J from 0.09 to 0.2 with Wilson's T_K
(`nrgthermo.kondo_temperature`, the slope of T chi_imp, which DMRG at T=0
cannot give since S^z_tot commutes with H). The paper reads T_K from the
crossing T chi_imp = 0.07 at Lambda=2 (that crossing is 1.03 T_K(W) at
Lambda=2 and 1.11 T_K(W) at Lambda=4), and read that way we find 3.83. So
the 4.4 against 3.7 is the definition of T_K, not the mesh. Always say which
T_K a ratio uses. The absolute Delta_c at Lambda=4 depends on the mesh (29%
above the paper on mean, 0.7% below on zp) while the ratio does not, so
compare ratios at Lambda=4, not absolute values. In `nrgthermo.py` the free
chain is grown site by site with relative precision and the (Q, 2Sz) and
(Q, -2Sz) spectra are symmetrized at B=0; without either, chains past about
n=46 at Lambda=4 gave a wrong, even negative, T_K. At Lambda=2 use
beta_bar=1.5 and keep >= 800.

On v3 a sector solve keeps the impurity charges of its start, since no term
of H moves charge on or off an impurity site: S > 1/2 (`ysr.py`,
`_seed_impurity_charge`) and several surrogate spins
(`surrogate._warm_sector_solve`) first solve with a weak hopping onto the
impurity, and both raise when a solved state is not singly occupied there.
Do not remove the warm stage. The z-averaged continuum
(`ysrtk/zaverage.py`) broadens in xi = sqrt(w^2 - Delta^2), not in
w - Delta (2e4 times worse), and is resolved near the gap edge by `nkry`,
not by `maxm`. The surrogate bath has no Kondo logarithm: it is
quantitative only for Delta >= 0.1 D at L=4.

`ysrtk/bcsperturbative.py` (second order with Dynes densities) must give
back `conductance.second_order_dIdV` at Delta=0 to roundoff; its substrate
relaxation rate carries the BCS coherence factor (the F F term), without
which it was 2.9x too large just above 2 Delta, and a reused `kernel=` raises outside its
range. `ysrtk/ysrcurrent.py` (Ruby et al. rate model) normalizes the
composite-vertex weights by their total per tip spin, and `ysr.bcs_dos` at
gamma=0 is exact now (zero inside the gap), while the rate model keeps its
own `tip_dos`.

Near the singlet-doublet transition a plain DMRG solve stops at random at a
mixture of the two parities (elastic weight up to 1.2e-2 of the total on 8
sites, poles split and moved, u2 halved); `parity=True` in
`zaverage.stm_poles`, `zaverage_stm_spectrum` and
`surrogate.get_stm_spectrum`, or `paritystate.parity_ground_state(fc)`,
solves inside the parity sector and promotes the state to dense indices on
its own clone (Index identity is per session, so never move it to another
chain). `parity="auto"` solves plainly and takes the sector route only when
the elastic weight is above 1e-6 of the sum rule (`paritystate.MIXTURE_TOL`).
Since 2026-11-07 "auto" is the default, chosen by the user; it never
trips at Delta=0, where there is no parity sector. On "python" a
sector solve stops its local Davidson solves on the scale of H without the
charge penalty (`dmrg(penalty=)`, passed by `Chain.gs_energy` in sector
mode only); with the penalty in the scale it stalled at 2 matvecs per bond
and needed twice the sweeps. `dmrg_excited`, `dmrg_generalized` and
`cbe.dmrg_cbe` take `penalty=` the same way, and the default penalty
strength follows each new Hamiltonian. A promoted sector state of a degenerate pair
has a smaller one-site tangent space than a plain mixture, and TASK at J=0
keeps 1/3 of the ED weight at U=0 and 0.405 at U=0.3; on short chains
`pad_bonds="full"` (TASK, `pyitensor/task.py`) makes both exact. Padding is
opt-in (None default). An int pad is not a long-chain cure (it slowed the
Krylov and moved poles on 21 sites) and it inflates the pole variances, so
they stop being an error bar; `task.MAX_PADDED_DIMENSION` = 4e6 caps the
padded space. `expand_bonds=k` (with `expand_rounds=2` at
J=0) adds at most k directions per bond from the part of B|GS> and H P B|GS>
the tangent space drops: exact on the short-chain cases at a fraction of the
"full" space, and the pole variances stay an error bar. It is not a
long-chain cure either: on 21 sites it made nkry=100 worse and nkry=300
1.5x slower, because there the missing piece is Krylov convergence, not the
space. Opt-in, None default.

`ysrtk/subgapscan.py` (warm-started scans over Delta, Dz and the Delta_c
bisection) needs `set_initial_wf_guess` on every backend, since "python"
drops its state when the terms change. The warm-sweep default
`warm_sweeps` = min(nsweeps, max(8, nsweeps//3)) has its floor 8 from
"python" (6 sweeps missed ED by 7e-9 on 4 sites); v3 does mix the chain's
noise into a warm start that is already at maxm. `nsweeps_warm="auto"`
(blocks of 2 sweeps until the energy moves by less than `warm_tol` |Delta|)
is opt-in and was not measured on the 26-site chain, so it is not the
default. A warm start does not see
a large jump in Delta, so only `cold_every` catches it.
