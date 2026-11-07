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
(T=0 only) reproduces its two T=0-reducible pieces via the dynamical
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
<psi_l(-s-tau)|Sk|psi_j(-s)> from checkpoints of the same trajectories (1.5x
the TDVP cost), and `kondo_term_from_two_time` takes (t, G, Gx) triples with
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
`itensor_version="python"` and 3 (v3 exports its site tensors to the same
Python code, `taskdmrg.session_spectrum`), the default on both; the near-zero Ritz energies are snapped
to 0 or the zero-bias point is rounding-dependent on a degenerate ground
state. It is in `_N_GS_SUBMODES`, checked on a Kramers doublet. A Ritz
energy below the snapping window means the state handed over is not the
ground state (DMRG stopped at a mixture of a near-degenerate doublet), and
`TaskKondoSpectrum` warns (2026-10 audit, finding 3).

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
symmetrize it. T=0 only, fixed Gamma0, no level renormalization; NaN plus a
warning when the sample coupling or the per-bias tip coupling
(`info["coupling"]`) reaches order one, a warning above 0.5.
`tests/test_kondo_leading_log.py`. Public as `get_kondo_spectrum(order="LL")`,
`mode="ED"` and `T=0` only; `mode="DMRG"` raises because the TASK space does
not close under the products of spin operators the flow needs.

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
