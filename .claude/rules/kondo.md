---
paths:
  - "src/dmrgpy/kondospectrumtk/**"
  - "tests/test_kondo_spectrum*.py"
  - "tests/test_kondo_screening.py"
  - "src/dmrgpy/wilsonchain.py"
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
