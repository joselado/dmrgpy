---
paths:
  - "src/dmrgpy/kondospectrumtk/**"
  - "tests/test_kondo_spectrum*.py"
  - "tests/test_atom_iets.py"
  - "examples/kondo/**"
---

# STM/Kondo tunneling spectra (`kondospectrumtk/`)

These are the notes CLAUDE.md used to carry as its "STM/Kondo tunneling
spectra" section; they load when a file under `kondospectrumtk/`, a Kondo
test or a Kondo example is read or edited.

`Spin_Chain.get_kondo_spectrum` is the third-order perturbation theory of
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
state. It is in `_N_GS_SUBMODES`, checked on a Kramers doublet.
