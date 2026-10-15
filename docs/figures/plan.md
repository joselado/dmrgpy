# Figures for the user guide: the plan

This file plans one figure per capability of the user guide
(`docs/guide/*.md`, in the order of the `nav` of `mkdocs.yml`). Each figure
shows physics a condensed matter physicist recognizes at a glance, and each
one carries a reference curve (an exact result, a free-fermion calculation,
or ED on a small chain), so that every plot is also a check of the code.
Nothing here has been run at full size yet. The parameters come from the
guide pages, from the scripts under `examples/` that compute something
close, and from a few probes of seconds each on the development machine
(an Intel Core Ultra 5 225U, one thread), scaled to the sizes below. The
runtimes are therefore estimates, and the figure scripts print their own.

## Conventions

### Folder layout

- One folder per figure, `docs/figures/<id>/`, with two scripts:
  `compute.py` runs the calculation and writes `data.npz`; `plot.py` reads
  only `data.npz` and writes the PNG files into `docs/guide/img/`.
  `plot.py` never imports `dmrgpy`, so a figure can be replotted on any
  machine without the compiled backend.
- `compute.py` prints every parameter it uses before the calculation
  starts, and stores them in `data.npz` under the key `params`, as one JSON
  string (`np.savez(..., params=json.dumps(params))`). The dictionary holds
  the model, the sizes, `maxm`, `nsweeps`, `delta`, the time step, the
  backend, every seed, the git commit of `dmrgpy` (`git rev-parse HEAD`)
  and the wall time of the run. Every array in `data.npz` is named in a
  comment block at the top of `compute.py`.
- A reference curve that needs no `dmrgpy` (an analytic formula, a
  free-fermion diagonalization) is computed in `compute.py` too and stored
  in `data.npz`, so that `plot.py` only draws.
- `docs/figures/style.py` is the shared style (below). `plot.py` loads it
  with `sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))`.
- `docs/figures/run.sh <id>` runs `compute.py` and then `plot.py` with the
  environment pinned:

```bash
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH=<repo>/src
cd docs/figures/<id> && taskset -c <P-core> python3 compute.py && python3 plot.py
```

  The i5-12600K has 6 P-cores (threads 0 to 11) and 4 E-cores (threads 12
  to 15), and DMRG runs about 2 times slower on an E-core, so each of the
  three concurrent jobs is pinned to one P-core with `taskset`. At most
  three `compute.py` run at the same time. `PYTHONPATH` must point to
  this checkout, since `site-packages/dmrgpy` is a symlink to one specific
  checkout (see `CLAUDE.md`, "Worktree/symlink pitfall").
- `data.npz` files are committed next to their scripts (they are small,
  below 1 MB each), so a figure can be replotted without recomputing.

### Shared matplotlib style (`docs/figures/style.py`)

- Two PNG files per figure: `<id>.png` for the light theme and
  `<id>-dark.png` for the dark theme, both from the same `data.npz`. The
  page embeds them as a pair, which MkDocs Material shows by theme:

```markdown
![caption](img/<id>.png#only-light)
![caption](img/<id>-dark.png#only-dark)
```

  `docs/build_guide_pdf.py` passes the pages to pandoc, which would look
  for a file with the `#only-light` suffix. When the first figure lands,
  `page_source()` gets one more regex: drop the `#only-dark` line and strip
  `#only-light` from the other one, so the PDF keeps the light figure.
- Transparent background (`savefig(..., transparent=True)`), so the page
  color shows through in both themes. The ink (axes, ticks, labels, legend
  text, reference curves) is `#0b0b0b` on light and `#f0efec` on dark;
  secondary ink (tick labels, minor text) `#52514e` and `#c3c2b7`. No grid.
  Top and right spines off. `legend.frameon=False`, and the colorbar
  outline and its tick labels in the ink color, since the matplotlib
  defaults (a white legend box, a black colorbar frame) show up as
  rectangles on the dark theme.
- Series colors, in a fixed order, from the validated palette of the
  `dataviz` skill (light theme, then dark theme): blue `#2a78d6`/`#3987e5`,
  orange `#eb6834`/`#d95926`, aqua `#1baf7a`/`#199e70`, yellow
  `#eda100`/`#c98500`, magenta `#e87ba4`/`#d55181`, green `#008300`,
  violet `#6250d6`/`#9085e9`, red `#e34948`/`#e66767`. Color follows the
  entity (the same chain length or the same backend keeps its color across
  panels). A reference curve (exact, free fermions, ED) is always drawn in
  the ink color, dashed, 1.2 pt, and labeled "exact" or "ED"; the computed
  data are colored lines of 1.6 pt or markers of 5 pt.
- Heat maps use `magma` on both themes, with a colorbar and its label.
  Lines drawn on a heat map (a dispersion, a continuum edge) are white
  dashed, 1.2 pt.
- Sizes: one panel `figsize=(6.4, 4.0)`, two panels `(9.6, 3.8)`, three
  panels `(11.0, 3.6)`, at `dpi=160`, `bbox_inches="tight"`. The content
  column of MkDocs Material is about 800 px wide, so these read at full
  width without zoom. Font size 10 pt, axis labels 11 pt, legend 9 pt,
  font family DejaVu Sans (the matplotlib default, no extra install).
- Axis labels carry the unit in brackets, for example `omega [J]`,
  `T [J]`, `eV [meV]`. Momentum axes are labeled in units of pi. No panel
  title; the caption under the figure in the page says what it shows. A
  legend whenever a panel has two or more curves.
- `plot.py` takes an optional `--theme light|dark` and by default writes
  both files.

### Sizes and references

- When DMRG is the point, the chain is longer than ED can reach (24 sites
  or more for spin 1/2, 16 or more for spin 1 and for spinful fermions),
  and the reference is exact: Bethe ansatz, free fermions, Calabrese-Cardy,
  or a closed form.
- When the figure checks a method against ED, the chain is 8 to 12 sites
  and the ED curve is the reference.
- The thermal ED reference uses the full spectrum, `eigvalsh` of the ED
  Hamiltonian as in `examples/finite_temperature/thermal_purification_VS_exact`,
  not `vev(mode="ED", T=...)`, which truncates to the lowest states above
  about 2000 states.
- Every KPM figure checks its sum rule and every SECTOR figure its
  `info["captured"]` before the full run, with a short calibration loop
  written in `compute.py` (the loop runs on the smallest size of the
  figure and prints the number it converged).

## Figures, page by page

### Home (`index.md`): no new figure

The page is an introduction and a table of contents. It can show the
$S(q,\omega)$ figure of the dynamical correlators page as a first image,
without a new calculation.

### Models and Hilbert spaces (`models.md`): no figure

A table of chain classes and operators. Nothing to plot.

### Hamiltonians and operators (`operators.md`): no figure

An API listing (how to write $H$, the state operations). Nothing to plot.

### Ground state (`ground-state.md`)

#### 1. `ground-state_convergence`

1. File `ground-state_convergence.png`, page `ground-state.md`, after the
   keyword list of "Ground-state energy", with a forward reference from
   "Energy fluctuation".
2. The energy error $E_0(\chi)-E_0^{\rm ED}$ and the energy fluctuation
   $\delta E(\chi)$ of the Heisenberg chain as functions of the bond
   dimension $\chi$. Both fall exponentially with $\chi$, and the energy
   error falls roughly as $\delta E^2$, which is the statement of the
   page. Reference: ED on the same chain.
3. Spin-1/2 Heisenberg chain, $n=16$ (ED has 65536 states, seconds),
   `maxm` in `[2, 4, 8, 12, 16, 24, 32, 48, 64, 96]`, `nsweeps=20`,
   `cutoff=1e-14` (so that `maxm` is the only truncation that acts),
   `itensor_version=3`.
   For each `maxm`: `gs_energy()` and `gs_energy_fluctuation()`. ED:
   `gs_energy(mode="ED")` once.
4. One panel, log scale on $y$: the energy error per site in blue, $\delta E$
   per site in orange, and $\delta E^2$ (times a constant fitted at the
   largest two points) as a dashed ink line, against $\chi$ on a log $x$
   axis. The reader sees two straight lines, and that the dashed $\delta E^2$
   line runs parallel to the energy error.
5. About 1 minute (ten solves of a 16-site chain, each under 2 s at
   `maxm=96`, plus ten variance calculations).

#### 2. `ground-state_correlators`

1. File `ground-state_correlators.png`, page `ground-state.md`, in
   "Static correlators".
2. The decay of $\langle S^z_iS^z_{i+r}\rangle$ with distance in a spin-1/2
   and in a spin-1 Heisenberg chain: a power law (a gapless Luttinger
   liquid) against an exponential (the Haldane gap, correlation length
   $\xi\approx 6.0$). This is the README figure of dmrgpy, now with the
   references on it: the Luttinger form $(-1)^r\sqrt{\ln r}/r$ for spin
   1/2, and $\xi=6.03$ for spin 1.
3. Open chains of $n=64$ sites, pairs measured from the site
   $i=n/2-16$ outwards ($r=1$ to 32, so that both sites stay in the bulk).
   Spin 1/2: `maxm=120`, `nsweeps=20`. Spin 1: `maxm=100`, `nsweeps=20`.
   `itensor_version=3`. Correlators with `vev(Sz[i]*Sz[i+r])`.
4. One panel, $|\langle S^z_iS^z_{i+r}\rangle|$ against $r$ on a log-log
   scale, spin 1/2 in blue, spin 1 in orange, with a dashed ink line
   $a\sqrt{\ln r}/r$ fitted at $r=8$ to 32 for spin 1/2 and a dashed
   $b\,e^{-r/6.03}$ for spin 1 (both prefactors fitted, the exponents and
   $\xi$ fixed). The reader sees the straight power law and the curve
   bending down for spin 1, and reads $\xi\approx6$ from the slope.
5. About 3 minutes (two 64-site ground states, 64 correlators).

#### 3. `ground-state_sectors`

1. File `ground-state_sectors.png`, page `ground-state.md`, in "Fixed
   particle number or fixed magnetization".
2. The magnetization curve $M(h)$ of a spin chain from the sector energies
   $E_0(S^z)$: the field at which the magnetization steps from $M$ to
   $M+1$ is $h_M=E_0(M+1)-E_0(M)$, so the sector energies give the whole
   staircase. For the XX chain the curve is free fermions,
   $m(h)=\frac{1}{\pi}\arcsin(h/J)$ for $|h|\le J$, which is overlaid.
   For the Heisenberg chain the saturation field $h_s=2J$ is exact and is
   marked.
3. $n=32$ spin-1/2 chains, `itensor_version=3`, `maxm=60`, `nsweeps=12`,
   `set_conserved_sector(Sz=2M)` for $M=0$ to 16 (17 solves per model).
   Models: XX ($S^xS^x+S^yS^y$) and Heisenberg, $J=1$.
4. One panel, $M/n$ against $h$ as two staircases (blue XX, orange
   Heisenberg), the exact XX curve dashed in ink, and a vertical dotted
   line at $h_s=2$. The reader sees the XX staircase sitting on the
   arcsine curve within $1/n$, and the Heisenberg curve reaching $m=1/2$
   exactly at $h=2$.
5. About 2 minutes (34 sector solves of a 32-site chain, about 2 s each).

Expectation values (`get_magnetization`, `get_density`) and `promote_to_dense`
get no figure of their own: the profiles appear in figures 5 and 13.

### Excited states and gaps (`excited-states.md`)

#### 4. `excited-states_haldane`

1. File `excited-states_haldane.png`, page `excited-states.md`, after the
   spin-1 example of "Low-lying spectrum". The degeneracy page links to it.
2. The Haldane chain: the two spin-1/2 edge states of an open spin-1
   chain, and the bulk gap $\Delta=0.4105\,J$. On an open chain the four
   lowest states are the singlet and the triplet of the two edge spins,
   split by an amount that decays as $e^{-n/\xi}$, and the bulk gap is the
   energy of one more magnon. Both come from ground states of sectors:
   $E_0(S^z=1)-E_0(S^z=0)$ is the edge splitting and
   $E_0(S^z=2)-E_0(S^z=1)$ is the gap. A probe at $n=16$ and $24$ gave
   splittings 0.0496 and 0.0134 ($\xi\approx6$) and gaps 0.599 and 0.512,
   which converge as $1/n^2$ (the magnon has momentum $\pi-\pi/n$ on an
   open chain) to the known 0.4105. Reference: the thermodynamic-limit gap
   0.41048 (White and Huse) and the extrapolation in $1/n^2$,
   `extrapolate.size_extrapolation(ns**2, gaps, n=1)` (a fit of
   $a+b/n^2$; a generic quadratic in $1/n$ would spend a parameter on a
   $1/n$ term that is zero here).
3. Open spin-1 Heisenberg chains, $n$ in `[16, 24, 32, 48, 64, 96]`,
   `itensor_version=3`, `maxm=100`, `nsweeps=14`, three sector solves per
   length (`set_conserved_sector(Sz=0)`, `Sz=2`, `Sz=4`; `Sz` is in units
   of 1/2). The edge profile $\langle S^z_i\rangle$ in the `Sz=2` sector at
   $n=48$ (`vev(Sz[i])`). This replaces `get_excited(n=6)` at these
   lengths, which is slower and can land on a wrong state; the guide's
   snippet stays, the figure says the gap was obtained from the sectors.
4. Two panels. Left: the gap against $1/n^2$ (blue points) with the
   $a+b/n^2$ extrapolation (dashed) and a horizontal dotted line at
   0.4105; in the same panel, the edge splitting against $n$ on a second
   small inset, log scale, with the $e^{-n/6.03}$ line. Right: the profile
   $\langle S^z_i\rangle$ at $n=48$ in the $S^z=1$ edge state: a spin-1/2
   worth of magnetization at each end, oscillating and decaying over
   about 6 sites, zero in the bulk. The reader sees the extrapolation hit
   0.41, the four-fold edge manifold collapse, and two edge spins 1/2.
5. About 8 minutes (18 sector solves; the three at $n=96$ take about 1
   minute each).

#### 5. `excited-states_nonhermitian`

1. File `excited-states_nonhermitian.png`, page `excited-states.md`, in
   "Non-Hermitian Hamiltonians".
2. The many-body skin effect of the interacting Hatano-Nelson chain: with
   asymmetric hopping $t\pm\gamma$ the right eigenstate with the smallest
   real part piles its particles against one edge. The density profile
   $\langle n_i\rangle$ in the right ground state for several $\gamma$,
   with ED on a short chain as the check of NH-DMRG, and the residual
   $\lVert H\psi_R-E\psi_R\rVert$ printed. For open boundaries and real
   $\gamma$ the spectrum is real, so the ground state is well defined.
3. Spinless `Fermionic_Chain`, hopping $(t+\gamma)c^\dagger_ic_{i+1}+(t-\gamma)c^\dagger_{i+1}c_i$
   with $t=1$, interaction $V\sum_i(n_i-\tfrac12)(n_{i+1}-\tfrac12)$ with
   $V=1$ (so that the ground state is at half filling), $\gamma$ in
   `[0, 0.2, 0.4]`. DMRG: $n=20$, `e, psil, psir = fc.nhdmrg(tol=1e-5)`,
   `maxm=80`, `itensor_version=3`. ED check: $n=8$,
   `get_excited_states(mode="ED",n=1)`. The density is the right-vector
   one, $\langle\psi_R|n_i|\psi_R\rangle/\langle\psi_R|\psi_R\rangle$,
   computed explicitly from `psir` with `fc.aMb(psir, N[i], psir)` and
   from the ED right eigenvector in the same way. The biorthogonal
   density $\langle\psi_L|n_i|\psi_R\rangle$, which is what `vev` gives
   after an NH solve, is flat: the Hatano-Nelson chain is a gauge
   transform of a Hermitian chain and $n_i$ commutes with the transform.
4. One panel: $\langle n_i\rangle$ against $i$ for the three $\gamma$ at
   $n=20$ (blue, orange, aqua lines), and the $n=8$ ED profiles as dashed
   ink lines on the same axis (plotted against $i$, so they occupy the
   first 8 sites). The reader sees a flat profile at $\gamma=0$ and the
   density moving to the right edge as $\gamma$ grows.
5. About 5 minutes (three NH-DMRG solves at 20 sites, with up to 5
   restarts each). Calibration: if the residual at $\gamma=0.4$ stays above
   `tol`, reduce $\gamma$ to 0.3; the skin state is a product-like state
   and is not hard for an MPS.

#### 6. `excited-states_heff`

1. File `excited-states_heff.png`, page `excited-states.md`, in
   "Effective low-energy Hamiltonians".
2. The exchange of the Hubbard dimer: `get_heff` fits $J$ in the four
   lowest states, with and without the singly occupied projector. With it
   the fit is the des Cloizeaux effective Hamiltonian and $J$ is the exact
   singlet-triplet splitting $E_T-E_S=\tfrac12\big(\sqrt{U^2+16t^2}-U\big)$,
   which tends to $4t^2/U$ at strong coupling. Without it the fit is
   $(E_T-E_S)/(1/4+3w/4)$, too large at small $U$.
3. `Spinful_Fermionic_Chain(2)`, $t=1$, $U$ in `np.geomspace(1, 20, 25)`,
   $20(N-2)^2$ added to pin half filling, `mode="ED"`,
   `get_heff(operators=[SxSx, SySy, SzSz], n=4, mode="ED", projector=Q)`
   with `Q = single_occupancy_projector(fc)`, and the same call without
   `projector`; the fitted $J$ is the common coefficient.
4. One panel, $J$ against $U$ on log axes: the fit with the projector
   (filled points), the fit without it (open squares), $E_T-E_S$ (dashed
   ink) and $4t^2/U$ (dash-dotted ink). The reader sees the filled points
   on the exact curve at every $U$, and the open squares above it at
   small $U$.
5. Under 1 minute.

"Gap to another particle number" and "Generalized eigenvalue problem"
get no figure: the charge gap is the same calculation as figure 3 (sector
energies), and the generalized problem has no standard physics picture
that a reader recognizes.

### Ground-state degeneracy (`excited-states.md`): no new figure

The section is one method, and its physics (the four-fold edge manifold of
the Haldane chain, with its finite-size splitting) is figure 4, which the
page links to.

### Entanglement (`entanglement.md`)

#### 7. `entanglement_central_charge`

1. File `entanglement_central_charge.png`, page `entanglement.md`, in
   "Central charge", with a back reference from "Entanglement entropy of a
   cut".
2. The Calabrese-Cardy formula on an open critical chain,
   $S(\ell)=\frac c6\ln\big[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}\big]+c'$,
   for three chains: XX ($c=1$), Heisenberg ($c=1$) and the critical
   transverse-field Ising chain ($c=1/2$), against the gapped spin-1
   chain, whose entropy saturates. For the XX chain the exact $S(\ell)$ at
   the same $L$ comes from free fermions (the eigenvalues $\nu_k$ of the
   correlation matrix $C_{ij}=\langle c^\dagger_ic_j\rangle$ restricted to
   the block, $S=-\sum_k[\nu_k\ln\nu_k+(1-\nu_k)\ln(1-\nu_k)]$), and is
   overlaid point by point.
3. $L=100$, `itensor_version=3`, `nsweeps=20`. XX and Heisenberg:
   `maxm=200`. Ising $H=\sum S^z_iS^z_{i+1}+\tfrac12\sum S^x_i$ (the
   cookbook's critical point): `maxm=120`. Spin-1 Heisenberg: `maxm=100`.
   Entropies with `wf.get_bond_entropy(i, i+1)` for every bond, and
   `wf.get_CFT_central_charge()` printed for each chain, plus a fit on even
   bonds only (the alternating term, see the page).
4. Two panels. Left: $S(\ell)$ against $\ell$ for the four chains (XX blue,
   Heisenberg orange, Ising aqua, spin-1 yellow), the exact XX entropies
   as ink dashes. Right: $S$ against $x=\ln[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}]$,
   even bonds only, with the lines of slope $1/6$ and $1/12$ dashed in
   ink; the legend gives the fitted $c$ of each chain. The reader sees
   two parallel straight lines for the $c=1$ chains, a line of half the
   slope for Ising, a flat curve for the Haldane chain, and the XX points
   on the exact ones.
5. About 12 minutes (three critical 100-site ground states at `maxm` 120
   to 200, 2 to 5 minutes each, plus the gapped one; 400 bond entropies
   are seconds).

#### 8. `entanglement_natural_orbitals`

1. File `entanglement_natural_orbitals.png`, page `entanglement.md`, in
   "Single-particle correlation matrix".
2. The natural-orbital occupations $n_\alpha$ of the Hubbard chain as $U$
   grows: 0 and 1 at $U=0$ (one Slater determinant), moving towards $1/2$
   in the Mott insulator, and the correlation entropy $S$ growing with
   $U$. Reference: at $U=0$ the occupations are exactly 0 and 1, so
   $S=0$; at $n=4$ orbitals the whole curve is checked with `mode="ED"`.
3. `Spinful_Fermionic_Chain(12)` (24 spinless sites), hopping $t=1$,
   interaction written as `U*sum(Nup[i]*Ndn[i])` (not `set_hubbard`,
   whose diagonal entry adds a $U/2$ chemical potential, harmless at
   fixed $N_f$ but better left out), $U$ in `[0, 1, 2, 4, 8, 16]`, half
   filling through
   `set_conserved_sector(Nf=12)`, `itensor_version=3`, `maxm=100`,
   `nsweeps=14`; `get_correlation_eigenvalues()` and
   `get_correlation_entropy()`. The ED check at `Spinful_Fermionic_Chain(4)`
   with the same $U$ values.
4. Two panels. Left: the sorted occupations $n_\alpha$ against the orbital
   index, one line per $U$ (palette order). Right: $S$ against $U$ (blue
   points) with the $n=4$ ED curve (dashed ink) and the DMRG $n=4$ values
   (orange markers) to show the two agree. The reader sees the step
   function at $U=0$ melting into a smooth curve around $1/2$.
5. About 4 minutes (six 24-site sector solves, each under 30 s, and the
   correlation matrices).

Site entropy, pair entropy, mutual information and the four-point tensor
get no figure: they are the same measurement on other reduced density
matrices, and the page's text says when to use them.

### Correlators and distributions (`correlators.md`)

#### 9. `correlators_distribution`

1. File `correlators_distribution.png`, page `correlators.md`,
   in "Operator distributions".
2. The full counting statistics of the magnetization in the
   transverse-field Ising chain: $P(m)$ of $m=\frac1n\sum_iS^z_i$ is one
   sharp peak at $m_0\approx0.45$ in the ordered phase (a symmetry-broken
   state, selected by a small pinning field), broad at the critical
   field, and one peak at $m=0$ in the paramagnet. The pinning field is
   needed because at $h=0.2$ the two ordered states are split by about
   $(2h)^n$, so without it DMRG returns an arbitrary member of the
   doublet and ED the symmetric cat state, which would disagree for a
   reason that is not an error. Reference: `mode="ED"` on a 12-site chain
   at the same fields and the same pinning, and the normalization
   $\int P\,dm=1$ printed.
3. $H=-\sum S^z_iS^z_{i+1}-h\sum S^x_i-h_z\sum S^z_i$ with $h_z=10^{-3}$,
   critical at $h=1/2$ in this normalization, $h$ in `[0.2, 0.5, 1.0]`.
   DMRG: $n=24$,
   `itensor_version=3`, `maxm=60`, `kpmmaxm=60`,
   `get_distribution(X=M, delta=0.02)`. ED: $n=12$, the same call with
   `mode="ED"`. (A probe at $n=8$, $h=0.2$ put the DMRG peak at $m=0.505$
   against 0.500 from ED, with the integral 0.9998.)
4. One panel, $P(m)$ against $m$ for the three fields at $n=24$ (blue,
   orange, aqua) with the $n=12$ ED curves dashed in ink. The reader sees
   a sharp peak near $m=0.45$ at $h=0.2$ broaden through the critical
   point and settle at $m=0$ in the paramagnet.
5. About 2 minutes.

The reduced density matrix `get_rdm` gets no figure (a $2\times2$ matrix).

### Topological invariants (`topology.md`)

#### 10. `topology_berry_phase`

1. File `topology_berry_phase.png`, page `topology.md`, one figure for
   the three sections.
2. Three invariants that jump. (a) The many-body Berry phase of the
   dimerized Heisenberg ring, $J_i=1+\delta(-1)^i$, with a twist
   $e^{i\phi}$ on one bond: $\gamma=\pi$ for one dimerization and 0 for
   the other, jumping at $\delta=0$. (b) The single-particle Zak phase of
   the SSH chain from `topology.berry_phase_matrix`, which jumps at the
   same $\delta$; this is the exact free-fermion reference of the same
   transition. (c) The fermion parity $P$ and $\langle N\rangle$ of a
   short spinful chain with onsite pairing and a magnetic impurity
   against the chemical potential: $P$ flips where a Yu-Shiba-Rusinov
   level crosses zero, and $\langle N\rangle$ jumps by one there.
3. (a) `Spin_Chain` of $n=12$ spins 1/2, ring closed by the twisted bond
   as in `examples/topological/zak_phase`, 16 twists
   `np.linspace(0, 2*np.pi, 16, endpoint=False)`, $\delta$ in
   `np.linspace(-0.6, 0.6, 12)` (an even number of points, so that no
   point sits at $\delta=0$, where the uniform ring has level crossings
   as $\phi$ winds and the phase is undefined), `itensor_version=3`, `maxm=60`,
   `nsweeps=12`. (b) `berry_phase_matrix(hk, nk=60)` with the two-band
   SSH Bloch matrix
   $h(k)=[(1+\delta)+(1-\delta)\cos k]\sigma_x+(1-\delta)\sin k\,\sigma_y$
   written in `compute.py` (`k` from 0 to 1 in units of $2\pi$, as the
   page says), same $\delta$ grid. (c) `Spinful_Fermionic_Chain(6)`
   with the Hamiltonian of `examples/topological/parity` ($U=0.5$,
   $B_z=1$, pairing 0.5), $\mu$ in `np.linspace(-1, 0, 41)`,
   `itensor_version=3`, `wf.get_fermionic_parity()` and `vev(sum(N))`.
4. Three panels. Left: $\gamma/\pi$ against $\delta$ for the spin ring
   (blue points) and the SSH Zak phase (dashed ink) on the same axis.
   Middle: $\langle N\rangle$ against $\mu$. Right: $P$ against $\mu$,
   with the crossing marked by a vertical dotted line. The reader sees
   two step functions at $\delta=0$ on top of each other, and a parity
   flip at the same $\mu$ where $\langle N\rangle$ jumps.
5. About 3 minutes (208 ground states of a 12-site ring, about 0.5 s each,
   and 41 of the 6-site chain).

### Fidelity susceptibility (`fidelity.md`)

#### 11. `fidelity_ising`

1. File `fidelity_ising.png`, page `fidelity.md`, after the keyword list.
2. The fidelity susceptibility per site of the transverse-field Ising
   chain across its critical point, $\lambda=1$ in the page's
   normalization ($H=\sum\sigma^z_i\sigma^z_{i+1}+\lambda\sum\sigma^x_i$):
   a peak at $\lambda=1$ that grows with the length. Reference: the same
   quantity from `mode="ED"` at $n=12$, and the position of the peak.
3. Open chains, $n$ in `[12, 16, 24, 32]`, $\lambda$ in
   `np.linspace(0.6, 1.4, 41)`, `fidelity.get_fidelity(sc, h0, h1, lam, dl=0.02)`
   (the default `dl=0.1` is coarse next to a peak of width $\sim1/n$),
   `fmode="derivative"`, `itensor_version=3`, `maxm=60`, `nsweeps=12`. ED
   at $n=12$ with the same `dl`.
4. One panel, $\chi/n$ against $\lambda$ for the four lengths (palette
   order) and the $n=12$ ED curve as dashed ink on top of the blue one.
   The reader sees the peak sharpen and grow at $\lambda=1$, and the ED
   curve covering the DMRG one.
5. About 10 minutes (41 points times 4 lengths, each point about 6 solves
   on a copy of the chain, 0.2 to 1.5 s per solve).

### Mean-field decoupling (`mean-field.md`)

#### 12. `mean-field_neel`

1. File `mean-field_neel.png`, page `mean-field.md`, after the keyword
   list.
2. How mean-field theory overestimates order: at $p=0$ the Heisenberg
   chain has a Neel state with $|\langle S^z_i\rangle|=1/2$, and as the
   fraction $p$ of the exchange treated exactly grows, the staggered
   magnetization melts to zero at $p=1$ (the exact chain), while the
   energy per site goes from the classical $-1/4$ to the Bethe-ansatz
   $1/4-\ln2=-0.4431$. Both limits are exact references.
3. Spin-1/2 Heisenberg chain of $n=24$ sites, `itensor_version=3`,
   `maxm=60`, `meanfield.spinchain_meanfield(sc, p=p, m0=staggered, mix=0.5)`
   for $p$ in `[0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]`,
   `mf.get_magnetization()` and `mf.gs_energy()`.
4. Two panels. Left: $\langle S^z_i\rangle$ against $i$ for $p=0$, 0.6,
   0.9 and 1 (palette order). Right: the staggered magnetization of the
   central 12 sites (blue, left of the panel) and the energy per site
   (orange points) against $p$, with dashed ink lines at $-1/4$ and
   $-0.4431$. The reader sees the Neel pattern fade and the energy drop
   to the Bethe-ansatz value.
5. About 4 minutes (eight self-consistent loops of up to 100 iterations,
   each iteration a 24-site DMRG solve of 0.3 s).

### Dynamical correlators (`dynamical-correlators.md`)

#### 13. `dynamical-correlators_sqw`

1. File `dynamical-correlators_sqw.png`, page `dynamical-correlators.md`,
   in "Dynamical spin structure factor". Also the first image of the
   Home page.
2. The two-spinon continuum of the spin-1/2 Heisenberg chain:
   $S^{zz}(q,\omega)$ with its weight concentrated just above the des
   Cloizeaux-Pearson lower edge $\frac\pi2J|\sin q|$, bounded above by
   $\pi J|\sin(q/2)|$, and piled up at $q=\pi$. Both edges are exact and
   are drawn on the map. The sum rule $\int S(q,\omega)d\omega=\langle S^z_{-q}S^z_q\rangle$
   is checked at every $q$ against `vev(Sq*Sq)`.
3. Open chain of $n=24$ spins 1/2, `itensor_version=3`, `maxm=80`,
   `nsweeps=14`, `kpmmaxm=80`. The route is one KPM run per momentum with
   the operator of the open chain's standing wave,
   $S^z_q=\sqrt{\tfrac{2}{n+1}}\sum_i\sin(q(i+1))S^z_i$,
   `get_dynamical_correlator(name=(Sq, Sq), es=es, delta=0.1)` for
   $q_k=\pi k/(n+1)$, $k=1$ to 24. The moments are computed once per
   $q$ with `get_dynamical_correlator_moments(name=(Sq, Sq), delta=0.1)`,
   and `kpmdmrg.dynamical_correlator_from_moments` builds from them both
   the plot grid `es=np.linspace(0, 4, 400)` and the sum-rule grid
   `es=np.linspace(-1, 20, 2000)`, so there is one moment run per $q$,
   not two. A probe at $n=16$,
   `delta=0.3`, `kpmmaxm=40` gave 92 moments in 0.8 s, with the sum rule
   at 0.858 against 0.914, so the calibration loop runs $q=\pi$ at
   `kpmmaxm` 40, 80 and 120 and stops when the sum rule holds to 1 per
   cent; the number of moments at $n=24$, `delta=0.1` is about 400.
4. Two panels. Left: the map $S^{zz}(q,\omega)$ (`magma`), $q/\pi$ on $x$,
   $\omega$ in units of $J$ on $y$, with the two edges dashed in white.
   Right: cuts at $q=\pi/2$ and $q=\pi$ (blue, orange) with the lower
   edge of each marked by a vertical dotted line, and in the legend the
   sum rule of each cut as a number. The reader sees the continuum with
   a sharp lower edge following $\frac\pi2|\sin q|$ and no weight below
   it.
5. About 20 minutes (24 runs of about 400 moments at `kpmmaxm=80`, 40 to
   50 s each, scaled from the probe). The heaviest figure of the plan; if
   the calibration asks for `kpmmaxm=120` the estimate doubles.

#### 14. `dynamical-correlators_spectral_function`

1. File `dynamical-correlators_spectral_function.png`, page
   `dynamical-correlators.md`, in "Single-particle spectral function".
2. The Mott gap of the half-filled Hubbard chain opening with $U$: the
   local spectral function $A_{ii}(\omega)$ of a bulk site, one band at
   $U=0$ splitting into a lower and an upper Hubbard band, with
   $\mu=U/2$ exactly by particle-hole symmetry. Reference: the Lieb-Wu
   charge gap of the infinite chain,
   $\Delta(U)=\frac{16t^2}{U}\int_1^\infty\frac{\sqrt{y^2-1}}{\sinh(2\pi ty/U)}dy$,
   drawn against the finite-chain gap from `info["gap"]`, and the exact
   $\mu=U/2$.
3. `Spinful_Fermionic_Chain_Native(12)`, hopping $-t$ with $t=1$, $U$ in
   `[0, 2, 4, 8]` for the spectra and `[0, 1, 2, 3, 4, 6, 8]` for the gap,
   `set_conserved_sector(Nf=12, Sz=0)`, `itensor_version=3`, `maxm=80`,
   `nsweeps=12`, `get_spectral_function(i=6, spin="up", nex=20, es=np.linspace(-8, 8, 800), delta=0.25, return_poles=True)`,
   as in `examples/dynamical_correlator/sector_spectral_function_hubbard`
   (which runs $n=6$). The calibration prints `info["captured"]` for each
   side and raises `nex` to 30 if it is below 0.9.
4. Two panels. Left: $A(\omega-\mu)$ for the four $U$, offset vertically,
   with a vertical line at $\omega=\mu$. Right: the gap against $U$ (blue
   points) with the Lieb-Wu curve dashed in ink, and $\mu$ (orange points)
   on the $U/2$ line. The reader sees the gap open and the finite-chain
   points lying slightly above the infinite-chain curve at small $U$.
5. About 8 minutes (seven values of $U$, three sector solves each on 12
   four-state sites plus `nex` excited states per side).

"Any pair of operators" gets no figure: figure 13 already uses an
arbitrary operator pair (a sum over sites).

### Choosing a correlator method (`correlator-methods.md`)

#### 15. `correlator-methods_comparison`

1. File `correlator-methods_comparison.png`, page
   `correlator-methods.md`, after the list of methods.
2. Every method computes the same quantity: the local
   $S^{zz}_{00}(\omega)$ of the Heisenberg chain from KPM, CVM, TD, TDZ,
   EX and SECTOR, each against the exact Lehmann sum
   (`mode="ED", submode="ED"`) with the same broadening. The reader sees
   where each method agrees and where it does not (KPM lines a little
   narrower and taller than the Lorentzian, see the page; EX missing the
   weight above its `nex` states).
3. Open spin-1/2 Heisenberg chain of $n=10$ sites, `itensor_version=3`,
   `maxm=60`, `kpmmaxm=60`, `cvm_maxm=60`, pair `(Sz[0], Sz[0])`,
   `es=np.linspace(-0.5, 3.5, 300)`, `delta=0.1` for all; TD and TDZ with
   `dt=0.05`; EX with `nex=20`; SECTOR with `nex=30`; CVM on a coarser
   `es` of 100 points (one linear solve per point).
4. One panel, the spectra stacked with a vertical offset, one method per
   color in palette order, the ED curve dashed in ink under each. The
   residual of each method against ED, as the maximum relative
   deviation, is printed in the legend. The reader sees six curves on top
   of their references.
5. About 8 minutes (CVM at 100 frequencies dominates).

### Real-time evolution (`real-time.md`)

#### 16. `real-time_light_cone`

1. File `real-time_light_cone.png`, page `real-time.md`, after the
   keyword list of the first section.
2. A light cone: the XX chain prepared in a domain wall
   $|\uparrow\cdots\uparrow\downarrow\cdots\downarrow\rangle$ and released.
   The magnetization front spreads at the maximal group velocity $v=J$,
   with the profile $\langle S^z_x(t)\rangle\to-\frac1\pi\arcsin(x/t)$
   inside the cone (Antal, Racz, Rakos and Schutz). The model is free
   fermions, so the exact profile at the same finite $n$ and $t$ is
   computed in `compute.py` by evolving the $n\times n$ correlation matrix
   with the hopping matrix, $C(t)=e^{iht}C(0)e^{-iht}$, and overlaid.
   The Jordan-Wigner map of $\sum_i(S^x_iS^x_{i+1}+S^y_iS^y_{i+1})$ gives
   the hopping amplitude $1/2$, so $h_{i,i+1}=h_{i+1,i}=1/2$, the
   dispersion is $\cos k$ and $v=1$; $C(0)$ is the diagonal matrix with
   ones on the left half, and $\langle S^z_i\rangle(t)=C_{ii}(t)-1/2$.
3. Spin-1/2 chain of $n=40$ sites, $H=\sum_i(S^x_iS^x_{i+1}+S^y_iS^y_{i+1})$
   ($J=1$, fermion dispersion $\cos k$, $v=1$), initial state the ground
   state of $-\sum_i\epsilon_iS^z_i$ with $\epsilon_i=+1$ on the left
   half and $-1$ on the right half (a product state, checked
   $\langle S^z_i\rangle=\pm1/2$ by a probe). Evolution with
   `timedependent.evolve_and_measure(sc, operator=Sz[n//2], nt=10, dt=0.05, wf=wf, return_wf=True)`
   in 32 chunks of 10 steps (so $t$ up to 16), measuring all
   $\langle S^z_i\rangle$ on the returned state with `sc.aMb(wf, Sz[i], wf)`
   after each chunk. `tevol_method="TDVP"` (the default), `maxm=100`,
   `itensor_version=3`. `evolve_and_measure` measures one operator, which
   is why the chunks are needed; `evolve_WF` would use the second-order
   MPO propagator instead of TDVP.
4. Two panels. Left: the map $\langle S^z_x(t)\rangle$ (diverging colors,
   `RdBu`, symmetric about 0), $x$ on $x$ and $t$ on $y$, with the lines
   $x-x_0=\pm t$ dashed in white. Right: the profile at $t=4$, 8 and 12
   (blue, orange, aqua) with the exact free-fermion profiles dashed in
   ink and the arcsine form dotted for $t=12$. The reader sees a cone of
   slope one and the TDVP profiles on the exact ones.
5. About 10 minutes (320 TDVP steps on 40 sites at `maxm=100`; the
   entanglement of this quench grows only logarithmically, so the bond
   dimension stays well below 100).

"Local quench" and "Evolved states" get no separate figure: the domain
wall is a local quench of the same kind, and the states route is the
same propagator.

### Local density of states (`dynamical-correlators.md`)

#### 17. `dynamical-correlators_ldos`

1. File `dynamical-correlators_ldos.png`, page `dynamical-correlators.md`,
   in "Local density of states".
2. The local density of states of a spinless fermion chain at a bulk site
   and at the edge site: at $V=0$ the exact free-fermion result on the
   same open chain, $A_i(\omega)=\sum_k|\varphi_k(i)|^2L_\delta(\omega-\varepsilon_k+\mu)$,
   is overlaid, and at $V=3t$ (inside the charge-density-wave phase,
   $V>2t$) a gap opens at $\mu$. The edge site shows the Friedel-like
   structure of the open chain.
3. `Fermionic_Chain(24)`, hopping $t=1$, `set_conserved_sector(Nf=12)`,
   `itensor_version=3`, `maxm=80`, `nsweeps=12`,
   `get_spectral_function(i=12, nex=30, es=np.linspace(-4, 4, 800), delta=0.1, return_poles=True)`
   and the same with `i=0`, for $V=0$ and $V=3$ with
   $V\sum_i(n_i-\tfrac12)(n_{i+1}-\tfrac12)$. The free-fermion reference
   from `numpy.linalg.eigh` of the $24\times24$ hopping matrix, with the
   same $\delta$ and the same $\mu$ convention (midpoint of the two
   levels around half filling).
4. Two panels, bulk site left and edge site right: $A_i(\omega)$ against
   $\omega-\mu$ for $V=0$ (blue) and $V=3$ (orange), the free-fermion
   curve dashed in ink over the blue one. The reader sees 24 levels
   merging into the $1/\sqrt{4t^2-\omega^2}$ shape at $V=0$, with the
   SECTOR poles on the exact levels, and a gap at $V=3$.
5. About 5 minutes (four sector solves and $2\times30$ excited states per
   case; `captured` printed, `nex` raised to 40 if below 0.9).

### STM and Kondo spectra (`kondo.md`)

#### 18. `kondo_dIdV`

1. File `kondo_dIdV.png`, page `kondo.md`, one figure after "Third
   order", with a reference from "ED and DMRG".
2. (a) Inelastic spin-flip steps: the second-order $dI/dV$ of a single
   spin with magnetic anisotropy, an $S=2$ atom with $D=-1.55$ meV and
   $E=0.31$ meV (the values of Fe on Cu$_2$N), at $B=0$: steps at the
   excitation energies of the anisotropic multiplet, which are known
   exactly from the $5\times5$ eigenvalues and are marked. (b) The Kondo
   peak of an $S=1/2$ at zero bias splitting in a field, third order,
   the paper's Fig. 3b: the $B=0$ peak height is $1.13$ in units of
   $e^2T_0^2/h$ and the tails at $\pm4$ mV are 0.886 (the values of
   `examples/kondo/kondo_spectrum_VS_paper`, which are pinned by tests).
   (c) The $T=0$ third-order spectrum of a 3-site chain from `mode="DMRG"`
   (TDVP two-time construction) on top of `mode="ED"`, as in
   `examples/kondo/kondo_third_order_dmrg_VS_ED`.
3. (a) `Spin_Chain(["S=2"])`, $H=D(S^z)^2+E[(S^x)^2-(S^y)^2]$ in eV,
   `get_kondo_spectrum(eV, order=2, T=0.5)`, `eV` from $-10$ to 10 meV,
   401 points. (b) `Spin_Chain(["S=1/2"])`, $H=g\mu_BBS^z$, $g=2$, $B$ in
   `[0, 2, 4, 6, 8, 10]` T, `Jrho_s=-0.05`, `U=0`, `T=1.0`, `order=3`,
   `omega0=20e-3`, `eV` from $-4$ to 4 meV, 161 points, divided by $2\pi$.
   (c) the 3-site chain, field and parameters of the example (`maxm=30`,
   `nsweeps=15`, `delta=2e-5`, `dt2=150`, `n_t2_half=20`, `dtau=150`,
   `n_tau_half=20`, `itensor_version=3`).
4. Three panels. Left: $dI/dV$ against $eV$ for the $S=2$ atom with
   vertical dotted lines at the exact excitation energies. Middle: the
   six field curves (palette order) with dotted ink lines at 1.13 and
   0.886. Right: ED (dashed ink) and DMRG (blue markers) for the 3-site
   chain. The reader sees steps at the right energies, a zero-bias peak
   splitting into two at $\pm g\mu_BB$, and DMRG on top of ED.
5. (a) and (b) under 1 minute; (c) about 10 minutes (the two-time grids
   are 41 TDVP evolutions on 3 sites); the example runs it today.

"Orbital-resolved IETS of a magnetic atom" gets no figure: it is the
same step spectrum on a specific atom model, and the example covers it.

### Finite temperature (`finite-temperature.md`)

#### 19. `finite-temperature_specific_heat`

1. File `finite-temperature_specific_heat.png`, page
   `finite-temperature.md`, after "Purification", with the METTS points
   referenced from "METTS".
2. The thermal energy and the specific heat of a Heisenberg chain against
   temperature, from purification and from METTS, on top of the exact
   result from the full spectrum. $C(T)=(\langle H^2\rangle_T-\langle H\rangle_T^2)/T^2$
   is computed directly in the purified state (no numerical derivative),
   and shows the broad peak near $T\approx0.5J$ of the antiferromagnetic
   chain. A probe at $n=8$ gave the purification energy equal to the exact
   one to 6 digits at $T=0.25$ and $T=1$.
3. Spin-1/2 Heisenberg chain of $n=10$ sites (the purified chain has 20).
   Purification: `thermal.Thermal_Spin_Chain(spins, T=T)`, `tc.MBChain.maxm=100`,
   `itensor_version=3`, $T$ in `np.geomspace(0.1, 5, 20)`,
   `E=tc.MBChain.vev(h)`, `E2=tc.MBChain.vev(h, npow=2)`. METTS:
   `sc.metts_vev(h, T, nsamples=300, nwarmup=30, seed=1)` at $T$ in
   `[0.25, 0.5, 1, 2, 4]`, `itensor_version=3`. Exact: `eigvalsh` of the
   $1024\times1024$ ED Hamiltonian, $E(T)$ and $C(T)$ from the Boltzmann
   sums.
4. Two panels. Left: $E(T)/n$ from purification (blue line), METTS (orange
   points with error bars) and exact (dashed ink). Right: $C(T)/n$ from
   purification (blue) and exact (dashed ink), log scale on $T$. The
   reader sees the purification line on the exact curve at every $T$, the
   METTS points within their error bars, and the specific-heat peak.
5. About 8 minutes (purification at $T=0.1$ is the slowest point, about
   `2*beta*W` factor applications, around 1 minute at $n=10$; METTS 1500
   samples of a 10-site chain, a few minutes).

`metts_dynamical_correlator` and the ED thermal correlator get no figure:
the example `dynamical_metts_VS_ED` covers them and a time-domain
Monte Carlo curve carries no physics the reader recognizes at a glance.

### Infinite chains (`infinite-chains.md`)

#### 20. `infinite-chains_energy`

1. File `infinite-chains_energy.png`, page `infinite-chains.md`, after
   "Two solvers", with the right panel referenced from "Static correlators".
2. The energy per site of the infinite Heisenberg chain converging to the
   Bethe-ansatz value $e_\infty=1/4-\ln2=-0.443147$ as the bond
   dimension grows (an algebraic convergence, since the chain is
   gapless), and the correlator $\langle S^z_0S^z_r\rangle$ at any
   distance with no edge, decaying as a power law. A probe with
   `gs_method="idmrg"` on `itensor_version=3` gave errors $1.7\times10^{-4}$,
   $3.1\times10^{-5}$ and $1.4\times10^{-5}$ at `maxm` 10, 20 and 40 in
   under a second each (`maxiter=60`, not converged to `etol`, as the
   page says happens on a gapless model). VUMPS on the same model at
   `maxm=30` took 4 minutes and did not converge either, so iDMRG is the
   solver of this figure.
3. `Infinite_Spin_Chain(["1/2"])`, `gs_method="idmrg"`, `itensor_version=3`,
   `maxm` in `[8, 12, 16, 24, 32, 48, 64, 96]`, `maxiter=400`,
   `etol=1e-10`, `ic.converged` stored for each. Correlators
   `ic.correlator("Sz", 0, "Sz", r)` for $r=1$ to 40 at `maxm=96`.
4. Two panels. Left: $|e(\chi)-e_\infty|$ against $\chi$ on a log-log
   scale (blue points). Right: $|\langle S^z_0S^z_r\rangle|$ against $r$,
   log-log, with the Luttinger form $a\sqrt{\ln r}/r$ dashed in ink
   fitted at large $r$. The reader sees the error fall by three orders
   of magnitude and a power law with no finite-size bending.
5. About 5 minutes (the `maxm=96` solve with 400 iterations dominates).

#### 21. `infinite-chains_excitations`

1. File `infinite-chains_excitations.png`, page `infinite-chains.md`,
   after "Dynamical structure factor".
2. Excitation bands from the quasiparticle ansatz. (a) The
   transverse-field Ising chain in the paramagnetic phase: the single
   mode $E(k)$ from VUMPS on the exact free-fermion band
   $2\sqrt{1+g^2-2g\cos k}$, drawn over the map $S^{xx}(k,\omega)$ from
   `dynamical_structure_factor`. A probe at $D=4$ reproduced the band to
   $6\times10^{-9}$ at 21 momenta in 0.3 s. (b) The magnon branch of the
   spin-1 Heisenberg chain over $S^{zz}(k,\omega)$, with the Haldane gap
   $0.4105$ at $k=\pi$ marked, as in `examples/idmrg/haldane_structure_factor`.
3. (a) `Infinite_Spin_Chain(["1/2"], itensor_version="python")`,
   $H=-4S^xS^x-2gS^z$ with $g=2.5$, `gs_method="vumps"`, `maxm=4`,
   `maxiter=800`, `vumps_nrestarts=6`, `ks=np.linspace(-np.pi, np.pi, 61)`,
   `excitation_energies(k, n=1)`,
   `dynamical_structure_factor("Sx", ks=ks, energies=np.linspace(1.5, 8.5, 400), delta=0.08, n=1)`.
   (b) `Infinite_Spin_Chain(["1"], itensor_version="python")`, Heisenberg,
   `gs_method="vumps"`, `maxm=10`, `maxiter=800`, `vumps_nrestarts=6`,
   `ks=np.linspace(0, 2*np.pi, 49)`, `spectral_weights("Sz", k, n=12, return_total=True)`,
   map broadened with `delta=0.08` on `np.linspace(0, 6, 400)` as in the
   example; the random VUMPS start is recorded by storing `e0` and
   `converged`. Excitations need `"python"`; see the page.
4. Two panels, both maps in `magma` with the branch as white dashes.
   Left: TFIM, with the exact band dotted in white on top of the VUMPS
   branch. Right: Haldane chain, log color scale (the weight spans five
   orders of magnitude across the zone), with a horizontal dotted line
   at $\omega=0.4105$ touching the branch at $k=\pi$. The reader sees the
   VUMPS band on the exact one, and the magnon minimum at $k=\pi$ at the
   Haldane gap.
5. About 10 minutes, almost all in (b) (a $D=10$ spin-1 VUMPS solve with
   6 restarts and 49 momenta with 12 branches each); (a) is seconds. If
   (b) exceeds 30 minutes, `maxm=8` is still within 1 per cent of the gap.

`kpm_finite`, `td_dynamical_correlator`, `local_excitation_gap` and the
lower-level tools get no figure: `td_dynamical_correlator` gives a second
$S(k,\omega)$ of the same kind as figure 21 at a higher cost, and the
others are checks rather than physics pictures.

### Cookbook (`cookbook.md`): no new figure

Every recipe points to a page that now has the figure (figures 4, 7, 11,
13, 22).

### Post-processing tools (`post-processing.md`)

#### 22. `post-processing_extrapolation`

1. File `post-processing_extrapolation.png`, page `post-processing.md`,
   after the `size_extrapolation` bullet.
2. The energy per site of open Heisenberg chains against $1/L$, with the
   polynomial extrapolation of `size_extrapolation` to $L\to\infty$ and
   the Bethe-ansatz value $1/4-\ln2$ as the reference. On an open chain
   the leading correction is the boundary energy, $E/L=e_\infty+a/L+b/L^2$,
   so the fit is quadratic in $1/L$ and the extrapolated value is
   printed next to the exact one.
3. Spin-1/2 Heisenberg chains, $L$ in `[16, 24, 32, 48, 64, 96, 128]`,
   `itensor_version=3`, `maxm=100`, `nsweeps=16`,
   `extrapolate.size_extrapolation(Ls, E/L, n=2)` and `n=1` for
   comparison.
4. One panel, $E/L$ against $1/L$ (blue points), the quadratic fit
   (dashed ink) extended to $1/L=0$, and the Bethe-ansatz value as a
   marker on the $y$ axis. The reader sees the points on a slightly
   curved line reaching $-0.4431$ at the axis.
5. About 5 minutes (the $L=128$ solve at `maxm=100` is about 2 minutes).

### Performance (`performance.md`): no figure

The page quotes measured speedups from `docs/gpu_cpu_performance.md`,
which needs a data-centre GPU; nothing to compute here. The BLAS-thread
and parallel-scan sections are one-command instructions; a timing plot
would be hardware-specific and is not physics.

## Summary

| id | page | runtime (one P-core) |
|---|---|---|
| 1 `ground-state_convergence` | ground-state | 1 min |
| 2 `ground-state_correlators` | ground-state | 3 min |
| 3 `ground-state_sectors` | ground-state | 2 min |
| 4 `excited-states_haldane` | excited-states (and degeneracy) | 8 min |
| 5 `excited-states_nonhermitian` | excited-states | 5 min |
| 6 `excited-states_heff` | excited-states | 1 min |
| 7 `entanglement_central_charge` | entanglement | 12 min |
| 8 `entanglement_natural_orbitals` | entanglement | 4 min |
| 9 `correlators_distribution` | correlators | 2 min |
| 10 `topology_berry_phase` | topology | 3 min |
| 11 `fidelity_ising` | fidelity | 10 min |
| 12 `mean-field_neel` | mean-field | 4 min |
| 13 `dynamical-correlators_sqw` | dynamical-correlators (and index) | 20 min |
| 14 `dynamical-correlators_spectral_function` | dynamical-correlators | 8 min |
| 15 `correlator-methods_comparison` | correlator-methods | 8 min |
| 16 `real-time_light_cone` | real-time | 10 min |
| 17 `dynamical-correlators_ldos` | dynamical-correlators | 5 min |
| 18 `kondo_dIdV` | kondo | 11 min |
| 19 `finite-temperature_specific_heat` | finite-temperature | 8 min |
| 20 `infinite-chains_energy` | infinite-chains | 5 min |
| 21 `infinite-chains_excitations` | infinite-chains | 10 min |
| 22 `post-processing_extrapolation` | post-processing | 5 min |

About 2.5 hours of single-core time in total, under one hour of wall
time with three jobs. No figure is planned above 30 minutes; figure 13
is the only one that can reach it, and its calibration step says so
before the full run.

Capabilities judged not worth a figure: the chain classes and operator
API (models, operators), `get_rdm`, the site and pair entropies and the
four-point tensor, the charge gap and the generalized eigenproblem, the
degeneracy counter on its own, `get_distribution_moments`, the
orbital-resolved atom IETS, the METTS and ED thermal correlators,
`kpm_finite`, `td_dynamical_correlator`, `local_excitation_gap`, the
infinite-chain lower-level tools, the cookbook, GPU and BLAS threads.

## Order of work

1. Write `style.py` and `run.sh`, and the pandoc regex in
   `build_guide_pdf.py`.
2. Run the cheap ED-referenced figures first (1, 3, 6, 9, 10, 18a-b),
   which also settle the style on real plots.
3. Run the calibration loops of 13, 14 and 17 and fix `kpmmaxm` and
   `nex` from them.
4. Run the rest, three at a time, longest first (13, 7, 18, 11, 16, 21).
5. Embed each figure in its page with a one-sentence caption that reads
   the physics off the plot, and run `mkdocs build` (strict) and
   `docs/build_guide_pdf.py`.
