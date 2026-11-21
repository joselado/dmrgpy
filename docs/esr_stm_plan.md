# ESR-STM spectra: plan for ED and MPS

This note says what the ESR-STM signal of a many-body spin system is,
which published method we follow and why, and how we implement it on the
ED and on the MPS side of dmrgpy. The weak-drive limit of the formula of
"The formula we implement" was checked numerically, to second order in
the drive, against an exact Floquet solution of the driven master
equation on two spins (`docs/esr_stm_plan_check.py`, see "The check
already made"). The plan was implemented on 2026-11-20; the next section
says what differs from it.

## What was built

`get_esr_spectrum` (`esrstm.py`, `esrtk/`, `tests/test_esr_stm.py`, the
guide page `docs/guide/esr-stm.md`, `docs/development.md` 4.8l) follows
the plan below, with the three picks of "Open choices". Saturation is now
checked too: one spin 1/2 gives the Bloch formula to 5e-14 at a
saturation parameter of 4.5, and the Floquet solution at finite drive
(`method="floquet"`) agrees to 6e-5 of the peak. The MPS route agrees with
ED to 2e-14 on a peak of 1.5e-6 on a 6-site chain, on `"python"` and v3.
Seven things differ from the plan:

- The degeneracy rule raises on an exact degeneracy only (1e-10 of the
  width of the kept spectrum) between states coupled by V, D or the
  relaxation, and not on two levels closer than a line width. The jumps of
  eq. 42 are between eigenstates, so the model is secular by construction,
  and a raise on levels closer than a few MHz would refuse any chain with
  two coupled levels that close, which a spectrum of hundreds of states
  within a meV can be expected to have (an estimate, not measured).
- The relaxation is block diagonal over the states the spin operators
  connect: S does not move charge, so on a fermion chain the empty and
  doubly occupied states form blocks of their own and keep their thermal
  weight. A drive that connects two blocks raises.
- The kept states end 30 free-spin widths above the window, not 1000
  (with 1000 the DMRG route asked for 160 GHz of states at $T_2$ = 2 ns),
  plus $kT\ln(1/p_{\rm cut})$, the populated states, so that $p_{\rm cut}$
  is also the accuracy of the DC part. The first commit kept a second
  $kT\ln(1/p_{\rm cut})$ for the thermal upward rates out of the upper
  state of a line: 6.6 times more states at the default $p_{\rm cut}$, for
  a window dependence that turned out to be the homodyne tail of the next
  point. Among the kept states, pairs are still taken up to 1000 widths
  away.
- The homodyne part of each line is cut smoothly between 15 and 30
  free-spin widths from its resonance (picked for the author on
  2026-11-21). Its dispersive tail falls only as one over the detuning, so
  without the cut it depended on which states are kept, by 1e-2 to 4e-2
  of the DC line under a zoom of the window. The background of the far
  lines needs every state up to the exchange scale and the counter-rotating
  terms, and experiments subtract a smooth background.
- `f_ref` is a required keyword (the author's pick of 2026-11-20): T1 and
  T2 are given at that frequency. The first version used the centre of the
  window, and zooming onto a line moved it by 2 to 8%.
- The interface takes `T1` and `T2` of a free spin 1/2 (the numbers an
  experiment quotes) instead of `Tphi`, with $1/T_\phi=1/T_2-1/(2T_1)$, and
  on the DMRG route the degeneracy tolerance adds the energy uncertainty
  $\sigma_n+\sigma_m$ of the two states.
- The third example is a spin chain with ED against DMRG
  (`examples/kondo/esr_spin_chain_dmrg_VS_ED`) instead of the Haldane
  chain, and the dimer example is `esr_dimer_singlet_triplet`.

One result was not in the plan. At zero tip field the S-T0 line of the
dimer keeps a small DC part, from the relaxation that carries the driven
population of S and T0 on to T- and T+: 1.2e-3 of the homodyne part at
rabi = 1e-4 and V_RF/V_DC = 0.2, and growing as rabi T1 over V_RF/V_DC.
Not done: the polarized-tip rates (the first choice, left for later), Sz
sectors on the MPS route, the digitized peak heights of Yang et al.
Fig. 3, the PyESRSTM cross-check, a run on `julia_live`, and section S5 of
Bae et al.

## The Krylov routes (2026-11-21)

The author asked whether the final states of the lines could come from a
Krylov space instead of from excited-state DMRG, since a long chain needs
more excited states than DMRG can give, and picked the regime of an
exchange well above kT. Fable took the decisions of the night. What was
measured and built:

- A Krylov space can only replace the final states. The populated states
  and every state below the window are still needed as eigenstates: the
  upper state of a line relaxes into them. Without the one empty state at
  10 GHz of a 14-spin chain, the widths were 5% off and the curve 8%. At
  the default `p_cut`, on 14 spins, 110 of the 139 needed states are
  populated at J = 240 GHz and T = 1.2 K, and 10 of 99 at J = 30 GHz and
  T = 0.05 K. So the Krylov routes pay in the second regime only. The ED
  route takes every state below the window since the 2026-11-21 hunt
  (finding 2).
- One common basis is needed. One two-level Bloch line per transition, as
  separate Krylov spaces give, is 46 to 61% of the peak off.
- A Krylov space of full MPS does not pay: it needs 3 to 4 vectors per
  eigenstate it replaces (273 vectors for 100 eigenstates at 2e-5 of the
  peak on 12 spins).
- ED, `submode="KRYLOV"`: every eigenstate below the window, one block
  Krylov space of V|m>, and a closing one from S^a_k|n> and V|n>
  (findings 1 to 3 of the 2026-11-21 hunt). Measured before the closing
  stage, it was 3.2 times faster than the eigenbasis on 16 spins at
  J = 30 GHz and T = 0.05 K, but 0.45 times at J = 240 GHz and T = 1.2 K,
  so it is opt-in.
- DMRG, `submode="KRYLOV"`: the Krylov space runs in the tangent space of
  each populated state (TASK), where its vectors are cheap. The TASK poles
  out of the ground state of 14 spins were off by a median of 1 to 1.7
  half widths at nkry=100, and by 2.4e-5 (maxm=32) and 3.4e-7 (maxm=64)
  at nkry=300: the error was the Krylov space, not the tangent
  projection. With the states below the window from DMRG, nkry=600 (300
  left 14 window lines unconverged) and one common basis, the curve agrees
  with ED to 3e-7 of the peak on 14 spins at J = 30 GHz, T = 0.05 K.
  It took 318 s on v3 on one E-core (328 s after the last changes, at the
  same accuracy), where the eigenstate route spent 1376 s on its 58
  eigenstates alone and then refused the result: two of its states were
  closer than their energy errors at maxm=64. On 12 spins the eigenstates
  of the common basis between 35.9 and 39.2 GHz, above the window and its
  margin, where the basis does not span every eigenstate, came out with
  energy errors of 0.1 to 2 GHz (1.4 to 25 free-spin half widths), the
  same for a cutoff of 1e-12 or 1e-14 and for an overlap cut of 1e-10,
  1e-7 or 1e-4: they are mixtures, not noise. Their lines from the second
  populated state end in the margin, so their errors are kept only for
  lines that reach the window give or take the error
  (`lineshape.READ_MARGIN`), in the degeneracy tolerance and in the
  position warning; the first version of this said 3 GHz, which was the
  error of a pair misread. A drive with terms on two MPS sites (the spin
  of a spinful fermion) is refused on this route, since its start in the
  tangent space would only be a projection.
- A read-only review of the night's diff found a silent gap: states that
  V never reaches (on a chain without a tip field in H, the singlets and
  quintets when the drive starts from a singlet) were missing as
  relaxation partners. On ED the Krylov route keeps every Ritz state
  below emax and agreed with ED to 1e-12 on 8 and 10 spins with or
  without the tip field (with a tilted tip; with the tip moment along the
  field it kept one Sz_tot sector, finding 1 of the 2026-11-21 hunt); on
  DMRG the window-line states alone gave 4.3e-5
  of the peak on 10 spins without the tip field, and with every converged
  Ritz state between the window bottom and emax as a partner 1.1e-6
  (1.2e-7 with the tip field), for 15% more time. On 14 spins with the tip
  field the partners cost 488 s against 328 s, at 1.1e-6 against 3.3e-7.
  The same review found the homodyne cut inside a power-broadened line;
  it now counts in the width of the line.
- Left open: the crossover between the routes, as a function of the
  populated fraction; chains above 14 spins; the MPS sums that form the
  eigenstates of the common basis, a third of the time; how nkry has to
  grow with the length; `julia_live`; and the regime of an exchange close
  to kT, where almost every state is populated and only a typicality
  method (finite-temperature Lanczos) would avoid eigenstates.

## Why it is not the spin spectral function

The d2I/dV2 of IETS at T=0 is, up to a constant,

$$\frac{d^2I}{dV^2}(eV)\propto\sum_{a=x,y,z}\sum_n
|\langle n|S^a_i|0\rangle|^2\,\delta(eV-\omega_n),$$

since dI/dV has a step at each pole. So d2I/dV2 is the local spin spectral
function of the tip site i itself, all three components, each pole
weighted only by its matrix element, and broadened by the temperature of
the leads (a width of about 5.4 kT, so 0.5 meV at 1 K).

The ESR-STM signal is a different object, even if its poles are the same
transitions. Let us take the simplest case, T=0, weak drive, and every
excited state decaying back to the ground state with a time T1. The DC
part of the signal is then

$$\Delta I_{\rm DC}(\omega)\propto\Omega^2T_1\sum_n
|\langle n|V|0\rangle|^2\,
\big(\langle n|D|n\rangle-\langle 0|D|0\rangle\big)\,
\frac{\gamma}{(\omega-\omega_n)^2+\gamma^2},$$

meaning that each pole of the local spectral function of the drive
operator V gets multiplied by a form factor, the change of the detection
operator D between the final and the initial state. Four things separate
it from d2I/dV2:

- V is one component, set by the driving mechanism, not the sum over a.
  For every mechanism with the tip on site i (piezoelectric modulation of
  the tip exchange or dipolar field, exchange field of the tip, modulated
  g-factor), V is a one-site operator $\mathbf b\cdot\mathbf S_i$, with
  $\mathbf b$ along the tip moment for the exchange-type mechanisms.
- The form factor $\Delta D_n$ with $D=\mathbf m\cdot\mathbf S_i$, $\mathbf m$
  the tip moment, comes from the tunnel magnetoresistance of the tip site.
  It can vanish or change sign. On a Ti dimer the singlet to $T_0$
  transition has $\Delta D=0$ at zero tip field, so it is invisible to the
  DC readout and is seen only through the homodyne part (Bae et al.,
  arXiv:1807.04322, Fig. 4). A local spectral function shows it as a full peak.
- At finite T the weight is the population difference $p_m-p_n$, and at
  1.2 K and 10 to 30 GHz the excited Zeeman states are populated
  ($e^{-hf/kT}=0.45$ at 20 GHz), so transitions between excited states
  are part of the spectrum.
- The line is a Lorentzian of width $1/T_2$, of order 1 to 10 MHz, that is
  4e-9 to 4e-8 eV. It is asymmetric when the homodyne part is present, and
  it saturates and broadens with the drive.

The assumption that the ESR-STM signal is the local dynamical structure
factor $-{\rm Im}\,G^R_{S^+_jS^-_j}$ (Maiellaro et al., arXiv:2312.07147) keeps
the poles and drops the form factor, the populations and the homodyne
part. Those authors say themselves that a better observable is an open
question.

## What the published theory says about the line shape

The recent theory review of Taherpour, Jankovic, Reina-Galvez, Le and Wolf
(arXiv:2606.27830) compares Heisenberg, Kondo and Anderson models against
the same experiments. Its section on CW spectra finds that the simulated
spectra of all three reproduce the line shape, the peak intensities and
the peak splittings of a Ti atom in a Ti-Ti dimer, and that the line shape
depends mostly on the angle between the spin and the tip moment. Its
Kondo panel there is taken from Ye et al., whose readout differs from ours
(see below), so the dimer spectrum alone does not test the readout; the
angle dependence of one spin does, and it is one of the tests. So the CW line shape does not tell the
driving mechanisms apart; what tells them apart is how the Rabi rate,
$T_2$ and the resonance frequency scale with the DC current and the DC
bias. The consequence for us is that the right method for CW spectra is
the master equation in the eigenbasis of H, with a local drive whose
strength is an input, and a readout made of two parts:

- The DC magnetoresistance part, the change of the time-averaged
  populations read through $\langle D\rangle$. Ast et al.
  (arXiv:2403.20247) derive it from a Kondo-type tunnelling Hamiltonian
  $c^\dagger_S(\delta_{\sigma\sigma'}+\lambda\,\boldsymbol\sigma\cdot\mathbf S)c_T$,
  the same form as the Ternes Hamiltonian of `get_kondo_spectrum`. The
  dominant term is the interference of the elastic and the spin-dependent
  channels, which is linear in $\langle\mathbf S\rangle$ and in the tip
  polarization. For one spin 1/2 they recover the Bloch formula
  $\Omega^2T_1T_2/(\delta^2T_2^2+1+\Omega^2T_1T_2)$, their eq. 41, with
  $T_1$ and $T_2$ set by the tunnelling rates.
- The homodyne part, the coherence of the driven transition precessing at
  the RF frequency, multiplied by the RF modulation of the conductance,
  which gives a DC current. It is odd in the detuning for a real phase, and
  it is what makes the line asymmetric (Ast et al. eq. 45, Delgado and
  Lorente review arXiv:2208.14105 eq. 12, Kim et al. arXiv:2103.09582 eq.
  S19, the review's eq. 34).

Ye, Zheng and Xu (arXiv:2402.01435) use a readout of the transverse spin
only, $\tilde I\propto\sin\theta\,\tilde{\mathbf S}\cdot\mathbf e_\perp$; the
review finds that this readout fails to reproduce the observed homodyne
response, since the current of a polarized tip reads $\mathbf m\cdot\mathbf S$,
so we do not follow it. Shavit et al. explain the asymmetry without homodyne detection, by
dropping the secular approximation of the dissipation; the review counts
it as the same physics in another form.

Two open-source codes exist from the group of Lorente: TimeESR
(github.com/qphensurf/TimeESR, Fortran, GPL-3.0, time-domain master
equation of one orbital coupled to a spin system) and PyESRSTM
(github.com/qphensurf/PyESRSTM, Python, AGPL-3.0, Floquet master equation,
Anderson impurity, ESR sweeps for one and two spins). Both are full ED of
a small impurity. We use PyESRSTM as an external reference for one and two
spins, run next to dmrgpy and not copied into it (AGPL).

The supplement of Yang et al. (PRL 119, 227206, open on the APS site; it
is not part of arXiv:1712.04082) models the peak amplitudes of a Ti dimer
in its Section 10 with linear response,

$$A(f)\propto\sum_nP_n\sum_{n'}|\langle n|S^a_1|n'\rangle|^2\,\delta\big(f-(f_n-f_{n'})\big),\qquad\rm{(S13)}$$

the thermal absorption spectrum of the tip-site spin, which is the
spectral-function picture with a Boltzmann prefactor. For the dimer of
its Section 7 this gives the weights $P_{00}\sin^2(\xi/2)$ and
$P_+\sin^2(\xi/2)$ for the lines I and II, and $P_{00}\cos^2(\xi/2)$ and
$P_-\cos^2(\xi/2)$ for III and IV (eq. S14), with $\tan\xi=(J-D)/\gamma\hbar
B_{\rm tip}$. The DC part of our formula keeps these matrix elements and
multiplies each line by its form factor, which for this dimer is exactly
$\Delta\langle S^z_1\rangle=\sin^2(\xi/2)$ for I and II and $\cos^2(\xi/2)$ for
III and IV, and it replaces $P_n$ by $P_n-P_{n'}$. So the weak lines I and
II fall as $\sin^4(\xi/2)$ instead of $\sin^2(\xi/2)$ when the tip field
grows. At J = 1 GHz, D = 0, B = 0.9 T and 1.2 K the ratio of line I to line III
is 0.67 in S14 and 0.44 in our formula at a tip field of 0.2 GHz, and 0.17
against 0.028 at 1 GHz (`docs/esr_stm_plan_yang_weights.py`). The peak
heights against the tip field of Yang et al. Fig. 3 decide between the
two, which makes them the reference for the dimer test. Section S5 of
Bae et al. (Science Advances, with the full dimer line shape) was not
read for this plan.

## The formula we implement

H is the Hamiltonian of the chain, including the static tip field on the
tip site, so every tip position is its own Hamiltonian. Its eigenstates
are $|n\rangle$ with energies $E_n$, and $\omega_{nm}=(E_n-E_m)/\hbar$. The
RF drive is $\Omega V\cos\omega t$, V Hermitian, and D is the detection
operator. Every transition is taken as resolved, meaning that two lines
closer than their width are the case where this formula is approximate
(see "The check already made"). Then the master equation reduces to three
pieces.

The driven rate between n and m, with $E_n>E_m$, is the golden rule with
a Lorentzian,

$$W_{nm}(\omega)=\frac{\Omega^2}{2}|V_{nm}|^2
\frac{\gamma_{nm}}{(\omega-\omega_{nm})^2+\gamma_{nm}^2},$$

$\gamma_{nm}=1/T_2$ of that transition, half the sum of the decay rates of
n and m plus the pure dephasing. The populations solve the rate equation
with the drive added to the relaxation matrix R,

$$\big(R+W(\omega)\big)\,p(\omega)=0,\qquad\sum_kp_k=1,$$

and $p^0$ is its solution without drive: the Boltzmann weights when R obeys
detailed balance at T, pumped populations when the tip current enters R.
Solved as written it includes saturation and power broadening, and for
one spin 1/2 it reduces to the Bloch formula, width
$(1/\pi T_2)\sqrt{1+\Omega^2T_1T_2}$ (derived, not run yet). With several
levels it holds while the Rabi rate $\Omega|V_{nm}|$ stays below the
separation between lines that share a level; at stronger drive the states
are dressed (the Autler-Townes splitting of the double-resonance
experiments), the rate equation fails, and only a Floquet solution at
finite drive describes it. The coherence of each transition is slaved to
its populations,

$$\rho_{nm}(t)=\frac{\Omega}{2}\,\frac{V_{nm}\,(p_m-p_n)}{\omega-\omega_{nm}+i\gamma_{nm}}\,e^{-i\omega t}.$$

The current is $I(t)=G_ja\,[V_{\rm DC}+V_{\rm RF}\cos(\omega t+\varphi)]\,
\langle D\rangle(t)$ plus terms that do not depend on the spin, so its DC
change in units of $G_jaV_{\rm DC}$ is

$$\Delta I(\omega)=\sum_k\big(p_k(\omega)-p^0_k\big)D_{kk}
+\frac{V_{\rm RF}}{V_{\rm DC}}\,{\rm Re}\Big[e^{i\varphi}\sum_{n,m}\rho_{nm}D_{mn}\Big],$$

the DC part reading the diagonal of D and the homodyne part its
off-diagonal elements. For one spin 1/2 this is eq. 12 of the Delgado and
Lorente review, with the asymmetry $T_2V_{\rm RF}/(2\Omega T_1V_{\rm DC})$.

Degenerate levels need a rule, because $p_k$ and $D_{kk}$ depend on the
basis chosen inside a degenerate level, and neither the eigensolver nor the
Ritz step of the MPS route fixes that basis. In the secular approximation
the coherences inside a level do not rotate, so the right object there is
a block of the density matrix, not a population. The first version raises
when two levels connected by V, D or R are closer than the smallest
$\gamma_{nm}$, and names the two usual causes, a zero field and a level
crossing in a field. In the experiments the field and the tip field lift
every degeneracy, so the raise costs nothing there; near a crossing the
Floquet reference of the ED route is the answer. A block treatment can
come later if a case needs it.

The relaxation matrix R is the one choice left. The first version takes
the collapse operators between eigenstates of the review's eq. 42 (from
Broekhoven et al., npj Quantum Inf. 10, 92), a Korringa rate through every
site,

$$\Gamma_{n\to m}=\sum_k g_k\sum_a|\langle m|S^a_k|n\rangle|^2\,
\frac{\epsilon_{nm}}{1-e^{-\epsilon_{nm}/kT}},$$

with one scale $g_k$ per site (equal on all sites by default, the tip site
optionally larger), plus a pure dephasing rate. This obeys detailed
balance, so $p^0$ is thermal. Microscopic rates with a polarized tip and
a DC bias, which give $1/T_2\propto I_{\rm DC}$ and spin pumping away from
the Boltzmann populations, need polarized leads, which `conductance.py`
does not have (its sums are unpolarized). That is the first open choice
below.

## The check already made

`docs/esr_stm_plan_check.py` takes two spins 1/2 with J = 1 GHz, B = 20
GHz, a tip field from 0.05 to 2 GHz on spin 0, the tip moment at 0.6 rad
from the field, V = D = m·S_0, eigenbasis jumps with detailed balance at
1.2 K, and pure dephasing. It compares the formula above, with the full
rate matrix, against the exact O(Omega^2) steady state of the driven
Lindblad equation, obtained from its first Floquet harmonics with no
secular and no rotating-wave approximation. At the peaks the DC part
agrees to 5e-5 to 6e-4 of its value on the strong lines. The worst case is
1.2% to 2% on a weak line 0.12 GHz from a strong one, with widths of 0.01
GHz, which is the expected error of the secular approximation. The
homodyne part agrees to 1e-5 in its absorptive part and to 4e-4 in its
dispersive part, the difference being the counter-rotating terms. The
singlet to $T_0$ line has $\Delta D$ = 0.041, 0.369 and 0.738 at tip
fields 0.05, 0.5 and 2 GHz, so its DC part grows from 0.17 to 1.3 while
its homodyne part stays of the same size, which is what Bae et al. see
(antisymmetric at small tip field, more symmetric at large tip field).
The check is at weak drive, Omega = 1e-4 GHz: it tests the populations and
coherences to second order in the drive, not the saturated rate equation.

## ED route

The ED route builds the eigenbasis once per Hamiltonian, as
`kondospectrumtk/edkondo.KondoSpectrum` does: full `eigh` up to 2000
states, with V, D and $S^a_k$ of every site taken into the eigenbasis by
`eigenbasis_operator`. Above 2000 states the lowest M eigenpairs of the
sparse H (eigsh) give the same manifold object the MPS route builds. The
rate equation is an M by M linear solve per frequency, so a sweep of a
few thousand frequencies costs nothing next to the diagonalization.

As a reference for small systems (up to about 64 states) the ED route
also gets the exact O(Omega^2) Floquet solution of the check above,
`method="floquet"`, for lines that overlap and for the comparison with
PyESRSTM.

## MPS route

The MPS route is the low-energy manifold: the ground state and the M-1
lowest excited states, from `get_excited_states` (penalty method, or
`method="state_averaged"` with `nguard=` on `"python"`), in a conserved
Sz sector when H has one and then promoted to dense indices, as
`sectordc.py` does for SECTOR. From those states we build $H_{ab}$, the
overlaps and the matrices of V, D and $S^a_k$, and diagonalize H in the
span, which removes the mixing of near-degenerate members that the
penalty method leaves. The ESR spectrum is then the same evaluator as on
ED. Three guards come with it:

- Completeness: the manifold must hold every state below
  $E_0+hf_{\max}+kT\ln(1/p_{\rm cut})$, checked by one guard state above
  it; otherwise the call raises.
- Eigenstates: the energy variance of every state (`vev.energy_variance`).
- Line positions: a transition energy is a difference of two DMRG
  energies, and it has to be right to below $1/T_2$, about 1e-8 eV, which
  is a relative 1e-4 of a 0.1 meV transition. The variance bounds it, and
  the output says when a line position is not resolved at that level.

What this route does not cover is a chain whose ESR window holds a large
number of states, for instance a gapless Ti chain with J of order hf at
kT of order J: M grows with the system and the route stops. The cases it
is made for are gapped chains (nanographene and porphyrin chains with J
of 1 to 10 meV, Haldane chains with their edge states) and clusters in a
tip field, where the window holds a few multiplets.

We do not use TASK, KPM or the correction vector from the ground state.
The populations of excited states are of order one at the temperatures of
the experiments, and the tangent space of the ground state does not hold
the excitations of an excited state (`.claude/rules/kondo.md` refuses a
thermal population of TASK Ritz states for this reason). The DC part also
needs $\langle n|D|n\rangle$ of the final states, which no two-point function
gives, and a width of 1e-8 eV rules out any broadened correlator: we need
the poles and an analytic line shape.

## Interface

The method goes on `Spin_Chain` and both spinful fermion chains through a
mixin, as `get_kondo_spectrum` does, with the implementation in a
`esrtk/` package (spectrum builders, rate matrices, line shape, Floquet
reference) behind a thin `esrstm.py`. Energies and temperatures follow
`get_kondo_spectrum` (eV and K, `kB=` to change), and frequencies are in
GHz. The relaxation keywords depend on the first open choice. With the
Korringa rates above, `T1` is the relaxation time a free spin 1/2 on the
tip site would have at the centre of the sweep, which fixes the common
scale $g$ through $1/T_1=(g/2)\,\epsilon\coth(\epsilon/2kT)$, and `Tphi`
is the pure dephasing time. A first sketch:

```python
f = np.linspace(17.5, 21.5, 2001)            # RF frequency in GHz
f, dI = sc.get_esr_spectrum(f, site=0,       # tip on spin 0
           tip=[0.56, 0., 0.83],             # tip moment m, D = m.S_0
           drive=None,                       # V, default m.S_0
           T=1.2, T1=200., Tphi=100.,        # K, ns, ns
           rabi=2e-3,                        # Rabi frequency of a free spin, GHz
           vrf_over_vdc=0.2, phase=0.,       # homodyne part
           mode="ED")                        # or "DMRG", nstates=M
```

with `return_parts=True` for the DC and homodyne parts separately, and
`return_poles=True` for the list of transitions (frequency, $|V_{nm}|^2$,
$\Delta D$, $p_m-p_n$, width), which is what one compares peak by peak
with an experiment.

## Tests, examples and documentation

- `tests/test_esr_stm.py`: one spin 1/2 against the Bloch formula with
  power broadening and against the asymmetry of the Delgado and Lorente
  eq. 12; one spin in the rotating frame against `lindblad.py`
  (`mode="ED"`), the exact steady state, which is also the test of
  saturation; the angle between field and tip moment for one spin, where
  the readout $\mathbf m\cdot\mathbf S$ gives a DC amplitude
  $\propto\sin^2\theta\cos\theta$, largest at a tilted angle and zero in
  plane, as the review reports, while the transverse readout of Ye et al.
  gives a large signal in plane (their Fig. 2b, which the review
  contrasts with the experiments), and Kim et al. eq. S19 for the
  asymmetry against the angle; the dimer
  of the check against the Floquet reference; the four lines of the Yang
  et al. dimer, with the weights of S14 times the form factors, and the
  peak heights of their Fig. 3 against the tip field, digitized, as the
  experimental reference between S13 and our formula; the singlet to $T_0$ DC part
  going to zero with the tip field; the raise on a degenerate level at
  zero field; ED against the MPS manifold on 8 to 12 sites; and the
  guards raising on an incomplete manifold and on a non-eigenstate.
- PyESRSTM on one and two spins as an external cross-check, recorded in
  the method's README and not run by the suite.
- `examples/kondo/esr_single_spin_bloch` (line shape and power
  broadening), `examples/kondo/esr_ti_dimer_singlet_triplet` (the four
  lines and the homodyne singlet line against the tip field), and
  `examples/kondo/esr_spin1_chain_edge` (ESR against the local structure
  factor on a Haldane chain, ED and DMRG), following the `new-example`
  rules.
- A guide page `docs/guide/esr-stm.md` next to "STM and Kondo spectra" in
  the `nav` of `mkdocs.yml`, written as the "Ground state" page, and a
  section in `docs/development.md` for the manifold route (`docs-sync`).

## Open choices

These are the decisions to take before writing code, one at a time:

- The relaxation model of the first version: Korringa rates of eq. 42 with
  thermal populations (as above), or also polarized-tip rates with bias
  and current (spin pumping, $1/T_2\propto I_{\rm DC}$), or a constant
  $T_1$ and $T_2$ for every transition. Taken by the author on
  2026-11-20: the Korringa rates of eq. 42 with thermal populations; the
  polarized-tip rates stay for later.
- Whether V and D are one operator by default (exchange-type driving) or
  V is always given separately. Taken by the author on 2026-11-20: V = D =
  m·S_tip by default, and `drive=` (a MultiOperator) overrides V.
- Where the guide text goes: a page of its own or a section of the STM
  and Kondo page. Taken by the author on 2026-11-20: a page of its own,
  `docs/guide/esr-stm.md`, next to "STM and Kondo spectra" in the nav.
- (Found while implementing.) What `f_ref` defaults to, the frequency at
  which T1 and T2 are given. Taken by the author on 2026-11-20: required,
  no default.

## Sources

- Taherpour et al., Theory of ESR-STM: the first decade, arXiv:2606.27830
- Ast et al., arXiv:2403.20247 (Kondo-type tunnelling and Lindblad, eqs. 41 and 45)
- Ye, Zheng, Xu, arXiv:2402.01435 (HEOM, Ti adatom and dimer)
- Delgado and Lorente, review, arXiv:2208.14105 (eqs. 11 to 13)
- Kim et al., arXiv:2103.09582 (DC and homodyne readout, eq. S19)
- Bae et al., arXiv:1807.04322 (singlet to T0 line, homodyne detection)
- Yang et al., arXiv:1712.04082 (eigenstates of coupled Ti spins)
- Reina-Galvez et al., arXiv:1904.08183, arXiv:2108.01011, arXiv:2303.09944, arXiv:2503.24046
- Maiellaro et al., arXiv:2312.07147 (local structure factor of spin-1 chains)
- TimeESR, github.com/qphensurf/TimeESR; PyESRSTM, github.com/qphensurf/PyESRSTM
