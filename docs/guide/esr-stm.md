# ESR-STM spectra

On this page we compute the electron spin resonance signal measured with
an STM (ESR-STM). A radio-frequency (RF) voltage on the tip drives
transitions between the spin states of the atom under the tip, and a
spin-polarized tip reads the change of populations through the tunnel
magnetoresistance. Scanning the RF frequency $f$ gives lines at the
transition frequencies, with widths of a few MHz. Here the sample is a
single spin-1/2 with a Zeeman splitting of 20 GHz, with $H$ in eV:

```python
import numpy as np
from dmrgpy import spinchain
h = 4.135667696e-6                     # Planck constant in eV/GHz
sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(h*20.0*sc.Sz[0])    # Zeeman splitting of 20 GHz
f = np.linspace(19.9, 20.1, 801)       # RF frequency in GHz
tip = [0.6, 0., 0.8]                   # direction of the tip moment
f, dI = sc.get_esr_spectrum(f, tip=tip, T=1.2, T1=100., T2=50., f_ref=20.,
                            rabi=1e-3)
```

The Hamiltonian holds everything static: the external field, and the
field of the tip on the site under it. So every tip position is its own
Hamiltonian.

## The signal

The RF voltage acts on the spin as a drive $\Omega V\cos(2\pi ft)$. For
the driving mechanisms with the tip on site $i$ (the exchange field of
the tip, or the tip field modulated by the motion of the atom) the drive
is $V=\mathbf m\cdot\mathbf S_i$, with $\mathbf m$ the tip moment. The
readout is the magnetoresistance of the same site, $D=\mathbf m\cdot\mathbf S_i$.
Between eigenstates $n$ (upper) and $m$ (lower) the drive gives the rate

$$W_{nm}=\frac{\Omega^2}{2}|V_{nm}|^2\frac{\gamma_{nm}}{(2\pi f-\omega_{nm})^2+\gamma_{nm}^2},$$

and the populations $p$ are the steady state of the rate equation with
the drive and the relaxation. The DC change of the current, in units of
$G_jaV_{\rm DC}$, is

$$\Delta I=\sum_k\big(p_k-p^0_k\big)D_{kk}
+\frac{V_{\rm RF}}{V_{\rm DC}}\,{\rm Re}\Big[e^{i\varphi}\sum_{n,m}\rho_{nm}D_{mn}\Big],$$

with $p^0$ the populations without drive and $\rho_{nm}$ the coherence of
each driven transition. At weak drive and $T=0$ the first term is

$$\Delta I\propto\sum_n|\langle n|V|0\rangle|^2\,\big(D_{nn}-D_{00}\big)\,
\frac{\gamma}{(2\pi f-\omega_n)^2+\gamma^2},$$

meaning that each line of the local spin spectral function gets a form
factor: the change of the spin under the tip between the two states. This
is what separates ESR-STM from IETS ([STM and Kondo
spectra](kondo.md#second-order)): one spin component instead of three,
the population difference $p_m-p_n$ instead of $p_m$, a form factor that
can vanish or change sign, and a Lorentzian of width $1/T_2$ instead of a
thermal step.

The relaxation between eigenstates is a Korringa rate through every site,
$\Gamma_{n\to m}\propto\sum_k\sum_a|\langle m|S^a_k|n\rangle|^2\,
\epsilon_{nm}/(1-e^{-\epsilon_{nm}/kT})$, which obeys detailed balance, so
$p^0$ is the Boltzmann distribution. Its scale is set by `T1` at the
frequency `f_ref`, and every coherence also decays by pure dephasing, set
by `T2`.

For the single spin, `dI` is a Lorentzian at 20 GHz of full width
$1/(\pi T_2)=6.4$ MHz and height 0.0101. It is positive here because the
tip moment has a positive component along the field.

## Saturation and the homodyne part

At stronger drive the line saturates and broadens. For one spin-1/2 this
is the Bloch formula, with the Rabi frequency $\Omega_R$ of the line,

$$\Delta I\propto\frac{\Omega_R^2T_1T_2}{1+(2\pi\,\delta f\,T_2)^2+\Omega_R^2T_1T_2},
\qquad {\rm FWHM}=\frac{\sqrt{1+\Omega_R^2T_1T_2}}{\pi T_2}.$$

With `rabi=1e-2` above, the width grows from 6.4 to 18 MHz and the height
from 0.0101 to 0.133. The second term of $\Delta I$, the homodyne part,
is the precessing spin read through the RF modulation of the conductance.
It is odd in the detuning for `phase=0` and makes the measured lines
asymmetric:

```python
f, dI, parts = sc.get_esr_spectrum(f, tip=tip, T=1.2, f_ref=20., rabi=1e-3,
                                   vrf_over_vdc=0.5, return_parts=True)
```

`parts["dc"]` is the symmetric line of height 0.0101 and
`parts["homodyne"]` an antisymmetric one between -0.0026 and 0.0026.
The dispersive tail of a line falls only as one over the detuning, so a
line enters the homodyne part only within 30 free-spin widths of its
resonance, $30/(\pi T_2)$ in frequency, with a smooth cut from 15 widths
on. The background of the lines further away is not computed: it needs
every state up to the exchange scale and the counter-rotating terms, and
experiments subtract a smooth background anyway. The keywords are:

- `tip`: the direction of the tip moment, which sets $D$ and the default
  $V$.
- `site`: the site under the tip.
- `drive`, `detection`: operators that replace $V$ and $D$, for instance
  the transverse spin `sc.Sx[0]` for a drive through a gradient of the tip
  field.
- `T`: the temperature in kelvin.
- `T1`, `T2`, `f_ref`: the relaxation and coherence times in ns of a free
  spin-1/2 on the tip site with Zeeman frequency `f_ref` (GHz), the
  frequency where they were measured. `f_ref` is required, since the
  relaxation rates grow with frequency.
- `rabi`: the Rabi frequency in GHz of a transition with
  $|V_{nm}|=1/2$. A small value gives the linear response.
- `vrf_over_vdc`, `phase`: the size and the phase of the homodyne part.
- `site_weights`: the weight of each site in the relaxation (equal by
  default).
- `h`, `kB`: the units. With $H$ in GHz pass `h=1` and `kB=20.8366`.

## Coupled spins

For several spins the lines are transitions between many-body states, and
`return_poles=True` lists them. As an example we take two spins with an
exchange of 15 GHz, a Zeeman energy of 25 GHz and a tip field of 4 GHz on
the first spin:

```python
sd = spinchain.Spin_Chain(["S=1/2", "S=1/2"])
J, B, Bt = 15.0, 25.0, 4.0     # exchange, Zeeman energy, tip field, in GHz
hd = J*(sd.Sx[0]*sd.Sx[1] + sd.Sy[0]*sd.Sy[1] + sd.Sz[0]*sd.Sz[1])
hd = hd + B*(sd.Sz[0] + sd.Sz[1]) + Bt*sd.Sz[0]
sd.set_hamiltonian(h*hd)
f = np.linspace(5., 45., 4001)
f, dI, poles = sd.get_esr_spectrum(f, tip=tip, T=1.2, f_ref=B, rabi=1e-3,
                                   return_poles=True)
```

`poles` holds, for each driven transition, its `frequency`, the matrix
element `weight` $|V_{nm}|^2$, the `form_factor` $D_{nn}-D_{mm}$, the
`population_difference` $p_m-p_n$ and the `width` (FWHM in GHz). Here
five transitions appear. The one between the singlet-like and the
$T_0$-like state, at $\sqrt{J^2+B_{\rm tip}^2}=15.52$ GHz, has the largest
weight (0.149) but a small form factor (0.21): the two states have almost
the same spin on the site under the tip. At zero tip field its form factor
is exactly zero, and the line is seen only through the homodyne part, as
an antisymmetric line, which is what Bae et al. observed
(arXiv:1807.04322). The run is in
`examples/kondo/esr_dimer_singlet_triplet`.

Two levels that are exactly degenerate and coupled by $V$, $D$ or the
relaxation raise an error, since their populations depend on the basis
inside the level. This happens at zero field.

## ED and DMRG

`mode="ED"` (the default) diagonalizes $H$, fully up to 2000 states and
for the lowest states above. `mode="DMRG"` reads the spectrum from the
`nstates` lowest DMRG eigenstates, with every matrix element taken between
them. On a 6-site chain with an exchange of 30 GHz:

```python
n = 6
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
hh = sum(30.0*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
         for i in range(n-1))
hh = hh + 20.0*sum(sc.Sz) + 2.0*sc.Sz[0] + 0.5*sc.Sx[0]   # field and tip field
sc.set_hamiltonian(h*hh)
sc.maxm, sc.nsweeps = 64, 20
f = np.linspace(15., 30., 1501)
kw = dict(tip=tip, T=0.05, T2=2., f_ref=20., rabi=1e-4, p_cut=1e-4)
f, ied = sc.get_esr_spectrum(f, mode="ED", **kw)
f, idm = sc.get_esr_spectrum(f, mode="DMRG", nstates=24, **kw)
```

The two agree to $7\times10^{-8}$ of the peak. The states needed are the
populated ones, up to $kT\ln(1/p_{\rm cut})$ above the ground state, and
those reached from them by a line in the window. When the `nstates`
lowest states do not reach that energy the call raises and asks for more
states; here 11 states are needed, and `nstates=8` is not enough. So DMRG is the method for long
chains whose window holds few states, a gapped chain or a cluster in a
field, and not for a gapless chain with an exchange of the order of $hf$,
where the number of states grows with the length. The line positions come
from DMRG energies, and a line is only placed within its width when the
energy variance of its two states is small enough; the call warns
otherwise. The keywords are:

- `mode`: `"ED"` or `"DMRG"`.
- `nstates`: the number of DMRG eigenstates.
- `p_cut`: the relative Boltzmann weight below which a state is empty.
  It is also roughly the relative accuracy of the result, since the
  states above are left out of the relaxation too. A larger value needs
  fewer states.
- `method="floquet"`: the driven Lindblad equation in the same
  eigenbasis, without the secular and rotating-wave approximations, on
  `mode="ED"` with up to 64 states. It is the reference for close lines
  and agrees with the default to $2\times10^{-4}$ of each line or better
  on the dimer above.

`examples/kondo/esr_single_spin_bloch` and
`examples/kondo/esr_spin_chain_dmrg_VS_ED` run the cases of this page.
