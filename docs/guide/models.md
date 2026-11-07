# Models and Hilbert spaces

Every model is a chain of $n$ local Hilbert spaces,
$\mathcal H=\bigotimes_{i=1}^n\mathcal H_i$. The chain class fixes the local
space and gives lists of local operators, indexed by site. The Hamiltonian
and the observables are then built from these operators ([Hamiltonians and operators](operators.md)).

| Chain class | Local Hilbert space | Operators |
|------------|----|-----|
| `spinchain.Spin_Chain(sites)` | spin $S=\tfrac12$ to $3$, one label per site (`"S=1/2"`, `"S=1"`, ...) | `Sx`, `Sy`, `Sz` |
| `fermionchain.Fermionic_Chain(n)` | spinless fermion, empty or occupied | `C`, `Cdag`, `N`, Jordan-Wigner string `F` |
| `fermionchain.Majorana_Chain(n)` | $n$ Majorana modes, built from $n/2$ fermions | `G`, with $\gamma_{2j}=(c_j+c_j^\dagger)/\sqrt2$, $\gamma_{2j+1}=i(c_j-c_j^\dagger)/\sqrt2$ |
| `fermionchain.Spinful_Fermionic_Chain(n)` | spin-1/2 fermion: $0,\uparrow,\downarrow,\uparrow\downarrow$ | `Cup`, `Cdagup`, `Cdn`, `Cdagdn`, `Nup`, `Ndn`, `Ntot`, `Sx`, `Sy`, `Sz`, `Delta` |
| `fermionchain.Spinful_Fermionic_Chain_Native(n)` | the same, one 4-state site per orbital | the same as above |
| `bosonchain.Bosonic_Chain(n, maxnb=...)` | boson with $0$ to $\text{maxnb}_i-1$ particles; `maxnb` is a list with one entry per site, 4 everywhere by default | `A`, `Adag`, `N`, projectors `D[i][k]` $=\lvert k\rangle\langle k\rvert$ |
| `bosonchain.SpinBoson_Chain(sites)` | boson and spin sites mixed, labels `"B"`, `"B6"`, `"S=1/2"`, ... | `A`, `Adag`, `N`, `D` on boson sites; `Sx`, `Sy`, `Sz` on spin sites |
| `parafermionchain.Parafermionic_Chain(n, Z=3)` | $\mathbb Z_N$ clock model, $N=2,3,4$ | clock and shift `Sig`, `Tau`, parafermions `Chi`, `Psi` |
| `mixedchain.Mixed_Spin_Fermion_Chain(sites)` | spin sites and spinful-fermion sites mixed, labels `"1/2"`, `"1"`, ..., `"F"` | `Sx`, `Sy`, `Sz` everywhere; `Cup`, `Cdagup`, `Cdn`, `Cdagdn`, `Nup`, `Ndn`, `Ntot`, `Delta` on fermion sites |

The spin operators of a spinful fermion are
$S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$,
$S^+_i=c^\dagger_{i\uparrow}c_{i\downarrow}$, and the onsite pairing is
$\Delta_i=\tfrac12c_{i\uparrow}c_{i\downarrow}$. In the mixed chains, an
operator that does not exist on a site (for example `Cup` on a spin site)
is zero.

`Spinful_Fermionic_Chain` stores each orbital as two spinless sites, and
`Spinful_Fermionic_Chain_Native` as one 4-state site. The physics and the
sign conventions are the same. In practice the native class is not faster,
so `Spinful_Fermionic_Chain` is the default choice.
`Mixed_Spin_Fermion_Chain` is for models with a true local moment next to
conduction electrons, such as the Kondo lattice
(`examples/fermion_models/mixed_spin_fermion_chain`).

`Mixed_Spin_Fermion_Chain` also runs with `mode="ED"`. The Hilbert space is
the product of the sites: a spin $S$ has $2S+1$ states, and a fermion site
has the four states $0,\uparrow,\downarrow,\uparrow\downarrow$. The fermion
operators carry the Jordan-Wigner string over the fermion sites only,

$$c_{j\sigma}=\Big(\prod_{\text{fermion modes } k<j\sigma}(1-2n_k)\Big)\,a_{j\sigma},$$

meaning that a spin site between two fermion sites gives no sign, and the
spins commute with the electrons. DMRG uses the same convention, so the two
solvers can be compared directly:

```python
from dmrgpy import mixedchain
mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"])
h = 0.5*mc.SS(0, 1) + 0.35*mc.SS(1, 2) + 0.15*mc.Ntot[0] - 0.2*mc.Ntot[2]
hop = 0.4*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
h = h + hop + hop.get_dagger() + 0.13*sum(mc.Sz)
mc.set_hamiltonian(h)
print(mc.gs_energy(mode="ED"), mc.gs_energy(mode="DMRG"))
```

Both numbers are $-0.992144$. The conserved sectors are the total $2S^z$,
of spins and electrons together, and the electron number. On a chain with
both spin and fermion sites they work only with `mode="ED"`: set
`mc.mode = "ED"` before `set_conserved_sector`. For the same reason
`submode="SECTOR"` is not available there; use `"KPM"`, `"TD"` or `"CVM"`.
The example `examples/fermion_models/mixed_spin_fermion_chain_VS_ED`
compares the magnetization profile of DMRG and ED on a Kondo chain.

!!! note "Backends"
    Some classes do not run on every backend:

    - `Mixed_Spin_Fermion_Chain`: `itensor_version=3` and `"python"` only.
    - `Spinful_Fermionic_Chain_Native`: not on `"julia_live"`.
    - `SpinBoson_Chain`: no `mode="ED"`. Use `Bosonic_Chain` for a purely
      bosonic model.
    - `Mixed_Spin_Fermion_Chain` with both spin and fermion sites: a
      conserved sector only with `mode="ED"`, and no `submode="SECTOR"`.
    - `Bosonic_Chain` and `SpinBoson_Chain` with a boson site that does not
      have 4 levels: not on `itensor_version=2`, it raises an error.

Every chain has the methods of this guide. A few methods exist only on
the class whose physics they need:

- `Spin_Chain`: `SS(i, j)`, `get_magnetization()`, `get_kondo_spectrum()`
  ([STM and Kondo spectra](kondo.md)), `get_full_SS_correlator()`.
- The fermionic chains: `get_density()`, `get_density_fluctuation()`,
  `get_charge_gap()`, the correlation matrix and its entropies ([Entanglement](entanglement.md)),
  and `get_fermionic_parity()` on their states ([Topological invariants](topology.md)).
- The spinful fermionic chains: `get_magnetization()`,
  `get_onsite_pairing()`, `get_spectral_function()` with `spin=`,
  `get_kondo_spectrum()` on the electron spin of a site
  ([STM and Kondo spectra](kondo.md#hubbard-chains)).
- `Bosonic_Chain`: `get_density()`, `get_density_fluctuation()`.

Every constructor takes `itensor_version=` and the chain settings as
keywords. A keyword has the same effect as an assignment after
construction, so

```python
from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["S=1/2"]*20, maxm=60, nsweeps=20)
```

is the same as setting `sc.maxm = 60` and `sc.nsweeps = 20` afterwards. A
keyword that is not a setting (for example `hubbard=` or `hamiltonian=`)
raises an error.
