# Audit, 2026-11-21: seven-lens hole hunt of the work since the Kondo hunt

On 2026-11-21 we ran a hole hunt over the 158 commits from `bbb92e4` to
`86ec994`, the work that landed after the 2026-11 hunt of the Kondo
machinery: the ESR-STM spectrum, the Kondo speedups and the NRG route, the
spin chains on metals, the `julia_live` parity ports, the new `"python"`
defaults, the imaginary-time and non-Hermitian routes and the C++ session
changes. The tree was clean at `86ec994`, and both compiled extensions were
current (`make -n pybind` had nothing to do in `mpscpp2` and `mpscpp3`).
This hunt differs from the earlier ones in one way. By the author's rule of
2026-11-15 no agent ran code: seven hunters read the code and wrote a repro
script for each candidate, the lead ran every script, one at a time on one
core (except that the first two ESR repros overlapped the Julia one for a
few seconds), and an independent reviewer briefed to refute the candidate then read
the script and its verbatim output. Where a reviewer needed a measurement,
it wrote a counter-repro and the rule that reads a verdict off it, and the
lead ran it and applied the rule. A candidate that a reviewer found while
reviewing another one went to a reviewer of its own. Twenty agents ran, all
on Opus: seven hunters and thirteen reviewers. The hunters returned 17
candidates and the reviewers 8 more. One was refuted and is not reproduced
here: the cost of the spin average of `ysr.get_stm_spectrum` on a
multiplet, which `docs/numerical_changes.md` records as the price of the
2026-11-20 change. The other 24 are below as 22 findings, since two hunters
found the sector correlation-matrix defect and two the per-bond cutoff one.
The last reviewer found one more candidate, which no agent was left to
review; it is under "Left open, and new leads" with its measured numbers.
No agent ran a dmrgpy calculation, a test or a build; the transcripts show
three small Python helpers that do not import dmrgpy (one printed where
`juliapkg` is installed, one added `flush=True` to its own scripts, which
were then rewritten from the text it returned, and one read its brief's
JSON file).

This file is the evidence, not a task list. Each finding keeps its repro
and gains a `**Status**` line when it is acted on.

## The seven lenses

| Lens | Brief |
|---|---|
| `esr` | `get_esr_spectrum`: `esrstm.py`, `esrtk/`, the ED and DMRG routes and both `submode="KRYLOV"` routes |
| `kondo-routes` | the Kondo speedups after the 2026-11 hunt: `kondospectrumtk/`, `kondocloud.py`, each faster route against the route before it |
| `nrg-ysr-metals` | the NRG route and the spin chains on metals, leads and superconductors: `nrgmps.py`, `nrgspectral.py`, `ysr.py`, `ysrtk/` |
| `julia` | `julia_live` after the parity night and ITensorMPS 0.4.2: `mpsjulialive/` |
| `pyitensor` | the new `"python"` defaults: the density-matrix sum and apply, the one-sweep correlation matrix, the per-bond MPO cutoff, the excited states |
| `solvers` | imaginary time, thermal states, `applyinverse` and the non-Hermitian ground state: `imaginarytime.py`, `thermal.py`, `nhlevels.py`, `nhdmrg.py`, `mpsalgebra.py` |
| `v3-cpp` | the C++ session code changed in the window, judged through its Python calls: `mpscpp3/`, `mpscpp2/` |

## Scope

The window `bbb92e4..86ec994`, and the behaviour next to it that the window
made reachable. Out of scope by construction, as in every record: the
vendored ITensor and TDVP folders, the vendored copy
`src/dmrgpy/mpsjulialive/ITensorNHDMRG/` (how dmrgpy calls it is in scope),
the legacy quirks of `.claude/rules/cpp-backend.md`, the open
`docs/known_issue_*.md` items, what earlier records hold in their findings,
"Left open" and "New leads" sections, the open findings of
`docs/audit_2026_11_15_hole_hunt.md`, and what `ROADMAP.md` marks as absent.
`julia_live` was in scope for reading on every lens and was run only by the
`julia` lens, in one script, since each Julia process pays about a minute of
JIT. Every repro ran as

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=<repo>/src python3 <scratch>/<script>
```

from its own folder, under a 6 GB memory cap, and every script also puts
`<repo>/src` first on `sys.path` before it imports dmrgpy. The scripts and
outputs below are spliced from disk. The repository and scratch paths are
written `<repo>` and `<scratch>`, and RuntimeWarning lines and ITensorMPS's
two-line bond-dimension warning are removed from the outputs; nothing else
is. A script whose name starts with the reviewer group (`reviews/gNN/`) is a
reviewer's counter-repro or a reviewer-found candidate.

## Findings

### 1. `get_esr_spectrum(submode="KRYLOV")` drops every Sz_tot sector the drive does not reach, so with the tip along the field the S-T0 line comes out -5.0 times the right value on a dimer, and exactly 0 at zero tip field

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `esr` &middot; candidate `esr-1`

**Where**: `src/dmrgpy/esrtk/spectrum.py:616-635 (ed_krylov_spectrum: Q = populated eigenstates only, then the block Krylov space of V P)`; `src/dmrgpy/esrtk/spectrum.py:458-473 (mps_krylov_spectrum: DMRG states only up to ecut = max(epop, window[0]))`; `src/dmrgpy/esrtk/spectrum.py:481-521 (TASK finals and partners from P V|m> only)`; `src/dmrgpy/esrstm.py:129-148 (the dispatch to both Krylov builders)`; `tests/test_esr_stm.py:225-241 (the test chain has no conserved quantity)`

The Krylov routes build their basis from the populated eigenstates and the
Krylov space of `V|m>`, m populated. The Korringa operators `S^a_k` that
relax the upper state of a line are not in that construction. When H and V
both conserve Sz_tot, which is the case for a tip moment along the field and
a tip field along z, and every populated state lies in one sector, the basis
is exactly that sector, deterministically, since `_orth` drops the 1e-16
leakage. The T- and T+ partners of an S-T0 line are then missing from the
rate matrix, from the widths and from the DC readout through `D_T-T- = -1/2`.
On the guide's dimer (J = 30 GHz, B = 5 GHz, tip field 3 GHz along z,
T = 0.03 K, `tip=[0,0,1]`) the ED Krylov route keeps 2 of the 4 states, the
DC peak is -5.0 times the right value, with the sign flipped, and the width
is 6.8% low. At zero tip field it returns a DC part of exactly 0, where full
ED gives -4.7e-5, the S-T0 DC part that the ESR rules in
`.claude/rules/kondo.md` call the model. On 4 spins the ED route is 5.5 of
the peak off. The DMRG route (`"python"` and v3) takes the states below the
window from DMRG and is 5.5 of the peak off when T- lies inside the window.
When T- lies below the window it is 1.3e-4 of the peak off, above the 5e-5
its own test asks, and the counter-repro shows that residual is exactly the
missing T+ between the window and `emax`. No route warns. With a transverse
tip field in H, the case of the test chain, the ED route agrees with full ED
to 2.3e-12. The counter-repro also shows the whole ED Krylov error is the
missing sectors: full ED cut by hand to the Sz=0 sector equals the Krylov
route to 2.5e-12.

**Expected**: the Krylov routes agreeing with the eigenstate route to their documented accuracy whatever the tip direction, or a refusal when the basis cannot reach the relaxation partners.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/esr/esr_1_krylov_partners.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: submode="KRYLOV" builds its basis from the populated
# eigenstates and the Krylov space of V|m> only. When V conserves a quantum
# number of H that the Korringa operators S^a_k do not (tip along the
# field: V = D = Sz_tip, H conserving Sz_tot), the states of the other Sz
# sectors are never in the basis, so the T- state below the upper state of
# the S-T0 line is missing as a relaxation partner and the line is wrong.
import warnings
from dmrgpy import spinchain
from dmrgpy.algebra import algebra
from dmrgpy.esrtk import spectrum, lineshape

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
T, T1, T2 = 0.03, 100., 50.


def chain(n, J, B, Bt, iv=None):
    kw = {} if iv is None else dict(itensor_version=iv)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    h = 0
    for i in range(n-1):
        h = h + H*J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                     + sc.Sz[i]*sc.Sz[i+1])
    h = h + H*B*sum(sc.Sz) + H*Bt*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 1e-14
    return sc


def levels(sc):
    """energies in GHz above the ground state, <Sz_tot>, |<n|Sz_0|0>|^2"""
    ed = sc.get_ED_obj()
    Hm = np.asarray(algebra.todense(ed.get_hamiltonian()))
    Szt = np.asarray(algebra.todense(ed.get_operator(sum(sc.Sz))))
    Sz0 = np.asarray(algebra.todense(ed.get_operator(sc.Sz[0])))
    w, v = np.linalg.eigh(Hm)
    e = (w - w[0])/H
    mz = np.real(np.einsum("in,ij,jn->n", v.conj(), Szt, v))
    v0 = np.abs(v.conj().T @ Sz0 @ v[:, 0])**2
    return e, mz, v0


def kw(tip):
    return dict(tip=tip, site=0, T=T, T1=T1, T2=T2, f_ref=20., rabi=1e-4,
                return_poles=True)


def compare(label, sc, f, tip, f0, extra_routes=()):
    _, a, pa = sc.get_esr_spectrum(f, **kw(tip))
    pk = np.max(np.abs(a))
    k = np.argmax(np.abs(a))
    ia = np.argmin(np.abs(pa["frequency"] - f0))
    print("%s f0_GHz = %.6f" % (label, f0))
    print("%s ED_full dI_peak = %.6e at f = %.6f width_MHz = %.4f"
          % (label, a[k], f[k], 1e3*pa["width"][ia]))
    routes = [("ED_KRYLOV", sc, dict(submode="KRYLOV"))] + list(extra_routes)
    for name, sc2, ex in routes:
        try:
            with warnings.catch_warnings(record=True) as wl:
                warnings.simplefilter("always")
                _, b, pb = sc2.get_esr_spectrum(f, **dict(kw(tip), **ex))
            ib = np.argmin(np.abs(pb["frequency"] - f0))
            print("%s %s dI_at_same_f = %.6e width_MHz = %.4f"
                  % (label, name, b[k], 1e3*pb["width"][ib]))
            print("%s %s max|diff|/peak = %.3e  ratio_at_peak = %.4f  "
                  "nwarn = %d" % (label, name, np.max(np.abs(a - b))/pk,
                                   b[k]/a[k], len(wl)))
        except Exception as ex_:
            print("%s %s EXCEPTION = %s: %s" % (label, name,
                                                type(ex_).__name__,
                                                str(ex_)[:200]))


def kept_states(label, sc, f, V):
    """number of states each ED route keeps (internal builders, same
    emax/epop/window as get_esr_spectrum)"""
    kT = KB*T
    p_cut = 1e-6
    margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
    emax = H*(np.max(f) + margin) + kT*np.log(1./p_cut)
    epop = kT*np.log(1./p_cut)
    window = (H*(np.min(f) - margin), H*(np.max(f) + margin))
    ritz_tol = 0.1*H/(2.*np.pi*T2)
    s1 = spectrum.ed_spectrum(sc, V, V, emax)
    s2 = spectrum.ed_krylov_spectrum(sc, V, V, emax, epop, window, kT, 100,
                                     ritz_tol)
    print("%s states_full = %d energies_GHz = %s" % (label, s1.dim,
                                                     np.round(s1.e/H, 4)))
    print("%s states_krylov = %d energies_GHz = %s" % (label, s2.dim,
                                                       np.round(s2.e/H, 4)))


# --- the dimer: J = 30 GHz, field 5 GHz, tip field Bt along z ------------
for Bt in (3.0, 0.0):
    sc = chain(2, 30., 5., Bt)
    e, mz, v0 = levels(sc)
    print("dimer_Bt%g levels_GHz = %s Sz_tot = %s" % (Bt, np.round(e, 4),
                                                       np.round(mz, 3)))
    f0 = e[np.argmax(v0[1:]) + 1]
    f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
    kept_states("dimer_Bt%g" % Bt, sc, f, sc.Sz[0])
    compare("dimer_Bt%g_tipz" % Bt, sc, f, [0., 0., 1.], f0)

# control: a tilted tip reaches every sector, so the Krylov basis is complete
sc = chain(2, 30., 5., 3.0)
e, mz, v0 = levels(sc)
f0 = e[np.argmax(v0[1:]) + 1]
f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
compare("dimer_Bt3_tilted_control", sc, f, [0.5, 0., 0.866], f0)

# --- 4 spins, same physics, with the DMRG Krylov route ("python") -------
sc = chain(4, 30., 5., 3.0)
e, mz, v0 = levels(sc)
print("chain4 levels_GHz = %s" % np.round(e[:8], 4))
print("chain4 Sz_tot = %s" % np.round(mz[:8], 3))
f0 = e[np.argmax(v0[1:]) + 1]
mz0 = mz[0]
etm = np.min(e[np.abs(mz - (mz0 - 1.)) < 1e-6])
print("chain4 line_S_T0_GHz = %.6f lowest_state_of_Sz0-1_GHz = %.6f"
      % (f0, etm))
# W1: the T- state lies below the window (the DMRG route takes it from DMRG)
f1 = np.linspace(f0 - 1.0, f0 + 0.5, 1501)
# W2: the T- state lies inside the window (the DMRG route needs it from TASK)
f2 = np.linspace(etm - 1.0, f0 + 0.5, 3001)
for lab, f in (("chain4_W1", f1), ("chain4_W2", f2)):
    kept_states(lab, sc, f, sc.Sz[0])
    scd = chain(4, 30., 5., 3.0, iv="python")
    compare(lab, sc, f, [0., 0., 1.], f0,
            extra_routes=[("DMRG_KRYLOV_python", scd,
                           dict(mode="DMRG", submode="KRYLOV", nstates=6))])
```

Observed:

```
dimer_Bt3 levels_GHz = [ 0.     23.5748 30.1496 36.5748] Sz_tot = [ 0. -1.  0.  1.]
dimer_Bt3 states_full = 4 energies_GHz = [ 0.     23.5748 30.1496 36.5748]
dimer_Bt3 states_krylov = 2 energies_GHz = [ 0.     30.1496]
dimer_Bt3_tipz f0_GHz = 30.149627
dimer_Bt3_tipz ED_full dI_peak = -1.159115e-05 at f = 30.149627 width_MHz = 7.6734
dimer_Bt3_tipz ED_KRYLOV dI_at_same_f = 5.793588e-05 width_MHz = 7.1501
dimer_Bt3_tipz ED_KRYLOV max|diff|/peak = 5.998e+00  ratio_at_peak = -4.9983  nwarn = 0
dimer_Bt0 levels_GHz = [ 0. 25. 30. 35.] Sz_tot = [ 0. -1.  0.  1.]
dimer_Bt0 states_full = 4 energies_GHz = [ 0. 25. 30. 35.]
dimer_Bt0 states_krylov = 2 energies_GHz = [ 0. 30.]
dimer_Bt0_tipz f0_GHz = 30.000000
dimer_Bt0_tipz ED_full dI_peak = -4.744056e-05 at f = 30.000000 width_MHz = 7.5601
dimer_Bt0_tipz ED_KRYLOV dI_at_same_f = 0.000000e+00 width_MHz = 7.1620
dimer_Bt0_tipz ED_KRYLOV max|diff|/peak = 1.000e+00  ratio_at_peak = -0.0000  nwarn = 0
dimer_Bt3_tilted_control f0_GHz = 30.149627
dimer_Bt3_tilted_control ED_full dI_peak = -7.530416e-06 at f = 30.149627 width_MHz = 7.6734
dimer_Bt3_tilted_control ED_KRYLOV dI_at_same_f = -7.530416e-06 width_MHz = 7.6734
dimer_Bt3_tilted_control ED_KRYLOV max|diff|/peak = 2.288e-12  ratio_at_peak = 1.0000  nwarn = 0
chain4 levels_GHz = [ 0.     13.5637 19.913  26.1218 35.2694 40.9996 46.7657 52.0461]
chain4 Sz_tot = [ 0. -1.  0.  1. -1.  0.  1.  0.]
chain4 line_S_T0_GHz = 19.912964 lowest_state_of_Sz0-1_GHz = 13.563745
chain4_W1 states_full = 4 energies_GHz = [ 0.     13.5637 19.913  26.1218]
chain4_W1 states_krylov = 2 energies_GHz = [ 0.    19.913]
chain4_W1 f0_GHz = 19.912964
chain4_W1 ED_full dI_peak = -1.256738e-05 at f = 19.912964 width_MHz = 6.8360
chain4_W1 ED_KRYLOV dI_at_same_f = 5.692995e-05 width_MHz = 6.4538
chain4_W1 ED_KRYLOV max|diff|/peak = 5.530e+00  ratio_at_peak = -4.5300  nwarn = 0
chain4_W1 DMRG_KRYLOV_python dI_at_same_f = -1.256899e-05 width_MHz = 6.8359
chain4_W1 DMRG_KRYLOV_python max|diff|/peak = 1.279e-04  ratio_at_peak = 1.0001  nwarn = 0
chain4_W2 states_full = 4 energies_GHz = [ 0.     13.5637 19.913  26.1218]
chain4_W2 states_krylov = 2 energies_GHz = [ 0.    19.913]
chain4_W2 f0_GHz = 19.912964
chain4_W2 ED_full dI_peak = -1.249158e-05 at f = 19.913231 width_MHz = 6.8360
chain4_W2 ED_KRYLOV dI_at_same_f = 5.654504e-05 width_MHz = 6.4538
chain4_W2 ED_KRYLOV max|diff|/peak = 5.527e+00  ratio_at_peak = -4.5267  nwarn = 0
chain4_W2 DMRG_KRYLOV_python dI_at_same_f = 5.654504e-05 width_MHz = 6.4538
chain4_W2 DMRG_KRYLOV_python max|diff|/peak = 5.527e+00  ratio_at_peak = -4.5267  nwarn = 0
[exit 0, 8 s]
```

Reviewer's counter-repro, the error as the missing sectors, and v3:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g00/esr_1_counter.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for esr-1. Question: is the whole error of submode="KRYLOV"
# the states of the other Sz_tot sectors missing from its basis, or is
# something else in the Krylov route wrong too? We take the FULL ED
# eigenbasis spectrum (spectrum.ed_spectrum, the reference), cut it down
# by hand to the states each Krylov route keeps, run the same rate
# equation (lineshape.Lines) on it, and compare with the route.
# (1) dimer, tip along z: ED cut to the Sz_tot=0 states against ED KRYLOV.
# (2) 4 spins, window W1 (T- below the window): ED cut to {GS, T-, T0}
#     (T+ dropped, the state between the window and emax in another
#     sector) against the DMRG KRYLOV route on "python".
import warnings
from dmrgpy import spinchain, esrstm
from dmrgpy.algebra import algebra
from dmrgpy.esrtk import spectrum, lineshape

np.random.seed(0)
H = esrstm.H_EV_GHZ
KB = esrstm.KB_EV_K
T, T1, T2, F_REF, RABI = 0.03, 100., 50., 20., 1e-4
P_CUT = esrstm.P_CUT
TIP = [0., 0., 1.]


def chain(n, J, B, Bt, iv=None):
    kw = {} if iv is None else dict(itensor_version=iv)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    h = 0
    for i in range(n-1):
        h = h + H*J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                     + sc.Sz[i]*sc.Sz[i+1])
    h = h + H*B*sum(sc.Sz) + H*Bt*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 1e-14
    return sc


def levels(sc):
    ed = sc.get_ED_obj()
    Hm = np.asarray(algebra.todense(ed.get_hamiltonian()))
    Szt = np.asarray(algebra.todense(ed.get_operator(sum(sc.Sz))))
    Sz0 = np.asarray(algebra.todense(ed.get_operator(sc.Sz[0])))
    w, v = np.linalg.eigh(Hm)
    e = (w - w[0])/H
    mz = np.real(np.einsum("in,ij,jn->n", v.conj(), Szt, v))
    v0 = np.abs(v.conj().T @ Sz0 @ v[:, 0])**2
    return e, mz, v0


def kw():
    return dict(tip=TIP, site=0, T=T, T1=T1, T2=T2, f_ref=F_REF, rabi=RABI)


def cut_curve(sc, f, keep_mask_fn):
    """the rate equation on the full ED spectrum cut to a subset of states"""
    kT = KB*T
    margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
    emax = H*(np.max(f) + margin) + kT*np.log(1./P_CUT)
    V = esrstm._tip_operator(sc, 0, TIP)
    full = spectrum.ed_spectrum(sc, V, V, emax)
    idx = np.flatnonzero(keep_mask_fn(full))
    sub = spectrum.ESRSpectrum(full.e[idx] - full.e[idx][0],
                               full.V[np.ix_(idx, idx)],
                               full.D[np.ix_(idx, idx)],
                               full.K[np.ix_(idx, idx)], route="ED")
    lines = lineshape.Lines(sub, f, H, kT, T1, T2, RABI, P_CUT, F_REF)
    dc, hom = lines.signal(f, 0., 0.)
    return full, idx, dc + hom


def report(label, ref, x, name):
    pk = np.max(np.abs(ref))
    print("%s %s max|diff|/peak_of_full = %.3e" % (label, name,
                                                     np.max(np.abs(ref - x))/pk))


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # (1) the dimer at tip field 3 GHz along z
    sc = chain(2, 30., 5., 3.0)
    e, mz, v0 = levels(sc)
    f0 = e[np.argmax(v0[1:]) + 1]
    f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
    _, a_full = sc.get_esr_spectrum(f, **kw())
    _, a_kry = sc.get_esr_spectrum(f, submode="KRYLOV", **kw())
    full, idx, a_cut = cut_curve(sc, f,
                                 lambda s: np.abs(mz[:s.dim]) < 1e-6)
    print("dimer kept_full_GHz = %s Sz_tot = %s" % (np.round(full.e/H, 4),
                                                     np.round(mz[:full.dim], 3)))
    print("dimer cut_states_GHz = %s" % np.round(full.e[idx]/H, 4))
    report("dimer", a_full, a_kry, "ED_KRYLOV_vs_full")
    report("dimer", a_full, a_cut, "ED_cut_to_Sz0_vs_full")
    print("dimer ED_KRYLOV_vs_ED_cut_to_Sz0 max|diff|/peak_of_cut = %.3e"
          % (np.max(np.abs(a_kry - a_cut))/np.max(np.abs(a_cut))))
    # (1b) the same cut with T- put back (only T+ missing)
    _, _, a_cutm = cut_curve(sc, f, lambda s: mz[:s.dim] < 0.5)
    report("dimer", a_full, a_cutm, "ED_cut_to_S_Tminus_T0_vs_full")

    # (2) 4 spins, window W1: T- below the window
    sc = chain(4, 30., 5., 3.0)
    e, mz, v0 = levels(sc)
    f0 = e[np.argmax(v0[1:]) + 1]
    f1 = np.linspace(f0 - 1.0, f0 + 0.5, 1501)
    _, b_full = sc.get_esr_spectrum(f1, **kw())
    full, idx, b_cut = cut_curve(sc, f1, lambda s: mz[:s.dim] < 0.5)
    print("chain4_W1 kept_full_GHz = %s Sz_tot = %s"
          % (np.round(full.e/H, 4), np.round(mz[:full.dim], 3)))
    print("chain4_W1 cut_states_GHz = %s" % np.round(full.e[idx]/H, 4))
    scd = chain(4, 30., 5., 3.0, iv="python")
    _, b_dmrg = scd.get_esr_spectrum(f1, mode="DMRG", submode="KRYLOV",
                                     nstates=6, **kw())
    report("chain4_W1", b_full, b_dmrg, "DMRG_KRYLOV_vs_full")
    report("chain4_W1", b_full, b_cut, "ED_cut_drop_Tplus_vs_full")
    print("chain4_W1 DMRG_KRYLOV_vs_ED_cut_drop_Tplus max|diff|/peak_of_cut "
          "= %.3e" % (np.max(np.abs(b_dmrg - b_cut))/np.max(np.abs(b_cut))))

    # (3) the same two windows on itensor_version=3, when it is compiled
    from dmrgpy import cppext
    etm = np.min(e[np.abs(mz - (mz[0] - 1.)) < 1e-6])
    f2 = np.linspace(etm - 1.0, f0 + 0.5, 3001)
    print("v3_available = %d" % int(cppext.available(3)))
    if cppext.available(3):
        for lab, ff in (("chain4_W1", f1), ("chain4_W2", f2)):
            _, c_full = sc.get_esr_spectrum(ff, **kw())
            sc3 = chain(4, 30., 5., 3.0, iv=3)
            try:
                _, c3 = sc3.get_esr_spectrum(ff, mode="DMRG", submode="KRYLOV",
                                             nstates=6, **kw())
                report(lab, c_full, c3, "DMRG_KRYLOV_v3_vs_full")
            except Exception as ex:
                print("%s DMRG_KRYLOV_v3 EXCEPTION = %s: %s"
                      % (lab, type(ex).__name__, str(ex)[:200]))
```

Observed:

```
dimer kept_full_GHz = [ 0.     23.5748 30.1496 36.5748] Sz_tot = [ 0. -1.  0.  1.]
dimer cut_states_GHz = [ 0.     30.1496]
dimer ED_KRYLOV_vs_full max|diff|/peak_of_full = 5.998e+00
dimer ED_cut_to_Sz0_vs_full max|diff|/peak_of_full = 5.998e+00
dimer ED_KRYLOV_vs_ED_cut_to_Sz0 max|diff|/peak_of_cut = 2.464e-12
dimer ED_cut_to_S_Tminus_T0_vs_full max|diff|/peak_of_full = 1.157e-04
chain4_W1 kept_full_GHz = [ 0.     13.5637 19.913  26.1218] Sz_tot = [ 0. -1.  0.  1.]
chain4_W1 cut_states_GHz = [ 0.     13.5637 19.913 ]
chain4_W1 DMRG_KRYLOV_vs_full max|diff|/peak_of_full = 1.279e-04
chain4_W1 ED_cut_drop_Tplus_vs_full max|diff|/peak_of_full = 1.279e-04
chain4_W1 DMRG_KRYLOV_vs_ED_cut_drop_Tplus max|diff|/peak_of_cut = 9.114e-13
v3_available = 1
chain4_W1 DMRG_KRYLOV_v3_vs_full max|diff|/peak_of_full = 1.279e-04
chain4_W2 DMRG_KRYLOV_v3_vs_full max|diff|/peak_of_full = 5.527e+00
[exit 0, 7 s]
```

**Reviewer (CONFIRMED)**: Nothing documents a symmetry restriction (the guide, the docstring,
development notes 4.8l and the plan's Krylov section), and the plan's own
first Krylov bullet asks for every state below the window as an
eigenstate, which the ED route does not take. The `2577c21` review recorded
the ED route as agreeing "with or without the tip field", measured with the
tilted V of the test chain. Struck: "the DC part changes by tens of %", since
what was measured is a sign flip and a factor of about 5; and the 4-spin W1
DMRG row as a passing control, since 1.3e-4 is above the route's own test
tolerance on a chain where DMRG is exact. HIGH: T = 0.03 K is the regime the
guide recommends the route for, and a tip moment along the tip field is the
physical case.

**Suggested fix**: Taking every eigenstate up to the top of the window removes the reason the
route exists, and a raise keyed on `[H,Q] = [V,Q] = 0` under-covers (the
reviewer-found extension below). On ED: take every eigenstate up to
`max(epop, window[0])`, as the DMRG route and the plan already require, and
add a second block Krylov stage seeded with `S^a_k|n>` for the final states
and the states below the window, then re-diagonalize. On both routes a cheap
guard: for each needed final state, measure the weight of `S^a_k|n>` outside
the basis and raise when part of it lies below `emax`. The regression test
is the guide's dimer with `tip=[0,0,1]` at 0.03 K on ED and DMRG.

### 2. Extension of the Krylov finding above: with a tilted tip the ED `submode="KRYLOV"` route still misses an empty level two Sz steps away, and a line comes out 1.6e-2 to 2.7e-2 of its peak off

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `esr` &middot; candidate found by the reviewer of `esr-1`

**Where**: `src/dmrgpy/esrtk/spectrum.py:616-635 (basis: populated eigenstates and the Krylov space of V P under H only)`; `src/dmrgpy/esrtk/spectrum.py:673-680 (K inside span(U) only)`; `src/dmrgpy/esrstm.py:129-134 (ED Krylov dispatch)`; `docs/esr_stm_plan.md:78-81 (every state below the window is needed as an eigenstate)`

When H conserves Sz_tot (field and tip field along z) and the tip is
tilted, V has an Sx part and reaches the sectors m0 and m0+-1 from the
populated sector m0, and H keeps each of them. An empty level in m0-2, which
the final state of a window line in m0-1 reaches only through `S^-_k`, is
then never in the basis. On 4 spins (J = 30 GHz, B = 8 to 9 GHz, a 3 GHz tip
field along z, `tip=[0.5, 0, 0.866]`, T = 0.03 K) the line at 53 to 54 GHz
into Sz = -1 comes out 1.6e-2, 2.1e-2 and 2.7e-2 of its peak off full ED,
with no warning. The empty Sz = -2 level sits 0.4 to 1.4 GHz below the
line's final state, with a Korringa weight 0.5 to it. A guard keyed on
`[V, Sz_tot] = 0`, the fallback the Krylov finding above first suggested,
would not catch it.

**Expected**: the eigenstate route of the same model, which the Krylov route is tested against to 1e-8.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g00/esr_1b_tilted_tip_two_steps.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis (a generalization of esr-1, ED route): with a TILTED tip the
# drive reaches Sz_tot = m0 and m0 +- 1 from a populated state in m0, but
# when H conserves Sz_tot (field and tip field along z) a final state in
# m0 - 1 relaxes through S^-_k into m0 - 2, which neither V nor H ever
# reach from m0. ed_krylov_spectrum takes as eigenstates only the
# populated ones, so an empty m0 - 2 level below the line's final state is
# missing as a relaxation partner, and the line is wrong, silently.
# 4 spins 1/2, J = 30 GHz (open chain: lowest triplet 0.66 J, highest
# triplet 2.07 J, quintet 2.37 J above the singlet; the 3 GHz tip field
# along z shifts them by about 1 GHz): a field B between about 7.3 and
# 9.9 GHz puts the Sz=-2 quintet below the highest Sz=-1 triplet while the
# lowest Sz=-1 state stays above epop = kT ln(1/p_cut) = 8.64 GHz, empty.
import warnings
from dmrgpy import spinchain, esrstm
from dmrgpy.algebra import algebra
from dmrgpy.esrtk import spectrum, lineshape

np.random.seed(0)
H = esrstm.H_EV_GHZ
KB = esrstm.KB_EV_K
T, T1, T2, F_REF, RABI = 0.03, 100., 50., 20., 1e-4
P_CUT = esrstm.P_CUT
TIP = [0.5, 0., 0.866]


def chain(n, J, B, Bt, Btx=0.):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + H*J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                     + sc.Sz[i]*sc.Sz[i+1])
    h = h + H*B*sum(sc.Sz) + H*Bt*sc.Sz[0]
    if Btx != 0.:
        h = h + H*Btx*sc.Sx[0]
    sc.set_hamiltonian(h)
    return sc


def levels(sc, V):
    ed = sc.get_ED_obj()
    Hm = np.asarray(algebra.todense(ed.get_hamiltonian()))
    Szt = np.asarray(algebra.todense(ed.get_operator(sum(sc.Sz))))
    Vm = np.asarray(algebra.todense(ed.get_operator(V)))
    w, v = np.linalg.eigh(Hm)
    e = (w - w[0])/H
    mz = np.real(np.einsum("in,ij,jn->n", v.conj(), Szt, v))
    vV = np.abs(v.conj().T @ Vm @ v[:, 0])**2
    return e, mz, vV


kw = dict(tip=TIP, site=0, T=T, T1=T1, T2=T2, f_ref=F_REF, rabi=RABI)
kT = KB*T
epop = kT*np.log(1./P_CUT)
margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
for label, B, Btx in (("B8", 8.0, 0.), ("B8.5", 8.5, 0.),
                      ("B9", 9.0, 0.), ("B8.5_control_Btx0.5", 8.5, 0.5)):
    sc = chain(4, 30., B, 3.0, Btx)
    V = esrstm._tip_operator(sc, 0, TIP)
    e, mz, vV = levels(sc, V)
    print("%s levels_GHz = %s" % (label, np.round(e[:12], 3)))
    print("%s Sz_tot = %s" % (label, np.round(mz[:12], 3)))
    popset = set(np.round(mz[e*H <= epop]).astype(int))
    reach = set()
    for m in popset:
        reach |= {m - 1, m, m + 1}
    print("%s populated_sectors = %s reached_sectors = %s"
          % (label, sorted(popset), sorted(reach)))
    i2 = np.flatnonzero(np.abs(mz + 2.) < 1e-6)
    if len(i2) == 0 or Btx != 0.:
        cand = np.flatnonzero(vV > 1e-6*np.max(vV[1:]))
        cand = cand[cand > 0]
    else:
        e2 = np.min(e[i2])
        cand = np.flatnonzero((np.abs(mz + 1.) < 1e-6) & (e > e2)
                              & (vV > 1e-6*np.max(vV[1:])))
    if len(cand) == 0:
        print("%s decisive = 0 (no driven Sz=-1 level above the Sz=-2 level)"
              % label)
        continue
    n = cand[np.argmax(vV[cand])]
    f0 = e[n]
    f = np.linspace(f0 - 0.15, f0 + 0.15, 601)
    emax = H*(np.max(f) + margin) + kT*np.log(1./P_CUT)
    full = spectrum.ed_spectrum(sc, V, V, emax)
    kept = np.arange(full.dim)
    unreached = [m for m in kept
                 if int(np.round(mz[m])) not in reach]
    kun = max([full.K[m, n] for m in unreached], default=0.)
    dec = int(len(unreached) > 0 and kun > 1e-6 and Btx == 0.)
    print("%s line_GHz = %.4f (Sz %.0f, |V_n0|^2 = %.3e) emax_GHz = %.3f"
          % (label, f0, mz[n], vV[n], emax/H))
    print("%s unreached_levels_below_emax_GHz = %s K_to_line_final = %.3e "
          "decisive = %d" % (label, np.round(e[unreached], 3), kun, dec))
    try:
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            _, a = sc.get_esr_spectrum(f, **kw)
            _, b = sc.get_esr_spectrum(f, submode="KRYLOV", **kw)
        pk = np.max(np.abs(a))
        k = np.argmax(np.abs(a))
        print("%s ED_full dI_peak = %.6e ED_KRYLOV dI_at_same_f = %.6e"
              % (label, a[k], b[k]))
        print("%s ED_KRYLOV max|diff|/peak = %.3e nwarn = %d"
              % (label, np.max(np.abs(a - b))/pk, len(wl)))
    except Exception as ex:
        print("%s EXCEPTION = %s: %s" % (label, type(ex).__name__,
                                         str(ex)[:200]))
```

Observed:

```
B8 levels_GHz = [ 0.    10.564 19.913 29.122 32.269 41.    49.766 52.046 53.574 53.994
 62.237 62.467]
B8 Sz_tot = [ 0. -1.  0.  1. -1.  0.  1.  0. -2. -1.  0. -1.]
B8 populated_sectors = [np.int64(0)] reached_sectors = [np.int64(-1), np.int64(0), np.int64(1)]
B8 line_GHz = 53.9942 (Sz -1, |V_n0|^2 = 3.099e-04) emax_GHz = 62.971
B8 unreached_levels_below_emax_GHz = [53.574] K_to_line_final = 5.000e-01 decisive = 1
B8 ED_full dI_peak = -2.133241e-07 ED_KRYLOV dI_at_same_f = -2.099486e-07
B8 ED_KRYLOV max|diff|/peak = 1.582e-02 nwarn = 0
B8.5 levels_GHz = [ 0.    10.064 19.913 29.622 31.769 41.    50.266 52.046 52.574 53.494
 61.967 62.237]
B8.5 Sz_tot = [ 0. -1.  0.  1. -1.  0.  1.  0. -2. -1. -1.  0.]
B8.5 populated_sectors = [np.int64(0)] reached_sectors = [np.int64(-1), np.int64(0), np.int64(1)]
B8.5 line_GHz = 53.4942 (Sz -1, |V_n0|^2 = 3.099e-04) emax_GHz = 62.471
B8.5 unreached_levels_below_emax_GHz = [52.574] K_to_line_final = 5.000e-01 decisive = 1
B8.5 ED_full dI_peak = -2.301978e-07 ED_KRYLOV dI_at_same_f = -2.253204e-07
B8.5 ED_KRYLOV max|diff|/peak = 2.119e-02 nwarn = 0
B9 levels_GHz = [ 0.     9.564 19.913 30.122 31.269 41.    50.766 51.574 52.046 52.994
 61.467 62.237]
B9 Sz_tot = [ 0. -1.  0.  1. -1.  0.  1. -2.  0. -1. -1.  0.]
B9 populated_sectors = [np.int64(0)] reached_sectors = [np.int64(-1), np.int64(0), np.int64(1)]
B9 line_GHz = 52.9942 (Sz -1, |V_n0|^2 = 3.099e-04) emax_GHz = 61.971
B9 unreached_levels_below_emax_GHz = [51.574] K_to_line_final = 5.000e-01 decisive = 1
B9 ED_full dI_peak = -2.486842e-07 ED_KRYLOV dI_at_same_f = -2.419755e-07
B9 ED_KRYLOV max|diff|/peak = 2.699e-02 nwarn = 0
B8.5_control_Btx0.5 levels_GHz = [ 0.    10.066 19.916 29.627 31.77  41.001 50.263 52.045 52.551 53.532
 61.941 62.271]
B8.5_control_Btx0.5 Sz_tot = [-0.    -1.    -0.     1.    -1.    -0.     0.997 -0.005 -1.973 -1.018
 -0.91  -0.09 ]
B8.5_control_Btx0.5 populated_sectors = [np.int64(0)] reached_sectors = [np.int64(-1), np.int64(0), np.int64(1)]
B8.5_control_Btx0.5 line_GHz = 19.9160 (Sz -0, |V_n0|^2 = 1.220e-01) emax_GHz = 28.893
B8.5_control_Btx0.5 unreached_levels_below_emax_GHz = [] K_to_line_final = 0.000e+00 decisive = 0
B8.5_control_Btx0.5 ED_full dI_peak = -4.033938e-05 ED_KRYLOV dI_at_same_f = -4.033938e-05
B8.5_control_Btx0.5 ED_KRYLOV max|diff|/peak = 2.933e-12 nwarn = 0
[exit 0, 2 s]
```

**Reviewer (CONFIRMED)**: The Hilbert space has 16 states and the sectors the drive reaches hold 14 of
them, which a block Krylov space at `nkry=100` exhausts; so the only states
missing are the two of Sz = +-2, and one of them lies below `emax`. The whole
difference from full ED is that level. Struck: the control at a transverse
tip field as evidence for this line, since it reads another line (19.9 GHz)
with no unreached level; the attribution rests on the count instead.
MEDIUM, below its parent, since it needs an empty level two Sz steps away
below a line's final state.

**Suggested fix**: The main fix of the Krylov finding covers it, its fallback does not: close
the basis under `S^a_k` applied to the converged final states of the window
lines, with one more block Krylov space orthogonalized against the basis
before the Rayleigh-Ritz step. Pin it with this repro's B8.5 case next to
the dimer with `tip=[0,0,1]`.

### 3. `get_esr_spectrum` keeps states only by temperature, so below about 70 mK at 20 GHz a second photon's state is dropped and a saturated spin-1 line comes out 50% low (4.8e-2 at the default drive)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `esr` &middot; candidate found by the reviewer of `esr-2`

**Where**: `src/dmrgpy/esrstm.py:107-124 (emax from the window, its margin and the thermal populations only)`; `src/dmrgpy/esrtk/spectrum.py:175 (keep = e <= emax) and 289-296 (DMRG route, same cut)`; `src/dmrgpy/esrtk/lineshape.py:138-147 (population grows along driven pairs among the kept states only)`; `src/dmrgpy/esrtk/floquet.py:28-33 (the Floquet reference uses the same basis)`

`emax` counts the window, its margin and the thermally populated states,
and its comment says `p_cut` is then the accuracy of the DC part. A state
that a second window photon reaches from the upper state of a driven line
is not counted, so it is dropped whenever `kT ln(1/p_cut)` is below about
h f, which is about 70 mK at 20 GHz and the default `p_cut`; the 50 mK of
the tests and examples is included. On a spin 1 in a field at 0 and 50 mK,
whose two rungs both lie in the window, the rung `|0> -> |+1>` is missing
(one pair kept against two). On the isotropic ladder the saturated line is
0.380 against 0.757 at `rabi=0.05` (50% low), 0.237 against 0.370 at
`rabi=5e-3` (36% low), and 4.8e-2 of the peak off at the default
`rabi=1e-3` (the counter-repro). The code's own Floquet reference uses the
same truncated basis and is low by the same amount, so it cannot catch it.
With an anisotropy of 5 MHz (rungs 10 MHz apart) keeping `|+1>` changes the
Floquet answer by 13% to 43% across the line at `rabi=5e-3` and by 2.0e-2
at `rabi=1e-3`. At 0.2 K the thermal term keeps `|+1>` and the two windows
agree exactly.

**Expected**: the answer of the window that keeps every state the drive can populate, which a hand-written three-level rate equation reproduces to 1e-13.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g01/esr_cascade_spin1.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: emax = h(f_max + margin) + kT ln(1/p_cut) counts only the
# thermally populated states. Under a saturating drive the upper state of a
# window line is populated too, and a line out of it whose upper state lies
# two window photons above the ground state is never kept when
# kT ln(1/p_cut) < h f (14.4 GHz at 50 mK and p_cut = 1e-6). On a spin 1 in a
# field the rung |0> -> |+1> is then missing. Window C is window A plus one
# point at 2 F0 + 0.5 GHz, which keeps |+1>. rate3/rate2 are an independent
# rate equation (Korringa eq. 42 as in esrtk/rates.py, written out here) with
# three and with the lowest two levels; flqA/flqC are the code's own Floquet
# reference (driven Lindblad, no secular, no rotating-wave approximation) on
# the same two windows. "iso": H = h F0 Sz, both rungs at F0. "aniso":
# + h d Sz^2, rungs at F0 -+ d, 2d = 10 MHz above the Rabi rate 4.6 MHz.
from dmrgpy import spinchain

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
F0, TH, T1, T2 = 20.0, 0.7, 100., 50.


def rate_eq(f, T, rabi, nlev, d):
    ms = np.array([-1., 0., 1.])[:nlev]
    E = H*(F0*ms + d*ms**2)
    E = E - np.min(E)
    kT = KB*T
    n = nlev
    Sp = np.zeros((n, n))
    for i in range(n - 1): # <m+1|S+|m> of a spin 1
        Sp[i + 1, i] = np.sqrt(2. - ms[i]*(ms[i] + 1.))
    Sx = 0.5*(Sp + Sp.T)
    Sz = np.diag(ms)
    K = 0.5*(Sp**2 + (Sp**2).T) # |Sx|^2 + |Sy|^2, Sz has no off-diagonal
    eps0 = H*F0
    c = 1. if kT == 0. else 1./np.tanh(eps0/(2.*kT))
    g = 2./(T1*eps0*c)
    def phi(x):
        if kT == 0.: return max(x, 0.)
        return x/(-np.expm1(-x/kT))
    Gam = np.zeros((n, n)) # Gam[m, k]: rate from k to m
    for m in range(n):
        for k in range(n):
            if m != k: Gam[m, k] = g*K[m, k]*phi(E[k] - E[m])
    out = Gam.sum(axis=0)
    gphi = 1./T2 - 0.5/T1
    V = np.sin(TH)*Sx + np.cos(TH)*Sz
    Om = 4.*np.pi*rabi
    w = 2.*np.pi*f
    R = Gam.copy()
    R[np.diag_indices(n)] = -out
    for l in range(n - 1):
        u = l + 1
        gam = 0.5*(out[u] + out[l]) + gphi
        wul = 2.*np.pi*(E[u] - E[l])/H
        W = 0.5*Om**2*V[u, l]**2*gam/((w - wul)**2 + gam**2)
        R[u, l] += W; R[l, u] += W; R[u, u] -= W; R[l, l] -= W
    A = np.vstack([R, np.ones(n)])
    b = np.zeros(n + 1); b[-1] = 1.
    p = np.linalg.lstsq(A, b, rcond=None)[0]
    if kT == 0.:
        p0 = np.zeros(n); p0[0] = 1.
    else:
        p0 = np.exp(-E/kT); p0 /= p0.sum()
    return float(np.sum((p - p0)*np.diag(V)))


fA = np.linspace(19.9, 20.1, 41)
fC = np.append(fA, 2*F0 + 0.5)
cases = [("iso", 0.0, T, rabi) for T in (0.0, 0.05, 0.2)
         for rabi in (1e-4, 5e-3, 5e-2)]
cases += [("aniso", 0.005, T, 5e-3) for T in (0.0, 0.05)]
for name, d, T, rabi in cases:
    sc = spinchain.Spin_Chain(["S=1"])
    Hm = H*F0*sc.Sz[0]
    if d != 0.:
        Hm = Hm + H*d*sc.Sz[0]*sc.Sz[0]
    sc.set_hamiltonian(Hm)
    tip = [np.sin(TH), 0., np.cos(TH)]
    kw = dict(tip=tip, T=T, T1=T1, T2=T2, f_ref=F0, rabi=rabi)
    _, a, pa = sc.get_esr_spectrum(fA, return_poles=True, **kw)
    _, c, pc = sc.get_esr_spectrum(fC, return_poles=True, **kw)
    _, qa = sc.get_esr_spectrum(fA, method="floquet", nharm=3, **kw)
    _, qc = sc.get_esr_spectrum(fC, method="floquet", nharm=3, **kw)
    lab = "%s_T%g_rabi%g" % (name, T, rabi)
    print("%s pairs_A = %d pairs_C = %d weights_C = %s" % (
        lab, len(pa["frequency"]), len(pc["frequency"]),
        np.array2string(pc["weight"], precision=6)))
    fs = [F0] if d == 0. else [F0 - d, F0, F0 + d]
    for fx in fs:
        i = int(np.argmin(np.abs(fA - fx)))
        r3 = rate_eq(fA[i], T, rabi, 3, d)
        r2 = rate_eq(fA[i], T, rabi, 2, d)
        lf = "%s_f%.3f" % (lab, fA[i])
        print("%s A = %.6e C = %.6e rate2 = %.6e rate3 = %.6e flqA = %.6e "
              "flqC = %.6e" % (lf, a[i], c[i], r2, r3, qa[i], qc[i]))
        print("%s (C-A)/C = %.3e |A-rate2|/|rate2| = %.3e |C-rate3|/|rate3| "
              "= %.3e (flqC-flqA)/flqC = %.3e |flqC-C|/|C| = %.3e" % (
                  lf, (c[i] - a[i])/c[i], abs(a[i] - r2)/abs(r2),
                  abs(c[i] - r3)/abs(r3), (qc[i] - qa[i])/qc[i],
                  abs(qc[i] - c[i])/abs(c[i])))
```

Observed:

```
iso_T0_rabi0.0001 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0_rabi0.0001_f20.000 A = 2.504623e-04 C = 2.505795e-04 rate2 = 2.504623e-04 rate3 = 2.505795e-04 flqA = 2.504623e-04 flqC = 2.505935e-04
iso_T0_rabi0.0001_f20.000 (C-A)/C = 4.677e-04 |A-rate2|/|rate2| = 3.876e-11 |C-rate3|/|rate3| = 1.020e-10 (flqC-flqA)/flqC = 5.238e-04 |flqC-C|/|C| = 5.607e-05
iso_T0_rabi0.005 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0_rabi0.005_f20.000 A = 2.374778e-01 C = 3.703562e-01 rate2 = 2.374778e-01 rate3 = 3.703562e-01 flqA = 2.374778e-01 flqC = 3.736114e-01
iso_T0_rabi0.005_f20.000 (C-A)/C = 3.588e-01 |A-rate2|/|rate2| = 5.961e-15 |C-rate3|/|rate3| = 1.319e-14 (flqC-flqA)/flqC = 3.644e-01 |flqC-C|/|C| = 8.789e-03
iso_T0_rabi0.05 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0_rabi0.05_f20.000 A = 3.801012e-01 C = 7.574350e-01 rate2 = 3.801012e-01 rate3 = 7.574350e-01 flqA = 3.801015e-01 flqC = 7.574381e-01
iso_T0_rabi0.05_f20.000 (C-A)/C = 4.982e-01 |A-rate2|/|rate2| = 0.000e+00 |C-rate3|/|rate3| = 2.932e-16 (flqC-flqA)/flqC = 4.982e-01 |flqC-C|/|C| = 4.071e-06
iso_T0.05_rabi0.0001 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.05_rabi0.0001_f20.000 A = 2.504623e-04 C = 2.505795e-04 rate2 = 2.504623e-04 rate3 = 2.505795e-04 flqA = 2.504623e-04 flqC = 2.505935e-04
iso_T0.05_rabi0.0001_f20.000 (C-A)/C = 4.677e-04 |A-rate2|/|rate2| = 1.144e-11 |C-rate3|/|rate3| = 6.006e-11 (flqC-flqA)/flqC = 5.238e-04 |flqC-C|/|C| = 5.607e-05
iso_T0.05_rabi0.005 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.05_rabi0.005_f20.000 A = 2.374778e-01 C = 3.703562e-01 rate2 = 2.374778e-01 rate3 = 3.703562e-01 flqA = 2.374778e-01 flqC = 3.736114e-01
iso_T0.05_rabi0.005_f20.000 (C-A)/C = 3.588e-01 |A-rate2|/|rate2| = 4.675e-15 |C-rate3|/|rate3| = 1.799e-15 (flqC-flqA)/flqC = 3.644e-01 |flqC-C|/|C| = 8.789e-03
iso_T0.05_rabi0.05 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.05_rabi0.05_f20.000 A = 3.801012e-01 C = 7.574350e-01 rate2 = 3.801012e-01 rate3 = 7.574350e-01 flqA = 3.801014e-01 flqC = 7.574381e-01
iso_T0.05_rabi0.05_f20.000 (C-A)/C = 4.982e-01 |A-rate2|/|rate2| = 4.381e-16 |C-rate3|/|rate3| = 4.397e-16 (flqC-flqA)/flqC = 4.982e-01 |flqC-C|/|C| = 4.071e-06
iso_T0.2_rabi0.0001 pairs_A = 2 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.2_rabi0.0001_f20.000 A = 2.526117e-04 C = 2.526117e-04 rate2 = 2.463704e-04 rate3 = 2.526117e-04 flqA = 2.526253e-04 flqC = 2.526253e-04
iso_T0.2_rabi0.0001_f20.000 (C-A)/C = 0.000e+00 |A-rate2|/|rate2| = 2.533e-02 |C-rate3|/|rate3| = 6.194e-11 (flqC-flqA)/flqC = 0.000e+00 |flqC-C|/|C| = 5.370e-05
iso_T0.2_rabi0.005 pairs_A = 2 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.2_rabi0.005_f20.000 A = 3.693852e-01 C = 3.693852e-01 rate2 = 2.335981e-01 rate3 = 3.693852e-01 flqA = 3.724991e-01 flqC = 3.724991e-01
iso_T0.2_rabi0.005_f20.000 (C-A)/C = 0.000e+00 |A-rate2|/|rate2| = 5.813e-01 |C-rate3|/|rate3| = 9.317e-15 (flqC-flqA)/flqC = 0.000e+00 |flqC-C|/|C| = 8.430e-03
iso_T0.2_rabi0.05 pairs_A = 2 pairs_C = 2 weights_C = [0.207508 0.207508]
iso_T0.2_rabi0.05_f20.000 A = 7.512055e-01 C = 7.512055e-01 rate2 = 3.738913e-01 rate3 = 7.512055e-01 flqA = 7.512084e-01 flqC = 7.512084e-01
iso_T0.2_rabi0.05_f20.000 (C-A)/C = 0.000e+00 |A-rate2|/|rate2| = 1.009e+00 |C-rate3|/|rate3| = 1.478e-16 (flqC-flqA)/flqC = 0.000e+00 |flqC-C|/|C| = 3.913e-06
aniso_T0_rabi0.005 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
aniso_T0_rabi0.005_f19.995 A = 2.375093e-01 C = 2.845573e-01 rate2 = 2.375093e-01 rate3 = 2.845573e-01 flqA = 2.375093e-01 flqC = 2.726398e-01
aniso_T0_rabi0.005_f19.995 (C-A)/C = 1.653e-01 |A-rate2|/|rate2| = 6.544e-15 |C-rate3|/|rate3| = 1.284e-13 (flqC-flqA)/flqC = 1.289e-01 |flqC-C|/|C| = 4.188e-02
aniso_T0_rabi0.005_f20.000 A = 1.485821e-01 C = 2.112560e-01 rate2 = 1.485821e-01 rate3 = 2.112560e-01 flqA = 1.485581e-01 flqC = 2.598401e-01
aniso_T0_rabi0.005_f20.000 (C-A)/C = 2.967e-01 |A-rate2|/|rate2| = 1.364e-14 |C-rate3|/|rate3| = 2.856e-13 (flqC-flqA)/flqC = 4.283e-01 |flqC-C|/|C| = 2.300e-01
aniso_T0_rabi0.005_f20.005 A = 6.997870e-02 C = 1.176617e-01 rate2 = 6.997870e-02 rate3 = 1.176617e-01 flqA = 6.995327e-02 flqC = 1.136966e-01
aniso_T0_rabi0.005_f20.005 (C-A)/C = 4.053e-01 |A-rate2|/|rate2| = 4.978e-14 |C-rate3|/|rate3| = 8.657e-14 (flqC-flqA)/flqC = 3.847e-01 |flqC-C|/|C| = 3.370e-02
aniso_T0.05_rabi0.005 pairs_A = 1 pairs_C = 2 weights_C = [0.207508 0.207508]
aniso_T0.05_rabi0.005_f19.995 A = 2.375093e-01 C = 2.845573e-01 rate2 = 2.375093e-01 rate3 = 2.845573e-01 flqA = 2.375093e-01 flqC = 2.726398e-01
aniso_T0.05_rabi0.005_f19.995 (C-A)/C = 1.653e-01 |A-rate2|/|rate2| = 2.688e-15 |C-rate3|/|rate3| = 1.407e-13 (flqC-flqA)/flqC = 1.289e-01 |flqC-C|/|C| = 4.188e-02
aniso_T0.05_rabi0.005_f20.000 A = 1.485821e-01 C = 2.112560e-01 rate2 = 1.485821e-01 rate3 = 2.112560e-01 flqA = 1.485581e-01 flqC = 2.598401e-01
aniso_T0.05_rabi0.005_f20.000 (C-A)/C = 2.967e-01 |A-rate2|/|rate2| = 1.027e-14 |C-rate3|/|rate3| = 2.537e-13 (flqC-flqA)/flqC = 4.283e-01 |flqC-C|/|C| = 2.300e-01
aniso_T0.05_rabi0.005_f20.005 A = 6.997870e-02 C = 1.176617e-01 rate2 = 6.997870e-02 rate3 = 1.176617e-01 flqA = 6.995327e-02 flqC = 1.136966e-01
aniso_T0.05_rabi0.005_f20.005 (C-A)/C = 4.053e-01 |A-rate2|/|rate2| = 2.241e-14 |C-rate3|/|rate3| = 1.452e-13 (flqC-flqA)/flqC = 3.847e-01 |flqC-C|/|C| = 3.370e-02
[exit 0, 7 s]
```

Reviewer's counter-repro, the default drive:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g12/g12_g01_default_rabi.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Severity probe for new-g01 (review group 12): the size of the missing second
# rung at the DEFAULT rabi=1e-3, where the candidate measured only 1e-4, 5e-3
# and 5e-2. Spin 1 in a field (iso: both rungs at 20 GHz; aniso: rungs
# 10 MHz apart), T = 0.05 K, default p_cut. Window A = 19.9..20.1 GHz;
# window C = A plus one point at 40.5 GHz, which keeps |+1>. The code's own
# Floquet reference in each window, and the rate route.
from dmrgpy import spinchain

np.random.seed(0)
H = 4.135667696e-6
F0, TH, T1, T2, T = 20.0, 0.7, 100., 50., 0.05
fA = np.linspace(19.9, 20.1, 41)
fC = np.append(fA, 2*F0 + 0.5)
for name, d in (("iso", 0.0), ("aniso", 0.005)):
    for rabi in (1e-3, 2e-3):
        sc = spinchain.Spin_Chain(["S=1"])
        Hm = H*F0*sc.Sz[0]
        if d != 0.:
            Hm = Hm + H*d*sc.Sz[0]*sc.Sz[0]
        sc.set_hamiltonian(Hm)
        kw = dict(tip=[np.sin(TH), 0., np.cos(TH)], T=T, T1=T1, T2=T2,
                  f_ref=F0, rabi=rabi)
        _, a = sc.get_esr_spectrum(fA, **kw)
        _, c = sc.get_esr_spectrum(fC, **kw)
        _, qa = sc.get_esr_spectrum(fA, method="floquet", nharm=3, **kw)
        _, qc = sc.get_esr_spectrum(fC, method="floquet", nharm=3, **kw)
        i = int(np.argmax(np.abs(qc[:len(fA)])))
        print("ROW %s rabi=%g f=%.3f rate (C-A)/C = %.3e floquet (flqC-flqA)/flqC "
              "= %.3e |flqC-C|/|C| = %.3e max_f |flqC-flqA|/max|flqC| = %.3e"
              % (name, rabi, fA[i], (c[i]-a[i])/c[i], (qc[i]-qa[i])/qc[i],
                 abs(qc[i]-c[i])/abs(c[i]),
                 np.max(np.abs(qc[:len(fA)]-qa))/np.max(np.abs(qc[:len(fA)]))),
              flush=True)
print("done")
```

Observed:

```
ROW iso rabi=0.001 f=20.000 rate (C-A)/C = 4.307e-02 floquet (flqC-flqA)/flqC = 4.755e-02 |flqC-C|/|C| = 4.709e-03 max_f |flqC-flqA|/max|flqC| = 4.755e-02
ROW iso rabi=0.002 f=20.000 rate (C-A)/C = 1.386e-01 floquet (flqC-flqA)/flqC = 1.486e-01 |flqC-C|/|C| = 1.172e-02 max_f |flqC-flqA|/max|flqC| = 1.486e-01
ROW aniso rabi=0.001 f=19.995 rate (C-A)/C = 1.073e-02 floquet (flqC-flqA)/flqC = 7.017e-03 |flqC-C|/|C| = 3.744e-03 max_f |flqC-flqA|/max|flqC| = 1.994e-02
ROW aniso rabi=0.002 f=19.995 rate (C-A)/C = 3.941e-02 floquet (flqC-flqA)/flqC = 2.656e-02 |flqC-C|/|C| = 1.320e-02 max_f |flqC-flqA|/max|flqC| = 8.077e-02
done
[exit 0, 3 s]
```

**Reviewer (CONFIRMED, NARROWED)**: The anchor is independent: an exact spin-1 Floquet solution with all three
levels kept, and the hand rate equations, which reproduce the narrow window
with two levels to 1e-14 and the wide one with three to 1e-13. Struck: the
anisotropic rows as inside the rate equation's valid regime, since the rate
route is 4.2e-2 to 0.23 off the full-basis Floquet answer there; the window
dependence on those rows stands on the Floquet columns alone. The
`submode="KRYLOV"` part stands by reading only: for `H = F0 Sz` the Krylov
space of `V|-1>` is `span{|-1>, |0>}` whatever `emax` is, the mechanism of
the two findings above. MEDIUM, since the default drive is reached.

**Suggested fix**: Let the kept set follow the drive. On the eigenbasis routes, after the first
cut, run the population-growth loop of `lineshape` with a saturation
estimate above `p_cut`, extend `emax` by `h(f_max + margin)` above the
highest state that can carry population, and repeat until it stops growing,
one extra pass for a ladder of 2S rungs. The margin finding below changes
the same `emax` line, so both should be one change. On `submode="KRYLOV"`
the basis needs V applied to the final states of the window lines, the
closure the two findings above ask for.

### 4. At T=0 `get_esr_spectrum` keeps only 30 free-spin widths above the window, so a saturated line just above it is dropped with its tail: 5.0e-2 of its peak at the window top at xi=14

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `esr` &middot; candidate `esr-2`

**Where**: `src/dmrgpy/esrstm.py:121-124 (margin = WINDOW_MARGIN/(pi T2), emax = h(f_max + margin))`; `src/dmrgpy/esrstm.py:132,141 (the Krylov routes' window, same margin; read, not run)`; `src/dmrgpy/esrtk/spectrum.py:175 (ED: keep = e <= emax)`; `src/dmrgpy/esrtk/lineshape.py:43-48 (the comment claiming the 3e-4 bound)`; `docs/development.md (4.8l, the same bound)`

The states are kept up to `emax = h(f_max + margin) + kT ln(1/p_cut)`, with
`margin = WINDOW_MARGIN/(pi T2)`, 30 full widths of a free spin. The
comment and the development notes bound the dropped tail by
`1/(4 WINDOW_MARGIN^2) = 3e-4` of a line. That holds only at weak drive. A
line just above the band leaves `(1+xi^2)` times that bound, xi the
saturation parameter, since its tail falls in its power-broadened width. At
T=0, xi=14 and a window ending 0.19 GHz below the line, the window returns
exactly 0 at its top, where the Bloch formula and a window that keeps the
line both give 5.0e-2 of the line's peak (3.0e-2 with `V_RF/V_DC = 0.3`).
At T>0 the band grows by `kT ln(1/p_cut)`, so the error falls fast with T:
8.9e-3 at 1 mK, 7.8e-4 at 5 mK and 5.9e-5 at 20 mK, for T2 = 50 ns. With a
short T2 = 2 ns it stays at 3.5e-3 at 50 mK and xi=14.

**Expected**: the dropped tail below the stated bound at every drive, or a warning when a kept band is narrower than the power-broadened width.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/esr/esr_2_t0_window_saturated.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: the kept states end WINDOW_MARGIN = 30 FREE-spin full widths
# above the window, plus kT ln(1/p_cut), which is 0 at T = 0. A saturated
# line (power-broadened to xi free-spin widths) just above that margin is
# dropped with its whole DC tail and homodyne part, which inside the window
# are a few % (DC) to tens of % (homodyne) of its own peak, far above the
# documented 1/(4 WINDOW_MARGIN^2) = 3e-4. A window that reaches the line
# keeps it, so the value at one frequency depends on where the window ends.
from dmrgpy import spinchain

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
F0, TH, T1, T2 = 20.0, 0.7, 100., 50.


def bloch(f, T, rabi, r=0., phase=0.):
    """one spin 1/2, field along z, tip at TH from it, V = D = m.S"""
    OR = 2*np.pi*rabi*np.sin(TH)
    xi2 = OR**2*T1*T2
    x = 2*np.pi*(f - F0)
    y0 = -np.tanh(H*F0/(2*KB*T)) if T > 0 else -1.
    y = y0*(1 + (x*T2)**2)/(1 + (x*T2)**2 + xi2)
    dc = np.cos(TH)*(y - y0)/2
    rho = -(2*np.pi*rabi)*(np.sin(TH)/2)*y/(x + 1j/T2)
    return dc + r*np.real(np.exp(1j*phase)*rho*np.sin(TH)/2), xi2


sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(H*F0*sc.Sz[0])
tip = [np.sin(TH), 0., np.cos(TH)]
fA = np.linspace(19.70, 19.80, 11)     # window ending 0.2 GHz below the line
fB = np.linspace(19.70, 20.30, 61)     # window reaching past the line
iB = int(np.argmin(np.abs(fB - fA[-1])))
print("free_full_width_MHz = %.4f margin_GHz = %.4f window_A_top_plus_margin"
      "_GHz = %.4f" % (1e3/(np.pi*T2), 30./(np.pi*T2), fA[-1] + 30./(np.pi*T2)))
for T in (0.0, 0.05):
    for rabi in (1e-4, 0.05):
        for r in (0.0, 0.3):
            kw = dict(tip=tip, T=T, T1=T1, T2=T2, f_ref=F0, rabi=rabi,
                      vrf_over_vdc=r, phase=0.)
            _, a = sc.get_esr_spectrum(fA, **kw)
            _, b = sc.get_esr_spectrum(fB, **kw)
            ref, xi2 = bloch(fA[-1], T, rabi, r)
            refB, _ = bloch(fB, T, rabi, r)
            pk = np.max(np.abs(refB))
            lab = "T%g_rabi%g_r%g" % (T, rabi, r)
            print("%s xi = %.3f peak_B = %.6e" % (lab, np.sqrt(xi2), pk))
            print("%s dI_A(19.80) = %.6e dI_B(19.80) = %.6e bloch(19.80) = "
                  "%.6e" % (lab, a[-1], b[iB], ref))
            print("%s (B-A)/peak = %.3e (bloch-A)/peak = %.3e "
                  "max|B-bloch|/peak = %.3e" % (lab, (b[iB] - a[-1])/pk,
                                                 (ref - a[-1])/pk,
                                                 np.max(np.abs(b - refB))/pk))
```

Observed:

```
free_full_width_MHz = 6.3662 margin_GHz = 0.1910 window_A_top_plus_margin_GHz = 19.9910
T0_rabi0.0001_r0 xi = 0.029 peak_B = 3.130266e-04
T0_rabi0.0001_r0 dI_A(19.80) = 0.000000e+00 dI_B(19.80) = 7.933542e-08 bloch(19.80) = 7.933541e-08
T0_rabi0.0001_r0 (B-A)/peak = 2.534e-04 (bloch-A)/peak = 2.534e-04 max|B-bloch|/peak = 8.399e-10
T0_rabi0.0001_r0.3 xi = 0.029 peak_B = 3.130266e-04
T0_rabi0.0001_r0.3 dI_A(19.80) = 0.000000e+00 dI_B(19.80) = 7.933542e-08 bloch(19.80) = -1.547984e-05
T0_rabi0.0001_r0.3 (B-A)/peak = 2.534e-04 (bloch-A)/peak = -4.945e-02 max|B-bloch|/peak = 5.317e-02
T0_rabi0.05_r0 xi = 14.311 peak_B = 3.805629e-01
T0_rabi0.05_r0 dI_A(19.80) = 0.000000e+00 dI_B(19.80) = 1.885591e-02 bloch(19.80) = 1.885591e-02
T0_rabi0.05_r0 (B-A)/peak = 4.955e-02 (bloch-A)/peak = 4.955e-02 max|B-bloch|/peak = 7.220e-14
T0_rabi0.05_r0.3 xi = 14.311 peak_B = 3.805629e-01
T0_rabi0.05_r0.3 dI_A(19.80) = 0.000000e+00 dI_B(19.80) = 1.145991e-02 bloch(19.80) = 1.145991e-02
T0_rabi0.05_r0.3 (B-A)/peak = 3.011e-02 (bloch-A)/peak = 3.011e-02 max|B-bloch|/peak = 7.818e-14
T0.05_rabi0.0001_r0 xi = 0.029 peak_B = 3.130266e-04
T0.05_rabi0.0001_r0 dI_A(19.80) = 7.933542e-08 dI_B(19.80) = 7.933542e-08 bloch(19.80) = 7.933541e-08
T0.05_rabi0.0001_r0 (B-A)/peak = 0.000e+00 (bloch-A)/peak = -5.266e-11 max|B-bloch|/peak = 8.400e-10
T0.05_rabi0.0001_r0.3 xi = 0.029 peak_B = 3.130266e-04
T0.05_rabi0.0001_r0.3 dI_A(19.80) = 7.933542e-08 dI_B(19.80) = 7.933542e-08 bloch(19.80) = -1.547984e-05
T0.05_rabi0.0001_r0.3 (B-A)/peak = 0.000e+00 (bloch-A)/peak = -4.971e-02 max|B-bloch|/peak = 5.317e-02
T0.05_rabi0.05_r0 xi = 14.311 peak_B = 3.805629e-01
T0.05_rabi0.05_r0 dI_A(19.80) = 1.885591e-02 dI_B(19.80) = 1.885591e-02 bloch(19.80) = 1.885591e-02
T0.05_rabi0.05_r0 (B-A)/peak = 0.000e+00 (bloch-A)/peak = -2.270e-15 max|B-bloch|/peak = 7.206e-14
T0.05_rabi0.05_r0.3 xi = 14.311 peak_B = 3.805629e-01
T0.05_rabi0.05_r0.3 dI_A(19.80) = 1.145991e-02 dI_B(19.80) = 1.145991e-02 bloch(19.80) = 1.145991e-02
T0.05_rabi0.05_r0.3 (B-A)/peak = 0.000e+00 (bloch-A)/peak = -1.837e-15 max|B-bloch|/peak = 7.804e-14
[exit 0, 3 s]
```

Reviewer's counter-repro, the size against T:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g01/esr2_low_T_boundary.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for esr-2: how fast does the dropped DC tail of a saturated
# line shrink with T, and does a short T2 bring it back at 50 mK? The line is
# put 2 MHz above the kept band of window A (f_top + margin + kT ln(1/p_cut)/h),
# the worst position, and window B is the same frequencies plus one far point,
# which keeps the line. pred is the Bloch tail (1+xi^2)/(1+(2 pi d T2)^2+xi^2).
from dmrgpy import spinchain, esrstm
from dmrgpy.esrtk import lineshape

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
F0, TH, T1 = 20.0, 0.7, 100.


def bloch_dc(f, T, T2, rabi):
    OR = 2*np.pi*rabi*np.sin(TH)
    xi2 = OR**2*T1*T2
    x = 2*np.pi*(f - F0)
    y0 = -np.tanh(H*F0/(2*KB*T)) if T > 0 else -1.
    y = y0*(1 + (x*T2)**2)/(1 + (x*T2)**2 + xi2)
    return np.cos(TH)*(y - y0)/2, xi2


sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(H*F0*sc.Sz[0])
tip = [np.sin(TH), 0., np.cos(TH)]
rows = [(0.0, 50., 0.05), (1e-3, 50., 0.05), (5e-3, 50., 0.05),
        (2e-2, 50., 0.05), (0.05, 2., 0.25), (0.05, 2., 0.05)]
for T, T2, rabi in rows:
    margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
    thermal = KB*T*np.log(1./esrstm.P_CUT)/H if T > 0 else 0.
    ftop = F0 - 0.002 - margin - thermal
    fA = np.linspace(ftop - 0.1, ftop, 11)
    fB = np.append(fA, F0 + 0.3)
    kw = dict(tip=tip, T=T, T1=T1, T2=T2, f_ref=F0, rabi=rabi)
    _, a = sc.get_esr_spectrum(fA, **kw)
    _, b = sc.get_esr_spectrum(fB, **kw)
    pk, xi2 = bloch_dc(F0, T, T2, rabi)
    pk = abs(pk)
    ref, _ = bloch_dc(ftop, T, T2, rabi)
    u2 = (2*np.pi*(F0 - ftop)*T2)**2
    lab = "T%g_T2%g_rabi%g" % (T, T2, rabi)
    print("%s xi = %.3f f_top = %.5f detuning_GHz = %.5f pred = %.3e"
          % (lab, np.sqrt(xi2), ftop, F0 - ftop, (1 + xi2)/(1 + u2 + xi2)))
    print("%s dI_A(top) = %.6e dI_B(top) = %.6e bloch(top) = %.6e"
          % (lab, a[-1], b[len(fA) - 1], ref))
    print("%s (B-A)/peak = %.3e (bloch-B)/peak = %.3e"
          % (lab, (b[len(fA) - 1] - a[-1])/pk, (ref - b[len(fA) - 1])/pk))
```

Observed:

```
T0_T250_rabi0.05 xi = 14.311 f_top = 19.80701 detuning_GHz = 0.19299 pred = 5.302e-02
T0_T250_rabi0.05 dI_A(top) = 0.000000e+00 dI_B(top) = 2.017747e-02 bloch(top) = 2.017747e-02
T0_T250_rabi0.05 (B-A)/peak = 5.302e-02 (bloch-B)/peak = -2.580e-15
T0.001_T250_rabi0.05 xi = 14.311 f_top = 19.51915 detuning_GHz = 0.48085 pred = 8.938e-03
T0.001_T250_rabi0.05 dI_A(top) = 0.000000e+00 dI_B(top) = 3.401346e-03 bloch(top) = 3.401346e-03
T0.001_T250_rabi0.05 (B-A)/peak = 8.938e-03 (bloch-B)/peak = -2.302e-16
T0.005_T250_rabi0.05 xi = 14.311 f_top = 18.36767 detuning_GHz = 1.63233 pred = 7.820e-04
T0.005_T250_rabi0.05 dI_A(top) = 0.000000e+00 dI_B(top) = 2.975927e-04 bloch(top) = 2.975927e-04
T0.005_T250_rabi0.05 (B-A)/peak = 7.820e-04 (bloch-B)/peak = -1.073e-16
T0.02_T250_rabi0.05 xi = 14.311 f_top = 14.04964 detuning_GHz = 5.95036 pred = 5.889e-05
T0.02_T250_rabi0.05 dI_A(top) = 0.000000e+00 dI_B(top) = 2.241122e-05 bloch(top) = 2.241122e-05
T0.02_T250_rabi0.05 (B-A)/peak = 5.889e-05 (bloch-B)/peak = 9.090e-18
T0.05_T22_rabi0.25 xi = 14.311 f_top = 0.82993 detuning_GHz = 19.17007 pred = 3.534e-03
T0.05_T22_rabi0.25 dI_A(top) = 0.000000e+00 dI_B(top) = 1.344845e-03 bloch(top) = 1.344845e-03
T0.05_T22_rabi0.25 (B-A)/peak = 3.534e-03 (bloch-B)/peak = -3.818e-17
T0.05_T22_rabi0.05 xi = 2.862 f_top = 0.82993 detuning_GHz = 19.17007 pred = 1.584e-04
T0.05_T22_rabi0.05 dI_A(top) = 0.000000e+00 dI_B(top) = 5.397600e-05 bloch(top) = 5.397600e-05
T0.05_T22_rabi0.05 (B-A)/peak = 1.584e-04 (bloch-B)/peak = -1.422e-16
[exit 0, 2 s]
```

**Reviewer (CONFIRMED, NARROWED)**: The homodyne part "up to 37% of the line's peak" was struck, since no line
measured it (the output gives 1.9e-2 of the DC peak), and the 5% holds only
at the top of the window; the tail falls to 2.3e-2 of the peak 0.1 GHz
further in. The counter-repro reproduces each size from the Bloch formula to
better than 1e-14 of the peak. The severity is LOW: it needs T of a few mK or
a very short T2 together with strong saturation and a line just outside
the window.

**Suggested fix**: Count the band kept above the window in the power-broadened width, the same
unit the homodyne taper got in `2577c21`, using the operator norm of V and
calling the factor an estimate, since a slow relaxation path makes xi
larger than the free-spin value. On ED the wider band is cheap and can be
applied directly; on DMRG and the Krylov routes warn when the estimated
broadened band exceeds the kept one. Correct the comment and the
development notes to `(1+xi^2)/(4 WINDOW_MARGIN^2)` of a line.

### 5. `site_weights` of `get_esr_spectrum` rescale the T1 and T2 of the tip spin itself: `[0.5]` on a free spin gives a width of 5.570 MHz for 1/(pi T2) = 6.366 MHz

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `esr` &middot; candidate `esr-3`

**Where**: `src/dmrgpy/esrtk/lineshape.py:109-110 (g = korringa_scale(T1, h f_ref, kT), no tip weight)`; `src/dmrgpy/esrtk/rates.py:39-56 (korringa_scale fixes g for K = 1/2)`; `src/dmrgpy/esrtk/spectrum.py:99-110 (_weights) and the K sums of every route`; `src/dmrgpy/esrstm.py:190-197 (docstring: T1, T2 of a free spin 1/2 on the tip site)`; `docs/guide/esr-stm.md (the T1, T2 and site_weights entries)`

The docstring and the guide define T1 and T2 as those of a free spin 1/2 on
the tip site, and `site_weights` as the weight of each site in the
relaxation. The Korringa scale g is fixed by T1 for a site of weight 1, and
the weights multiply it as absolute factors. So whenever
`site_weights[site] != 1` the tip spin relaxes in `T1/w[site]` and its
coherence time is `1/(1/T2 + (w[site]-1)/(2 T1))`, with no warning. On one
free spin `[0.5]` gives a weak-drive width of 5.570 MHz and `[2.0]` 7.958
MHz, against 1/(pi T2) = 6.366 MHz (0.244 and 0.330 of the Bloch peak off).
On the dimer of the guide `[0.5, 0.5]` moves the spectrum by 1.20 times the
peak of the default, while `[1, 1]` equals the default exactly. All four
routes share it. The counter-repro shows the weights are a pure rescaling:
weights `c*w` at (T1, T2) equal `w` at `T1/c` and
`1/T2' = 1/T2 - 1/(2T1) + c/(2T1)` to roundoff, which is also the
workaround today.

**Expected**: the T1 and T2 passed are those of the tip spin whatever its weight, or a docstring that says T1 refers to a site of weight 1.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/esr/esr_3_site_weights_scale.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: T1 and T2 are documented as those of a free spin 1/2 on the
# tip site, and site_weights as "the weight of each site in the relaxation
# (equal by default)". The Korringa scale g is fixed from T1 assuming the
# tip site has weight 1 (rates.korringa_scale) and is never divided by
# site_weights[site], so equal weights other than 1 (for instance
# normalized ones) change T1 and T2 of the tip spin silently: a free spin
# with site_weights=[0.5] is not the Bloch line of its own T1 and T2, and a
# dimer with [0.5, 0.5] differs from the default, although both are equal.
from dmrgpy import spinchain

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
F0, TH, T, T1, T2 = 20.0, 0.7, 1.2, 100., 50.


def bloch(f, rabi):
    OR = 2*np.pi*rabi*np.sin(TH)
    xi2 = OR**2*T1*T2
    x = 2*np.pi*(f - F0)
    y0 = -np.tanh(H*F0/(2*KB*T))
    y = y0*(1 + (x*T2)**2)/(1 + (x*T2)**2 + xi2)
    return np.cos(TH)*(y - y0)/2


tip = [np.sin(TH), 0., np.cos(TH)]
sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(H*F0*sc.Sz[0])
f = np.linspace(19.95, 20.05, 2001)
rabi = 5e-3
ref = bloch(f, rabi)
pk = np.max(np.abs(ref))
for w in (None, [1.0], [0.5], [2.0]):
    _, dI, poles = sc.get_esr_spectrum(f, tip=tip, T=T, T1=T1, T2=T2,
                                       f_ref=F0, rabi=rabi, site_weights=w,
                                       return_poles=True)
    print("free_spin site_weights=%s width_weak_MHz = %.4f (free spin "
          "1/(pi T2) = %.4f) peak = %.6e max|dI-bloch|/peak = %.3e"
          % (w, 1e3*poles["width"][0], 1e3/(np.pi*T2), np.max(dI),
             np.max(np.abs(dI - ref))/pk))

# dimer with equal weights, 1 and 0.5 on each site
sd = spinchain.Spin_Chain(["S=1/2"]*2)
hd = H*15.*(sd.Sx[0]*sd.Sx[1] + sd.Sy[0]*sd.Sy[1] + sd.Sz[0]*sd.Sz[1])
hd = hd + H*25.*(sd.Sz[0] + sd.Sz[1]) + H*4.*sd.Sz[0]
sd.set_hamiltonian(hd)
f = np.linspace(5., 45., 8001)
out = {}
for w in (None, [1.0, 1.0], [0.5, 0.5]):
    _, dI, poles = sd.get_esr_spectrum(f, tip=tip, T=T, f_ref=25., rabi=1e-3,
                                       site_weights=w, return_poles=True)
    out[str(w)] = dI
    print("dimer site_weights=%s widths_MHz = %s peak = %.6e"
          % (w, np.round(1e3*poles["width"], 4), np.max(np.abs(dI))))
pk = np.max(np.abs(out["None"]))
print("dimer max|[1,1]-None|/peak = %.3e" % (np.max(np.abs(out["[1.0, 1.0]"]
                                                            - out["None"]))/pk))
print("dimer max|[0.5,0.5]-None|/peak = %.3e" % (np.max(np.abs(
    out["[0.5, 0.5]"] - out["None"]))/pk))
```

Observed:

```
free_spin site_weights=None width_weak_MHz = 6.3662 (free spin 1/(pi T2) = 6.3662) peak = 9.761624e-02 max|dI-bloch|/peak = 6.571e-13
free_spin site_weights=[1.0] width_weak_MHz = 6.3662 (free spin 1/(pi T2) = 6.3662) peak = 9.761624e-02 max|dI-bloch|/peak = 6.571e-13
free_spin site_weights=[0.5] width_weak_MHz = 5.5704 (free spin 1/(pi T2) = 6.3662) peak = 1.197078e-01 max|dI-bloch|/peak = 2.441e-01
free_spin site_weights=[2.0] width_weak_MHz = 7.9577 (free spin 1/(pi T2) = 6.3662) peak = 6.542108e-02 max|dI-bloch|/peak = 3.298e-01
dimer site_weights=None widths_MHz = [ 7.4552  8.2977  9.0136  9.1948 10.0373] peak = 1.909163e-03
dimer site_weights=[1.0, 1.0] widths_MHz = [ 7.4552  8.2977  9.0136  9.1948 10.0373] peak = 1.909163e-03
dimer site_weights=[0.5, 0.5] widths_MHz = [6.1149 6.5362 6.8941 6.9847 7.406 ] peak = 4.202449e-03
dimer max|[1,1]-None|/peak = 0.000e+00
dimer max|[0.5,0.5]-None|/peak = 1.201e+00
[exit 0, 3 s]
```

Reviewer's counter-repro, the weights as a rescaling:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g02/esr_3_weights_identity.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for esr-3: if site_weights are absolute multipliers of the
# scale g that T1 fixes (and enter nowhere else), then weights c*w with
# (T1, T2) give exactly the spectrum of weights w with T1' = T1/c and
# 1/T2' = 1/T2 - 1/(2 T1) + c/(2 T1) (same rates, same pure dephasing).
# vrf_over_vdc = 0, so T2 enters only through the kept-state margin, and
# every state of these small systems is kept either way.
from dmrgpy import spinchain

np.random.seed(0)
H = 4.135667696e-6 # eV/GHz
TH, T, T1, T2 = 0.7, 1.2, 100., 50.
tip = [np.sin(TH), 0., np.cos(TH)]


def t2prime(c):
    return 1./(1./T2 - 0.5/T1 + 0.5*c/T1)


# one free spin
sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(H*20.*sc.Sz[0])
f = np.linspace(19.95, 20.05, 2001)
for c in (0.5, 2.0):
    _, a = sc.get_esr_spectrum(f, tip=tip, T=T, T1=T1, T2=T2, f_ref=20.,
                               rabi=5e-3, site_weights=[c])
    _, b = sc.get_esr_spectrum(f, tip=tip, T=T, T1=T1/c, T2=t2prime(c),
                               f_ref=20., rabi=5e-3)
    print("free_spin c=%.1f max|w=c vs (T1/c, T2')|/peak = %.3e"
          % (c, np.max(np.abs(a - b))/np.max(np.abs(b))))

# the dimer of the repro, equal weights 0.5 and unequal weights [1, 0.25]
sd = spinchain.Spin_Chain(["S=1/2"]*2)
hd = H*15.*(sd.Sx[0]*sd.Sx[1] + sd.Sy[0]*sd.Sy[1] + sd.Sz[0]*sd.Sz[1])
hd = hd + H*25.*(sd.Sz[0] + sd.Sz[1]) + H*4.*sd.Sz[0]
sd.set_hamiltonian(hd)
f = np.linspace(5., 45., 8001)
for w0, c in (([1., 1.], 0.5), ([1., 0.25], 2.0)):
    w = [c*x for x in w0]
    _, a = sd.get_esr_spectrum(f, tip=tip, T=T, T1=T1, T2=T2, f_ref=25.,
                               rabi=1e-3, site_weights=w)
    _, b = sd.get_esr_spectrum(f, tip=tip, T=T, T1=T1/c, T2=t2prime(c),
                               f_ref=25., rabi=1e-3, site_weights=w0)
    print("dimer w=%s max|w vs w0=%s at (T1/c, T2')|/peak = %.3e"
          % (w, w0, np.max(np.abs(a - b))/np.max(np.abs(b))))
```

Observed:

```
free_spin c=0.5 max|w=c vs (T1/c, T2')|/peak = 0.000e+00
free_spin c=2.0 max|w=c vs (T1/c, T2')|/peak = 3.182e-16
dimer w=[0.5, 0.5] max|w vs w0=[1.0, 1.0] at (T1/c, T2')|/peak = 0.000e+00
dimer w=[2.0, 0.5] max|w vs w0=[1.0, 0.25] at (T1/c, T2')|/peak = 4.047e-17
[exit 0, 5 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: "equal weights of any size equal the default" as an expected
behaviour (nothing promises it), and the trigger "equal weights other than
1", replaced by `site_weights[site] != 1`, since `[1, 0.3]` leaves the tip
spin as documented. A raise at `w[site] == 0`, as the hunter suggests,
would refuse a legal input (relaxation only through the other sites).

**Suggested fix**: This is the author's pick between two conventions. (a) Normalize to the
tip: divide g by `w[site]` when it is positive, so the tip spin has exactly
the T1 and T2 passed, as an experiment quotes them. This changes numbers
for every call with `w[site] != 1` and needs a ledger entry. (b) Document
it: T1 and T2 are those of a free spin of weight 1, and the weights multiply
each site's Korringa scale. Either way, add a test with `site_weights`,
which no test uses today.

### 6. The METTS collapse-basis guard runs only in `metts_vev` on `"python"` and v3, so `julia_live` and `metts_dynamical_correlator` on every backend freeze in one charge sector: `<Sztot^2>` 1, 4 and 0 for an exact 0.370

`bug` &middot; severity **HIGH** &middot; CONFIRMED, NARROWED &middot; lens `julia` &middot; candidate `julia-1`

**Where**: `src/dmrgpy/vevtk/mettsvev.py:194-195 (the guard, python and v3 only)`; `src/dmrgpy/vevtk/mettsdynamicalcorrelator.py:98-160 (no basis check)`; `src/dmrgpy/mpsjulialive/metts.jl:41-56 (eigen(Hermitian(M)) reads the upper triangle)`; `src/dmrgpy/mpsjulialive/metts.py:33-104`; `src/dmrgpy/mpscpp3/chain_session.h:2437-2508 (the v3 dynamical correlator, no check)`; `src/dmrgpy/pyitensor/metts.py:89-95 (np.linalg.eigh reads the lower triangle)`; `tests/test_audit_2026_10_23.py, tests/test_audit_2026_10_23_leads.py (python and v3 only)`

The guard of the 2026-10-23 hunt (its findings 1 and 10,
`vevtk/mettsvev.py::_check_basis`) refuses a collapse basis METTS cannot
sample with: a non-Hermitian name, or an all-diagonal basis when H conserves
the matching charge. It runs only inside `metts_vev`, and only on
`"python"` and v3. `metts_dynamical_correlator` calls no basis check on any
backend, and on `julia_live` `metts_vev` has none either. An all-diagonal
basis on a charge-conserving H then returns the thermal average inside the
charge sector of the first random sample, and the seed picks the sector. On
`julia_live` `metts_vev`, `("Sz",)` on a 4-site Heisenberg chain at T=0.5
gives `<Sztot^2>` = 1, 4 and 0 for an exact 0.3702, and `("N",)` on a 4-site
hopping chain at T=1 gives `<Ntot^2>` = 9 and 0 for 4.583. On
`metts_dynamical_correlator`, `C(0) = <(3+Sztot)^2>` comes out 16, 4, 4 on
`"python"`, 16, 1 on `julia_live` and 25, 1, 9 on v3, for an exact 9.370. A
non-Hermitian name is not refused either: `julia_live`'s `metts.jl` reads
the upper triangle of M, which is zero for `S-`, and froze at
`<Sztot^2> = 4`; on `"python"` the dynamical correlator froze with `("S+",)`
and scattered around ED with `("S-",)`, since `eigh` reads the other
triangle. The error bar is 0 only for a function of the frozen charge; for
`<H>` it is the scatter inside the sector (-0.8308 +- 3e-6 and
0.7500 +- 0 for an exact -1.1197), so the error bar does not flag the
problem. `CLAUDE.md` and the METTS section of the guide state the refusal
with no backend restriction.

**Expected**: the refusal of the 2026-10-23 guard on every backend and in both METTS calls.

Repro (python and julia_live; its p3 lines are the julia_live half of the sector correlation-matrix finding):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/julia/julia_1_metts_basis_julia.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Candidate julia-1: the METTS collapse-basis guard of the 2026-10-23 hunt
# (vevtk/mettsvev.py::_check_basis) runs for itensor_version "python" and 3
# only, and metts_dynamical_correlator never calls it on any backend. An
# all-diagonal basis on an H that conserves the charge then never leaves
# the charge sector of the first sample: the value is that of one sector,
# chosen by the seed, with stderr 0.
# Part 1: python, metts_dynamical_correlator (no Julia needed).
# Part 2: julia_live, metts_vev.
# Part 3: julia_live, candidate julia-2 (sector correlation matrix, operators=).
# Part 4: julia_live, metts_dynamical_correlator (last, optional).
# The dynamical pairs are (3+Sztot, 3+Sztot) and (1+Ntot, 1+Ntot), so that
# B|phi> is never the zero vector; a stuck chain gives (3+M)^2 or (1+M)^2,
# an exact integer square, with stderr 0. nt=1: only C(0) is measured.
import time

np.random.seed(0)

from dmrgpy import spinchain, fermionchain

T0 = time.time()
NAN = float("nan")


def say(name, value):
    print("%s = %s   [%.0fs]" % (name, value, time.time() - T0), flush=True)


def err(e):
    return "ERROR %s: %s" % (type(e).__name__, str(e)[:160])


def dense(ed, m):
    m = ed.MO2matrix(m)
    if hasattr(m, "todense"):
        m = m.todense()
    return np.asarray(m, dtype=complex)


def thermal(chain, ops, T):
    """exact Boltzmann averages from the dense ED matrices"""
    ed = chain.get_ED_obj()
    H = dense(ed, chain.hamiltonian)
    e, U = np.linalg.eigh(H)
    w = np.exp(-(e - e.min()) / T)
    w = w / w.sum()
    out = []
    for O in ops:
        Od = U.conj().T @ dense(ed, O) @ U
        out.append(float(np.real(np.sum(w * np.diag(Od)))))
    return out


def spin_chain(backend, n=4):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=backend,
                              maxm=16, nsweeps=6, cutoff=1e-12)
    h = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    sc.set_hamiltonian(h)
    return sc, h, sum(sc.Sz)


def fermion_chain(backend, n=4):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=backend, maxm=16,
                                      nsweeps=6, cutoff=1e-12)
    h = sum(-1.0 * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
            for i in range(n - 1))
    h = h + 0.3 * fc.N[0] - 0.2 * fc.N[n - 1] + 0.1 * fc.N[1]
    fc.set_hamiltonian(h)
    return fc, h, sum(fc.N)


TS, TF = 0.5, 1.0

# ---------------------------------------------------------------- Part 1
ex = [NAN, NAN, NAN]     # <Sztot^2>, <H>, <(3+Sztot)^2> at TS
exf = [NAN, NAN, NAN]    # <Ntot^2>, <Ntot>, <(1+Ntot)^2> at TF
sc = fc = None
try:
    sc, h, Sz = spin_chain("python")
    Bs = 3.0 + Sz
    ex = thermal(sc, [Sz * Sz, h, Bs * Bs], TS)
    say("ED spin4 T=0.5 <Sztot^2> <H> <(3+Sztot)^2>",
        "%.6f %.6f %.6f" % tuple(ex))
except Exception as e:
    say("p1 spin setup", err(e))
if sc is not None:
    try:
        sc.metts_vev(Sz * Sz, TS, nsamples=4, nwarmup=1, basis_ops=("Sz",),
                     seed=1)
        say("p1 python metts_vev spin basis=(Sz,)", "RETURNED (no guard)")
    except ValueError as e:
        say("p1 python metts_vev spin basis=(Sz,)", "ValueError (guard fires)")
    except Exception as e:
        say("p1 python metts_vev spin basis=(Sz,)", err(e))
    for seed in (1, 2, 3):
        try:
            ts, m, s = sc.metts_dynamical_correlator(
                (Bs, Bs), TS, nt=1, dt=0.1, nsamples=8, nwarmup=2,
                basis_ops=("Sz",), seed=seed)
            say("p1 python dynMETTS spin basis=(Sz,) seed=%d" % seed,
                "C(0)=<(3+Sztot)^2> %.6f stderr %.2e | ED %.6f | diff %.3e"
                % (m[0].real, s[0], ex[2], m[0].real - ex[2]))
        except Exception as e:
            say("p1 python dynMETTS spin basis=(Sz,) seed=%d" % seed, err(e))
    try:
        ts, m, s = sc.metts_dynamical_correlator((Bs, Bs), TS, nt=1, dt=0.1,
                                                 nsamples=8, nwarmup=2, seed=1)
        say("p1 python dynMETTS spin basis=(Sz,Sx) control seed=1",
            "C(0)=<(3+Sztot)^2> %.6f stderr %.2e | ED %.6f | diff %.3e"
            % (m[0].real, s[0], ex[2], m[0].real - ex[2]))
    except Exception as e:
        say("p1 python dynMETTS spin control", err(e))

try:
    fc, hf, Nt = fermion_chain("python")
    Bf = 1.0 + Nt
    exf = thermal(fc, [Nt * Nt, Nt, Bf * Bf], TF)
    say("ED fermion4 T=1 <Ntot^2> <Ntot> <(1+Ntot)^2>",
        "%.6f %.6f %.6f" % tuple(exf))
except Exception as e:
    say("p1 fermion setup", err(e))
if fc is not None:
    try:
        fc.metts_vev(Nt * Nt, TF, nsamples=4, nwarmup=1, basis_ops=("N",),
                     seed=1)
        say("p1 python metts_vev fermion basis=(N,)", "RETURNED (no guard)")
    except ValueError as e:
        say("p1 python metts_vev fermion basis=(N,)", "ValueError (guard fires)")
    except Exception as e:
        say("p1 python metts_vev fermion basis=(N,)", err(e))
    for seed in (1, 2, 3):
        try:
            ts, m, s = fc.metts_dynamical_correlator(
                (Bf, Bf), TF, nt=1, dt=0.1, nsamples=8, nwarmup=2,
                basis_ops=("N",), seed=seed)
            say("p1 python dynMETTS fermion basis=(N,) seed=%d" % seed,
                "C(0)=<(1+Ntot)^2> %.6f stderr %.2e | ED %.6f | diff %.3e"
                % (m[0].real, s[0], exf[2], m[0].real - exf[2]))
        except Exception as e:
            say("p1 python dynMETTS fermion basis=(N,) seed=%d" % seed, err(e))

# ---------------------------------------------------------------- Part 2
sj = None
try:
    sj, hj, Szj = spin_chain("julia_live")
    from dmrgpy.mpsjulialive.juliasession import Main
    try:
        say("julia ITensorMPS ITensors versions",
            "%s %s" % (Main.seval("string(pkgversion(ITensorMPS))"),
                       Main.seval("string(pkgversion(ITensors))")))
    except Exception as e:
        say("julia versions", err(e))
except Exception as e:
    say("p2 julia spin setup", err(e))
if sj is not None:
    for seed in (1, 2, 3):
        try:
            (m1, s1), (m2, s2) = sj.metts_vev([Szj * Szj, hj], TS, nsamples=10,
                                              nwarmup=2, basis_ops=("Sz",),
                                              seed=seed)
            say("p2 julia metts_vev spin basis=(Sz,) seed=%d" % seed,
                "<Sztot^2> %.6f stderr %.2e diff %.3e | <H> %.6f stderr %.2e diff %.3e"
                % (np.real(m1), s1, np.real(m1) - ex[0],
                   np.real(m2), s2, np.real(m2) - ex[1]))
        except Exception as e:
            say("p2 julia metts_vev spin basis=(Sz,) seed=%d" % seed, err(e))
    try:
        (m1, s1), (m2, s2) = sj.metts_vev([Szj * Szj, hj], TS, nsamples=10,
                                          nwarmup=2, seed=1)
        say("p2 julia metts_vev spin basis=(Sz,Sx) control seed=1",
            "<Sztot^2> %.6f stderr %.2e diff %.3e | <H> %.6f stderr %.2e diff %.3e"
            % (np.real(m1), s1, np.real(m1) - ex[0],
               np.real(m2), s2, np.real(m2) - ex[1]))
    except Exception as e:
        say("p2 julia metts_vev spin control", err(e))
    # a non-Hermitian basis name: python refuses it; julia_live takes
    # eigen(Hermitian(M)), the upper triangle, which is zero for S-
    if sc is not None:
        try:
            sc.metts_vev(Sz * Sz, TS, nsamples=4, nwarmup=1,
                         basis_ops=("S-",), seed=1)
            say("p2 python metts_vev basis=(S-,)", "RETURNED (no guard)")
        except ValueError as e:
            say("p2 python metts_vev basis=(S-,)", "ValueError (guard fires)")
        except Exception as e:
            say("p2 python metts_vev basis=(S-,)", err(e))
    try:
        m1, s1 = sj.metts_vev(Szj * Szj, TS, nsamples=10, nwarmup=2,
                              basis_ops=("S-",), seed=2)
        say("p2 julia metts_vev spin basis=(S-,) seed=2",
            "<Sztot^2> %.6f stderr %.2e diff %.3e"
            % (np.real(m1), s1, np.real(m1) - ex[0]))
    except Exception as e:
        say("p2 julia metts_vev spin basis=(S-,)", err(e))
try:
    fj, hfj, Ntj = fermion_chain("julia_live")
    for seed in (1, 2):
        try:
            (m1, s1), (m2, s2) = fj.metts_vev([Ntj * Ntj, Ntj], TF, nsamples=10,
                                              nwarmup=2, basis_ops=("N",),
                                              seed=seed)
            say("p2 julia metts_vev fermion basis=(N,) seed=%d" % seed,
                "<Ntot^2> %.6f stderr %.2e diff %.3e | <Ntot> %.6f stderr %.2e diff %.3e"
                % (np.real(m1), s1, np.real(m1) - exf[0],
                   np.real(m2), s2, np.real(m2) - exf[1]))
        except Exception as e:
            say("p2 julia metts_vev fermion basis=(N,) seed=%d" % seed, err(e))
except Exception as e:
    say("p2 julia fermion setup", err(e))

# ---------------------------------------------------------------- Part 3
# candidate julia-2 on julia_live: in an Nf sector the default dmmode is
# "full" (correlations.jl), which reads the chain's own C, not operators=
try:
    fs, _, _ = fermion_chain("julia_live")
    fs.set_conserved_sector(Nf=2)
    r = float(2 ** -0.5)
    ops = [r * (fs.C[0] + fs.C[1]), r * (fs.C[0] - fs.C[1]), fs.C[2], fs.C[3]]
    U = np.array([[r, r, 0, 0], [r, -r, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    M0 = np.array(fs.get_correlation_matrix())
    Md = np.array(fs.get_correlation_matrix(operators=ops))
    Me = np.array(fs.get_correlation_matrix(operators=ops, dmmode="explicit"))
    Mr = U @ M0 @ U.T
    say("p3 julia sector default(ops) shape", Md.shape)
    say("p3 julia sector max|default(ops) - explicit(ops)|",
        "%.3e" % np.max(np.abs(Md - Me)))
    say("p3 julia sector max|default(ops) - default(no ops)|",
        "%.3e" % np.max(np.abs(Md - M0)))
    say("p3 julia sector max|explicit(ops) - U M0 U^T|",
        "%.3e" % np.max(np.abs(Me - Mr)))
except Exception as e:
    say("p3 julia sector correlation matrix", err(e))

# ---------------------------------------------------------------- Part 4
if sj is not None:
    Bj = 3.0 + Szj
    for seed in (1, 2):
        try:
            ts, m, s = sj.metts_dynamical_correlator(
                (Bj, Bj), TS, nt=1, dt=0.1, nsamples=6, nwarmup=2,
                basis_ops=("Sz",), seed=seed)
            say("p4 julia dynMETTS spin basis=(Sz,) seed=%d" % seed,
                "C(0)=<(3+Sztot)^2> %.6f stderr %.2e | ED %.6f | diff %.3e"
                % (m[0].real, s[0], ex[2], m[0].real - ex[2]))
        except Exception as e:
            say("p4 julia dynMETTS seed=%d" % seed, err(e))
say("done", "ok")
```

Observed:

```
ED spin4 T=0.5 <Sztot^2> <H> <(3+Sztot)^2> = 0.370227 -1.119723 9.370227   [1s]
p1 python metts_vev spin basis=(Sz,) = ValueError (guard fires)   [1s]
p1 python dynMETTS spin basis=(Sz,) seed=1 = C(0)=<(3+Sztot)^2> 16.000000 stderr 2.30e-15 | ED 9.370227 | diff 6.630e+00   [3s]
p1 python dynMETTS spin basis=(Sz,) seed=2 = C(0)=<(3+Sztot)^2> 4.000000 stderr 5.81e-16 | ED 9.370227 | diff -5.370e+00   [4s]
p1 python dynMETTS spin basis=(Sz,) seed=3 = C(0)=<(3+Sztot)^2> 4.000000 stderr 4.85e-16 | ED 9.370227 | diff -5.370e+00   [6s]
p1 python dynMETTS spin basis=(Sz,Sx) control seed=1 = C(0)=<(3+Sztot)^2> 10.002268 stderr 8.58e-01 | ED 9.370227 | diff 6.320e-01   [8s]
ED fermion4 T=1 <Ntot^2> <Ntot> <(1+Ntot)^2> = 4.583062 1.963635 9.510333   [8s]
p1 python metts_vev fermion basis=(N,) = ValueError (guard fires)   [8s]
p1 python dynMETTS fermion basis=(N,) seed=1 = C(0)=<(1+Ntot)^2> 16.000000 stderr 4.89e-15 | ED 9.510333 | diff 6.490e+00   [9s]
p1 python dynMETTS fermion basis=(N,) seed=2 = C(0)=<(1+Ntot)^2> 4.000000 stderr 2.45e-16 | ED 9.510333 | diff -5.510e+00   [10s]
p1 python dynMETTS fermion basis=(N,) seed=3 = C(0)=<(1+Ntot)^2> 4.000000 stderr 5.93e-17 | ED 9.510333 | diff -5.510e+00   [11s]
julia ITensorMPS ITensors versions = 0.4.2 0.9.32   [29s]
p2 julia metts_vev spin basis=(Sz,) seed=1 = <Sztot^2> 1.000000 stderr 2.34e-17 diff 6.298e-01 | <H> -0.830847 stderr 2.95e-06 diff 2.889e-01   [122s]
p2 julia metts_vev spin basis=(Sz,) seed=2 = <Sztot^2> 4.000000 stderr 0.00e+00 diff 3.630e+00 | <H> 0.750000 stderr 0.00e+00 diff 1.870e+00   [123s]
p2 julia metts_vev spin basis=(Sz,) seed=3 = <Sztot^2> 0.000000 stderr 6.23e-18 diff -3.702e-01 | <H> -1.373139 stderr 1.38e-02 diff -2.534e-01   [127s]
p2 julia metts_vev spin basis=(Sz,Sx) control seed=1 = <Sztot^2> 0.422315 stderr 1.20e-01 diff 5.209e-02 | <H> -1.076149 stderr 1.18e-01 diff 4.357e-02   [131s]
p2 python metts_vev basis=(S-,) = ValueError (guard fires)   [131s]
p2 julia metts_vev spin basis=(S-,) seed=2 = <Sztot^2> 4.000000 stderr 0.00e+00 diff 3.630e+00   [132s]
p2 julia metts_vev fermion basis=(N,) seed=1 = <Ntot^2> 9.000000 stderr 3.74e-16 diff 4.417e+00 | <Ntot> 3.000000 stderr 1.15e-16 diff 1.036e+00   [135s]
p2 julia metts_vev fermion basis=(N,) seed=2 = <Ntot^2> 0.000000 stderr 0.00e+00 diff -4.583e+00 | <Ntot> 0.000000 stderr 0.00e+00 diff -1.964e+00   [135s]
p3 julia sector default(ops) shape = (4, 4)   [283s]
p3 julia sector max|default(ops) - explicit(ops)| = 5.067e-01   [283s]
p3 julia sector max|default(ops) - default(no ops)| = 0.000e+00   [283s]
p3 julia sector max|explicit(ops) - U M0 U^T| = 2.220e-16   [283s]
p4 julia dynMETTS spin basis=(Sz,) seed=1 = C(0)=<(3+Sztot)^2> 16.000000 stderr 4.59e-16 | ED 9.370227 | diff 6.630e+00   [292s]
p4 julia dynMETTS spin basis=(Sz,) seed=2 = C(0)=<(3+Sztot)^2> 1.000000 stderr 0.00e+00 | ED 9.370227 | diff -8.370e+00   [293s]
done = ok   [293s]
[exit 0, 296 s]
```

Reviewer's counter-repro, v3 and a non-Hermitian name on python:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g05/counter_v3_dyn_and_nonhermitian.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Reviewer counter-repro for candidate julia-1 (group 05).
# Part A: itensor_version=3, metts_dynamical_correlator with basis ("Sz",) on a
#   4-site Heisenberg chain (H conserves Sztot). The hunter measured python
#   and julia_live only; this measures v3, which shares the C++ Markov chain
#   of metts_vev (random_cps/collapse_to_cps) but has no basis check.
# Part B: "python", metts_dynamical_correlator with a NON-Hermitian basis
#   name, ("S+",) and ("S-",). metts_vev refuses both on python; the dynamical
#   correlator has no check, and pyitensor's eigh reads one triangle only.
# C(0) = <(3+Sztot)^2>; a chain stuck in one Sz sector gives (3+M)^2 exactly,
# M in {-2,...,2}, with stderr ~1e-15.
np.random.seed(0)
from dmrgpy import spinchain

T = 0.5
n = 4


def dense(ed, m):
    m = ed.MO2matrix(m)
    if hasattr(m, "todense"):
        m = m.todense()
    return np.asarray(m, dtype=complex)


def thermal(chain, ops, T):
    ed = chain.get_ED_obj()
    H = dense(ed, chain.hamiltonian)
    e, U = np.linalg.eigh(H)
    w = np.exp(-(e - e.min()) / T)
    w = w / w.sum()
    return [float(np.real(np.sum(w * np.diag(U.conj().T @ dense(ed, O) @ U))))
            for O in ops]


def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=v, maxm=16,
                              nsweeps=6, cutoff=1e-12)
    h = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    sc.set_hamiltonian(h)
    return sc


def err(e):
    return "ERROR %s: %s" % (type(e).__name__, str(e)[:140])


ref = None
for v, bases in ((3, [("Sz",), ("Sz", "Sx")]),
                 ("python", [("S+",), ("S-",)])):
    try:
        sc = chain(v)
        B = 3.0 + sum(sc.Sz)
        if ref is None:
            ref = thermal(sc, [B * B], T)[0]
            print("ED <(3+Sztot)^2> T=0.5 = %.6f" % ref, flush=True)
    except Exception as e:
        print("setup v=%s %s" % (v, err(e)), flush=True)
        continue
    for basis in bases:
        for seed in (1, 2, 3):
            tag = "v=%s dynMETTS basis=%s seed=%d" % (v, basis, seed)
            try:
                ts, m, s = sc.metts_dynamical_correlator(
                    (B, B), T, nt=1, dt=0.1, nsamples=8, nwarmup=2,
                    dbeta_half_step=0.05, basis_ops=basis, seed=seed)
                print("%s C(0) %.6f stderr %.2e | ED %.6f | diff %.3e"
                      % (tag, m[0].real, s[0], ref, m[0].real - ref), flush=True)
            except ValueError as e:
                print("%s ValueError (guard fires): %s" % (tag, str(e)[:100]), flush=True)
            except Exception as e:
                print("%s %s" % (tag, err(e)), flush=True)
print("done")
```

Observed:

```
ED <(3+Sztot)^2> T=0.5 = 9.370227
v=3 dynMETTS basis=('Sz',) seed=1 C(0) 25.000000 stderr 0.00e+00 | ED 9.370227 | diff 1.563e+01
v=3 dynMETTS basis=('Sz',) seed=2 C(0) 1.000000 stderr 0.00e+00 | ED 9.370227 | diff -8.370e+00
v=3 dynMETTS basis=('Sz',) seed=3 C(0) 9.000000 stderr 3.36e-16 | ED 9.370227 | diff -3.702e-01
v=3 dynMETTS basis=('Sz', 'Sx') seed=1 C(0) 8.623988 stderr 6.72e-01 | ED 9.370227 | diff -7.462e-01
v=3 dynMETTS basis=('Sz', 'Sx') seed=2 C(0) 11.219452 stderr 1.54e+00 | ED 9.370227 | diff 1.849e+00
v=3 dynMETTS basis=('Sz', 'Sx') seed=3 C(0) 10.904589 stderr 1.11e+00 | ED 9.370227 | diff 1.534e+00
v=python dynMETTS basis=('S+',) seed=1 C(0) 4.000000 stderr 6.22e-16 | ED 9.370227 | diff -5.370e+00
v=python dynMETTS basis=('S+',) seed=2 C(0) 16.000000 stderr 1.20e-15 | ED 9.370227 | diff 6.630e+00
v=python dynMETTS basis=('S+',) seed=3 C(0) 16.000000 stderr 1.14e-15 | ED 9.370227 | diff 6.630e+00
v=python dynMETTS basis=('S-',) seed=1 C(0) 9.541062 stderr 8.17e-03 | ED 9.370227 | diff 1.708e-01
v=python dynMETTS basis=('S-',) seed=2 C(0) 9.532892 stderr 6.67e-07 | ED 9.370227 | diff 1.627e-01
v=python dynMETTS basis=('S-',) seed=3 C(0) 9.590075 stderr 8.17e-03 | ED 9.370227 | diff 2.198e-01
done
[exit 0, 21 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: "every such call returns with stderr 0" as a general statement,
since it holds only for a function of the frozen charge; the expected size
is replaced by the measured 0.37 to 15.6 on the charge observables and 0.25
to 1.87 on `<H>`. The v3 part of the dynamical correlator, first confirmed by
reading, is measured by the counter-repro (25, 1 and 9 at stderr 0).

**Suggested fix**: Drop the version filter at `mettsvev.py:194`; everything `_check_basis`
needs works on `julia_live`. Call `_check_basis` from
`mettsdynamicalcorrelator.py` on every backend, after the cheap argument
checks and before `str2MO`, with a caller argument so the message names the
call. As defence in depth, each backend's eigenbasis should raise when
`||M - M^dag|| > 1e-12 ||M||`, since the three backends read three different
triangles; on v3 as a catchable C++ throw. Pin it by adding `julia_live`
and the dynamical correlator to the 2026-10-23 tests. The guard itself is
also incomplete, which is the reviewer-found finding below.

### 7. Where the METTS guard runs, it misses a single non-diagonal name and every Z2 parity: `("Sx",)` on a Heisenberg chain and `("N",)` on a Kitaev chain freeze at stderr 0

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `julia` &middot; candidate found by the reviewer of `julia-1`

**Where**: `src/dmrgpy/vevtk/mettsvev.py:57 (_check_basis returns at the first non-diagonal name)`; `src/dmrgpy/vevtk/mettsvev.py (_check_basis tests U(1) charges only, never a parity)`

`_check_basis` misses two kinds of frozen Markov chain on `"python"` and
v3. First, it tests a conserved charge only when every basis name is
diagonal, so a single non-diagonal name O whose sum H conserves is never
tested. With `("Sx",)` on the isotropic 4-site Heisenberg chain at T = 0.5,
`<Sxtot^2>` is 1 on every seed on `"python"` and 4, 4, 0 by seed on v3, at
stderr below 1e-14, against 0.370 from ED; two v3 seeds sit on the fully
polarized eigenstate, `<H> = 0.750` at stderr 0 against -1.120. Second, it
never tests a Z2 parity. With `("Sx",)` on the transverse-field Ising chain
with its field along x, `<Px>` is -1 on `"python"` and +1 on v3 on all
three seeds at stderr 0, against 0.038. With `("N",)` on a 4-site Kitaev
chain, which does not conserve N, `<(-1)^N>` is -1 on `"python"` and +1 on
v3 at stderr 0, against 0.043. On a spinless fermion site every Hermitian
name the sites define is diagonal, so every METTS run on a spinless fermion
chain without a conserved N is frozen in one fermion parity, and nothing
warns.

**Expected**: a refusal, as for an all-diagonal basis on a conserved U(1) charge, or for the fermion parity a stated choice.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g05/new_single_offdiagonal_basis.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# New candidate from the review of julia-1: vevtk/mettsvev.py::_check_basis
# tests a conserved charge only when EVERY basis name is diagonal, and only
# U(1) charges (sums of one-site operators), never a Z2 parity (a product of
# one-site operators). Three bases it lets through, each of which freezes
# the Markov chain in one symmetry sector of its first sample:
#   P1  ("Sx",) on the isotropic Heisenberg chain: H conserves Sxtot, the Sx
#       collapse gives states of definite Sxtot. Stuck: <Sxtot^2> = M^2.
#   P2  ("Sx",) on the transverse-field Ising chain (field along x): H
#       conserves Px = prod_i 2Sx_i, an Sx product state has Px = +-1.
#   P3  ("N",) on a spinless Kitaev (pairing) chain: H does not conserve N,
#       so the guard passes, but every fermionic H conserves (-1)^N, and an
#       occupation state has definite parity. On a spinless fermion site
#       every Hermitian parity-even one-site operator is diagonal, so no
#       basis exists that leaves the parity sector.
# Controls: ("Sz","Sx") on P1 and P2. Anchor: dense ED Boltzmann average.
np.random.seed(0)
from dmrgpy import spinchain, fermionchain

T = 0.5
n = 4


def dense(ed, m):
    m = ed.MO2matrix(m)
    if hasattr(m, "todense"):
        m = m.todense()
    return np.asarray(m, dtype=complex)


def thermal(chain, ops, T):
    ed = chain.get_ED_obj()
    H = dense(ed, chain.hamiltonian)
    e, U = np.linalg.eigh(H)
    w = np.exp(-(e - e.min()) / T)
    w = w / w.sum()
    return [float(np.real(np.sum(w * np.diag(U.conj().T @ dense(ed, O) @ U))))
            for O in ops]


def err(e):
    return "ERROR %s: %s" % (type(e).__name__, str(e)[:140])


def build(v, model):
    if model == "kitaev":
        c = fermionchain.Fermionic_Chain(n, itensor_version=v, maxm=16,
                                         nsweeps=6, cutoff=1e-12)
        h = sum(-1.0 * (c.Cdag[i] * c.C[i + 1] + c.Cdag[i + 1] * c.C[i])
                + 0.6 * (c.C[i] * c.C[i + 1] + c.Cdag[i + 1] * c.Cdag[i])
                for i in range(n - 1))
        h = h + 0.3 * sum(c.N)
        P = 1.0
        for i in range(n):
            P = P * (1.0 - 2.0 * c.N[i])
        obs = [("parity", P), ("Ntot", sum(c.N)), ("H", h)]
        bases = [("N",)]
    else:
        c = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=v, maxm=16,
                                 nsweeps=6, cutoff=1e-12)
        if model == "heis":
            h = sum(c.Sx[i] * c.Sx[i + 1] + c.Sy[i] * c.Sy[i + 1]
                    + c.Sz[i] * c.Sz[i + 1] for i in range(n - 1))
            obs = [("Sxtot^2", sum(c.Sx) * sum(c.Sx)), ("H", h)]
        else:  # transverse-field Ising, field along x
            h = sum(c.Sz[i] * c.Sz[i + 1] for i in range(n - 1)) + 0.5 * sum(c.Sx)
            Px = 1.0
            for i in range(n):
                Px = Px * (2.0 * c.Sx[i])
            obs = [("Px", Px), ("H", h)]
        bases = [("Sx",), ("Sz", "Sx")]
    c.set_hamiltonian(h)
    return c, obs, bases


for v in ("python", 3):
    for model in ("heis", "tfim", "kitaev"):
        try:
            c, obs, bases = build(v, model)
            ex = thermal(c, [o for _, o in obs], T)
            print("ED %s T=0.5 %s" % (model, " ".join("<%s> %.6f" % (nm, x)
                  for (nm, _), x in zip(obs, ex))), flush=True)
        except Exception as e:
            print("setup v=%s %s %s" % (v, model, err(e)), flush=True)
            continue
        for basis in bases:
            seeds = (1, 2, 3) if len(basis) == 1 else (1,)
            for seed in seeds:
                tag = "v=%s %s metts_vev basis=%s seed=%d" % (v, model, basis, seed)
                try:
                    res = c.metts_vev([o for _, o in obs], T, nsamples=10,
                                      nwarmup=2, dbeta_half_step=0.05,
                                      basis_ops=basis, seed=seed)
                    print(tag + " " + " | ".join(
                        "<%s> %.6f stderr %.2e diff %.3e"
                        % (nm, np.real(m), s, np.real(m) - x)
                        for (nm, _), (m, s), x in zip(obs, res, ex)), flush=True)
                except ValueError as e:
                    print("%s ValueError (guard fires): %s" % (tag, str(e)[:100]), flush=True)
                except Exception as e:
                    print("%s %s" % (tag, err(e)), flush=True)
print("done")
```

Observed:

```
ED heis T=0.5 <Sxtot^2> 0.370227 <H> -1.119723
v=python heis metts_vev basis=('Sx',) seed=1 <Sxtot^2> 1.000000 stderr 1.13e-15 diff 6.298e-01 | <H> -0.788112 stderr 4.27e-02 diff 3.316e-01
v=python heis metts_vev basis=('Sx',) seed=2 <Sxtot^2> 1.000000 stderr 1.41e-15 diff 6.298e-01 | <H> -0.830844 stderr 2.95e-06 diff 2.889e-01
v=python heis metts_vev basis=('Sx',) seed=3 <Sxtot^2> 1.000000 stderr 1.36e-16 diff 6.298e-01 | <H> -0.531750 stderr 6.53e-02 diff 5.880e-01
v=python heis metts_vev basis=('Sz', 'Sx') seed=1 <Sxtot^2> 0.145092 stderr 5.72e-02 diff -2.251e-01 | <H> -1.279493 stderr 7.61e-02 diff -1.598e-01
ED tfim T=0.5 <Px> 0.037647 <H> -0.722737
v=python tfim metts_vev basis=('Sx',) seed=1 <Px> -1.000000 stderr 0.00e+00 diff -1.038e+00 | <H> -0.665928 stderr 7.73e-02 diff 5.681e-02
v=python tfim metts_vev basis=('Sx',) seed=2 <Px> -1.000000 stderr 3.70e-17 diff -1.038e+00 | <H> -0.726147 stderr 3.70e-17 diff -3.410e-03
v=python tfim metts_vev basis=('Sx',) seed=3 <Px> -1.000000 stderr 3.70e-17 diff -1.038e+00 | <H> -0.726147 stderr 3.70e-17 diff -3.410e-03
v=python tfim metts_vev basis=('Sz', 'Sx') seed=1 <Px> -0.283238 stderr 2.16e-01 diff -3.209e-01 | <H> -0.718969 stderr 7.91e-02 diff 3.768e-03
ED kitaev T=0.5 <parity> 0.043112 <Ntot> 1.756881 <H> -1.718767
v=python kitaev metts_vev basis=('N',) seed=1 <parity> -1.000000 stderr 7.12e-17 diff -1.043e+00 | <Ntot> 1.440755 stderr 7.13e-02 diff -3.161e-01 | <H> -1.711103 stderr 1.44e-02 diff 7.663e-03
v=python kitaev metts_vev basis=('N',) seed=2 <parity> -1.000000 stderr 1.33e-16 diff -1.043e+00 | <Ntot> 1.876269 stderr 1.58e-01 diff 1.194e-01 | <H> -1.661249 stderr 2.54e-02 diff 5.752e-02
v=python kitaev metts_vev basis=('N',) seed=3 <parity> -1.000000 stderr 1.05e-16 diff -1.043e+00 | <Ntot> 1.545203 stderr 1.14e-01 diff -2.117e-01 | <H> -1.694090 stderr 1.90e-02 diff 2.468e-02
ED heis T=0.5 <Sxtot^2> 0.370227 <H> -1.119723
v=3 heis metts_vev basis=('Sx',) seed=1 <Sxtot^2> 4.000000 stderr 0.00e+00 diff 3.630e+00 | <H> 0.750000 stderr 0.00e+00 diff 1.870e+00
v=3 heis metts_vev basis=('Sx',) seed=2 <Sxtot^2> 4.000000 stderr 0.00e+00 diff 3.630e+00 | <H> 0.750000 stderr 0.00e+00 diff 1.870e+00
v=3 heis metts_vev basis=('Sx',) seed=3 <Sxtot^2> 0.000000 stderr 1.20e-18 diff -3.702e-01 | <H> -1.373139 stderr 1.38e-02 diff -2.534e-01
v=3 heis metts_vev basis=('Sz', 'Sx') seed=1 <Sxtot^2> 0.686971 stderr 3.82e-01 diff 3.167e-01 | <H> -0.968916 stderr 2.21e-01 diff 1.508e-01
ED tfim T=0.5 <Px> 0.037647 <H> -0.722737
v=3 tfim metts_vev basis=('Sx',) seed=1 <Px> 1.000000 stderr 7.40e-17 diff 9.624e-01 | <H> -0.310573 stderr 9.26e-03 diff 4.122e-01
v=3 tfim metts_vev basis=('Sx',) seed=2 <Px> 1.000000 stderr 3.31e-17 diff 9.624e-01 | <H> -1.023916 stderr 9.08e-02 diff -3.012e-01
v=3 tfim metts_vev basis=('Sx',) seed=3 <Px> 1.000000 stderr 1.02e-16 diff 9.624e-01 | <H> -0.471957 stderr 2.60e-02 diff 2.508e-01
v=3 tfim metts_vev basis=('Sz', 'Sx') seed=1 <Px> -0.286346 stderr 2.16e-01 diff -3.240e-01 | <H> -0.713816 stderr 8.96e-02 diff 8.921e-03
ED kitaev T=0.5 <parity> 0.043112 <Ntot> 1.756881 <H> -1.718767
v=3 kitaev metts_vev basis=('N',) seed=1 <parity> 1.000000 stderr 1.02e-16 diff 9.569e-01 | <Ntot> 1.789965 stderr 1.15e-01 diff 3.308e-02 | <H> -1.719779 stderr 3.40e-02 diff -1.013e-03
v=3 kitaev metts_vev basis=('N',) seed=2 <parity> 1.000000 stderr 7.02e-17 diff 9.569e-01 | <Ntot> 1.900023 stderr 2.84e-02 diff 1.431e-01 | <H> -1.811124 stderr 1.88e-02 diff -9.236e-02
v=3 kitaev metts_vev basis=('N',) seed=3 <parity> 1.000000 stderr 7.76e-17 diff 9.569e-01 | <Ntot> 1.813408 stderr 7.69e-02 diff 5.653e-02 | <H> -1.786279 stderr 2.95e-02 diff -6.751e-02
done
[exit 0, 46 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: "chosen by the seed" for the two parity parts, since `<Px>` and
`<(-1)^N>` came out the same on all three seeds of each backend; and "the
in-sector `<H>` with a small, finite error bar", since some runs froze on a
single state. The `(Sz, Sx)` controls scatter around ED with a finite error
bar, so the Markov chain is right when the basis is.

**Suggested fix**: Drop the early return and test the charges whatever the names are: freezing
happens when every basis name commutes, on one site, with some one-site q
whose sum H conserves, q being each name itself or a sector charge; for
spin-1/2 names also test the parity `prod_i 2O_i`. The fermion parity is
the author's pick, since every fermionic H conserves it and no basis leaves
it on a spinless site: raise on every spinless fermion chain without a
conserved N, warn once naming the parity sector, or document it as an
ensemble at fixed parity whose bias falls off as `Tr(P e^-bH)/Z`.

### 8. In a conserved sector, `get_correlation_matrix(operators=X)` ignores X on every DMRG backend: a Hadamard rotation is off by 0.507 and `C[:2]` returns the (n, n) matrix

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `julia, v3-cpp` &middot; candidate `julia-2` and `v3-cpp-1`, one defect

**Where**: `src/dmrgpy/entropytk/correlationentropy.py:89-95 (the operators= guard, skipped in a sector)`; `src/dmrgpy/entropytk/correlationentropy.py:334-339 (_default_dmmode's sector branch returns "full")`; `src/dmrgpy/entropytk/correlationentropy.py:116-117 and 265-279 ("full" never sees operators)`

Two hunters found the same defect, one on `julia_live` and one on v3, and
each measured it on `"python"` as well. In a sector `_default_dmmode`
returns `"full"`, the one-sweep matrix of the chain's own `Cdag_i C_j`, and
the guard that sends a caller's `operators=` to a route that reads them is
skipped when a sector is set. On chains of spinless fermion sites the
default then ignores X, on `"python"`, v3 and `julia_live`, with no warning.
Against `dmmode="explicit"`, which reads X and agrees with ED to 1e-13 in
the sector: a Hadamard rotation of modes 0 and 1 on 4 sites at Nf=2 is off
by 0.507 on all three backends, and the default there equals the matrix
without `operators=` exactly; on 6 sites at Nf=3 a reversed list is off by
0.123 and a momentum basis by 0.541; `C[:2]` comes back as the (4, 4) or
(6, 6) matrix of the whole chain instead of (2, 2). The readers that forward
`operators=` inherit it: the counter-repro shows the correlation entropy of
a 3-mode subsystem in a sector is 0.000, the entropy of the whole pure
chain, against 0.720 from ED. Outside a sector the default reads X to 1e-12.
An explicit `dmmode="full"` with `operators=` ignores X in or out of a
sector. A `Spin_Chain` in an Sz sector called with `operators=[Sz_i]` raises
instead of answering, while `dmmode="explicit"` answers. The defect predates
the window: `1b87543` (2026-08-30, the fix of 2026-08 finding 11) made
`"full"` the sector default without looking at `operators=`, and `225be82`
added the guard for the no-sector case only.

**Expected**: the matrix of the operators passed, as `dmmode="explicit"` and the no-sector default return.

Repro (python and v3, 4 sites):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/julia/julia_2_sector_corrmatrix.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Candidate julia-2: get_correlation_matrix(operators=ops) in a conserved
# sector, with the default dmmode, answers with the chain's own
# <Cdag_i C_j> and never reads ops. _default_dmmode returns "full" in a
# sector (julia_live spinless, or any session with correlation_matrix), and
# the guard that sends "full" back to "fast" when operators= is given has
# `not conserved_sector` in it, so cpp_correlation_matrix(wf) answers.
# dmmode="explicit" reads ops and is the reference; U M0 U^T is the exact
# rotation of the chain's own matrix M0.
import time

np.random.seed(0)

from dmrgpy import fermionchain

T0 = time.time()


def say(name, value):
    print("%s = %s   [%.0fs]" % (name, value, time.time() - T0), flush=True)


def chain(backend, n=4):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=backend, maxm=16,
                                      nsweeps=8, cutoff=1e-12)
    h = sum(-1.0 * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
            for i in range(n - 1))
    h = h + 0.3 * fc.N[0] - 0.2 * fc.N[n - 1] + 0.1 * fc.N[1]
    fc.set_hamiltonian(h)
    return fc


r = float(2 ** -0.5)
U = np.array([[r, r, 0, 0], [r, -r, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

for backend in ("python", 3):
    tag = "backend=%r" % (backend,)
    for sector in (True, False):
        stag = tag + (" sector Nf=2" if sector else " no sector")
        try:
            fc = chain(backend)
            if sector:
                fc.set_conserved_sector(Nf=2)
            ops = [r * (fc.C[0] + fc.C[1]), r * (fc.C[0] - fc.C[1]),
                   fc.C[2], fc.C[3]]
            M0 = np.array(fc.get_correlation_matrix())
            Md = np.array(fc.get_correlation_matrix(operators=ops))
            Me = np.array(fc.get_correlation_matrix(operators=ops,
                                                    dmmode="explicit"))
            Mr = U @ M0 @ U.T
            say(stag + " mode of the solve", fc.get_mode())
            say(stag + " max|default(ops) - explicit(ops)|",
                "%.3e" % np.max(np.abs(Md - Me)))
            say(stag + " max|default(ops) - default(no ops)|",
                "%.3e" % np.max(np.abs(Md - M0)))
            say(stag + " max|explicit(ops) - U M0 U^T|",
                "%.3e" % np.max(np.abs(Me - Mr)))
            say(stag + " [0,0] default(ops) explicit(ops) rotated",
                "%.6f %.6f %.6f" % (Md[0, 0].real, Me[0, 0].real, Mr[0, 0].real))
            sub = [fc.C[0], fc.C[1]]
            Ms = np.array(fc.get_correlation_matrix(operators=sub))
            Mse = np.array(fc.get_correlation_matrix(operators=sub,
                                                     dmmode="explicit"))
            say(stag + " shape default(C0,C1) explicit(C0,C1)",
                "%s %s" % (Ms.shape, Mse.shape))
            # explicit dmmode="full" with operators=, outside the default
            Mf = np.array(fc.get_correlation_matrix(operators=ops,
                                                    dmmode="full"))
            say(stag + " max|dmmode='full'(ops) - explicit(ops)|",
                "%.3e" % np.max(np.abs(Mf - Me)))
        except Exception as e:
            say(stag, "ERROR %s: %s" % (type(e).__name__, str(e)[:200]))
say("done", "ok")
```

Observed:

```
backend='python' sector Nf=2 mode of the solve = DMRG   [1s]
backend='python' sector Nf=2 max|default(ops) - explicit(ops)| = 5.067e-01   [1s]
backend='python' sector Nf=2 max|default(ops) - default(no ops)| = 0.000e+00   [1s]
backend='python' sector Nf=2 max|explicit(ops) - U M0 U^T| = 1.111e-15   [1s]
backend='python' sector Nf=2 [0,0] default(ops) explicit(ops) rotated = 0.401066 0.907759 0.907759   [1s]
backend='python' sector Nf=2 shape default(C0,C1) explicit(C0,C1) = (4, 4) (2, 2)   [1s]
backend='python' sector Nf=2 max|dmmode='full'(ops) - explicit(ops)| = 5.067e-01   [1s]
backend='python' no sector mode of the solve = DMRG   [1s]
backend='python' no sector max|default(ops) - explicit(ops)| = 5.968e-16   [1s]
backend='python' no sector max|default(ops) - default(no ops)| = 5.067e-01   [1s]
backend='python' no sector max|explicit(ops) - U M0 U^T| = 7.425e-16   [1s]
backend='python' no sector [0,0] default(ops) explicit(ops) rotated = 0.907759 0.907759 0.907759   [1s]
backend='python' no sector shape default(C0,C1) explicit(C0,C1) = (2, 2) (2, 2)   [1s]
backend='python' no sector max|dmmode='full'(ops) - explicit(ops)| = 5.067e-01   [1s]
backend=3 sector Nf=2 mode of the solve = DMRG   [1s]
backend=3 sector Nf=2 max|default(ops) - explicit(ops)| = 5.067e-01   [1s]
backend=3 sector Nf=2 max|default(ops) - default(no ops)| = 0.000e+00   [1s]
backend=3 sector Nf=2 max|explicit(ops) - U M0 U^T| = 4.441e-16   [1s]
backend=3 sector Nf=2 [0,0] default(ops) explicit(ops) rotated = 0.401066 0.907759 0.907759   [1s]
backend=3 sector Nf=2 shape default(C0,C1) explicit(C0,C1) = (4, 4) (2, 2)   [1s]
backend=3 sector Nf=2 max|dmmode='full'(ops) - explicit(ops)| = 5.067e-01   [1s]
backend=3 no sector mode of the solve = DMRG   [1s]
backend=3 no sector max|default(ops) - explicit(ops)| = 2.220e-16   [1s]
backend=3 no sector max|default(ops) - default(no ops)| = 5.067e-01   [1s]
backend=3 no sector max|explicit(ops) - U M0 U^T| = 5.551e-16   [1s]
backend=3 no sector [0,0] default(ops) explicit(ops) rotated = 0.907759 0.907759 0.907759   [1s]
backend=3 no sector shape default(C0,C1) explicit(C0,C1) = (2, 2) (2, 2)   [1s]
backend=3 no sector max|dmmode='full'(ops) - explicit(ops)| = 5.067e-01   [1s]
done = ok   [1s]
[exit 0, 3 s]
```

Repro (v3 and python, 6 sites, reversed, momentum and subset):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/v3-cpp/v3cpp_2_corrmat_operators_sector.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# get_correlation_matrix(operators=X) on a chain with a conserved sector:
# _default_dmmode returns "full" in the sector branch, and the override
# back to an operator-reading route (correlationentropy.py, 225be82) is
# skipped "not conserved_sector", so the session matrix <Cdag_i C_j> of the
# chain's own site operators comes back whatever X is (on v3 through
# Chain::correlation_matrix_sweep). Expected: <X_i^dagger X_j>, which
# dmmode="explicit" and mode="ED" both compute.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import fermionchain, cppext

n = 6


def build(version, sector):
    np.random.seed(3)
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version, maxm=40,
                                      nsweeps=10)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    for i in range(n-1):
        h = h + 0.8*fc.N[i]*fc.N[i+1]
    h = h + 0.4*fc.N[0] - 0.2*fc.N[2] # no reflection symmetry
    fc.set_hamiltonian(h)
    if sector:
        fc.set_conserved_sector(Nf=3)
    return fc


def momentum_ops(fc):
    out = []
    for k in range(n):
        o = 0
        for j in range(n):
            o = o + (np.exp(2j*np.pi*k*j/n)/np.sqrt(n))*fc.C[j]
        out.append(o)
    return out


def mx(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape: return float("nan")
    return float(np.max(np.abs(a-b)))


versions = ["python"]
if cppext.available(3): versions.insert(0, 3)
else: print("backend v3 = not compiled, skipped")

for v in versions:
    for sector in (False, True):
        fc = build(v, sector)
        tag = "v=%s sector=%s" % (v, "Nf=3" if sector else "none")
        site = np.asarray(fc.get_correlation_matrix())
        # 1. reversed site order
        rev = fc.C[::-1]
        got = np.asarray(fc.get_correlation_matrix(operators=rev))
        exp = np.asarray(fc.get_correlation_matrix(operators=rev,
                                                   dmmode="explicit"))
        print("%s reversed: |default-explicit| = %.3e  |default-site| = %.3e"
              "  |explicit-site[::-1,::-1]| = %.3e" % (tag, mx(got, exp),
              mx(got, site), mx(exp, site[::-1, ::-1])))
        # 2. momentum basis c_k = sum_j e^{ikj} c_j / sqrt(n)
        ks = momentum_ops(fc)
        got = np.asarray(fc.get_correlation_matrix(operators=ks))
        exp = np.asarray(fc.get_correlation_matrix(operators=ks,
                                                   dmmode="explicit"))
        ed = np.asarray(fc.get_correlation_matrix(operators=ks, mode="ED"))
        print("%s momentum: |default-explicit| = %.3e  |default-ED| = %.3e"
              "  |explicit-ED| = %.3e  |default-site| = %.3e"
              % (tag, mx(got, exp), mx(got, ed), mx(exp, ed), mx(got, site)))
        print("%s momentum: diag default = %s" % (tag,
              np.round(np.real(np.diag(got)), 6)))
        print("%s momentum: diag ED      = %s" % (tag,
              np.round(np.real(np.diag(ed)), 6)))
        # 3. a subset of the operators
        sub = np.asarray(fc.get_correlation_matrix(operators=fc.C[:2]))
        print("%s subset operators=C[:2]: shape = %s (expected (2, 2))"
              % (tag, sub.shape))
```

Observed:

```
v=3 sector=none reversed: |default-explicit| = 1.443e-15  |default-site| = 1.231e-01  |explicit-site[::-1,::-1]| = 5.551e-16
v=3 sector=none momentum: |default-explicit| = 7.677e-13  |default-ED| = 4.996e-06  |explicit-ED| = 4.996e-06  |default-site| = 5.413e-01
v=3 sector=none momentum: diag default = [0.021632 0.060325 0.933633 0.990443 0.933633 0.060325]
v=3 sector=none momentum: diag ED      = [0.021632 0.060325 0.933637 0.990443 0.933637 0.060325]
v=3 sector=none subset operators=C[:2]: shape = (2, 2) (expected (2, 2))
v=3 sector=Nf=3 reversed: |default-explicit| = 1.231e-01  |default-site| = 0.000e+00  |explicit-site[::-1,::-1]| = 1.221e-15
v=3 sector=Nf=3 momentum: |default-explicit| = 5.413e-01  |default-ED| = 5.413e-01  |explicit-ED| = 1.436e-13  |default-site| = 0.000e+00
v=3 sector=Nf=3 momentum: diag default = [0.478481 0.501596 0.497682 0.496691 0.423959 0.60159 ]
v=3 sector=Nf=3 momentum: diag ED      = [0.021632 0.060325 0.933637 0.990443 0.933637 0.060325]
v=3 sector=Nf=3 subset operators=C[:2]: shape = (6, 6) (expected (2, 2))
v=python sector=none reversed: |default-explicit| = 1.388e-15  |default-site| = 1.230e-01  |explicit-site[::-1,::-1]| = 1.555e-15
v=python sector=none momentum: |default-explicit| = 1.145e-12  |default-ED| = 6.135e-05  |explicit-ED| = 6.135e-05  |default-site| = 5.412e-01
v=python sector=none momentum: diag default = [0.021634 0.060322 0.933579 0.990436 0.933576 0.060318]
v=python sector=none momentum: diag ED      = [0.021632 0.060325 0.933637 0.990443 0.933637 0.060325]
v=python sector=none subset operators=C[:2]: shape = (2, 2) (expected (2, 2))
v=python sector=Nf=3 reversed: |default-explicit| = 1.231e-01  |default-site| = 0.000e+00  |explicit-site[::-1,::-1]| = 1.110e-15
v=python sector=Nf=3 momentum: |default-explicit| = 5.413e-01  |default-ED| = 5.413e-01  |explicit-ED| = 1.070e-12  |default-site| = 0.000e+00
v=python sector=Nf=3 momentum: diag default = [0.478481 0.501596 0.497682 0.496691 0.423959 0.60159 ]
v=python sector=Nf=3 momentum: diag ED      = [0.021632 0.060325 0.933637 0.990443 0.933637 0.060325]
v=python sector=Nf=3 subset operators=C[:2]: shape = (6, 6) (expected (2, 2))
[exit 0, 4 s]
```

Reviewer's counter-repro, the subsystem entropy and a spin chain:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g06/g06_subsystem_entropy_and_spin.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Reviewer counter-repro for julia-2 / v3-cpp-1 (one defect): the size of
# the defect on the user-facing scalar, and its reach to a spin chain.
# Part A: the correlation entropy of a subsystem, operators=fc.C[:3], on a
# free-fermion chain (one Slater determinant, so the whole-chain entropy is
# zero and the subsystem one is not). Anchor: ED's full matrix, its 3x3
# block, Peschel's formula.
# Part B: a Spin_Chain in an Sz sector, operators=[Sz_i] (each conserves
# Sz). Anchor: mode="ED" on the same chain.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import fermionchain, spinchain, cppext

n = 6
L = 3


def ent(vs):
    vs = np.clip(np.real(vs), 1e-12, 1. - 1e-12)
    return float(-np.sum(vs*np.log(vs) + (1. - vs)*np.log(1. - vs)))


def free_chain(v, sector):
    np.random.seed(1)
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v, maxm=40,
                                      nsweeps=10, cutoff=1e-12)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    h = h + 0.4*fc.N[0] - 0.2*fc.N[2]
    fc.set_hamiltonian(h)
    if sector:
        fc.set_conserved_sector(Nf=3)
    return fc


def spin_chain(v, sector):
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v, maxm=40,
                              nsweeps=10, cutoff=1e-12)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] \
            + 0.7*sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.3*sc.Sz[0]
    sc.set_hamiltonian(h)
    if sector:
        sc.set_conserved_sector(Sz=0)
    return sc


def mx(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return float("nan")
    return float(np.max(np.abs(a - b)))


versions = ["python"]
if cppext.available(3):
    versions.append(3)  # after python, in case a C++ error aborts
else:
    print("v3 = not compiled, skipped")

for v in versions:
    for sector in (False, True):
        tag = "A v=%s sector=%s" % (v, "Nf=3" if sector else "none")
        try:
            fc = free_chain(v, sector)
            Med = np.asarray(fc.get_correlation_matrix(mode="ED"))
            s_ref = ent(np.linalg.eigvalsh(Med[:L, :L]))
            s_def = fc.get_correlation_entropy(operators=fc.C[:L])
            s_exp = fc.get_correlation_entropy(operators=fc.C[:L],
                                               dmmode="explicit")
            s_all = fc.get_correlation_entropy()
            print("%s S_sub default = %.6f  explicit = %.6f  ED block = %.6f"
                  "  S_whole default = %.6f" % (tag, s_def, s_exp, s_ref,
                                                s_all), flush=True)
            print("%s |S_sub default - ED| = %.3e  |S_sub explicit - ED| = %.3e"
                  % (tag, abs(s_def - s_ref), abs(s_exp - s_ref)), flush=True)
        except Exception as e:
            print("%s ERROR %s: %s" % (tag, type(e).__name__, str(e)[:200]),
                  flush=True)

for v in versions:
    for sector in (False, True):
        tag = "B v=%s sector=%s" % (v, "Sz=0" if sector else "none")
        try:
            sc = spin_chain(v, sector)
            ops = [sc.Sz[i] for i in range(n)]
            ed = np.asarray(sc.get_correlation_matrix(operators=ops,
                                                      mode="ED"))
            ex = np.asarray(sc.get_correlation_matrix(operators=ops,
                                                      dmmode="explicit"))
            print("%s |explicit-ED| = %.3e  ED[0,1] = %.6f" %
                  (tag, mx(ex, ed), ed[0, 1].real), flush=True)
            try:
                d = np.asarray(sc.get_correlation_matrix(operators=ops))
                print("%s default returned shape %s  |default-ED| = %.3e" %
                      (tag, d.shape, mx(d, ed)), flush=True)
            except Exception as e:
                print("%s default RAISED %s: %s" % (tag, type(e).__name__,
                      str(e)[:200]), flush=True)
        except Exception as e:
            print("%s ERROR %s: %s" % (tag, type(e).__name__, str(e)[:200]),
                  flush=True)
print("done")
```

Observed:

```
A v=python sector=none S_sub default = 0.720480  explicit = 0.720480  ED block = 0.720480  S_whole default = 0.000000
A v=python sector=none |S_sub default - ED| = 1.410e-14  |S_sub explicit - ED| = 6.661e-15
A v=python sector=Nf=3 S_sub default = 0.000000  explicit = 0.720480  ED block = 0.720480  S_whole default = 0.000000
A v=python sector=Nf=3 |S_sub default - ED| = 7.205e-01  |S_sub explicit - ED| = 4.774e-15
A v=3 sector=none S_sub default = 0.720480  explicit = 0.720480  ED block = 0.720480  S_whole default = 0.000000
A v=3 sector=none |S_sub default - ED| = 2.771e-12  |S_sub explicit - ED| = 2.769e-12
A v=3 sector=Nf=3 S_sub default = 0.000000  explicit = 0.720480  ED block = 0.720480  S_whole default = 0.000000
A v=3 sector=Nf=3 |S_sub default - ED| = 7.205e-01  |S_sub explicit - ED| = 6.106e-14
B v=python sector=none |explicit-ED| = 1.101e-13  ED[0,1] = -0.212588
B v=python sector=none default returned shape (6, 6)  |default-ED| = 1.101e-13
B v=python sector=Sz=0 |explicit-ED| = 1.538e-13  ED[0,1] = -0.212588
B v=python sector=Sz=0 default RAISED ValueError: Operator 'Cdag' not recognized for site type SpinHalfSite
B v=3 sector=none |explicit-ED| = 1.378e-13  ED[0,1] = -0.212588
B v=3 sector=none default returned shape (6, 6)  |default-ED| = 1.379e-13
B v=3 sector=Sz=0 |explicit-ED| = 7.498e-14  ED[0,1] = -0.212588
B v=3 sector=Sz=0 default RAISED ValueError: Chain::correlation_matrix_sweep: only spinless fermion sites (type 0) are covered
done
[exit 0, 6 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Both reviewers confirmed it and placed the defect in the dispatch of
`correlationentropy.py` only: `correlation_matrix_cdc` on `julia_live` and
`Chain::correlation_matrix_sweep` on v3 return what they are asked for, so
they are struck as defect sites. The commit that brought it in is `1b87543`,
not `921c0be`, which only moved the choice. The subset shape on
`julia_live` holds by reading only.

**Suggested fix**: One condition at lines 93 to 95: when `operators` is given and the mode was
picked automatically, take `"explicit"` in a sector and `"fast"` outside
one. An explicit `dmmode="full"` with `operators=` that are not the chain's
own C should raise naming `dmmode="explicit"`. A faster route for both
branches: when every X_a is a linear combination of the chain's own C_j,
compute the one-sweep matrix once and return `conj(U) M0 U^T`. The numbers
change for every sector call with `operators=`, so it needs a ledger entry,
and since the file dispatches into `julia_live` its tests should run there.

### 9. On a native spinful chain in a conserved sector, `get_correlation_matrix()` and `get_correlation_entropy()` with no keywords raise on `"python"` and v3, where the same calls answer without a sector

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `julia, v3-cpp` &middot; candidate found by the reviewer of `julia-2`

**Where**: `src/dmrgpy/entropytk/correlationentropy.py:334-339 (_default_dmmode's sector branch returns "full" with no site-type test)`

A second missing condition on the branch of the sector correlation-matrix
finding above. On `Spinful_Fermionic_Chain_Native` in an Nf sector the
default calls raise: on `"python"` "Operator 'Cdag' not recognized for site
type ElectronSite", on v3 "Chain::correlation_matrix_sweep: only spinless
fermion sites (type 0) are covered". In the same Nf=4 sector
`dmmode="explicit"` agrees with ED restricted to the sector (trace 4) to
3e-12 on `"python"` and 3e-11 on v3, and without a sector the same default
call answers to 5e-7. So turning on a sector breaks a default call that
works without one. It is the shape of 2026-08 finding 11 and 2026-09
finding 17, for a site type those fixes did not test.

**Expected**: the default taking a route that answers, as without a sector.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g06/g06_native_spinful_sector_default.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# New candidate: on a native spinful chain (Electron sites, site code 1) in
# a conserved sector, get_correlation_matrix() with no keywords takes
# dmmode="full" (_default_dmmode's sector branch has no spinless check,
# while its no-sector branch has one), and "full" is the session's
# <Cdag_i C_j> walk of spinless sites: v3's correlation_matrix_sweep throws
# "only spinless fermion sites", pyitensor's correlation_sweep asks an
# Electron site for "Cdag". Without the sector the same call takes "fast"
# and answers. Anchor: mode="ED" on the same chain (ED restricts to the
# sector), and dmmode="explicit".
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import fermionchain, cppext


def native(v, sector):
    # sector is None or a dict of set_conserved_sector keywords
    np.random.seed(4)
    fc = fermionchain.Spinful_Fermionic_Chain_Native(4, itensor_version=v,
                                                     maxm=40, nsweeps=10,
                                                     cutoff=1e-12)
    h = 0
    for i in range(3):
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + h.get_dagger()
    h = h + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(4))
    h = h + 0.3*fc.Nup[0] - 0.1*fc.Ndn[2]
    fc.set_hamiltonian(h)
    if sector:
        fc.set_conserved_sector(**sector)
    return fc


def mx(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return float("nan")
    return float(np.max(np.abs(a - b)))


versions = ["python"]
if cppext.available(3):
    versions.append(3)  # after python, in case a C++ error aborts
else:
    print("v3 = not compiled, skipped")

for v in versions:
    # the decision rests on Nf=4: there every entry <Cdag_a C_b> conserves
    # the sector, so dmmode="explicit" is a valid anchor; Nf=4,Sz=0 is
    # printed for information only (its up-down entries break Sz)
    for sector in (None, {"Nf": 4}, {"Nf": 4, "Sz": 0}):
        tag = "v=%s sector=%s" % (v, ",".join("%s=%d" % kv for kv in
                                      sorted(sector.items()))
                                  if sector else "none")
        try:
            fc = native(v, sector)
        except Exception as e:
            print("%s setup RAISED %s: %s" % (tag, type(e).__name__,
                  str(e)[:200]), flush=True)
            continue

        def get(label, **kw):
            try:
                return np.asarray(fc.get_correlation_matrix(**kw))
            except Exception as e:
                print("%s %s RAISED %s: %s" % (tag, label, type(e).__name__,
                      str(e)[:200]), flush=True)
                return None
        d = get("default")
        ed = get("mode=ED", mode="ED")
        ex = get("explicit", dmmode="explicit")
        if ed is not None:
            print("%s trace ED = %.6f" % (tag, np.trace(ed).real), flush=True)
        if ed is not None and ex is not None:
            print("%s |explicit-ED| = %.3e" % (tag, mx(ex, ed)), flush=True)
        if d is not None and ed is not None:
            print("%s default returned shape %s  |default-ED| = %.3e" %
                  (tag, d.shape, mx(d, ed)), flush=True)
        try:
            s = fc.get_correlation_entropy()
            print("%s get_correlation_entropy() default = %.6f" % (tag, s),
                  flush=True)
        except Exception as e:
            print("%s get_correlation_entropy() RAISED %s: %s" %
                  (tag, type(e).__name__, str(e)[:200]), flush=True)
print("done")
```

Observed:

```
v=python sector=none trace ED = 3.000000
v=python sector=none |explicit-ED| = 5.164e-07
v=python sector=none default returned shape (8, 8)  |default-ED| = 4.753e-07
v=python sector=none get_correlation_entropy() default = 0.690103
v=python sector=Nf=4 default RAISED ValueError: Operator 'Cdag' not recognized for site type ElectronSite
v=python sector=Nf=4 trace ED = 4.000000
v=python sector=Nf=4 |explicit-ED| = 2.889e-12
v=python sector=Nf=4 get_correlation_entropy() RAISED ValueError: Operator 'Cdag' not recognized for site type ElectronSite
v=python sector=Nf=4,Sz=0 default RAISED ValueError: Operator 'Cdag' not recognized for site type ElectronSite
v=python sector=Nf=4,Sz=0 mode=ED RAISED ValueError: this operator (Cdagup[0]*Cdn[0]) does not conserve ['Sz'], which the conserved sector {'Nf': 4, 'Sz': 0} fixes. Measure it after clearing the sector with set_conserved_sector().
v=python sector=Nf=4,Sz=0 explicit RAISED ValueError: conserved-sector mode (Nf=4, Sz=0): the term Adagup(1)*Fup(1)*Adn(1) changes the conserved quantum numbers by Sz=2 -- every operator built on a chain with a sector set must conserve them. (Terms are n
v=python sector=Nf=4,Sz=0 get_correlation_entropy() RAISED ValueError: Operator 'Cdag' not recognized for site type ElectronSite
v=3 sector=none trace ED = 3.000000
v=3 sector=none |explicit-ED| = 9.469e-07
v=3 sector=none default returned shape (8, 8)  |default-ED| = 8.797e-07
v=3 sector=none get_correlation_entropy() default = 0.690103
v=3 sector=Nf=4 default RAISED ValueError: Chain::correlation_matrix_sweep: only spinless fermion sites (type 0) are covered
v=3 sector=Nf=4 trace ED = 4.000000
v=3 sector=Nf=4 |explicit-ED| = 3.074e-11
v=3 sector=Nf=4 get_correlation_entropy() RAISED ValueError: Chain::correlation_matrix_sweep: only spinless fermion sites (type 0) are covered
v=3 sector=Nf=4,Sz=0 default RAISED ValueError: Chain::correlation_matrix_sweep: only spinless fermion sites (type 0) are covered
v=3 sector=Nf=4,Sz=0 mode=ED RAISED ValueError: this operator (Cdagup[0]*Cdn[0]) does not conserve ['Sz'], which the conserved sector {'Nf': 4, 'Sz': 0} fixes. Measure it after clearing the sector with set_conserved_sector().
v=3 sector=Nf=4,Sz=0 explicit RAISED ValueError: Conserved-sector mode (Nf=4, Sz=0): the term Adagup(1)*Fup(1)*Adn(1) changes the conserved quantum numbers by Nf=0, Sz=2 -- every operator built on a chain with a sector set must conserve them. (Terms
v=3 sector=Nf=4,Sz=0 get_correlation_entropy() RAISED ValueError: Chain::correlation_matrix_sweep: only spinless fermion sites (type 0) are covered
done
[exit 0, 6 s]
```

**Reviewer (CONFIRMED)**: Confirmed as stated, LOW, since it raises. The Nf=4, Sz=0 rows decide
nothing: there every route refuses an operator that changes Sz, which is
right.

**Suggested fix**: Write one predicate, "full applies to this state and these operators"
(spinless sites, no `operators=`, a session with the method), and use it in
both branches of `_default_dmmode`; the fix of the finding above then lands
in the same place, and the two branches cannot drift a fourth time.

### 10. `Thermal_Spin_Chain.get_gs()` ignores a temperature set after the first solve, and returns the T=2.0 state at T=0.2 (energy 0.785 off)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `solvers` &middot; candidate `solvers-3`

**Where**: `src/dmrgpy/thermal.py:52-53 (get_gs returns self.wf0 whenever computed_gs)`; `src/dmrgpy/thermal.py:65-68 (the T check a stored state skips)`; `src/dmrgpy/thermal.py:73-80 (anneal_method and anneal_options validation, skipped the same way)`; `src/dmrgpy/thermal.py:141-143 (set_hamiltonian, the only reset)`; `src/dmrgpy/thermal.py:121-140 (class attributes overridden on the instance)`

`get_gs()` returns its stored purified state whenever `computed_gs` is set,
and only `set_hamiltonian()` clears that flag. The class sets its
temperature as an attribute (`tc.T = 0.2`), and its own comments, tests and
examples do it this way. So a temperature sweep on one object returns the
state of the first temperature at every later one, silently. On a 3-site
Heisenberg chain the same object gives -0.206813 at T=0.2, the energy at
T=2.0, against the exact Boltzmann value -0.991668, identically on
`mode="ED"` and on `"python"` DMRG; a fresh object at T=0.2 gives the right
number. The same short circuit ignores `anneal_step`, `anneal_order`,
`anneal_method`, `anneal_options`, `tc.mode` and `tc.MBChain.maxm` set after
the first solve, and it skips the validation, so an invalid `anneal_method`
set after a solve returns the old state where a fresh object raises (the
counter-repro). The defect is older than the window; `c146d47` added
`anneal_method` and `anneal_options` inside the same short circuit.

**Expected**: the purified state at the temperature and the settings the object holds when `get_gs()` is called.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/solvers/solvers_3_thermal_T_cache.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: Thermal_Spin_Chain.get_gs() returns its stored purified state
# whenever computed_gs is set, and only set_hamiltonian() clears that flag,
# so a temperature assigned after the first get_gs() (tc.T = T, the way the
# class documents T is set) is ignored: every later get_gs() returns the
# first temperature's state. The same holds for anneal_method/
# anneal_options, added in c146d47.
import io, contextlib
from dmrgpy import spinchain, thermal

n = 3


def heis(c):
    h = 0
    for i in range(n-1):
        h = h + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
    return h


ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
ref.set_hamiltonian(heis(ref))
hm = ref.get_ED_obj().get_hamiltonian()
hm = hm.toarray() if hasattr(hm, "toarray") else np.array(hm)
ev = np.linalg.eigvalsh(hm)


def boltzmann(T):
    w = np.exp(-(ev-ev.min())/T)
    return float(np.sum(ev*w)/np.sum(w))


def energy(tc, h):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):  # anneal() prints every step
        wf = tc.get_gs()
    return float(np.real(wf.dot(h*wf)/wf.dot(wf)))


for mode in ("ED", "DMRG"):
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=2.0, mode=mode,
                                    itensor_version="python")
    tc.MBChain.maxm, tc.MBChain.nsweeps = 16, 10
    h = heis(tc)
    tc.set_hamiltonian(h)
    e_first = energy(tc, h)
    tc.T = 0.2  # a new temperature on the same object
    e_second = energy(tc, h)
    tc2 = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=0.2, mode=mode,
                                     itensor_version="python")
    tc2.MBChain.maxm, tc2.MBChain.nsweeps = 16, 10
    h2 = heis(tc2)
    tc2.set_hamiltonian(h2)
    e_fresh = energy(tc2, h2)
    print("mode=%s: exact <H>(T=2.0) = %.8f   first get_gs at T=2.0 = %.8f"
          % (mode, boltzmann(2.0), e_first))
    print("mode=%s: exact <H>(T=0.2) = %.8f   same object after tc.T=0.2 "
          "= %.8f   fresh object at T=0.2 = %.8f" %
          (mode, boltzmann(0.2), e_second, e_fresh))
    print("mode=%s: stale_after_T_change = %s   |second - exact(0.2)| = %.3e"
          % (mode, abs(e_second-e_first) < 1e-10,
             abs(e_second-boltzmann(0.2))))
```

Observed:

```
mode=ED: exact <H>(T=2.0) = -0.20681258   first get_gs at T=2.0 = -0.20681258
mode=ED: exact <H>(T=0.2) = -0.99166816   same object after tc.T=0.2 = -0.20681258   fresh object at T=0.2 = -0.99166811
mode=ED: stale_after_T_change = True   |second - exact(0.2)| = 7.849e-01
mode=DMRG: exact <H>(T=2.0) = -0.20681258   first get_gs at T=2.0 = -0.20681258
mode=DMRG: exact <H>(T=0.2) = -0.99166816   same object after tc.T=0.2 = -0.20681258   fresh object at T=0.2 = -0.99166811
mode=DMRG: stale_after_T_change = True   |second - exact(0.2)| = 7.849e-01
[exit 0, 5 s]
```

Reviewer's counter-repro, a setting other than T:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g11/g11_anneal_setting_ignored.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for solvers-3, the anneal_* face: a setting assigned after
# the first get_gs() is never read. An invalid anneal_method is the cleanest
# witness, since get_gs() validates it (thermal.py:73-77) only on a solve.
# ED only, 3 physical sites, deterministic.
import io, contextlib
from dmrgpy import thermal

np.random.seed(1)
n = 3


def heis(c):
    h = 0
    for i in range(n-1):
        h = h + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
    return h


def build(T):
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=T, mode="ED",
                                    itensor_version="python")
    h = heis(tc)
    tc.set_hamiltonian(h)
    return tc, h


def try_gs(tc, h):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            wf = tc.get_gs()
    except ValueError as err:
        return True, None, str(err)
    return False, float(np.real(wf.dot(h*wf)/wf.dot(wf))), ""


tc, h = build(0.5)
r0, e0, _ = try_gs(tc, h)
tc.anneal_method = "not_a_method"
r1, e1, _ = try_gs(tc, h)
tf, hf = build(0.5)
tf.anneal_method = "not_a_method"
r2, e2, msg = try_gs(tf, hf)
print("first get_gs (default method): raised = %s  E = %.10f" % (r0, e0))
print("same object, invalid anneal_method after first get_gs: raised = %s"
      "  E = %s" % (r1, "%.10f" % e1 if e1 is not None else None))
print("fresh object, invalid anneal_method: raised = %s  (%s)" % (r2, msg))
print("anneal_setting_ignored = %s" % ((not r1) and r2 and (e1 == e0)))
```

Observed:

```
first get_gs (default method): raised = False  E = -0.7694594961
same object, invalid anneal_method after first get_gs: raised = False  E = -0.7694594961
fresh object, invalid anneal_method: raised = True  (Thermal_Spin_Chain: anneal_method must be None, 'taylor' or 'implicit', got 'not_a_method')
anneal_setting_ignored = True
[exit 0, 3 s]
```

**Reviewer (CONFIRMED)**: Nothing documents the stored state as frozen, and the class comments
describe setting T on the instance. The anchor is the Boltzmann average from
ED. The counter-repro measures the sub-claim on the `anneal_*` settings: an
invalid method set after a solve is not refused, while a fresh object
raises.

**Suggested fix**: Store, next to `wf0`, the key the state was purified under, and return the
stored state only when the current key equals it, the repository's own
idiom (`groundstate.solver_key`, `stored_answer_holds`). The key holds
`float(self.T)`, `anneal_step`, `anneal_order`, `anneal_method`, a frozen
copy of `anneal_options` (the dict is mutable), `self.mode`,
`MBChain.itensor_version` and `groundstate.solver_key(self.MBChain)`. Run
the T check and the `anneal_*` validation before the key comparison. A
property setter on T alone is not enough, since it misses the other
settings. Pin it with a T sweep on one object against fresh objects and the
Boltzmann value, on `mode="ED"` and `"python"`.

### 11. The default non-Hermitian `gs_energy()` on DMRG returns a level up to 1e-4 above the smallest real part, with Im E off by 0.20 and `<Sz_total>` off by 1 or 2, in 3 of 4 seeded runs

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `solvers` &middot; candidate `solvers-1`

**Where**: `src/dmrgpy/nhlevels.py:216-221 (the tie window 1e-4*(c+|E_min-e_id|) and the pick by smallest residual)`; `src/dmrgpy/nhlevels.py:157 (_accepted, sorted by (Re, Im))`; `src/dmrgpy/nhlevels.py:275 (the certificate checks the eigen-residual only)`; `src/dmrgpy/groundstate.py:970-973 (default dispatch to gs_energy_imaginary_time)`; `src/dmrgpy/manybodychain.py:180 (nh_gs_method default)`

Since `94367c8` a non-Hermitian `gs_energy()` on v3 and `"python"` takes the
imaginary-time route by default. Its `ground_pair` treats every accepted
level within `1e-4*(c+|E_min-e_id|)` of the smallest real part as tied, and
returns the one of smallest Ritz residual (`ef9a885`). The window is wider
than the level spacing of a chain with a weak field. On the 6-site
PT-symmetric XX chain with a field `5e-5*Sz_total` the exact lowest level is
`E0-eps` with `<Sz_total> = -1`, followed by the pair `E0 +- 0.2012i`
(`<Sz_total> = 0`) and by `E0+eps` (`+1`). The route returned `E0 +- 0.2012i`
or `E0+eps` in 3 of 4 seeded runs on `"python"` and in 3 of 3 on v3: the
real part is high by 5e-5 to 1e-4, the imaginary part off by 0.20 and the
magnetization off by 1 or 2. Each returned pair is an exact eigenpair of
another level, so the certificate passes (3e-7 to 6e-13) and nothing warns.
The counter-repro shows the minimum was accepted first in every run: the
window alone picks the wrong level. On v3 the start is not seeded from
numpy, and the counter-repro's three v3 runs returned the right level, so on
v3 it depends on the draw. NH-DMRG, the old default, returned `E0 +- 0.2012i`
in every run of the counter-repro too, and so does `mode="ED"` (the
finding on the ED sort below), so neither is an independent anchor
here; the anchor is `np.linalg.eig` of the dense matrix.

**Expected**: the eigenvalue of smallest real part, `E0-eps` with `<Sz_total> = -1`, whenever the Ritz values resolve the spacing (here they are 7.8e-7 from the exact levels).

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/solvers/solvers_1_nh_tie_window.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: the default non-Hermitian gs_energy() on DMRG
# (nh_gs_method="imaginary_time", nhlevels.ground_pair) treats every
# accepted Ritz value whose real part is within 1e-4*(c+|E-e_id|) of the
# smallest as "tied" and returns the one of smallest residual, so when two
# or more DISTINCT levels lie inside that window it can return a level that
# is not the one of smallest real part (different Im E, different sector),
# and the certificate passes because it is an exact eigenpair. NH-DMRG
# (nh_gs_method="nhdmrg", the old default) and the exact spectrum give the
# level of smallest real part.
import random, time
from dmrgpy import spinchain, cppext

G = 0.3


def pt_chain(n, eps, version):
    """PT-symmetric XX chain, i G (-1)^i Sz_i, plus a tiny uniform real
    field eps*Sz_total: the four levels sharing Re E0 at eps=0 (Sz=-1, 0,
    0, +1) split to E0-eps, E0+-i*gamma, E0+eps."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 16
    sc.nsweeps = 12
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        h = h + 1j*G*(-1)**i*sc.Sz[i]
    obs = 0
    for i in range(n): obs = obs + sc.Sz[i]
    h = h + eps*obs
    sc.set_hamiltonian(h)
    return sc, obs


def impurity_chain(n, eps, version):
    """Decoupled impurity: -sum_{i>=1} Sz_i + (eps + 0.5i) Sz_0, exact
    levels -(n-1)/2 -+ (eps+0.5i)/2; smallest real part has Sz_0=-1/2."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 8
    sc.nsweeps = 12
    h = 0
    for i in range(1, n):
        h = h - 1.0*sc.Sz[i]
    h = h + (eps + 0.5j)*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc, sc.Sz[0]


def exact(sc, obs):
    ed = sc.get_ED_obj()
    hm = ed.get_hamiltonian()
    hm = hm.toarray() if hasattr(hm, "toarray") else np.array(hm)
    om = ed.get_operator(obs)
    om = om.toarray() if hasattr(om, "toarray") else np.array(om)
    ev, vr = np.linalg.eig(hm)
    idx = np.argsort(ev.real, kind="stable")
    out = []
    for k in idx[:5]:
        v = vr[:, k]
        out.append((ev[k], float(np.real(np.vdot(v, om@v)/np.vdot(v, v)))))
    return out


def report(tag, build, n, eps, version, seeds, run_nhdmrg=True, levels=False):
    sc, obs = build(n, eps, version)
    ex = exact(sc, obs)
    emin, omin = ex[0]
    window = 1e-4*(1.0 + abs(emin))
    print("== %s n=%d eps=%g backend=%s" % (tag, n, eps, version))
    for e, o in ex[:4]:
        print("   exact level E = %.8f%+.8fj  <obs> = %+.4f  Re-ReMin = %.2e"
              % (e.real, e.imag, o, e.real-emin.real))
    print("   tie window of ground_pair = %.2e" % window)
    print("   ED gs_energy(mode='ED') = %r" % (complex(sc.gs_energy(mode="ED")),))
    nwrong = 0
    for s in seeds:
        if s is not None:
            np.random.seed(s); random.seed(s)
        sc, obs = build(n, eps, version)
        t0 = time.time()
        e = complex(sc.gs_energy())
        o = float(np.real(sc.vev(obs)))
        res = getattr(sc, "_nh_last_residual", None)
        wrong = abs(e - emin) > 1e-6
        nwrong += int(wrong)
        print("   imaginary_time seed=%s E = %.8f%+.8fj  |E-Emin| = %.2e  "
              "Re-ReMin = %+.2e  <obs> = %+.4f  certificate = %s  "
              "wrong_level = %s  (%.1f s)"
              % (s, e.real, e.imag, abs(e-emin), e.real-emin.real, o,
                 None if res is None else "%.1e (tol %.0e)" % res, wrong,
                 time.time()-t0))
    print("   SUMMARY %s eps=%g backend=%s: wrong_level in %d of %d runs"
          % (tag, eps, version, nwrong, len(seeds)))
    if levels:
        np.random.seed(0); random.seed(0)
        sc, obs = build(n, eps, version)
        es, rs = sc.get_nh_levels()
        for e, r in list(zip(es, rs))[:6]:
            print("   get_nh_levels: E = %.8f%+.8fj  rel.residual = %.1e"
                  % (e.real, e.imag, r))
    if run_nhdmrg:
        np.random.seed(1); random.seed(1)
        sc, obs = build(n, eps, version)
        sc.nh_gs_method = "nhdmrg"
        t0 = time.time()
        e = complex(sc.gs_energy())
        o = float(np.real(sc.vev(obs)))
        print("   nhdmrg (old default) E = %.8f%+.8fj  |E-Emin| = %.2e  "
              "<obs> = %+.4f  (%.1f s)" % (e.real, e.imag, abs(e-emin), o,
                                            time.time()-t0))


# PT chain, near-tied (inside the window) and control (outside it)
report("PT_XX_field", pt_chain, 6, 5e-5, "python", [0, 1, 2, 3], levels=True)
report("PT_XX_field_control", pt_chain, 6, 1e-3, "python", [0, 1],
       run_nhdmrg=False)
# decoupled impurity: two levels, Re split by eps, Im split by 0.5
report("impurity", impurity_chain, 4, 5e-5, "python", [0, 1, 2])
if cppext.available(3):
    report("PT_XX_field", pt_chain, 6, 5e-5, 3, [None, None, None],
           run_nhdmrg=False)
else:
    print("v3 extension not available: skipped")
```

Observed:

```
== PT_XX_field n=6 eps=5e-05 backend=python
   exact level E = -1.39617641-0.00000000j  <obs> = -1.0000  Re-ReMin = 0.00e+00
   exact level E = -1.39612641-0.20120744j  <obs> = -0.0000  Re-ReMin = 5.00e-05
   exact level E = -1.39612641+0.20120744j  <obs> = -0.0000  Re-ReMin = 5.00e-05
   exact level E = -1.39607641-0.00000000j  <obs> = +1.0000  Re-ReMin = 1.00e-04
   tie window of ground_pair = 2.40e-04
   ED gs_energy(mode='ED') = (-1.3961264106354192-0.20120744009904132j)
   imaginary_time seed=0 E = -1.39612641+0.20120744j  |E-Emin| = 2.01e-01  Re-ReMin = +5.00e-05  <obs> = -0.0000  certificate = 1.3e-07 (tol 1e-04)  wrong_level = True  (5.4 s)
   imaginary_time seed=1 E = -1.39607641+0.00000000j  |E-Emin| = 1.00e-04  Re-ReMin = +1.00e-04  <obs> = +1.0000  certificate = 6.5e-10 (tol 1e-04)  wrong_level = True  (7.2 s)
   imaginary_time seed=2 E = -1.39612641-0.20120744j  |E-Emin| = 2.01e-01  Re-ReMin = +5.00e-05  <obs> = -0.0000  certificate = 1.4e-12 (tol 1e-04)  wrong_level = True  (5.8 s)
   imaginary_time seed=3 E = -1.39617641+0.00000000j  |E-Emin| = 6.49e-13  Re-ReMin = -6.36e-13  <obs> = -1.0000  certificate = 4.4e-10 (tol 1e-04)  wrong_level = False  (5.5 s)
   SUMMARY PT_XX_field eps=5e-05 backend=python: wrong_level in 3 of 4 runs
   get_nh_levels: E = -1.39617563+0.00000193j  rel.residual = 1.1e-05
   get_nh_levels: E = -1.39612642+0.20120745j  rel.residual = 2.9e-07
   get_nh_levels: E = -1.39612640-0.20120744j  rel.residual = 3.5e-07
   get_nh_levels: E = -1.39607803-0.00000160j  rel.residual = 1.2e-05
   nhdmrg (old default) E = -1.39612641-0.20120744j  |E-Emin| = 2.01e-01  <obs> = -0.0000  (0.2 s)
== PT_XX_field_control n=6 eps=0.001 backend=python
   exact level E = -1.39712641+0.00000000j  <obs> = -1.0000  Re-ReMin = 0.00e+00
   exact level E = -1.39612641+0.20120744j  <obs> = -0.0000  Re-ReMin = 1.00e-03
   exact level E = -1.39612641-0.20120744j  <obs> = -0.0000  Re-ReMin = 1.00e-03
   exact level E = -1.39512641+0.00000000j  <obs> = +1.0000  Re-ReMin = 2.00e-03
   tie window of ground_pair = 2.40e-04
   ED gs_energy(mode='ED') = (-1.3971264106354193+4.621749220398127e-16j)
   imaginary_time seed=0 E = -1.39712643-0.00000000j  |E-Emin| = 2.26e-08  Re-ReMin = -2.25e-08  <obs> = -1.0000  certificate = 1.6e-06 (tol 1e-04)  wrong_level = False  (6.4 s)
   imaginary_time seed=1 E = -1.39712641-0.00000000j  |E-Emin| = 2.05e-14  Re-ReMin = +2.00e-14  <obs> = -1.0000  certificate = 1.1e-12 (tol 1e-04)  wrong_level = False  (5.8 s)
   SUMMARY PT_XX_field_control eps=0.001 backend=python: wrong_level in 0 of 2 runs
== impurity n=4 eps=5e-05 backend=python
   exact level E = -1.50002500-0.25000000j  <obs> = -0.5000  Re-ReMin = 0.00e+00
   exact level E = -1.49997500+0.25000000j  <obs> = +0.5000  Re-ReMin = 5.00e-05
   exact level E = -0.50002500-0.25000000j  <obs> = -0.5000  Re-ReMin = 1.00e+00
   exact level E = -0.50002500-0.25000000j  <obs> = -0.5000  Re-ReMin = 1.00e+00
   tie window of ground_pair = 2.52e-04
   ED gs_energy(mode='ED') = (-1.500025-0.25j)
   imaginary_time seed=0 E = -1.50002500-0.25000000j  |E-Emin| = 8.95e-16  Re-ReMin = +8.88e-16  <obs> = -0.5000  certificate = 2.4e-12 (tol 1e-04)  wrong_level = False  (2.8 s)
   imaginary_time seed=1 E = -1.50002500-0.25000000j  |E-Emin| = 9.09e-16  Re-ReMin = +8.88e-16  <obs> = -0.5000  certificate = 6.0e-13 (tol 1e-04)  wrong_level = False  (2.7 s)
   imaginary_time seed=2 E = -1.50002500-0.25000000j  |E-Emin| = 7.45e-16  Re-ReMin = +6.66e-16  <obs> = -0.5000  certificate = 8.5e-13 (tol 1e-04)  wrong_level = False  (2.7 s)
   SUMMARY impurity eps=5e-05 backend=python: wrong_level in 0 of 3 runs
   nhdmrg (old default) E = -1.50002500-0.25000000j  |E-Emin| = 9.39e-16  <obs> = -0.5000  (0.1 s)
== PT_XX_field n=6 eps=5e-05 backend=3
   exact level E = -1.39617641-0.00000000j  <obs> = -1.0000  Re-ReMin = 0.00e+00
   exact level E = -1.39612641-0.20120744j  <obs> = -0.0000  Re-ReMin = 5.00e-05
   exact level E = -1.39612641+0.20120744j  <obs> = -0.0000  Re-ReMin = 5.00e-05
   exact level E = -1.39607641-0.00000000j  <obs> = +1.0000  Re-ReMin = 1.00e-04
   tie window of ground_pair = 2.40e-04
   ED gs_energy(mode='ED') = (-1.3961264106354192-0.20120744009904132j)
   imaginary_time seed=None E = -1.39612641-0.20120744j  |E-Emin| = 2.01e-01  Re-ReMin = +5.00e-05  <obs> = -0.0000  certificate = 6.2e-13 (tol 1e-04)  wrong_level = True  (2.5 s)
   imaginary_time seed=None E = -1.39612641+0.20120744j  |E-Emin| = 2.01e-01  Re-ReMin = +5.00e-05  <obs> = -0.0000  certificate = 3.1e-07 (tol 1e-04)  wrong_level = True  (2.6 s)
   imaginary_time seed=None E = -1.39612641-0.20120744j  |E-Emin| = 2.01e-01  Re-ReMin = +5.00e-05  <obs> = +0.0000  certificate = 3.0e-07 (tol 1e-04)  wrong_level = True  (2.7 s)
   SUMMARY PT_XX_field eps=5e-05 backend=3: wrong_level in 3 of 3 runs
[exit 0, 58 s]
```

Reviewer's counter-repro, the mechanism and the candidate rules:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g10/g10_solvers1_mechanism.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Reviewer counter-repro for solvers-1. For every default gs_energy() run on
# the 6-site PT XX chain in a field eps=5e-5, it records the accepted Ritz
# pairs that ground_pair() selects from (by wrapping nhlevels._accepted) and
# prints, per run:
#   min_in_acc    : whether the exact level of smallest real part was accepted
#   acc0_is_min   : whether the first accepted pair (smallest Re) is that level
#   returned_wrong: whether gs_energy() returned another level
# and what four selection rules pick from the same accepted list:
#   current   : the code's rule (window 1e-4*(c+|E-e_id|), then min residual)
#   strict    : acc[0], the smallest real part
#   hunter_fN : tie if dRe <= N*(r_i + r_0) (absolute residuals), then min residual
#   proposed  : smallest Re among pairs with relative residual < 1e-4 (the
#               certificate tol), ties only within 1e-6*(c+|E-e_id|), min residual
# It also runs NH-DMRG (nh_gs_method="nhdmrg") from several seeds.
import random
from dmrgpy import spinchain, cppext, nhlevels

G = 0.3
EPS = 5e-5


def pt_chain(n, eps, version):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 16
    sc.nsweeps = 12
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        h = h + 1j*G*(-1)**i*sc.Sz[i]
    obs = 0
    for i in range(n): obs = obs + sc.Sz[i]
    h = h + eps*obs
    sc.set_hamiltonian(h)
    return sc, obs


def exact_levels(sc, obs):
    ed = sc.get_ED_obj()
    hm = ed.get_hamiltonian()
    hm = hm.toarray() if hasattr(hm, "toarray") else np.array(hm)
    om = ed.get_operator(obs)
    om = om.toarray() if hasattr(om, "toarray") else np.array(om)
    ev, vr = np.linalg.eig(hm)
    idx = np.argsort(ev.real, kind="stable")
    out = []
    for k in idx:
        v = vr[:, k]
        out.append((complex(ev[k]),
                    float(np.real(np.vdot(v, om@v)/np.vdot(v, v)))))
    return out


REC = []
_orig = nhlevels._accepted


def _spy(pairs, H, tol):
    out = _orig(pairs, H, tol)
    REC.append([(complex(p[0]), float(p[1])) for p in out])
    return out


nhlevels._accepted = _spy


def nearest(ex, z):
    k = int(np.argmin([abs(e - z) for e, _ in ex]))
    return k


def select(acc, rscale, e_id, rule):
    a0 = acc[0]
    def rel(p): return p[1]/(rscale + abs(p[0] - e_id))
    if rule == "current":
        w = 1e-4*(rscale + abs(a0[0] - e_id))
        tied = [p for p in acc if p[0].real <= a0[0].real + w]
        return min(tied, key=lambda p: p[1])
    if rule == "strict":
        return a0
    if rule.startswith("hunter_f"):
        f = float(rule[len("hunter_f"):])
        tied = [p for p in acc if p[0].real - a0[0].real <= f*(p[1] + a0[1])]
        return min(tied, key=lambda p: p[1])
    if rule == "proposed":
        cands = [p for p in acc if rel(p) < 1e-4]
        if not cands: return select(acc, rscale, e_id, "current")
        f0 = cands[0]
        w = 1e-6*(rscale + abs(f0[0] - e_id))
        tied = [p for p in cands if p[0].real - f0[0].real <= w]
        return min(tied, key=lambda p: p[1])
    raise ValueError(rule)


RULES = ["current", "strict", "hunter_f1", "hunter_f2", "hunter_f3", "proposed"]


def run_backend(version, seeds):
    sc, obs = pt_chain(6, EPS, version)
    ex = exact_levels(sc, obs)
    emin = ex[0][0]
    rscale, e_id = nhlevels._scale(sc.hamiltonian)
    print("== backend=%s eps=%g  exact Emin = %.8f%+.8fj  <obs> = %+.3f"
          % (version, EPS, emin.real, emin.imag, ex[0][1]))
    # the hunter's script made this call before its seeded runs; kept so
    # that the seeded runs start from the same module state
    print("  ED gs_energy(mode='ED') = %r" % (complex(sc.gs_energy(mode="ED")),))
    nwrong = 0
    counts = dict((r, 0) for r in RULES)
    mech_window = 0
    mech_lost = 0
    for s in seeds:
        if s is not None:
            np.random.seed(s); random.seed(s)
        sc, obs = pt_chain(6, EPS, version)
        del REC[:]
        e = complex(sc.gs_energy())
        o = float(np.real(sc.vev(obs)))
        acc = REC[0] if REC else []
        wrong = abs(e - emin) > 1e-6
        nwrong += int(wrong)
        ks = [nearest(ex, p[0]) for p in acc]
        min_in_acc = 0 in ks
        acc0_is_min = bool(acc) and ks[0] == 0
        if wrong and min_in_acc: mech_window += 1
        if wrong and not min_in_acc: mech_lost += 1
        print("  run seed=%s returned E = %.8f%+.8fj <obs> = %+.3f "
              "returned_wrong = %s min_in_acc = %s acc0_is_min = %s n_acc = %d"
              % (s, e.real, e.imag, o, wrong, min_in_acc, acc0_is_min, len(acc)))
        for p, k in list(zip(acc, ks))[:6]:
            print("     accepted: Re-ReMin = %+.3e Im = %+.5f rel.res = %.2e "
                  "exact level #%d <obs> = %+.3f"
                  % (p[0].real - emin.real, p[0].imag,
                     p[1]/(rscale + abs(p[0] - e_id)), k, ex[k][1]))
        if acc:
            line = []
            for r in RULES:
                q = select(acc, rscale, e_id, r)
                w = nearest(ex, q[0]) != 0
                counts[r] += int(w)
                line.append("%s=%s" % (r, "wrong" if w else "ok"))
            print("     rules: " + " ".join(line))
    print("  SUMMARY backend=%s default gs_energy wrong in %d of %d runs; "
          "wrong with the minimum accepted (window) %d, wrong with the "
          "minimum not accepted (lost) %d"
          % (version, nwrong, len(seeds), mech_window, mech_lost))
    print("  RULES backend=%s " % version
          + " ".join("%s_wrong=%d/%d" % (r, counts[r], len(seeds)) for r in RULES))
    return emin


def run_nhdmrg(version, seeds, emin):
    nw = 0
    for s in seeds:
        if s is not None:
            np.random.seed(s); random.seed(s)
        sc, obs = pt_chain(6, EPS, version)
        sc.nh_gs_method = "nhdmrg"
        e = complex(sc.gs_energy())
        o = float(np.real(sc.vev(obs)))
        w = abs(e - emin) > 1e-6
        nw += int(w)
        print("  nhdmrg seed=%s E = %.8f%+.8fj <obs> = %+.3f wrong = %s"
              % (s, e.real, e.imag, o, w))
    print("  NHDMRG backend=%s wrong in %d of %d runs" % (version, nw, len(seeds)))


emin = run_backend("python", [0, 1, 2, 3, 4, 5])
run_nhdmrg("python", [1, 2, 3, 4], emin)
if cppext.available(3):
    emin3 = run_backend(3, [None, None, None])
    run_nhdmrg(3, [None, None], emin3)
else:
    print("v3 extension not available: skipped")
```

Observed:

```
== backend=python eps=5e-05  exact Emin = -1.39617641-0.00000000j  <obs> = -1.000
  ED gs_energy(mode='ED') = (-1.3961264106354192-0.20120744009904132j)
  run seed=0 returned E = -1.39612641+0.20120744j <obs> = -0.000 returned_wrong = True min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = -1.371e-06 Im = -0.00000 rel.res = 2.72e-06 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.002e-05 Im = -0.20121 rel.res = 2.18e-07 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +5.002e-05 Im = +0.20121 rel.res = 9.93e-08 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +1.005e-04 Im = +0.00000 rel.res = 1.75e-05 exact level #3 <obs> = +1.000
     rules: current=wrong strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=1 returned E = -1.39607641+0.00000000j <obs> = +1.000 returned_wrong = True min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = -2.210e-12 Im = +0.00000 rel.res = 1.95e-07 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 6.18e-09 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 1.24e-08 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = -0.00000 rel.res = 0.00e+00 exact level #3 <obs> = +1.000
     rules: current=wrong strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=2 returned E = -1.39612641-0.20120744j <obs> = -0.000 returned_wrong = True min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = -6.082e-12 Im = -0.00000 rel.res = 3.22e-07 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 0.00e+00 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 0.00e+00 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = +0.00000 rel.res = 0.00e+00 exact level #3 <obs> = +1.000
     rules: current=wrong strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=3 returned E = -1.39617641+0.00000000j <obs> = -1.000 returned_wrong = False min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = +2.175e-12 Im = -0.00000 rel.res = 0.00e+00 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 0.00e+00 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 0.00e+00 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = +0.00000 rel.res = 4.40e-07 exact level #3 <obs> = +1.000
     rules: current=ok strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=4 returned E = -1.39617641+0.00000000j <obs> = -1.000 returned_wrong = False min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = +1.919e-12 Im = -0.00000 rel.res = 0.00e+00 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 0.00e+00 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 8.74e-09 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = +0.00000 rel.res = 1.31e-07 exact level #3 <obs> = +1.000
     rules: current=ok strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=5 returned E = -1.39612641-0.20120744j <obs> = -0.000 returned_wrong = True min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = -1.732e-07 Im = +0.00000 rel.res = 2.40e-06 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 5.73e-08 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 5.17e-08 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +1.002e-04 Im = -0.00000 rel.res = 2.51e-06 exact level #3 <obs> = +1.000
     rules: current=wrong strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  SUMMARY backend=python default gs_energy wrong in 4 of 6 runs; wrong with the minimum accepted (window) 4, wrong with the minimum not accepted (lost) 0
  RULES backend=python current_wrong=4/6 strict_wrong=0/6 hunter_f1_wrong=0/6 hunter_f2_wrong=0/6 hunter_f3_wrong=0/6 proposed_wrong=0/6
  nhdmrg seed=1 E = -1.39612641-0.20120744j <obs> = -0.000 wrong = True
  nhdmrg seed=2 E = -1.39612641-0.20120744j <obs> = -0.000 wrong = True
  nhdmrg seed=3 E = -1.39612641-0.20120744j <obs> = -0.000 wrong = True
  nhdmrg seed=4 E = -1.39612641+0.20120744j <obs> = -0.000 wrong = True
  NHDMRG backend=python wrong in 4 of 4 runs
== backend=3 eps=5e-05  exact Emin = -1.39617641-0.00000000j  <obs> = -1.000
  ED gs_energy(mode='ED') = (-1.3961264106354192-0.20120744009904132j)
  run seed=None returned E = -1.39617641+0.00000002j <obs> = -1.000 returned_wrong = False min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = +3.533e-10 Im = +0.00000 rel.res = 0.00e+00 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 1.07e-08 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 1.24e-08 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = -0.00000 rel.res = 0.00e+00 exact level #3 <obs> = +1.000
     rules: current=ok strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=None returned E = -1.39617641-0.00000000j <obs> = -1.000 returned_wrong = False min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = -1.354e-10 Im = -0.00000 rel.res = 0.00e+00 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 6.18e-09 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 0.00e+00 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = +0.00000 rel.res = 6.98e-08 exact level #3 <obs> = +1.000
     rules: current=ok strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  run seed=None returned E = -1.39617641+0.00000000j <obs> = -1.000 returned_wrong = False min_in_acc = True acc0_is_min = True n_acc = 4
     accepted: Re-ReMin = +2.485e-12 Im = -0.00000 rel.res = 0.00e+00 exact level #0 <obs> = -1.000
     accepted: Re-ReMin = +5.000e-05 Im = +0.20121 rel.res = 1.07e-08 exact level #2 <obs> = -0.000
     accepted: Re-ReMin = +5.000e-05 Im = -0.20121 rel.res = 0.00e+00 exact level #1 <obs> = -0.000
     accepted: Re-ReMin = +1.000e-04 Im = +0.00000 rel.res = 1.65e-06 exact level #3 <obs> = +1.000
     rules: current=ok strict=ok hunter_f1=ok hunter_f2=ok hunter_f3=ok proposed=ok
  SUMMARY backend=3 default gs_energy wrong in 0 of 3 runs; wrong with the minimum accepted (window) 0, wrong with the minimum not accepted (lost) 0
  RULES backend=3 current_wrong=0/3 strict_wrong=0/3 hunter_f1_wrong=0/3 hunter_f2_wrong=0/3 hunter_f3_wrong=0/3 proposed_wrong=0/3
  nhdmrg seed=None E = -1.39612641-0.20120744j <obs> = -0.000 wrong = True
  nhdmrg seed=None E = -1.39612641-0.20120744j <obs> = -0.000 wrong = True
  NHDMRG backend=3 wrong in 2 of 2 runs
[exit 0, 78 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Two sub-claims were struck. The impurity chain shows no defect (0 of 3 runs
wrong), and the window came in with `ef9a885`, after `94367c8` made the
route the default, not with `4ed4891`. The pick is by residual and seed
dependent, so it returns another level often, not always. The 2026-11-15
ledger entry still says the smallest real part is taken, which has been
stale since `ef9a885`.

**Suggested fix**: The counter-repro scores five rules on the accepted lists of every run:
the plain smallest real part, the hunter's three residual-based windows and
the reviewer's rule (keep the pairs whose relative residual is below the
certificate tolerance, take the smallest real part among them, and break
only true ties by residual, with a window at the accuracy of the Ritz
values, NH-DMRG's `1e-6*(Re spread)+100*eps*max|E|`). All five are right in
every run, 6 of 6 on `"python"` and 3 of 3 on v3, where the current window
is wrong in 4 of 6. The reviewer's objection to the hunter's windows, that
with the printed residuals a factor of 2 already ties the 5e-5 split again,
is arithmetic on one run's residuals, not a measurement. By reading, the
reviewer's rule also still skips the poorly converged member of the 16-site
exact tie that motivated `ef9a885` (its residual of 5.7e-4 is above 1e-4).
When a skipped pair lies below the chosen one by more than the tie window,
take it and let the certificate decide, or warn. Fix it in one cluster with
the ED sort below, and anchor the test on `np.linalg.eig`, not on
`get_excited(mode="ED")`.

### 12. ED orders non-Hermitian eigenvalues by Re + 1e-3 Im, so `gs_energy(mode="ED")` returns a level 1e-4 above the smallest real part, with `<Sz_0>` = +0.5 for -0.5

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `solvers` &middot; candidate `solvers-2`

**Where**: `src/dmrgpy/algebra/algebra.py:381-388 (sorteigen: Re(E - min Re) + 1e-3*Im E)`; `src/dmrgpy/algebra/algebra.py:323, 347, 356 (lowest_eigenvalues, lowest_states)`; `src/dmrgpy/edtk/edchain.py:204-206 and 214-225 (ED gs_energy and get_gs_array)`; `src/dmrgpy/excited.py:250-251 (the same key on DMRG excited states; by reading)`

`algebra.sorteigen` orders non-Hermitian eigenvalues by
`Re(E - min Re) + 1e-3*Im E`, so `gs_energy`, `get_gs` and `vev` with
`mode="ED"` return a level whose real part lies above the minimum by less
than 1e-3 times its deficit in Im. On a 4-spin decoupled impurity
(eps = 1e-4, g = 0.5) ED returns -1.49995-0.25i with `<Sz_0> = +0.5`, for the
true -1.50005+0.25i with -0.5. On the 6-site PT-symmetric XX chain in a field
5e-5 it returns `E0-0.2012i` with `<Sz_total> = 0`, for `E0-5e-5` with -1. The
controls at eps = 1e-2 and 1e-3 are right. The defect is older than the
window; it matters here because ED is the reference the non-Hermitian DMRG
tests anchor on.

**Expected**: the eigenvalue of smallest real part, with Im only breaking an exact tie.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/solvers/solvers_2_ed_nh_sort_tilt.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: the ED "ground state" of a non-Hermitian Hamiltonian
# (algebra.sorteigen, used by lowest_eigenvalues/lowest_states, hence by
# gs_energy(mode="ED"), get_gs(mode="ED"), vev(mode="ED") and
# biorthogonal_ground_state) sorts by Re(E-Emin) + 1e-3*Im(E), not by the
# real part, so a level whose real part is HIGHER than the minimum by less
# than 1e-3 times the imaginary-part difference is returned as the ground
# state. Every DMRG route documents "the eigenvalue of smallest real part".
import random
from dmrgpy import spinchain
from dmrgpy.algebra import algebra


def impurity(n, eps, g):
    """-sum_{i>=1} Sz_i + (eps - i g) Sz_0: levels -(n-1)/2 -+ (eps - i g)/2.
    Smallest real part: Sz_0 = -1/2, E = -(n-1)/2 - eps/2 + i g/2."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(1, n): h = h - 1.0*sc.Sz[i]
    h = h + (eps - 1j*g)*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc, sc.Sz[0]


def pt_field(n, eps, g=0.3):
    """PT XX chain + eps*Sz_total: levels E0-eps (Sz=-1), E0+-i gamma
    (Sz=0), E0+eps (Sz=+1)."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(n-1): h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n): h = h + 1j*g*(-1)**i*sc.Sz[i]
    obs = 0
    for i in range(n): obs = obs + sc.Sz[i]
    h = h + eps*obs
    sc.set_hamiltonian(h)
    return sc, obs


def check(tag, sc, obs):
    ed = sc.get_ED_obj()
    hm = ed.get_hamiltonian()
    hm = hm.toarray() if hasattr(hm, "toarray") else np.array(hm)
    om = ed.get_operator(obs)
    om = om.toarray() if hasattr(om, "toarray") else np.array(om)
    ev, vr = np.linalg.eig(hm)
    k = int(np.argmin(ev.real))
    v = vr[:, k]
    o_true = float(np.real(np.vdot(v, om@v)/np.vdot(v, v)))
    e_ed = complex(sc.gs_energy(mode="ED"))
    o_ed = float(np.real(sc.vev(obs, mode="ED")))
    from scipy.sparse import csr_matrix
    try:
        e_bi = complex(algebra.biorthogonal_ground_state(csr_matrix(hm))[0])
    except Exception as err:  # context only, not the verdict
        print("%s: biorthogonal_ground_state raised %r" % (tag, err))
        e_bi = complex(np.nan)
    wrong = e_ed.real > ev[k].real + 1e-12
    print("%s: true smallest-Re E = %.8f%+.8fj <obs> = %+.4f" %
          (tag, ev[k].real, ev[k].imag, o_true))
    print("%s: ED gs_energy = %.8f%+.8fj <obs>(ED vev) = %+.4f  "
          "Re(ED)-ReMin = %+.2e  biorthogonal_ground_state E = %.8f%+.8fj  "
          "ED_not_smallest_real_part = %s"
          % (tag, e_ed.real, e_ed.imag, o_ed, e_ed.real-ev[k].real,
             e_bi.real, e_bi.imag, wrong))


for eps in (1e-4, 1e-2):  # 1e-2 is the control, outside 1e-3*g
    sc, obs = impurity(4, eps, 0.5)
    check("impurity eps=%g g=0.5" % eps, sc, obs)
for eps in (5e-5, 1e-3):  # 1e-3 is the control
    sc, obs = pt_field(6, eps)
    check("PT_XX_field eps=%g" % eps, sc, obs)
```

Observed:

```
impurity eps=0.0001 g=0.5: biorthogonal_ground_state raised ValueError('Right and left ground states are (near) orthogonal, <vl|vr> = 0j - biorthogonal normalization failed')
impurity eps=0.0001 g=0.5: true smallest-Re E = -1.50005000+0.25000000j <obs> = -0.5000
impurity eps=0.0001 g=0.5: ED gs_energy = -1.49995000-0.25000000j <obs>(ED vev) = +0.5000  Re(ED)-ReMin = +1.00e-04  biorthogonal_ground_state E = nan+0.00000000j  ED_not_smallest_real_part = True
impurity eps=0.01 g=0.5: true smallest-Re E = -1.50500000+0.25000000j <obs> = -0.5000
impurity eps=0.01 g=0.5: ED gs_energy = -1.50500000+0.25000000j <obs>(ED vev) = -0.5000  Re(ED)-ReMin = +0.00e+00  biorthogonal_ground_state E = -1.50500000+0.25000000j  ED_not_smallest_real_part = False
PT_XX_field eps=5e-05: biorthogonal_ground_state raised ValueError('Right and left ground states are (near) orthogonal, <vl|vr> = (-9.103828801926284e-15-1.5265566588595902e-15j) - biorthogonal normalization failed')
PT_XX_field eps=5e-05: true smallest-Re E = -1.39617641-0.00000000j <obs> = -1.0000
PT_XX_field eps=5e-05: ED gs_energy = -1.39612641-0.20120744j <obs>(ED vev) = -0.0000  Re(ED)-ReMin = +5.00e-05  biorthogonal_ground_state E = nan+0.00000000j  ED_not_smallest_real_part = True
PT_XX_field eps=0.001: true smallest-Re E = -1.39712641+0.00000000j <obs> = -1.0000
PT_XX_field eps=0.001: ED gs_energy = -1.39712641+0.00000000j <obs>(ED vev) = -1.0000  Re(ED)-ReMin = +0.00e+00  biorthogonal_ground_state E = -1.39712641+0.00000000j  ED_not_smallest_real_part = False
[exit 0, 3 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: `biorthogonal_ground_state` as a silently wrong consumer, since it
raised `ValueError` in both cases (that is the reviewer-found finding
below); and the DMRG excited states, read only. The condition is narrowed to
its direction: the level returned must have the smaller Im, `Re A - Re M <
1e-3 (Im M - Im A)`. LOW, since it needs a near tie of real parts.

**Suggested fix**: A lexicographic key: real part first, grouped at a tolerance set by the
eigenvalues (for dense ED `1e-10 max|E|` is enough), and Im ascending inside
a group, in `sorteigen` and so in both `lowest_states` branches and in
`excited.py`. Record it in the ledger, since ED non-Hermitian references move
in this band. Fix it with the imaginary-time tie window above, in one cluster,
and pin both against `np.linalg.eig`.

### 13. `algebra.biorthogonal_ground_state` raises at any scale when the smallest real part is shared by a conjugate pair, and the ED non-Hermitian KPM correlator raises with it

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `solvers` &middot; candidate found by the reviewer of `solvers-2`

**Where**: `src/dmrgpy/algebra/algebra.py (biorthogonal_ground_state sorts h and dagger(h) with the same key)`; `src/dmrgpy/algebra/algebra.py:381-388 (sorteigen)`; `src/dmrgpy/nonhermitian/kpm.py:111 (the ED NH-KPM route)`

`biorthogonal_ground_state` sorts h and `dagger(h)` with the same key,
`Re + 1e-3 Im`. On `dagger(h)` that key acts on `conj(E)`, so when the
smallest real part is shared by levels of different imaginary part, the
right state belongs to one member of the pair and the left state to the
other. The two are biorthogonal, and the call raises `ValueError`. On the
6-site PT-symmetric XX chain at g = 0.3, four levels tie at Re = -1.396126,
among them the pair `-1.396126 +- 0.201207i`; on the 4-site Heisenberg chain
with `0.3 Sz_tot + 1j Sz_0` the pair is `-1.242697 +- 0.291419i`. Both raise
at unit scale, and `get_dynamical_correlator(mode="ED")` with the default
non-Hermitian KPM submode raises with them, while scipy's left and right
eigenvectors of the same level have `|<vl|vr>|` = 0.55 and 0.62. This
sharpens a lead of `docs/audit_2026_09_24c_hole_hunt.md` ("New leads, not
reviewed"), which saw the raise at s = 1e-3 and reported s = 1 as passing;
its s = 1 calls had stopped earlier, at the missing `E_max`.

**Expected**: the biorthogonal pair of one level, as `scipy.linalg.eig(h, left=True)` gives it.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g10/g10_new_biorth_conj_tie.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis (sharpens the unreviewed 2026-09-24c lead at
# docs/audit_2026_09_24c_hole_hunt.md:10102): algebra.biorthogonal_ground_state
# sorts the eigenvalues of h and of dagger(h) with the same tilted key
# Re + 1e-3*Im (sorteigen). The eigenvalues of dagger(h) are conj(E), so the
# tilt has the opposite sign for the left state: on any tie of the smallest
# real part between levels of different Im (a PT conjugate pair), the right
# state is the member with Im<0 and the left state belongs to the member with
# Im>0, the two are biorthogonal and the call raises ValueError, at unit
# scale. The ED NH-KPM correlator (get_dynamical_correlator(mode="ED") on a
# non-Hermitian chain, nonhermitian/kpm.py:111) raises with it. Deterministic
# (dense ED, no random start).
from scipy.sparse import csr_matrix
import scipy.linalg as sla
from dmrgpy import spinchain
from dmrgpy.algebra import algebra


def pt_chain(n, g=0.3):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(n-1): h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n): h = h + 1j*g*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def heis_pair(L=4, gamma=1.0):
    """The 2026-09-24c lead's model at unit scale."""
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L): h = h + 0.3*sc.Sz[i]
    h = h + 1j*gamma*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc


def check(tag, sc):
    hm = sc.get_ED_obj().get_hamiltonian()
    hm = hm.toarray() if hasattr(hm, "toarray") else np.array(hm)
    ev = np.linalg.eigvals(hm)
    ev = ev[np.argsort(ev.real, kind="stable")]
    ntie = int(np.sum(np.abs(ev.real - ev[0].real) < 1e-10*(1+abs(ev[0]))))
    print("%s: lowest levels %s  tied at the smallest Re: %d"
          % (tag, np.round(ev[:4], 6), ntie))
    er, _ = algebra.lowest_states(csr_matrix(hm), n=1)
    el, _ = algebra.lowest_states(algebra.dagger(csr_matrix(hm)), n=1)
    print("%s: right state picked at E = %s, left state picked for E = conj(%s) = %s"
          % (tag, np.round(complex(er[0]), 6), np.round(complex(el[0]), 6),
             np.round(np.conj(complex(el[0])), 6)))
    print("%s: left_belongs_to_another_level = %s"
          % (tag, abs(np.conj(complex(el[0])) - complex(er[0])) > 1e-8))
    try:
        e0 = algebra.biorthogonal_ground_state(csr_matrix(hm))[0]
        print("%s: biorthogonal_ground_state OK E = %s" % (tag, np.round(complex(e0), 6)))
        braised = False
    except Exception as err:
        print("%s: biorthogonal_ground_state raised %s: %s"
              % (tag, type(err).__name__, str(err)[:60]))
        braised = isinstance(err, ValueError)
    try:
        x, y = sc.get_dynamical_correlator(mode="ED", name=(sc.Sz[1], sc.Sz[1]),
                                           delta=0.1, es=np.linspace(-0.5, 3.0, 20),
                                           E_max=20.0, n=50)
        print("%s: ED NH-KPM correlator OK max|C| = %.4f" % (tag, np.max(np.abs(y))))
        kraised = False
    except Exception as err:
        print("%s: ED NH-KPM correlator raised %s: %s"
              % (tag, type(err).__name__, str(err)[:60]))
        kraised = "biorthogonal" in str(err) or "orthogonal" in str(err)
    # anchor: the pair exists. scipy's left eigenvectors of the same index
    w, vl, vr = sla.eig(hm, left=True, right=True)
    k = int(np.lexsort((w.imag, np.round(w.real, 10)))[0])
    ov = abs(np.conj(vl[:, k]) @ vr[:, k])
    print("%s: anchored pair at E = %s has |<vl|vr>| = %.3e (exists = %s)"
          % (tag, np.round(w[k], 6), ov, ov > 1e-8))
    print("SUMMARY %s: tie=%d biorthogonal_raised=%s nhkpm_ed_raised=%s anchored_pair_exists=%s"
          % (tag, ntie, braised, kraised, ov > 1e-8))


check("PT_XX_6sites_eps0", pt_chain(6))
check("heis4_gamma1_s1", heis_pair())
```

Observed:

```
PT_XX_6sites_eps0: lowest levels [-1.396126+0.j       -1.396126-0.201207j -1.396126+0.j
 -1.396126+0.201207j]  tied at the smallest Re: 4
PT_XX_6sites_eps0: right state picked at E = (-1.396126-0.201207j), left state picked for E = conj((-1.396126-0.201207j)) = (-1.396126+0.201207j)
PT_XX_6sites_eps0: left_belongs_to_another_level = True
PT_XX_6sites_eps0: biorthogonal_ground_state raised ValueError: Right and left ground states are (near) orthogonal, <vl|vr> 
PT_XX_6sites_eps0: ED NH-KPM correlator raised ValueError: Right and left ground states are (near) orthogonal, <vl|vr> 
PT_XX_6sites_eps0: anchored pair at E = (-1.396126-0.201207j) has |<vl|vr>| = 5.544e-01 (exists = True)
SUMMARY PT_XX_6sites_eps0: tie=4 biorthogonal_raised=True nhkpm_ed_raised=True anchored_pair_exists=True
heis4_gamma1_s1: lowest levels [-1.242697+0.291419j -1.242697-0.291419j -1.219455-0.456421j
 -0.619455+0.456421j]  tied at the smallest Re: 2
heis4_gamma1_s1: right state picked at E = (-1.242697-0.291419j), left state picked for E = conj((-1.242697-0.291419j)) = (-1.242697+0.291419j)
heis4_gamma1_s1: left_belongs_to_another_level = True
heis4_gamma1_s1: biorthogonal_ground_state raised ValueError: Right and left ground states are (near) orthogonal, <vl|vr> 
heis4_gamma1_s1: ED NH-KPM correlator raised ValueError: Right and left ground states are (near) orthogonal, <vl|vr> 
heis4_gamma1_s1: anchored pair at E = (-1.242697-0.291419j) has |<vl|vr>| = 6.234e-01 (exists = True)
SUMMARY heis4_gamma1_s1: tie=2 biorthogonal_raised=True nhkpm_ed_raised=True anchored_pair_exists=True
[exit 0, 3 s]
```

**Reviewer (CONFIRMED)**: Confirmed as stated. The lexicographic sort suggested for the ED sort
finding above has the same sign problem, so it does not fix this one. LOW,
since it raises.

**Suggested fix**: Pick the left state by value, not by its own sort: sort h only and take the
eigenvector of `dagger(h)` whose eigenvalue is nearest `conj(e0)`, or take
both from one `scipy.linalg.eig(h, left=True, right=True)` at the same
index. Inside a degenerate level choose `vl` in the left eigenspace so that
`<vl|vr>` is nonzero. Fix it in the cluster of the two findings above.

### 14. The opt-in `tol=` stop of the ED Krylov Kondo route counts an empty window as converged, and drops a whole inelastic step: 14.3% of the peak on a 12-spin chain, 64% on the 8-spin test chain at one bias

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kondo-routes` &middot; candidate `kondo-routes-1`

**Where**: `src/dmrgpy/kondospectrumtk/krylovkondo.py:254-264 (the check: an empty window counts as converged)`; `src/dmrgpy/kondospectrumtk/krylovkondo.py:476-501 (window() and state_spectrum)`; `src/dmrgpy/kondospectrumtk/krylovkondo.py:525-549 (krylov_kondo_dIdV passes emax = max|eV|)`; `src/dmrgpy/kondospectrumtk/kondomethod.py:377-387 (the public tol= keyword)`; `docs/guide/kondo.md:291-295 (the promise that every possible final state is converged)`

With `tol=` set, the residual stop of `krylovkondo.block_krylov` tests only
the Ritz states already inside the window of open final states, and an
empty window counts as converged (`a9d8d06`). Ritz values lie above the
eigenvalues they approach, so at T=0 the stop can fire at its first check,
10 blocks in, while a final state below `max|eV|` has no Ritz value in the
window yet. Its whole step is then missing from dI/dV, and nothing warns. On
the uniform 12-spin ferromagnetic chain (J = -1 meV, B = 0.1 meV, tip on site
0, T = 0, `nkry=100`) with `max|eV|` = 0.1056 meV, the only inelastic step
of the scan is missing: the k=0 magnon at exactly B, with `M2 = 1/24`. The
error is 14.3% of the peak at order 2 and 14.4% at order 3, and the absolute
value at order 2 is 0.262 = 2pi/24. The counter-repro shows this holds for
every `max|eV|` from B up to the first-check Ritz value at B + 0.0112 meV,
whatever B is (measured at B = 0.01 meV for `max|eV|/B` = 1.2, 1.5 and 2.0),
at 2 mK as at T=0, and not at 10 mK. On the 8-spin AF chain of the tests the
stop at the first check also drops weight, 64% of the peak at order 2, but
only in a band of width 2e-6 of `emax` at 0.31282 meV.

**Expected**: the stop firing only once every final state below `max|eV|` is converged, as the guide promises.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/kondo-routes/kr_1_tol_stop_empty_window.py
```
```python
"""kondo-routes 1: the opt-in residual stop of submode="KRYLOV" (tol=, a9d8d06)
tests that the Ritz states it HAS in the window are converged, and counts a
window with no Ritz state at all as converged (krylovkondo.py:258). Ritz
values are upper bounds, so at the first check (10 blocks) a true final
state inside the window can still sit in a Ritz value above it: the stop
fires and the step of that state is missing from dI/dV.

Ferromagnetic S=1/2 chains in a field along z: the ground state is the
polarized one, S^-_0|GS> is one column, so each block is one state and the
first check is at 10 states of the 12 the one-magnon space holds. The
window top emax is put between a true one-magnon energy and the Ritz values
of the 10-state space (computed here with an independent Lanczos), and the
public call with tol=1e-8 is compared with the default (no stop, exhausted
space, exact).
"""
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
import warnings
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import krylovkondo, conductance

np.random.seed(1) # eigsh draws its own start vector from fresh entropy; the
# ground state here is unique (the polarized one), so that does not matter
meV = 1e-3
N = 12


def chain(Jd, Jp, Jc, B0, Bc):
    """Sites 0-1 coupled by Jd, 1-2 by Jp, 2..N-1 by Jc (meV, negative is
    ferromagnetic), field -B Sz with B0 on sites 0,1 and Bc on the rest."""
    sc = spinchain.Spin_Chain(["S=1/2"]*N)

    def dot(i, j):
        return sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
    h = Jd*meV*dot(0, 1) + Jp*meV*dot(1, 2)
    for i in range(2, N-1):
        h = h + Jc*meV*dot(i, i+1)
    for i in range(N):
        h = h - (B0 if i < 2 else Bc)*meV*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def lanczos(H, v0, psi, m):
    """m Lanczos vectors from v0, orthogonal to psi, full reorthogonalization:
    Ritz values, vectors in Lanczos coordinates, residuals beta*|y_last|."""
    q = v0 - psi*(psi @ v0)
    q = q/np.linalg.norm(q)
    Q = [q]
    for _ in range(m - 1):
        w = H @ Q[-1]
        for _ in range(2):
            w = w - psi*(psi @ w)
            for qq in Q:
                w = w - qq*(qq @ w)
        Q.append(w/np.linalg.norm(w))
    Qm = np.stack(Q, axis=1)
    T = Qm.T @ (H @ Qm)
    theta, y = np.linalg.eigh((T + T.T)/2.)
    w = H @ Qm[:, -1]
    for _ in range(2):
        w = w - psi*(psi @ w)
        w = w - Qm @ (Qm.T @ w)
    return theta, y, np.linalg.norm(w)*np.abs(y[-1, :])


cases = [("uniform FM chain, B=0.1 meV", (-1.0, -1.0, -1.0, 0.1, 0.1)),
         ("dimer in 0.5 meV, chain in 0.1 meV", (-1.0, -0.3, -1.0, 0.5, 0.1)),
         ("dimer mode in window (non-empty)", (-2.0, -0.2, -0.1, 0.1, 1.0))]

for label, pars in cases:
    print("==== case:", label, "pars(Jd,Jp,Jc,B0,Bc) meV =", pars)
    sc = chain(*pars)
    ref = krylovkondo.KrylovKondoSpectrum(sc, 0, 0., nkry=100)
    print("su2 taken =", ref.su2, " real =", ref.real, " n_initial =", len(ref.weights))
    psi = ref.V[:, 0]/np.linalg.norm(ref.V[:, 0])
    Ei = float(psi @ (ref.H @ psi))
    Sp, Sm, Sz = ref._ladder
    v0 = Sm @ psi
    print("|S+ psi| =", np.linalg.norm(Sp @ psi), " |S- psi| =", np.linalg.norm(v0),
          " |Sz psi - <Sz>psi| =", np.linalg.norm(Sz @ psi - (psi @ (Sz @ psi))*psi))
    ss = ref.state_spectrum(0) # exhausted: exact
    e_ex = ss.e[1:]
    M2 = conductance._spin_matrix_elements_squared(ss)[1:, 0]
    o = np.argsort(e_ex)
    print("reference dim =", ss.dim, " (13 = initial + 12 one-magnon states)")
    print("exact e (meV) =", np.array2string(e_ex[o]/meV, precision=5, max_line_width=200))
    print("exact M2      =", np.array2string(M2[o], precision=3, max_line_width=200))
    theta, y, res = lanczos(ref.H, v0, psi, 10)
    r10 = theta - Ei
    w10 = 0.5*(np.linalg.norm(v0)**2)*y[0, :]**2
    tolabs = 1e-8*ref.scale
    print("ritz10 e (meV)=", np.array2string(r10/meV, precision=5, max_line_width=200))
    print("ritz10 resid/tolabs =", np.array2string(res/tolabs, precision=2, max_line_width=200))
    # a window top with more exact states than converged Ritz states below
    # it, half way between two neighbouring energies (exact or Ritz), so
    # that it sits as far as it can from both
    pts = np.unique(np.concatenate([e_ex, r10]))
    valid = []
    for lo, hi in zip(pts[:-1], pts[1:]):
        c = 0.5*(lo + hi)
        if c <= 0.:
            continue
        nin_ex, inr = int(np.sum(e_ex <= c)), r10 <= c
        if nin_ex > int(np.sum(inr)) and np.all(res[inr] < tolabs):
            missed = float(np.sum(M2[e_ex <= c]) - np.sum(w10[inr]))
            valid.append((c, missed, nin_ex, int(np.sum(inr)), 0.5*(hi - lo)))
    if not valid:
        print("design_ok = False (no window top with a missing state)")
        continue
    top_missed = max(v[1] for v in valid)
    # among the windows that miss at least half the largest weight, the one
    # farthest from every exact and Ritz energy
    emax, missed, nin_ex, nin_r, margin = max(
        (v for v in valid if v[1] >= 0.5*top_missed), key=lambda v: v[4])
    above = r10[r10 > emax]
    below = e_ex[e_ex <= emax]
    print("design_ok = True  emax (meV) =", emax/meV, " exact states in window =", nin_ex,
          " converged Ritz states in window =", nin_r, " missed M2 weight =", missed)
    print("margin: lowest Ritz above window - emax = %.4g meV, emax - highest exact in "
          "window = %.4g meV, margin/emax = %.3g" % (((np.min(above) if above.size
                                                       else np.inf) - emax)/meV,
                                                     (emax - np.max(below))/meV,
                                                     margin/emax))
    stop = krylovkondo.KrylovKondoSpectrum(sc, 0, 0., nkry=100, tol=1e-8)
    a = stop.state_spectrum(0, emax)
    print("stopped dim =", a.dim, " reference dim =", ss.dim,
          " stopped e (meV) =", np.array2string(np.sort(a.e[1:])/meV, precision=5,
                                                max_line_width=200))
    inwin = e_ex[e_ex <= emax]
    mids = 0.5*(inwin + emax)
    eV = np.unique(np.concatenate([np.linspace(-emax, emax, 201), mids, -mids]))
    for order in (2, 3):
        kw = dict(site=0, Jrho_s=-0.05, T=0., order=order, mode="ED", submode="KRYLOV",
                  nkry=100)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, d_stop = sc.get_kondo_spectrum(eV, tol=1e-8, **kw)
            _, d_ref = sc.get_kondo_spectrum(eV, **kw)
        diff = np.max(np.abs(d_stop - d_ref))
        print("order=%d  max|dIdV(tol=1e-8) - dIdV(default)| / peak = %.3e" %
              (order, diff/np.max(np.abs(d_ref))),
              "  abs = %.4e" % diff, "  2*pi*missed = %.4e" % (2*np.pi*missed))
```

Observed:

```
==== case: uniform FM chain, B=0.1 meV pars(Jd,Jp,Jc,B0,Bc) meV = (-1.0, -1.0, -1.0, 0.1, 0.1)
su2 taken = False  real = True  n_initial = 1
|S+ psi| = 8.393410165696443e-15  |S- psi| = 1.0  |Sz psi - <Sz>psi| = 8.393410165696443e-15
reference dim = 301  (13 = initial + 12 one-magnon states)
exact e (meV) = [0.1     0.13407 0.2     0.23397 0.23407 0.27079 0.3     0.33397 0.33407 0.36571 0.37079 0.39289 0.4     0.4107  0.43397 0.43407 0.46571 0.46985 0.47079 0.49289 0.49996 0.5     0.5107  0.52852
 0.53399 0.5341  0.55456 0.56575 0.56989 0.58015 0.5929  0.6     0.60633 0.62661 0.63161 0.63702 0.65605 0.66675 0.68222 0.69574 0.70002 0.72738 0.73735 0.75235 0.7677  0.78465 0.80001 0.81105
 0.84118 0.84756 0.86392 0.88135 0.9003  0.92446 0.94118 0.96449 0.98053 0.99914 1.02433 1.04117 1.06676 1.0815  1.1     1.125   1.14138 1.17719 1.2     1.22696 1.24089 1.27608 1.29986 1.31272
 1.34204 1.35882 1.39206 1.4083  1.44112 1.45883 1.49976 1.50637 1.55138 1.5591  1.6     1.63053 1.65741 1.69732 1.70007 1.74433 1.76799 1.80352 1.80711 1.84997 1.88306 1.9071  1.95042 1.96603
 2.00708 2.0364  2.06593 2.09975 2.11561 2.16584 2.18484 2.21773 2.25681 2.2713  2.31397 2.35147 2.37841 2.40959 2.44638 2.47904 2.50936 2.54916 2.57828 2.61108 2.64942 2.68231 2.71414 2.7506
 2.78551 2.81521 2.85444 2.88731 2.92409 2.96271 2.99487 3.02752 3.06587 3.09872 3.13778 3.17462 3.20643 3.24509 3.27988 3.31553 3.35151 3.39417 3.42319 3.45791 3.49655 3.53236 3.57643 3.60448
 3.64345 3.68235 3.71692 3.75449 3.78924 3.82655 3.86111 3.8961  3.93965 3.97117 4.00717 4.04421 4.08036 4.11544 4.15548 4.19226 4.22501 4.26105 4.30072 4.32862 4.37344 4.41452 4.44356 4.47803
 4.51387 4.54597 4.58546 4.62465 4.65319 4.69664 4.72805 4.76548 4.79907 4.83599 4.87215 4.90363 4.93864 4.96908 5.00962 5.03746 5.07524 5.11169 5.1362  5.17843 5.20531 5.23781 5.27831 5.30505
 5.34045 5.37161 5.40298 5.4317  5.47044 5.50414 5.53097 5.55793 5.59559 5.63018 5.65563 5.68693 5.7184  5.74233 5.77283 5.80121 5.83194 5.8601  5.89286 5.92182 5.95048 5.9783  6.00579 6.03205
 6.0575  6.08006 6.11233 6.14453 6.15964 6.18658 6.21982 6.25063 6.26589 6.29828 6.32108 6.35723 6.36763 6.39343 6.41997 6.44011 6.46149 6.47274 6.52016 6.52596 6.55202 6.57808 6.59513 6.62365
 6.62641 6.65408 6.67143 6.69463 6.7246  6.7298  6.75106 6.76012 6.76414 6.79051 6.79535 6.83075 6.85012 6.85392 6.86074 6.8896  6.8935  6.92741 6.92926 6.93079 6.95381 6.95555 6.96064 6.9896
 6.9935  7.03079 7.05134 7.05381 7.05555 7.06064 7.12312 7.14175 7.15381 7.15991 7.16064 7.16842 7.22312 7.25381 7.25991 7.26842 7.32312 7.35991 7.36842 7.42406 7.43826 7.44163 7.45991 7.53826
 7.54163 7.55991 7.63826 7.64163 7.75783 7.76329 7.86329 7.96329 8.11115 8.21115 8.31115 8.49209]
exact M2      = [4.167e-02 8.191e-02 1.188e-29 7.775e-02 2.469e-30 1.268e-30 6.151e-31 1.522e-28 2.407e-30 1.339e-28 7.351e-29 7.113e-02 3.729e-30 2.723e-30 2.474e-29 4.017e-32 2.001e-31 1.124e-31 2.217e-32
 6.406e-31 3.124e-32 3.672e-32 2.418e-31 1.965e-29 1.232e-30 4.284e-31 1.009e-31 1.618e-30 1.824e-29 1.473e-32 7.415e-29 6.250e-02 5.553e-27 1.467e-30 4.723e-30 1.517e-30 3.307e-31 1.681e-30
 4.244e-32 9.583e-31 1.663e-28 4.371e-33 4.030e-33 9.692e-33 3.218e-33 2.691e-33 9.539e-31 1.101e-32 5.245e-02 3.808e-32 2.076e-32 2.619e-32 5.460e-32 4.855e-33 5.690e-29 3.190e-32 6.316e-31
 1.787e-32 4.823e-33 9.320e-31 2.074e-33 4.071e-33 4.167e-02 2.182e-32 4.330e-32 6.765e-33 1.094e-29 1.036e-30 2.584e-29 1.672e-30 4.975e-30 2.201e-30 2.238e-28 3.088e-02 7.444e-30 1.671e-30
 7.601e-33 8.095e-29 3.215e-31 3.219e-31 1.357e-30 2.482e-29 2.083e-02 2.559e-31 1.019e-30 4.782e-31 1.189e-29 1.739e-33 5.454e-33 1.407e-32 1.220e-02 1.496e-33 2.504e-32 7.914e-30 1.920e-33
 5.582e-03 1.010e-31 5.120e-33 1.420e-03 2.762e-33 1.649e-33 4.233e-30 5.402e-32 4.081e-32 7.826e-32 1.436e-31 6.053e-32 2.030e-32 8.343e-33 1.035e-32 2.617e-32 1.252e-32 1.478e-33 1.093e-33
 4.488e-33 2.645e-33 7.641e-35 6.738e-34 3.671e-34 3.034e-33 7.384e-34 3.623e-35 8.297e-34 3.943e-33 1.479e-33 1.769e-34 1.123e-33 3.746e-33 5.372e-33 1.875e-33 1.655e-33 5.723e-33 6.545e-33
 2.494e-33 2.374e-34 4.436e-34 1.348e-33 7.677e-34 8.031e-35 8.881e-35 1.201e-33 8.087e-34 2.967e-36 8.196e-35 1.046e-34 4.669e-34 9.222e-35 9.430e-35 9.645e-34 2.423e-33 2.157e-33 5.488e-34
 7.465e-34 1.691e-33 1.311e-33 3.441e-34 1.422e-34 1.012e-33 1.394e-33 9.593e-34 1.141e-34 1.066e-34 4.376e-34 7.153e-35 7.120e-35 3.033e-35 4.039e-34 1.005e-33 5.410e-34 2.951e-34 5.877e-34
 1.078e-33 4.433e-34 1.117e-34 4.980e-34 1.242e-33 1.355e-33 4.627e-34 1.076e-33 8.266e-34 1.373e-33 1.273e-34 2.273e-35 1.064e-34 6.047e-34 6.642e-34 3.429e-34 2.164e-34 8.861e-34 8.733e-34
 2.642e-34 4.061e-34 9.594e-34 9.391e-34 1.221e-33 2.685e-34 1.446e-33 1.003e-33 2.194e-33 2.661e-34 3.699e-35 4.369e-34 6.396e-34 2.165e-34 1.794e-34 4.887e-34 7.176e-34 6.189e-34 1.617e-34
 1.037e-34 3.489e-34 3.922e-34 1.994e-34 5.257e-35 4.660e-34 3.776e-34 2.338e-33 4.756e-34 2.282e-34 5.238e-34 5.670e-34 7.802e-34 9.178e-34 1.776e-35 2.764e-34 5.536e-34 2.825e-34 9.856e-35
 5.673e-34 7.060e-35 4.724e-34 3.666e-34 1.823e-34 4.095e-34 4.605e-34 2.453e-34 1.777e-34 6.514e-34 3.052e-34 6.741e-34 1.704e-33 3.490e-34 1.278e-33 1.600e-34 3.144e-34 5.988e-34 5.046e-34
 6.889e-34 1.069e-33 3.845e-35 2.312e-33 5.552e-34 3.065e-34 1.441e-34 4.638e-34 7.105e-34 4.111e-34 5.625e-35 1.890e-34 4.695e-34 5.650e-34 8.658e-34 1.644e-34 1.117e-34 2.033e-33 3.048e-34
 1.494e-33 4.756e-34 6.098e-34 4.122e-34 4.606e-34 9.606e-34 2.164e-33 2.205e-35 7.292e-34 6.250e-35 2.093e-33 3.885e-33 2.732e-34 4.246e-34 8.944e-34 2.337e-34 6.537e-34 8.022e-34 1.354e-34
 3.150e-34 7.479e-34 1.778e-33 3.523e-34 9.748e-34 7.163e-33 1.625e-33 3.787e-33 3.313e-34 4.927e-34 6.801e-33 7.289e-34 8.102e-34 5.294e-33 1.408e-33]
ritz10 e (meV)= [0.11117 0.19903 0.36695 0.6     0.87748 1.17473 1.46534 1.72349 1.92624 2.05557]
ritz10 resid/tolabs = [ 394226.53 1147650.84 1799101.33 2290693.78 2578748.03 2637669.2  2462221.89 2067995.35 1490018.36  779646.69]
design_ok = True  emax (meV) = 0.10558458688747496  exact states in window = 1  converged Ritz states in window = 0  missed M2 weight = 0.041666666666665624
margin: lowest Ritz above window - emax = 0.005585 meV, emax - highest exact in window = 0.005585 meV, margin/emax = 0.0529
stopped dim = 11  reference dim = 301  stopped e (meV) = [0.11117 0.19903 0.36695 0.6     0.87748 1.17473 1.46534 1.72349 1.92624 2.05557]
order=2  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 1.429e-01   abs = 2.6180e-01   2*pi*missed = 2.6180e-01
order=3  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 1.440e-01   abs = 3.7275e-01   2*pi*missed = 2.6180e-01
==== case: dimer in 0.5 meV, chain in 0.1 meV pars(Jd,Jp,Jc,B0,Bc) meV = (-1.0, -0.3, -1.0, 0.5, 0.1)
su2 taken = False  real = True  n_initial = 1
|S+ psi| = 1.0028757704647817e-15  |S- psi| = 1.0  |Sz psi - <Sz>psi| = 1.0028757704647817e-15
reference dim = 301  (13 = initial + 12 one-magnon states)
exact e (meV) = [0.10681 0.16952 0.21414 0.27579 0.31037 0.32208 0.34286 0.38261 0.41688 0.43079 0.44774 0.47042 0.49012 0.50899 0.51956 0.52397 0.54045 0.55315 0.57721 0.59725 0.59853 0.60997 0.61592 0.62087
 0.62982 0.63176 0.65137 0.65929 0.68189 0.68971 0.69565 0.70065 0.7176  0.72339 0.73672 0.73813 0.75702 0.77879 0.78343 0.80462 0.81684 0.83189 0.84665 0.85854 0.86076 0.88384 0.91187 0.92458
 0.93639 0.96154 0.97821 0.9955  1.00648 1.03266 1.04826 1.05959 1.07878 1.11585 1.12657 1.15018 1.16592 1.18561 1.21574 1.22484 1.2482  1.27275 1.29908 1.31256 1.33305 1.36571 1.38905 1.41523
 1.45434 1.46841 1.48879 1.5236  1.55234 1.59137 1.60194 1.63249 1.65551 1.69726 1.70206 1.74557 1.75723 1.79706 1.80662 1.83825 1.8597  1.91356 1.91358 1.95432 1.98145 2.01664 2.05214 2.06311
 2.09904 2.12667 2.16246 2.19072 2.22544 2.25121 2.28165 2.31834 2.34794 2.38091 2.4167  2.45632 2.48352 2.52391 2.55571 2.59084 2.61774 2.65252 2.68813 2.71593 2.75574 2.78844 2.81662 2.85243
 2.88696 2.91947 2.95681 2.9881  3.02667 3.05852 3.09167 3.13    3.16176 3.20365 3.23327 3.27194 3.30626 3.34027 3.37696 3.41599 3.44959 3.4902  3.52182 3.56058 3.59987 3.63722 3.6722  3.70751
 3.74479 3.78142 3.8224  3.85913 3.895   3.92554 3.96935 4.00685 4.04397 4.08415 4.11695 4.15637 4.19121 4.22969 4.26808 4.3021  4.33194 4.36983 4.40513 4.43944 4.4773  4.5097  4.54675 4.57863
 4.61887 4.65173 4.69029 4.72837 4.76123 4.794   4.83481 4.87198 4.90654 4.94389 4.97595 5.0104  5.0514  5.08368 5.12577 5.1527  5.19319 5.2258  5.26004 5.28689 5.3238  5.36415 5.38903 5.42431
 5.45981 5.49436 5.52405 5.55842 5.58889 5.62179 5.65331 5.68042 5.72361 5.74899 5.7774  5.81272 5.84892 5.87628 5.90363 5.937   5.96833 5.9951  6.02136 6.04859 6.08731 6.11206 6.13922 6.16698
 6.19488 6.2322  6.25456 6.27566 6.31543 6.33657 6.36484 6.38783 6.41594 6.44068 6.46197 6.48242 6.5187  6.54156 6.55424 6.59052 6.61117 6.63823 6.65647 6.67941 6.70063 6.71391 6.72974 6.75873
 6.7733  6.79619 6.82516 6.83966 6.85164 6.8722  6.89133 6.8985  6.91651 6.93226 6.93686 6.95104 6.97168 6.99534 7.00529 7.02687 7.04492 7.05764 7.06579 7.0697  7.07468 7.09209 7.12439 7.14467
 7.14566 7.17376 7.17475 7.21    7.21538 7.22058 7.24005 7.24239 7.27335 7.31972 7.32058 7.3408  7.35277 7.3687  7.39549 7.40434 7.41119 7.43113 7.48441 7.49359 7.51615 7.59264 7.59756 7.60151
 7.69012 7.76598 7.76612 7.77171 7.82787 7.87578 7.94937 8.11066 8.16672 8.26993 8.37529 8.59381]
exact M2      = [1.830e-04 1.291e-03 5.239e-32 1.881e-31 4.814e-03 3.022e-32 1.503e-31 4.578e-34 2.153e-30 1.962e-33 7.832e-34 6.631e-31 7.439e-33 7.324e-02 1.912e-32 9.331e-30 6.084e-33 1.124e-30 7.328e-28
 1.959e-01 3.373e-33 5.381e-30 7.399e-28 1.446e-32 1.833e-32 2.825e-33 3.843e-32 2.770e-32 2.345e-31 7.071e-34 2.084e-32 9.383e-29 6.629e-31 4.466e-30 2.177e-30 2.171e-33 1.268e-33 1.689e-29
 4.945e-29 1.809e-29 1.293e-02 1.209e-30 4.175e-33 7.654e-32 8.170e-31 1.395e-31 3.140e-33 6.706e-31 1.593e-31 7.000e-32 9.293e-32 1.441e-31 1.856e-31 2.245e-30 5.900e-31 6.471e-32 7.263e-33
 5.503e-03 7.059e-33 8.872e-32 9.067e-32 7.278e-32 1.454e-31 1.858e-30 6.095e-32 2.952e-32 3.677e-33 3.083e-33 2.955e-32 1.724e-33 3.383e-34 1.115e-02 1.524e-32 4.765e-33 3.268e-33 2.724e-32
 5.038e-33 1.848e-01 3.663e-32 1.734e-30 3.820e-30 3.789e-30 9.770e-03 1.060e-29 2.219e-30 6.051e-30 1.433e-29 2.469e-31 2.732e-31 2.446e-30 3.957e-04 4.919e-31 3.868e-31 1.550e-31 4.080e-05
 5.056e-33 1.086e-33 5.001e-32 8.598e-32 9.580e-32 8.555e-32 5.262e-32 5.259e-32 3.021e-32 7.651e-33 2.600e-33 1.111e-33 5.281e-33 7.071e-33 1.391e-32 3.241e-32 3.785e-32 3.771e-32 3.288e-32
 3.187e-32 1.771e-32 1.492e-32 2.470e-32 2.530e-32 3.945e-32 2.249e-32 2.056e-32 8.645e-33 9.469e-33 4.182e-33 7.270e-33 1.132e-32 9.758e-33 1.476e-32 1.390e-32 6.379e-33 1.813e-33 9.131e-34
 7.444e-34 4.467e-34 7.824e-34 2.503e-33 8.176e-33 2.616e-33 1.456e-33 1.721e-34 8.348e-35 3.337e-35 4.494e-35 9.302e-35 9.286e-34 2.191e-33 2.377e-33 1.445e-33 1.488e-33 2.595e-34 1.840e-35
 2.206e-34 5.884e-35 4.547e-34 1.797e-33 1.980e-33 2.027e-33 2.844e-33 5.542e-34 9.050e-35 3.774e-35 3.891e-35 6.965e-34 1.683e-33 1.559e-33 2.088e-33 1.108e-33 5.961e-35 7.423e-35 2.629e-34
 2.349e-36 3.565e-34 7.103e-34 6.128e-34 1.392e-33 8.461e-34 1.271e-34 1.120e-34 2.856e-34 2.772e-34 8.776e-35 7.520e-35 1.608e-34 5.101e-34 1.682e-33 5.052e-34 2.129e-34 2.669e-35 1.399e-34
 1.515e-34 1.218e-34 3.629e-35 6.778e-34 1.072e-33 1.143e-33 2.524e-34 2.818e-35 3.011e-35 2.921e-34 2.555e-34 1.157e-34 8.295e-35 3.378e-34 9.309e-34 7.663e-34 3.046e-34 6.370e-35 5.697e-35
 3.313e-34 3.972e-34 1.106e-34 9.031e-35 1.755e-34 3.251e-34 4.106e-34 7.177e-34 6.592e-34 2.097e-34 2.988e-35 5.662e-34 3.931e-34 2.556e-34 1.912e-34 2.531e-35 4.137e-34 9.015e-34 4.264e-34
 7.268e-35 4.002e-35 7.162e-35 1.968e-34 1.083e-34 4.904e-34 3.178e-34 9.019e-35 4.106e-34 6.233e-34 2.122e-34 1.580e-34 2.961e-33 1.593e-34 1.994e-35 1.116e-34 8.896e-35 1.744e-34 4.293e-34
 1.817e-34 2.751e-34 1.310e-33 1.656e-33 2.754e-34 1.390e-33 1.560e-33 6.960e-34 1.068e-33 8.945e-36 1.199e-34 5.397e-35 2.208e-34 2.284e-33 1.778e-34 6.634e-34 5.280e-34 9.830e-35 3.654e-33
 1.834e-33 1.528e-33 1.717e-33 1.010e-35 7.571e-34 6.355e-33 3.638e-34 1.881e-34 1.538e-34 4.744e-34 8.829e-34 2.836e-33 1.671e-33 3.537e-33 1.004e-33 1.020e-33 9.427e-34 1.220e-33 2.677e-33
 8.149e-34 3.000e-33 6.578e-34 1.302e-33 2.895e-34 7.705e-34 2.014e-32 3.828e-33 1.113e-33 3.898e-33 9.290e-33 2.710e-33 8.098e-33 7.201e-34 5.638e-33]
ritz10 e (meV)= [0.13773 0.2737  0.49807 0.5971  0.85473 1.20792 1.53761 1.61213 1.84677 2.03424]
ritz10 resid/tolabs = [ 771881.93 1687338.35 2149284.76 1485969.02 2960252.15 3073235.28 2435111.27 1414209.61 2061376.68 1111671.66]
design_ok = True  emax (meV) = 0.1222698309458663  exact states in window = 1  converged Ritz states in window = 0  missed M2 weight = 0.00018295237655642258
margin: lowest Ritz above window - emax = 0.01546 meV, emax - highest exact in window = 0.01546 meV, margin/emax = 0.126
stopped dim = 11  reference dim = 301  stopped e (meV) = [0.13773 0.2737  0.49807 0.5971  0.85473 1.20792 1.53761 1.61213 1.84677 2.03424]
order=2  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 7.313e-04   abs = 1.1495e-03   2*pi*missed = 1.1495e-03
order=3  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 7.409e-04   abs = 1.5364e-03   2*pi*missed = 1.1495e-03
==== case: dimer mode in window (non-empty) pars(Jd,Jp,Jc,B0,Bc) meV = (-2.0, -0.2, -0.1, 0.1, 1.0)
su2 taken = False  real = True  n_initial = 1
|S+ psi| = 1.273391489193693e-14  |S- psi| = 1.0  |Sz psi - <Sz>psi| = 1.273391489193693e-14
reference dim = 301  (13 = initial + 12 one-magnon states)
exact e (meV) = [ 0.144    0.29293  1.00123  1.01091  1.02931  1.05463  1.08437  1.11563  1.14462  1.14533  1.15072  1.16503  1.17059  1.18696  1.18893  1.19858  1.2145   1.23254  1.24499  1.27549  1.29446  1.30307
  1.30623  1.32501  1.32776  1.33913  1.35613  1.38811  1.42034  1.44955  1.4728   1.48777  2.0026   2.01207  2.02244  2.03058  2.03837  2.05506  2.05631  2.06745  2.07174  2.08618  2.09779  2.1191
  2.13088  2.13946  2.15144  2.15649  2.18546  2.19785  2.21542  2.23006  2.2395   2.26516  2.28142  2.30145  2.31697  2.33698  2.35576  2.38487  2.39667  2.42117  2.43484  2.4679   2.48822  2.50067
  2.51927  2.53713  2.57086  2.59053  2.60401  2.62251  2.63914  2.66488  2.67206  3.00652  3.01895  3.03258  3.05769  3.08205  3.11844  3.15572  3.17469  3.18196  3.1978   3.22542  3.2823   3.29066
  3.31858  3.34892  3.38754  3.42366  3.45879  3.4941   3.52284  3.56076  3.60365  3.62981  3.66724  3.70507  3.72749  3.77265  3.80831  3.95101  4.02632  4.08071  4.1524   4.15974  4.19627  4.22066
  4.24838  4.28753  4.31757  4.35181  4.39025  4.42646  4.46933  4.50524  4.53998  4.58022  4.61718  4.65767  4.70092  4.73896  4.7683   4.81854  4.85665  4.8997   5.04433  5.08942  5.14579  5.17388
  5.21511  5.25238  5.29307  5.33358  5.37298  5.4152   5.45944  5.50117  5.54552  5.58988  5.63507  5.66858  5.7133   5.75746  5.7971   5.83569  5.87772  5.92013  6.04333  6.11407  6.16079  6.20292
  6.22752  6.2714   6.30507  6.34715  6.38662  6.42649  6.46598  6.50702  6.54169  6.58782  6.62439  6.66406  6.69988  6.73431  6.77391  6.82266  6.86231  7.04162  7.11151  7.15648  7.179    7.21715
  7.24604  7.28643  7.32286  7.36133  7.39996  7.43845  7.47428  7.51149  7.54372  7.58438  7.61691  7.64973  7.68577  7.70698  7.7401   7.76387  7.80523  8.02032  8.10464  8.15074  8.16776  8.18702
  8.21688  8.24734  8.2778   8.31338  8.34379  8.38136  8.40564  8.44244  8.47115  8.5056   8.53595  8.56435  8.58973  8.6184   8.65248  8.66648  8.69575  8.73155  8.77168  8.9833   9.03344  9.10202
  9.11403  9.15053  9.16556  9.1858   9.20972  9.22691  9.25197  9.2832   9.31187  9.31855  9.34754  9.36975  9.38804  9.40818  9.44053  9.45795  9.47474  9.50027  9.52074  9.52914  9.54549  9.56738
  9.57691  9.5966   9.603    9.61338  9.63976  9.67341  9.96023 10.1005  10.15368 10.15466 10.15887 10.16458 10.17305 10.17769 10.19349 10.2     10.20813 10.22276 10.22773 10.2428  10.25083 10.25733
 10.26525 10.2786  10.28835 10.29427 10.30556 10.31233 10.33276 10.33781 10.34176 10.34861 10.35657 10.36182 10.36682 10.37645 10.39009 10.39933 10.40506 10.42293 10.4325  10.43811 10.45636 10.47172
 10.48677 10.50227 10.52663 11.1534  11.15864 11.1728  11.19483 11.22259 11.25338 11.28418 11.31198 11.33406 11.34824 12.15312]
exact M2      = [2.600e-01 1.281e-28 6.287e-06 5.544e-05 1.480e-04 2.742e-04 4.218e-04 5.776e-04 3.061e-32 7.280e-04 3.851e-32 7.852e-31 8.577e-04 8.618e-31 9.490e-04 9.324e-04 7.366e-32 1.312e-31 3.097e-33
 6.728e-32 4.668e-32 2.507e-32 9.327e-32 1.466e-32 4.648e-32 1.734e-33 4.948e-33 5.706e-33 1.211e-32 1.399e-32 1.714e-32 2.401e-33 2.299e-33 5.499e-33 4.953e-34 1.089e-34 2.599e-33 2.565e-33
 1.190e-33 1.961e-33 5.029e-34 3.744e-33 1.226e-34 1.683e-32 6.448e-34 4.952e-33 5.554e-33 2.350e-01 4.061e-34 3.669e-33 1.098e-32 7.804e-33 2.460e-32 5.980e-33 4.818e-33 2.608e-32 1.241e-32
 2.016e-32 1.716e-32 1.479e-32 4.656e-33 5.334e-33 9.649e-33 1.899e-32 8.047e-33 2.706e-32 7.807e-33 5.193e-33 1.825e-33 2.564e-33 4.143e-33 4.140e-33 7.670e-34 4.568e-34 1.861e-33 2.285e-35
 8.844e-35 7.522e-35 3.962e-35 2.004e-35 2.685e-35 1.065e-31 5.355e-31 3.537e-32 1.520e-31 6.871e-32 3.399e-32 2.194e-32 1.633e-32 8.067e-33 5.379e-34 3.619e-34 2.027e-34 2.859e-34 3.539e-34
 5.842e-35 6.948e-35 1.061e-34 8.089e-35 5.547e-35 3.478e-34 9.559e-34 2.297e-34 4.687e-34 2.459e-34 5.937e-34 5.412e-33 9.771e-33 7.479e-33 1.623e-32 1.482e-32 1.711e-32 2.056e-32 2.423e-32
 9.620e-33 1.399e-32 9.347e-33 1.300e-32 8.293e-33 4.274e-33 4.073e-33 1.631e-33 1.048e-33 9.083e-34 6.710e-34 6.954e-34 4.225e-34 2.134e-34 1.981e-35 6.100e-35 2.784e-34 4.285e-34 4.587e-34
 9.319e-34 1.081e-33 1.183e-33 8.518e-34 8.314e-34 9.834e-34 9.255e-34 7.271e-34 6.594e-34 9.451e-34 9.084e-34 1.521e-34 3.361e-34 1.852e-34 1.417e-34 5.622e-35 3.100e-35 4.206e-36 4.545e-36
 1.264e-33 1.834e-33 2.640e-33 1.901e-33 3.058e-33 1.340e-33 1.833e-33 2.308e-33 8.579e-34 1.110e-33 5.941e-34 4.022e-34 4.229e-34 2.653e-34 6.425e-36 1.632e-34 1.154e-34 4.567e-36 9.872e-36
 9.372e-36 4.438e-35 4.511e-35 1.175e-35 3.581e-35 5.875e-35 4.614e-36 7.996e-35 2.728e-35 5.557e-35 1.594e-34 1.097e-35 8.124e-35 3.601e-35 2.796e-35 2.214e-34 6.101e-36 5.145e-36 4.438e-35
 2.750e-35 8.813e-35 6.098e-35 5.300e-36 5.473e-36 4.048e-35 9.038e-35 6.681e-35 1.028e-34 2.372e-34 2.835e-34 2.732e-34 2.053e-34 4.035e-34 5.801e-34 4.484e-34 4.309e-34 4.136e-34 3.094e-34
 9.949e-34 1.949e-34 1.576e-35 1.027e-34 1.011e-34 1.576e-34 2.961e-35 4.675e-35 4.327e-35 2.251e-34 1.573e-34 1.602e-34 7.865e-35 1.719e-34 9.021e-35 4.592e-34 9.901e-35 1.645e-34 7.582e-34
 9.340e-34 6.849e-34 2.948e-34 2.458e-34 9.544e-34 7.325e-35 5.048e-35 9.694e-35 2.164e-34 1.128e-34 3.576e-35 1.784e-34 3.674e-35 1.296e-34 1.622e-34 1.159e-34 8.583e-35 2.131e-34 4.807e-35
 1.034e-35 3.287e-35 3.454e-32 1.931e-35 3.299e-32 3.080e-35 1.195e-35 6.843e-35 2.609e-35 4.659e-35 7.382e-32 1.097e-35 2.602e-35 4.820e-35 7.959e-35 1.536e-34 6.731e-35 8.008e-36 5.226e-35
 1.420e-35 4.436e-35 4.671e-35 1.976e-34 2.306e-35 3.402e-35 4.337e-35 6.267e-35 1.257e-35 1.087e-35 5.772e-35 2.374e-35 4.296e-35 2.879e-35 1.054e-34 4.794e-35 8.633e-35 3.647e-36 1.794e-34
 1.242e-34 2.998e-35 4.819e-35 1.934e-34 1.546e-33 5.392e-34 5.441e-34 1.028e-33 5.991e-35 2.768e-33 4.703e-34 9.400e-35 7.470e-34 4.426e-34 5.585e-32]
ritz10 e (meV)= [ 0.144    1.01132  1.04276  1.08727  1.13477  1.17433  1.19676  2.15649  2.76719 10.39536]
ritz10 resid/tolabs = [3.07e-02 4.89e+05 8.75e+05 1.07e+06 1.03e+06 7.39e+05 2.75e+05 1.04e-01 5.96e+07 1.08e+07]
design_ok = True  emax (meV) = 1.0111158369538131  exact states in window = 4  converged Ritz states in window = 1  missed M2 weight = 6.172893572981009e-05
margin: lowest Ritz above window - emax = 0.0002067 meV, emax - highest exact in window = 0.0002067 meV, margin/emax = 0.000204
stopped dim = 11  reference dim = 301  stopped e (meV) = [ 0.144    1.01126  1.04258  1.08699  1.1345   1.17419  1.19674  2.13313  2.15649 10.74753]
order=2  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 1.210e-04   abs = 3.8785e-04   2*pi*missed = 3.8785e-04
order=3  max|dIdV(tol=1e-8) - dIdV(default)| / peak = 1.120e-04   abs = 5.1391e-04   2*pi*missed = 3.8785e-04
[exit 0, 12 s]
```

Reviewer's counter-repro, the band, the test chain and a better rule:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g03/g03_kr1_regime.py
```
```python
"""Reviewer counter-repro for kondo-routes 1 (the tol= stop of
get_kondo_spectrum(mode="ED", submode="KRYLOV")). The verdict does not
depend on this script. Part A checks the band of max|eV| argued by reading
(the uniform field shifts the one-magnon band by B, so the first-check Ritz
error of the k=0 magnon stays about 0.0112 meV) and records whether the tol
call warns. Part B asks whether the three-column AF test chain reaches the
regime, and Part C whether the fix the review proposes (a Ritz state is in
the window when theta - r <= top) still stops on the pinned cases.

Part A: uniform 12-spin FM chain, J=-1 meV, B=0.01 meV, site 0. The
first-check space is the code's own: state_spectrum(0, emax=-1.0) with tol
set has an empty window at the first check and returns that space.
Part B: heisenberg(8, "xz") of tests/test_kondo_krylov.py, site 1, T=0:
every max|eV| between consecutive exact final-state energies and
first-check Ritz values, those where the first check stops with weight
missing, and the public call at the worst one.
Part C: the stop schedule of block_krylov (CHECK_EVERY, CHECK_GROWTH)
replayed on spaces of fixed size (kks.nkry = blocks), with the current rule
and with the interval rule, on the two pinned cases of
tests/test_kondo_speedups_krylov.py (nkry=60, p_cut=1e-6: 0.1 K and 0.5
meV, documented stop at 121 states; T=0 and 0.3 meV, at 31).
"""
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
import warnings
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import krylovkondo, conductance

np.random.seed(1)
meV = 1e-3
TOL = 1e-8


def fm_chain(N, J, B):
    sc = spinchain.Spin_Chain(["S=1/2"]*N)
    h = 0
    for i in range(N-1):
        h = h + J*meV*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(N):
        h = h - B*meV*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def heisenberg_xz(n):
    """copy of tests/test_kondo_krylov.py::heisenberg(n, "xz")"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        J = 1e-3*(1.0 + 0.3*np.cos(i))
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + 2e-4*(np.cos(2*np.pi*i/3)*sc.Sx[i] + np.sin(2*np.pi*i/3)*sc.Sz[i])
    sc.set_hamiltonian(h)
    return sc


def compare(sc, site, T, emax, nkry, orders=(2,), extra=()):
    """max|dIdV(tol) - dIdV(default)| over a grid up to emax, relative to
    the peak of the default curve, and the warnings of the tol call."""
    pts = np.asarray([x for x in extra if 0. < x < emax])
    mids = 0.5*(pts + emax)
    eV = np.unique(np.concatenate([np.linspace(-emax, emax, 201), mids, -mids]))
    out = []
    for order in orders:
        kw = dict(site=site, Jrho_s=-0.05, T=T, order=order, mode="ED",
                  submode="KRYLOV", nkry=nkry)
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            _, a = sc.get_kondo_spectrum(eV, tol=TOL, **kw)
        msgs = [str(w.message)[:100] for w in wl]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, b = sc.get_kondo_spectrum(eV, **kw)
        d = np.max(np.abs(a - b))
        out.append((order, d/np.max(np.abs(b)), d, len(msgs), msgs[:1]))
    return out


print("==== Part A: uniform FM chain N=12, J=-1 meV, B=0.01 meV")
sc = fm_chain(12, -1.0, 0.01)
first = krylovkondo.KrylovKondoSpectrum(sc, 0, 0., nkry=100, tol=TOL)
a1 = first.state_spectrum(0, -1.0) # the code's own first-check space
o = np.argsort(a1.e[1:])
th, res = a1.e[1:][o], a1.residual[1:][o]
print("A first-check dim =", a1.dim, " lowest Ritz e (meV) = %.5f" % (th[0]/meV),
      " its residual = %.4g meV = %.3g tolabs" % (res[0]/meV, res[0]/(TOL*first.scale)))
print("A lowest Ritz - B (meV) = %.5f  (0.01117 at B=0.1 meV in the repro)" % ((th[0] - 0.01*meV)/meV))
for T, emax in ((0., 0.012), (0., 0.015), (0., 0.020), (0., 0.025),
                (0.002, 0.012), (0.010, 0.012)):
    stop = krylovkondo.KrylovKondoSpectrum(sc, 0, T, nkry=100, tol=TOL)
    dims = [stop.state_spectrum(n, emax*meV).dim for n in range(len(stop.weights))]
    top = stop.window(emax*meV)
    interval_in = bool(th[0] - res[0] <= top)
    for order, rel, d, nw, m in compare(sc, 0, T, emax*meV, 100, extra=(0.01*meV,)):
        print("A T=%.3f K emax=%.3f meV (emax/B=%.1f): n_initial=%d stopped dims=%s "
              "order=%d rel=%.3e abs=%.4e warnings=%d %s | theta1-top=%.4f meV, "
              "theta1-r1<=top: %s" % (T, emax, emax/0.01, len(stop.weights), dims, order,
                                      rel, d, nw, m, (th[0] - top)/meV, interval_in))

print("==== Part B: heisenberg(8, \"xz\") of the tests, site 1, T=0")
sc = heisenberg_xz(8)
ref = krylovkondo.KrylovKondoSpectrum(sc, 1, 0., nkry=100)
first = krylovkondo.KrylovKondoSpectrum(sc, 1, 0., nkry=100, tol=TOL)
tolabs = TOL*first.scale
print("B n_initial =", len(ref.weights), " dim =", ref.dim)
best = None
for n in range(len(ref.weights)):
    ss = ref.state_spectrum(n)
    e_ex = ss.e[1:]
    M2 = conductance._spin_matrix_elements_squared(ss)[1:, 0]
    a1 = first.state_spectrum(n, -1.0)
    th, res = a1.e[1:], a1.residual[1:]
    M2r = conductance._spin_matrix_elements_squared(a1)[1:, 0]
    total = float(np.sum(M2))
    print("B state %d: reference dim = %d (Hilbert dim %d), first-check dim = %d, "
          "converged first-check Ritz states = %d" % (n, ss.dim, ref.dim, a1.dim,
                                                      int(np.sum(res < tolabs))))
    wex = e_ex[M2 > 1e-12*total]
    pts = np.unique(np.concatenate([wex, th]))
    nfire = nmiss = 0
    for lo, hi in zip(pts[:-1], pts[1:]):
        c = 0.5*(lo + hi)
        if c <= 0.:
            continue
        inr = th <= c
        if not np.all(res[inr] < tolabs):
            continue
        nfire += 1
        missed = float(np.sum(M2[e_ex <= c]) - np.sum(M2r[inr]))
        if missed > 1e-6*total:
            nmiss += 1
            above = th[th > c]
            below = wex[wex <= c]
            band = ((np.min(above) if above.size else np.inf) - (np.max(below) if below.size else 0.))
            interval = bool(np.all(res[(th - res) <= c] < tolabs))
            if best is None or missed > best[1]:
                best = (c, missed, band, n, interval, total)
    print("B state %d: emax bands where the first check stops = %d, of them missing "
          "weight > 1e-6 of the total = %d" % (n, nfire, nmiss))
if best is None:
    print("B no emax with a first-check stop and a missing state")
else:
    c, missed, band, n, interval, total = best
    print("B worst: emax = %.5f meV, missed M2 = %.4e (total %.4e), band width = %.5f meV "
          "(band/emax = %.3g), state %d, interval rule also stops there: %s"
          % (c/meV, missed, total, band/meV, band/c, n, interval))
    for order, rel, d, nw, m in compare(sc, 1, 0., c, 100, orders=(2, 3)):
        print("B order=%d rel=%.3e abs=%.4e 2*pi*missed=%.4e warnings=%d"
              % (order, rel, d, 2*np.pi*missed, nw))

print("==== Part C: the stop schedule replayed, current rule and interval rule")


def replay(kks, k, top, rule, bmax):
    """the number of states at which block_krylov's schedule stops under a
    rule, or None if it reaches bmax blocks without stopping"""
    tolabs = TOL*kks.scale
    b = krylovkondo.CHECK_EVERY
    while b <= bmax:
        kks.nkry = b
        s = kks.state_spectrum(k) # kks.tol is None: a space of b blocks
        th, res = s.e[1:], s.residual[1:]
        win = (th <= top) if rule == "current" else ((th - res) <= top)
        conv = 1. if not np.any(win) else float(np.mean(res[win] < tolabs))
        if conv == 1.:
            return s.dim
        growth = krylovkondo.CHECK_GROWTH if conv >= 0.5 else 2.
        b = max(b + krylovkondo.CHECK_EVERY, int(np.ceil(growth*b)))
    return None


for T, emax, doc in ((0.1, 5e-4, 121), (0., 3e-4, 31)):
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, T, nkry=60, p_cut=1e-6)
    stopk = krylovkondo.KrylovKondoSpectrum(sc, 1, T, nkry=60, p_cut=1e-6, tol=TOL)
    for k in range(len(kks.weights)):
        Ei_top = kks.window(emax)
        actual = stopk.state_spectrum(k, emax).dim
        cur = replay(kks, k, Ei_top, "current", 59)
        itv = replay(kks, k, Ei_top, "interval", 59)
        print("C T=%.1f K emax=%.1e eV state %d: code stops at dim %d (documented %d), "
              "replayed current rule %s, interval rule %s (None: no stop before 181)"
              % (T, emax, k, actual, doc, cur, itv))
```

Observed:

```
==== Part A: uniform FM chain N=12, J=-1 meV, B=0.01 meV
A first-check dim = 11  lowest Ritz e (meV) = 0.02117  its residual = 0.03252 meV = 3.94e+05 tolabs
A lowest Ritz - B (meV) = 0.01117  (0.01117 at B=0.1 meV in the repro)
A T=0.000 K emax=0.012 meV (emax/B=1.2): n_initial=1 stopped dims=[11] order=2 rel=1.429e-01 abs=2.6180e-01 warnings=0 [] | theta1-top=0.0092 meV, theta1-r1<=top: True
A T=0.000 K emax=0.015 meV (emax/B=1.5): n_initial=1 stopped dims=[11] order=2 rel=1.429e-01 abs=2.6180e-01 warnings=0 [] | theta1-top=0.0062 meV, theta1-r1<=top: True
A T=0.000 K emax=0.020 meV (emax/B=2.0): n_initial=1 stopped dims=[11] order=2 rel=1.429e-01 abs=2.6180e-01 warnings=0 [] | theta1-top=0.0012 meV, theta1-r1<=top: True
A T=0.000 K emax=0.025 meV (emax/B=2.5): n_initial=1 stopped dims=[21] order=2 rel=2.157e-14 abs=3.9524e-14 warnings=0 [] | theta1-top=-0.0038 meV, theta1-r1<=top: True
A T=0.002 K emax=0.012 meV (emax/B=1.2): n_initial=1 stopped dims=[11] order=2 rel=1.428e-01 abs=2.6177e-01 warnings=0 [] | theta1-top=0.0023 meV, theta1-r1<=top: True
A T=0.010 K emax=0.012 meV (emax/B=1.2): n_initial=2 stopped dims=[21, 42] order=2 rel=2.993e-14 abs=5.3513e-14 warnings=0 [] | theta1-top=-0.0253 meV, theta1-r1<=top: True
==== Part B: heisenberg(8, "xz") of the tests, site 1, T=0
B n_initial = 1  dim = 256
B state 0: reference dim = 256 (Hilbert dim 256), first-check dim = 31, converged first-check Ritz states = 0
B state 0: emax bands where the first check stops = 1, of them missing weight > 1e-6 of the total = 1
B worst: emax = 0.31282 meV, missed M2 = 2.0596e-02 (total 7.3844e-01), band width = 0.00000 meV (band/emax = 2.24e-06), state 0, interval rule also stops there: False
B order=2 rel=6.405e-01 abs=1.2941e-01 2*pi*missed=1.2941e-01 warnings=0
B order=3 rel=6.212e-01 abs=1.5305e-01 2*pi*missed=1.2941e-01 warnings=0
==== Part C: the stop schedule replayed, current rule and interval rule
C T=0.1 K emax=5.0e-04 eV state 0: code stops at dim 121 (documented 121), replayed current rule 121, interval rule 121 (None: no stop before 181)
C T=0.0 K emax=3.0e-04 eV state 0: code stops at dim 31 (documented 31), replayed current rule 31, interval rule 31 (None: no stop before 181)
[exit 0, 14 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: the reference dimension "13", since the printed reference is the
301-state cap whose extra states carry no weight (that is the
reviewer-found performance finding below); "the documentation says
`tol=1e-8` moves dI/dV by about 1e-11", which is a measurement on the test
chain, not a promise; the promise that fails is the guide's. The
non-empty-window case is kept only as a knife edge of 1.4e-4 meV missing
1.2e-4 of the peak. MEDIUM, since the keyword is opt-in but the failure
drops a whole step silently.

**Suggested fix**: Not counting an empty window as converged is not enough, since a stop with a
non-empty window also misses weight. Count as part of the window every Ritz
state whose residual interval reaches it, `theta_j - r_j <= top`, and stop
only when all of those are converged; `r_j` is the residual the check
already forms. The counter-repro replays the stop schedule on the two pinned
cases of `tests/test_kondo_speedups_krylov.py` and finds this rule stops at
the documented sizes (121 and 31 states). Pin the fix with the one-magnon
chain, whose answer is analytic.

### 15. The ED Krylov Kondo route does not see an exhausted Krylov space when its start comes from eigsh, and builds 301 states where 13 suffice (64x the cost, same numbers)

`performance` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `kondo-routes` &middot; candidate found by the reviewer of `kondo-routes-1`

**Where**: `src/dmrgpy/kondospectrumtk/krylovkondo.py (block_krylov and its start states)`

On the polarized 12-spin ferromagnetic chain in a field, `S|GS>` reaches
only the 12 one-magnon states, yet the default `nkry=100` builds the full
300 Krylov states: `state_spectrum` has dimension 301 against 13 with the
exact basis vector as start, and 288 of the Ritz states carry `M2` below
1e-27. dI/dV is the same to 1e-15 of its peak, and `state_spectrum` plus
the third-order sums take 64 times as long (a ratio measured while other
sessions held the CPU). The eigsh start carries 1.8e-14
of roundoff outside its single component, and the exact start stops at 13,
so that roundoff keeps the recursion alive. The README and the Kondo rules
say the space stops when a block is empty.

**Expected**: the recursion stopping at the reachable space.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g03/g03_new_exhaustion_noise.py
```
```python
"""New candidate from the review of kondo-routes 1: block_krylov does not
see that the space is exhausted when the reachable space is small. On the
polarized 12-spin FM chain S|GS> reaches only the 12 one-magnon states, yet
the default call builds 300 Krylov states (the repro printed reference
dim = 301). Hypothesis: the eigsh ground state carries roundoff in other
Sz sectors, each Lanczos step multiplies it by about ||H||/beta, and after
12 steps it is far above DROP_TOL (1e-9), so the recursion goes on into
directions of weight about 1e-28. Test: the same call with the ground state
replaced by the exact basis vector it approximates.
"""
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
import time
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import krylovkondo, conductance

np.random.seed(1)
meV = 1e-3
N = 12
sc = spinchain.Spin_Chain(["S=1/2"]*N)
h = 0
for i in range(N-1):
    h = h - 1.0*meV*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
for i in range(N):
    h = h - 0.1*meV*sc.Sz[i]
sc.set_hamiltonian(h)
eV = np.linspace(-3e-3, 3e-3, 81) + 1.3e-5

kks = krylovkondo.KrylovKondoSpectrum(sc, 0, 0., nkry=100)
psi = kks.V[:, 0]/np.linalg.norm(kks.V[:, 0])
idx = int(np.argmax(np.abs(psi)))
rest = psi.copy(); rest[idx] = 0.
print("n_initial =", len(kks.weights), " dim =", kks.dim,
      " norm of psi outside its largest component = %.3e" % np.linalg.norm(rest))

kks.state_spectrum(0) # untimed warm-up
results = {}
for label in ("eigsh psi", "exact psi"):
    if label == "exact psi":
        v = np.zeros_like(kks.V[:, 0]); v[idx] = np.sign(psi[idx])
        kks.V[:, 0] = v
    t0 = time.perf_counter()
    ss = kks.state_spectrum(0)
    t1 = time.perf_counter()
    d = krylovkondo.krylov_kondo_dIdV(kks, eV, -0.05, order=3)
    t2 = time.perf_counter()
    M2 = conductance._spin_matrix_elements_squared(ss)[1:, 0]
    print("%s: state_spectrum dim = %d, Ritz states with M2 > 1e-20 = %d, "
          "largest M2 among the others = %.3e, state_spectrum %.3f s, dIdV(order=3) %.3f s"
          % (label, ss.dim, int(np.sum(M2 > 1e-20)),
             (np.max(M2[M2 <= 1e-20]) if np.any(M2 <= 1e-20) else 0.), t1 - t0, t2 - t1))
    results[label] = (ss.dim, d, t2 - t0)
da, db = results["eigsh psi"][1], results["exact psi"][1]
print("max|dIdV(eigsh psi) - dIdV(exact psi)| / peak = %.3e" % (np.max(np.abs(da - db))/np.max(np.abs(db))))
print("time ratio eigsh/exact = %.2f" % (results["eigsh psi"][2]/results["exact psi"][2]))
```

Observed:

```
n_initial = 1  dim = 4096  norm of psi outside its largest component = 1.841e-14
eigsh psi: state_spectrum dim = 301, Ritz states with M2 > 1e-20 = 12, largest M2 among the others = 5.501e-28, state_spectrum 0.602 s, dIdV(order=3) 1.126 s
exact psi: state_spectrum dim = 13, Ritz states with M2 > 1e-20 = 12, largest M2 among the others = 0.000e+00, state_spectrum 0.011 s, dIdV(order=3) 0.016 s
max|dIdV(eigsh psi) - dIdV(exact psi)| / peak = 1.043e-15
time ratio eigsh/exact = 63.64
[exit 0, 5 s]
```

**Reviewer (CONFIRMED)**: Struck: "the roundoff part grows by about `||H||` at each step", which was
not measured; the repro shows only that the start's roundoff causes it.

**Suggested fix**: When `[H, Sz_tot] = 0` (checked on a random state, as `su2_initial_states`
does) and the start has a definite Sz_tot, project the start and each
`S^a|i>` onto its own sector before `block_krylov`. That moves numbers only
at roundoff. Zeroing small entries instead depends on the basis and can cut
genuine small amplitudes. Pin it with this chain and assert a dimension of
at most 13.

### 16. The BiCGSTAB residual check forms A x and A x - b by density-matrix routes, so on an ill-conditioned converged solve it reports up to 1.0 for a true residual of 1.4e-7 and warns falsely, on `"python"`

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `pyitensor` &middot; candidate `pyitensor-2`

**Where**: `src/dmrgpy/mpsalgebra.py:430-468 (check_bicgstab; 447-448 the apply and the sum)`; `src/dmrgpy/mpsalgebra.py:277-298 (summps_dmrg, applyoperator_dmrg)`; `src/dmrgpy/pyitensor/chain.py:880-883, 907-908, 2448-2454 (apply_operator and sum_mps take APPLY_METHOD and SUM_METHOD)`; `src/dmrgpy/pyitensor/mpsalgebra.py:396, 967, 975-995 (both default to "densitymatrix")`; `docs/development.md and docs/numerical_changes.md (both say "summed exactly")`

`check_bicgstab` (`e7e1608`) says it forms A x and the difference exactly,
but it calls the session's `apply_operator` and `sum_mps`, which since
`c0c4d2d` take the density-matrix routes on `"python"`. For an
ill-conditioned `A = H - e` with `||x|| >> ||b||`, the shift-invert next to
a level, the residual it reports then has a floor unrelated to the true
residual. On 8 sites held exactly at cutoff 0, with x formed exactly and
`||x||/||b||` = 9e3, 9e5, 9e7, the dense residual is 1e-11, 1.1e-9 and
1.4e-7 and the check reads 2.1e-9, 4.4e-5 and 1.00. The result is a false
RuntimeWarning on a converged solve. On a real BiCGSTAB solve at `delta=1e-5`
and shift 1e-7 the dense residual is 5.0e-6 and the check reads 2.8e-2 and
warns; at the public `applyinverse` with `delta=1e-7` and shift 1e-5 the
dense residual is 5.0e-8 and the check reads 1.15e-6. The counter-repro
splits the routes: the density-matrix apply carries the floor, and the sum
does not add to it.

**Expected**: a residual within a small factor of the dense one, as the exact apply and sum give (2e-15 to 5e-15 on the same rows).

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/pyitensor/pyit_2_bicgstab_check_floor.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: mpsalgebra.check_bicgstab (e7e1608) says it measures
# ||A x - wf|| "with A applied and the difference summed exactly", but on
# itensor_version="python" both steps go to the session's default routes,
# the density-matrix apply and sum (c0c4d2d). Their floor is set by the
# channel magnitudes, about sqrt(eps) * (||H x|| + |e| ||x||), so for a
# well-solved but ill-conditioned A = H - e (shift-invert next to a level,
# ||x|| >> ||wf||) the check reports a residual many orders above the true
# one: a false RuntimeWarning, and a RuntimeError ("x is worse than x = 0")
# once the floor passes 1. Part A feeds the check an x whose true residual
# is roundoff (b = A x formed exactly); part B is the public applyinverse.
import warnings
from dmrgpy import spinchain, multioperator
from dmrgpy import mps as dmps
from dmrgpy import mpsalgebra as dma
from dmrgpy.pyitensor import mpsalgebra as ma

np.random.seed(2)
n = 8
sc = spinchain.Spin_Chain(["S=1/2"] * n)
sc.setup_python()
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
sc.set_hamiltonian(h)
sc.maxm = 32      # exact for 8 sites (bond dimension 16)
sc.nsweeps = 12
sc.cutoff = 0.0   # no truncation anywhere in the solves
sc.gs_energy()
gs = sc.get_gs()
gs = gs * (1. / np.sqrt(abs(gs.dot(gs))))
s = sc._session
s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
sites = s.sites
E0 = sc.vev(h, wf=gs).real
print("E0 (Rayleigh quotient) =", E0, flush=True)
print("APPLY_METHOD =", ma.APPLY_METHOD, " SUM_METHOD =", ma.SUM_METHOD, flush=True)


def vec(h_):
    T = h_.A(1)
    for i in range(2, h_.length() + 1):
        T = T * h_.A(i)
    arr = T.transpose_to([sites.si(i) for i in range(1, h_.length() + 1)])
    return np.asarray(arr).reshape(-1)


def mpo_dense(mpo):
    T = mpo.A(1)
    for i in range(2, mpo.length() + 1):
        T = T * mpo.A(i)
    order = ([sites.si(i).prime(1) for i in range(1, n + 1)]
             + [sites.si(i) for i in range(1, n + 1)])
    return np.asarray(T.transpose_to(order)).reshape(2**n, 2**n)


def with_routes(method, fn):
    old = (ma.SUM_METHOD, ma.APPLY_METHOD)
    ma.SUM_METHOD, ma.APPLY_METHOD = method, method
    try:
        return fn()
    finally:
        ma.SUM_METHOD, ma.APPLY_METHOD = old


def run_check(A, x, b, nb, delta, route):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            res = with_routes(route, lambda: dma.check_bicgstab(sc, A, x, b, nb, delta))
            verdict = "warned" if any("BiCGSTAB" in str(x_.message) for x_ in w) else "passed"
        except RuntimeError as err:
            res, verdict = float("nan"), "RAISED: " + str(err)[:90]
    return res, verdict


w1 = sc.Sx[0] * gs
delta = 1e-5   # cvm_tol, the default delta of applyinverse
print("Part A: x = gs/shift + Sx_0 gs, b = (H - e) x formed exactly, e = E0 - shift", flush=True)
for shift in (1e-2, 1e-4, 1e-6, 1e-8):
    e = E0 - shift
    A = h - e * multioperator.identity()
    Ampo = s._mpo(A.to_terms())
    x_h = ma.sum_many([gs.cpp_handle * complex(1. / shift), w1.cpp_handle * 1.0],
                      cutoff=0.0, maxdim=10**6, method="exact")
    b_h = ma.applyMPO(Ampo, x_h, cutoff=0.0, maxdim=None, method="exact")
    b_h.noPrime("Site")
    x = dmps.MPS(MBO=sc, cpp_handle=x_h)
    b = dmps.MPS(MBO=sc, cpp_handle=b_h)
    nb = np.sqrt(abs(b.dot(b)))
    Ad = mpo_dense(Ampo)
    vx, vb = vec(x_h), vec(b_h)
    res_dense = np.linalg.norm(Ad @ vx - vb) / np.linalg.norm(vb)
    r_dm, v_dm = run_check(A, x, b, nb, delta, "densitymatrix")
    r_ex, v_ex = run_check(A, x, b, nb, delta, "exact")
    print("  shift = %.0e  ||x||/||b|| = %.2e  res_dense = %.2e  "
          "check(default) = %.2e [%s]  check(exact) = %.2e [%s]"
          % (shift, np.linalg.norm(vx) / np.linalg.norm(vb), res_dense,
             r_dm, v_dm, r_ex, v_ex), flush=True)

print("Part B: BiCGSTAB solve of (H - e) x = b, delta=1e-7, maxn=400: the session solve, "
      "the check on both routes, and the verdict of the public applyinverse "
      "(inconclusive where res_dense > 1e-6: the solve itself did not converge)", flush=True)
sc.inverse_solver = "bicgstab"
dpub = 1e-7
for shift in (1e-3, 1e-5):
    e = E0 - shift
    A = h - e * multioperator.identity()
    b = gs + w1
    nb = np.sqrt(abs(b.dot(b)))
    s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    x_h = s.apply_inverse(A.to_terms(), b.cpp_handle, dpub * nb, 400)
    x = dmps.MPS(MBO=sc, cpp_handle=x_h)
    Ad = mpo_dense(s._mpo(A.to_terms()))
    vb = vec(b.cpp_handle)
    res_dense = np.linalg.norm(Ad @ vec(x_h) - vb) / np.linalg.norm(vb)
    r_dm, v_dm = run_check(A, x, b, nb, dpub, "densitymatrix")
    r_ex, v_ex = run_check(A, x, b, nb, dpub, "exact")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            sc.applyinverse(A, b, delta=dpub, maxn=400)
            pub = "warned" if any("BiCGSTAB" in str(x_.message) for x_ in w) else "passed"
        except RuntimeError as err:
            pub = "RAISED: " + str(err)[:90]
    print("  shift = %.0e  ||x||/||b|| = %.2e  res_dense = %.2e  check(default) = %.2e [%s]  "
          "check(exact) = %.2e [%s]  public applyinverse: %s"
          % (shift, np.linalg.norm(vec(x_h)) / np.linalg.norm(vb), res_dense,
             r_dm, v_dm, r_ex, v_ex, pub), flush=True)
```

Observed:

```
E0 (Rayleigh quotient) = -3.3749325986878915
APPLY_METHOD = densitymatrix  SUM_METHOD = densitymatrix
Part A: x = gs/shift + Sx_0 gs, b = (H - e) x formed exactly, e = E0 - shift
  shift = 1e-02  ||x||/||b|| = 8.98e+01  res_dense = 1.23e-13  check(default) = 1.02e-14 [passed]  check(exact) = 4.45e-15 [passed]
  shift = 1e-04  ||x||/||b|| = 9.00e+03  res_dense = 1.00e-11  check(default) = 2.12e-09 [passed]  check(exact) = 2.28e-15 [passed]
  shift = 1e-06  ||x||/||b|| = 9.00e+05  res_dense = 1.13e-09  check(default) = 4.37e-05 [passed]  check(exact) = 2.79e-15 [passed]
  shift = 1e-08  ||x||/||b|| = 9.00e+07  res_dense = 1.42e-07  check(default) = 1.00e+00 [warned]  check(exact) = 3.49e-15 [passed]
Part B: BiCGSTAB solve of (H - e) x = b, delta=1e-7, maxn=400: the session solve, the check on both routes, and the verdict of the public applyinverse (inconclusive where res_dense > 1e-6: the solve itself did not converge)
  shift = 1e-03  ||x||/||b|| = 8.94e+02  res_dense = 5.21e-08  check(default) = 5.20e-08 [passed]  check(exact) = 5.21e-08 [passed]  public applyinverse: passed
  shift = 1e-05  ||x||/||b|| = 8.94e+04  res_dense = 4.99e-08  check(default) = 1.15e-06 [warned]  check(exact) = 4.99e-08 [passed]  public applyinverse: warned
[exit 0, 6 s]
```

Reviewer's counter-repro, the routes split and real solves:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g08/g08_routes_and_raise.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for pyitensor-2. Two questions the hunter's output does not
# settle: (1) which step of check_bicgstab carries the floor on "python",
# the density-matrix APPLY of A x or the density-matrix SUM of A x - wf
# (this decides the fix); (2) does the check ever RAISE ("x is worse than
# x = 0") on a correct x, as the claim says, or only warn (at shift 1e-8
# it read 1.00 and warned). Same construction as the hunter's Part A: x
# = gs/shift + Sx_0 gs, b = (H - e) x formed by the exact route, so the
# true residual is roundoff, bounded by res_dense (a dense numpy product).
import warnings
from dmrgpy import spinchain, multioperator
from dmrgpy import mps as dmps
from dmrgpy import mpsalgebra as dma
from dmrgpy.pyitensor import mpsalgebra as ma

np.random.seed(2)
n = 8
sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
sc.set_hamiltonian(h)
sc.maxm = 32
sc.nsweeps = 12
sc.cutoff = 0.0
sc.gs_energy()
gs = sc.get_gs()
gs = gs * (1. / np.sqrt(abs(gs.dot(gs))))
s = sc._session
s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
sites = s.sites
E0 = sc.vev(h, wf=gs).real
print("E0 =", E0, flush=True)


def vec(h_):
    T = h_.A(1)
    for i in range(2, h_.length() + 1):
        T = T * h_.A(i)
    return np.asarray(T.transpose_to([sites.si(i) for i in range(1, n + 1)])).reshape(-1)


def mpo_dense(mpo):
    T = mpo.A(1)
    for i in range(2, mpo.length() + 1):
        T = T * mpo.A(i)
    order = ([sites.si(i).prime(1) for i in range(1, n + 1)]
             + [sites.si(i) for i in range(1, n + 1)])
    return np.asarray(T.transpose_to(order)).reshape(2**n, 2**n)


def run_check(A, x, b, nb, delta, apply_m, sum_m):
    # BICGSTAB_RAISE lifted so the value is returned even where the check
    # would raise; the verdict the real check gives is then read off res:
    # res > 1 is "would RAISE", res > 10*delta "warned", else "passed"
    old = (ma.APPLY_METHOD, ma.SUM_METHOD, dma.BICGSTAB_RAISE)
    ma.APPLY_METHOD, ma.SUM_METHOD, dma.BICGSTAB_RAISE = apply_m, sum_m, np.inf
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            res = dma.check_bicgstab(sc, A, x, b, nb, delta)
            verdict = ("would RAISE" if res > 1. else
                       "warned" if any("BiCGSTAB" in str(x_.message) for x_ in w) else "passed")
        except RuntimeError as err:
            res, verdict = float("nan"), "RAISED (not finite)"
        finally:
            ma.APPLY_METHOD, ma.SUM_METHOD, dma.BICGSTAB_RAISE = old
    return res, verdict


w1 = sc.Sx[0] * gs
delta = 1e-5
print("Part A: check by route (apply/sum: de = densitymatrix, ex = exact), delta=1e-5", flush=True)
for shift in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10):
    e = E0 - shift
    A = h - e * multioperator.identity()
    Ampo = s._mpo(A.to_terms())
    x_h = ma.sum_many([gs.cpp_handle * complex(1. / shift), w1.cpp_handle * 1.0],
                      cutoff=0.0, maxdim=10**6, method="exact")
    b_h = ma.applyMPO(Ampo, x_h, cutoff=0.0, maxdim=None, method="exact")
    b_h.noPrime("Site")
    x = dmps.MPS(MBO=sc, cpp_handle=x_h)
    b = dmps.MPS(MBO=sc, cpp_handle=b_h)
    nb = np.sqrt(abs(b.dot(b)))
    vx, vb = vec(x_h), vec(b_h)
    res_dense = np.linalg.norm(mpo_dense(Ampo) @ vx - vb) / np.linalg.norm(vb)
    ratio = np.linalg.norm(vx) / np.linalg.norm(vb)
    row = []
    for apply_m, sum_m in (("densitymatrix", "densitymatrix"), ("densitymatrix", "exact"),
                           ("exact", "densitymatrix"), ("exact", "exact")):
        r, v = run_check(A, x, b, nb, delta, apply_m, sum_m)
        if apply_m == sum_m == "densitymatrix":
            law = r / ratio**2
        row.append("%s/%s=%.3e(res-1=%+.1e)[%s]" % (apply_m[:2], sum_m[:2], r, r - 1., v))
    print("  shift=%.0e ||x||/||b||=%.2e res_dense=%.2e de/de/(||x||/||b||)^2=%.2e  "
          % (shift, ratio, res_dense, law) + "  ".join(row), flush=True)

# Part C: the same session solve as the hunter's Part B, but at the
# library default cutoff=1e-12 instead of cutoff=0. maxm=30 still holds
# 8 sites exactly (bond dimension 16), so this tests the default cutoff
# only, not truncation by maxm. Question: at the default cutoff, does the
# density-matrix floor change a verdict, or is the true (dense) residual
# of a near-level solve already above 10*delta wherever the check warns?
print("Part C: session BiCGSTAB of (H - e) x = gs + Sx_0 gs at cutoff=1e-12, maxm=30 (exact for 8 sites), maxn=400", flush=True)
sc.cutoff, sc.maxm = 1e-12, 30
b = gs + w1
nb = np.sqrt(abs(b.dot(b)))
vb = vec(b.cpp_handle)
for dl in (1e-5, 1e-7):
    for shift in (1e-3, 1e-5, 1e-7):
        e = E0 - shift
        A = h - e * multioperator.identity()
        s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
        x_h = s.apply_inverse(A.to_terms(), b.cpp_handle, dl * nb, 400)
        x = dmps.MPS(MBO=sc, cpp_handle=x_h)
        res_dense = np.linalg.norm(mpo_dense(s._mpo(A.to_terms())) @ vec(x_h) - vb) / np.linalg.norm(vb)
        r_dm, v_dm = run_check(A, x, b, nb, dl, "densitymatrix", "densitymatrix")
        r_ex, v_ex = run_check(A, x, b, nb, dl, "exact", "exact")
        print("  delta=%.0e shift=%.0e ||x||/||b||=%.2e res_dense=%.2e check(default)=%.2e[%s] check(exact)=%.2e[%s]"
              % (dl, shift, np.linalg.norm(vec(x_h)) / np.linalg.norm(vb), res_dense, r_dm, v_dm, r_ex, v_ex),
              flush=True)
```

Observed:

```
E0 = -3.3749325986878915
Part A: check by route (apply/sum: de = densitymatrix, ex = exact), delta=1e-5
  shift=1e-03 ||x||/||b||=9.00e+02 res_dense=1.25e-12 de/de/(||x||/||b||)^2=9.45e-18  de/de=7.647e-12(res-1=-1.0e+00)[passed]  de/ex=7.642e-11(res-1=-1.0e+00)[passed]  ex/de=3.237e-15(res-1=-1.0e+00)[passed]  ex/ex=8.125e-15(res-1=-1.0e+00)[passed]
  shift=1e-04 ||x||/||b||=9.00e+03 res_dense=1.00e-11 de/de/(||x||/||b||)^2=2.62e-17  de/de=2.122e-09(res-1=-1.0e+00)[passed]  de/ex=6.460e-09(res-1=-1.0e+00)[passed]  ex/de=2.625e-16(res-1=-1.0e+00)[passed]  ex/ex=2.277e-15(res-1=-1.0e+00)[passed]
  shift=1e-05 ||x||/||b||=9.00e+04 res_dense=9.75e-11 de/de/(||x||/||b||)^2=5.38e-17  de/de=4.358e-07(res-1=-1.0e+00)[passed]  de/ex=4.358e-07(res-1=-1.0e+00)[passed]  ex/de=4.329e-16(res-1=-1.0e+00)[passed]  ex/ex=2.210e-15(res-1=-1.0e+00)[passed]
  shift=1e-06 ||x||/||b||=9.00e+05 res_dense=1.13e-09 de/de/(||x||/||b||)^2=5.40e-17  de/de=4.374e-05(res-1=-1.0e+00)[passed]  de/ex=4.374e-05(res-1=-1.0e+00)[passed]  ex/de=9.127e-17(res-1=-1.0e+00)[passed]  ex/ex=2.789e-15(res-1=-1.0e+00)[passed]
  shift=1e-07 ||x||/||b||=9.00e+06 res_dense=1.11e-08 de/de/(||x||/||b||)^2=3.33e-16  de/de=2.696e-02(res-1=-9.7e-01)[warned]  de/ex=2.696e-02(res-1=-9.7e-01)[warned]  ex/de=9.544e-16(res-1=-1.0e+00)[passed]  ex/ex=3.446e-15(res-1=-1.0e+00)[passed]
  shift=1e-08 ||x||/||b||=9.00e+07 res_dense=1.42e-07 de/de/(||x||/||b||)^2=1.24e-16  de/de=9.999e-01(res-1=-1.0e-04)[warned]  de/ex=9.999e-01(res-1=-1.0e-04)[warned]  ex/de=1.295e-17(res-1=-1.0e+00)[passed]  ex/ex=3.487e-15(res-1=-1.0e+00)[passed]
  shift=1e-09 ||x||/||b||=9.00e+08 res_dense=1.30e-06 de/de/(||x||/||b||)^2=1.23e-18  de/de=1.000e+00(res-1=-4.4e-05)[warned]  de/ex=1.000e+00(res-1=-4.4e-05)[warned]  ex/de=2.337e-16(res-1=-1.0e+00)[passed]  ex/ex=5.504e-15(res-1=-1.0e+00)[passed]
  shift=1e-10 ||x||/||b||=9.01e+09 res_dense=1.21e-05 de/de/(||x||/||b||)^2=1.23e-20  de/de=9.999e-01(res-1=-1.3e-04)[warned]  de/ex=9.999e-01(res-1=-1.3e-04)[warned]  ex/de=1.194e-16(res-1=-1.0e+00)[passed]  ex/ex=2.184e-15(res-1=-1.0e+00)[passed]
Part C: session BiCGSTAB of (H - e) x = gs + Sx_0 gs at cutoff=1e-12, maxm=30 (exact for 8 sites), maxn=400
  delta=1e-05 shift=1e-03 ||x||/||b||=8.94e+02 res_dense=1.16e-05 check(default)=1.16e-05[passed] check(exact)=1.16e-05[passed]
  delta=1e-05 shift=1e-05 ||x||/||b||=8.94e+04 res_dense=6.70e-06 check(default)=6.70e-06[passed] check(exact)=6.70e-06[passed]
  delta=1e-05 shift=1e-07 ||x||/||b||=8.94e+06 res_dense=5.01e-06 check(default)=2.80e-02[warned] check(exact)=5.01e-06[passed]
  delta=1e-07 shift=1e-03 ||x||/||b||=8.94e+02 res_dense=3.14e-07 check(default)=3.14e-07[passed] check(exact)=3.14e-07[passed]
  delta=1e-07 shift=1e-05 ||x||/||b||=8.94e+04 res_dense=6.10e-07 check(default)=1.17e-06[warned] check(exact)=6.10e-07[passed]
  delta=1e-07 shift=1e-07 ||x||/||b||=8.94e+06 res_dense=3.50e-07 check(default)=1.66e-02[warned] check(exact)=3.54e-07[passed]
[exit 0, 9 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: the size law `sqrt(eps)|e| ||x||/||b||`, since the output does not
follow it; the measured numbers stand in its place, and they do not follow
`eps (||x||/||b||)^2` within a factor of 10 either. Struck: "near a level a
RuntimeError", since the saturated rows read 1.0 and warn, and a raise needs
more than 1; and "the check can in principle pass a bad solve", which stays a
lead. MEDIUM by the counter-repro's rule: at `delta=1e-5` a real solve gets a
different verdict from the default check than from the exact one.

**Suggested fix**: On `"python"` reuse what `17017a9` already wrote: `pyitensor/implicit.py`'s
`_residual(A, x, b)` applies and sums with `method="exact"` at cutoff 0 and
reads the norm off the orthogonality centre. Wrap it as a session method
that `check_bicgstab` uses when the session has it. The apply must move;
the sum may. On `julia_live` the check should take
`apply(A, x; alg="naive", truncate=false)` and an exact direct sum. The v3
side is the reviewer-found extension below. Correct the "summed exactly" in
the two docs.

### 17. Extension of the BiCGSTAB check finding above to v3, the default backend: the check overstates the true residual up to 2e6 times, and the public `applyinverse` warns at 2.6e-3 on a solve whose residual is 1.2e-4 at delta=1e-4

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `pyitensor` &middot; candidate found by the reviewer of `pyitensor-2`

**Where**: `src/dmrgpy/mpsalgebra.py:430-468 (check_bicgstab)`; `src/dmrgpy/mpscpp3/chain_session.h (Chain::apply_operator passes no Method to applyMPO, whose default is DensityMatrix)`

On v3 the check's A x also comes from a density-matrix apply. On an exactly
formed pair (x, b = A x) on 8 spins it overstates the true residual about
5e3 times at `||x||/||b||` = 9e3 (3.6e-7 against 6.9e-11), 5e5 times at 9e5
(1.2e-2 against 2.5e-8, a false warning at `delta=1e-5`) and 2e6 times at
9e7 (0.54 against 2.4e-7, half way to the RuntimeError). The counter-repro
reaches the public API: at shift 1e-6 and `delta=1e-4` the v3 BiCGSTAB
solve has a dense residual of 1.19e-4, below the warning threshold of 10
delta, and the check reads 2.58e-3 and warns. It also splits the steps: the
apply alone is 2.58e-3 off, the sum 1e-12, so the apply is the floor.

**Expected**: a residual within a small factor of the dense one.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g08/g08_v3_check_floor.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# New candidate from the review of pyitensor-2: on itensor_version=3, the
# default backend, check_bicgstab applies A with Chain::apply_operator,
# which calls ITensor v3's applyMPO(K,x,{Cutoff,MaxDim}) with no Method,
# and mpoalgs.cc:135 takes "DensityMatrix" then. The same floor as on
# "python" is therefore expected on v3. Part A: the hunter's exactly
# formed pair (x, b = A x), built on "python" by the exact route and moved
# into a v3 session through mps_site_arrays/mps_from_site_arrays; the
# anchor is a dense numpy product with v3's OWN MPO (mpo_site_arrays) and
# v3's own copies of x and b. Part B: v3's own BiCGSTAB solve of the same
# (H - e) x = gs + Sx_0 gs, its dense residual, the check, and the public
# applyinverse verdict.
import warnings
from dmrgpy import spinchain, multioperator
from dmrgpy import mps as dmps
from dmrgpy import mpsalgebra as dma
from dmrgpy.pyitensor import mpsalgebra as ma

np.random.seed(2)
n = 8


def heis(c):
    h = 0
    for i in range(n - 1):
        h = h + c.Sx[i] * c.Sx[i + 1] + c.Sy[i] * c.Sy[i + 1] + c.Sz[i] * c.Sz[i + 1]
    return h


sp = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
hp = heis(sp)
sp.set_hamiltonian(hp)
sp.maxm, sp.nsweeps, sp.cutoff = 32, 12, 0.0
sp.gs_energy()
gs = sp.get_gs()
gs = gs * (1. / np.sqrt(abs(gs.dot(gs))))
spy = sp._session
spy.set_sweep_params(sp.maxm, sp.nsweeps, sp.cutoff, sp.noise)
E0 = sp.vev(hp, wf=gs).real
print("E0 (python, Rayleigh quotient) =", E0, flush=True)

s3c = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=3)
print("v3 chain itensor_version =", s3c.itensor_version, flush=True)
assert s3c.itensor_version == 3 and getattr(s3c, "_session", None) is not None, \
    "no compiled v3 session: the script measures nothing"
print("v3 session type =", type(s3c._session).__module__, flush=True)
h3 = heis(s3c)
s3c.set_hamiltonian(h3)
s3c.maxm, s3c.nsweeps, s3c.cutoff = 32, 12, 0.0
s3 = s3c._session
s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)


def to_v3(handle_py):
    arrs = [np.ascontiguousarray(np.asarray(a), dtype=complex) for a in spy.mps_site_arrays(handle_py)]
    return dmps.MPS(MBO=s3c, cpp_handle=s3.mps_from_site_arrays(arrs))


def vec_arrays(arrs):
    arrs = [np.asarray(a) for a in arrs]
    v = arrs[0].reshape(arrs[0].shape[1], arrs[0].shape[2])
    for a in arrs[1:]:
        v = np.tensordot(v, a, axes=([1], [0])).reshape(-1, a.shape[2])
    return v.reshape(-1)


def mpo_dense_arrays(Ws):
    Ws = [np.asarray(W) for W in Ws]
    M = Ws[0][0]  # (s', s, r)
    for W in Ws[1:]:
        T = np.tensordot(M, W, axes=([2], [0]))  # (R, C, s', s, r)
        R, C = M.shape[0], M.shape[1]
        M = T.transpose(0, 2, 1, 3, 4).reshape(R * W.shape[1], C * W.shape[2], W.shape[3])
    return M[:, :, 0]


def dense_residual3(A, x3, b3):
    Ad = mpo_dense_arrays(s3.mpo_site_arrays(s3.build_operator(A.to_terms())))
    vx = vec_arrays(s3.mps_site_arrays(x3.cpp_handle))
    vb = vec_arrays(s3.mps_site_arrays(b3.cpp_handle))
    return np.linalg.norm(Ad @ vx - vb) / np.linalg.norm(vb), np.linalg.norm(vx) / np.linalg.norm(vb)


def run_check(chain, A, x, b, nb, delta):
    # BICGSTAB_RAISE lifted so the value is returned where the check would
    # raise: res > 1 is "would RAISE", res > 10*delta "warned"
    old = dma.BICGSTAB_RAISE
    dma.BICGSTAB_RAISE = np.inf
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            res = dma.check_bicgstab(chain, A, x, b, nb, delta)
            verdict = ("would RAISE" if res > 1. else
                       "warned" if any("BiCGSTAB" in str(x_.message) for x_ in w) else "passed")
        except RuntimeError:
            res, verdict = float("nan"), "RAISED (not finite)"
        finally:
            dma.BICGSTAB_RAISE = old
    return res, verdict


w1p = sp.Sx[0] * gs
delta = 1e-5
print("Part A (v3 check on an exactly formed pair, delta=1e-5)", flush=True)
for shift in (1e-2, 1e-4, 1e-6, 1e-8):
    e = E0 - shift
    Ap = hp - e * multioperator.identity()
    A3 = h3 - e * multioperator.identity()
    x_h = ma.sum_many([gs.cpp_handle * complex(1. / shift), w1p.cpp_handle * 1.0],
                      cutoff=0.0, maxdim=10**6, method="exact")
    b_h = ma.applyMPO(spy._mpo(Ap.to_terms()), x_h, cutoff=0.0, maxdim=None, method="exact")
    b_h.noPrime("Site")
    x3, b3 = to_v3(x_h), to_v3(b_h)
    nb3 = np.sqrt(abs(b3.dot(b3)))
    rd3, ratio = dense_residual3(A3, x3, b3)
    r3, v3 = run_check(s3c, A3, x3, b3, nb3, delta)
    print("  shift=%.0e ||x||/||b||=%.2e res_dense_v3=%.2e check_v3=%.2e [%s]"
          % (shift, ratio, rd3, r3, v3), flush=True)

print("Part B (v3 BiCGSTAB of (H-e) x = gs + Sx_0 gs, delta=1e-7, maxn=400)", flush=True)
dpub = 1e-7
gs3 = to_v3(gs.cpp_handle)
b3 = to_v3((gs + w1p).cpp_handle)
nb3 = np.sqrt(abs(b3.dot(b3)))
for shift in (1e-3, 1e-5):
    e = E0 - shift
    A3 = h3 - e * multioperator.identity()
    s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
    x3 = dmps.MPS(MBO=s3c, cpp_handle=s3.apply_inverse(A3.to_terms(), b3.cpp_handle, dpub * nb3, 400))
    rd3, ratio = dense_residual3(A3, x3, b3)
    r3, v3 = run_check(s3c, A3, x3, b3, nb3, dpub)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            s3c.applyinverse(A3, b3, delta=dpub, maxn=400)
            pub = "warned" if any("BiCGSTAB" in str(x_.message) for x_ in w) else "passed"
        except RuntimeError:
            pub = "RAISED"
    print("  shift=%.0e ||x||/||b||=%.2e res_dense_v3=%.2e check_v3=%.2e [%s] public applyinverse: %s"
          % (shift, ratio, rd3, r3, v3, pub), flush=True)
```

Observed:

```
E0 (python, Rayleigh quotient) = -3.3749325986878915
v3 chain itensor_version = 3
v3 session type = dmrgpy.mpscpp3._dmrgcpp
Part A (v3 check on an exactly formed pair, delta=1e-5)
  shift=1e-02 ||x||/||b||=8.98e+01 res_dense_v3=2.48e-12 check_v3=3.39e-12 [passed]
  shift=1e-04 ||x||/||b||=9.00e+03 res_dense_v3=6.85e-11 check_v3=3.57e-07 [passed]
  shift=1e-06 ||x||/||b||=9.00e+05 res_dense_v3=2.51e-08 check_v3=1.20e-02 [warned]
  shift=1e-08 ||x||/||b||=9.00e+07 res_dense_v3=2.44e-07 check_v3=5.42e-01 [warned]
Part B (v3 BiCGSTAB of (H-e) x = gs + Sx_0 gs, delta=1e-7, maxn=400)
  shift=1e-03 ||x||/||b||=8.94e+02 res_dense_v3=5.21e-08 check_v3=4.95e-08 [passed] public applyinverse: passed
  shift=1e-05 ||x||/||b||=8.94e+04 res_dense_v3=1.02e-06 check_v3=7.57e-06 [warned] public applyinverse: warned
[exit 0, 7 s]
```

Reviewer's counter-repro, the public applyinverse and the steps split:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g12/g12_g08_public_false_warning.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for new-g08 (review group 12). Part A of the candidate shows
# the v3 check overstating the residual of an exactly formed pair; its Part B
# (the public applyinverse on v3) showed only an honest warning. This asks
# whether the PUBLIC v3 route ever gives a false warning: v3's own BiCGSTAB
# solve, its dense residual against v3's own MPO, the check, and the public
# applyinverse verdict, at shifts and deltas where a converged solve would sit
# below 10*delta while the check's floor (Part A: 3.6e-7 at ||x||/||b|| = 9e3,
# 1.2e-2 at 9e5) would sit above it. 8 spins, maxm=32: no truncation at all.
import warnings
from dmrgpy import spinchain, multioperator
from dmrgpy import mps as dmps
from dmrgpy import mpsalgebra as dma

np.random.seed(2)
n = 8


def heis(c):
    h = 0
    for i in range(n - 1):
        h = h + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1] + c.Sz[i]*c.Sz[i+1]
    return h


s3c = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=3)
assert s3c.itensor_version == 3 and getattr(s3c, "_session", None) is not None, \
    "no compiled v3 session: the script measures nothing"
h3 = heis(s3c)
s3c.set_hamiltonian(h3)
s3c.maxm, s3c.nsweeps, s3c.cutoff = 32, 12, 0.0
s3 = s3c._session
s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
s3c.gs_energy()
gs = s3c.get_gs()
gs = gs*(1./np.sqrt(abs(gs.dot(gs))))
E0 = s3c.vev(h3, wf=gs).real
print("E0 (v3, Rayleigh quotient) = %.15f" % E0, flush=True)
b3 = gs + s3c.Sx[0]*gs
nb3 = np.sqrt(abs(b3.dot(b3)))


def vec_arrays(arrs):
    arrs = [np.asarray(a) for a in arrs]
    v = arrs[0].reshape(arrs[0].shape[1], arrs[0].shape[2])
    for a in arrs[1:]:
        v = np.tensordot(v, a, axes=([1], [0])).reshape(-1, a.shape[2])
    return v.reshape(-1)


def mpo_dense_arrays(Ws):
    Ws = [np.asarray(W) for W in Ws]
    M = Ws[0][0]
    for W in Ws[1:]:
        T = np.tensordot(M, W, axes=([2], [0]))
        R, C = M.shape[0], M.shape[1]
        M = T.transpose(0, 2, 1, 3, 4).reshape(R*W.shape[1], C*W.shape[2], W.shape[3])
    return M[:, :, 0]


def dense_residual3(A, x3, b3):
    Ad = mpo_dense_arrays(s3.mpo_site_arrays(s3.build_operator(A.to_terms())))
    vx = vec_arrays(s3.mps_site_arrays(x3.cpp_handle))
    vb = vec_arrays(s3.mps_site_arrays(b3.cpp_handle))
    return (np.linalg.norm(Ad @ vx - vb)/np.linalg.norm(vb),
            np.linalg.norm(vx)/np.linalg.norm(vb))


def run_check(A, x, delta):
    old = dma.BICGSTAB_RAISE
    dma.BICGSTAB_RAISE = np.inf
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            res = dma.check_bicgstab(s3c, A, x, b3, nb3, delta)
            verdict = ("would RAISE" if res > 1. else
                       "warned" if any("BiCGSTAB" in str(m.message) for m in w)
                       else "passed")
        except RuntimeError:
            res, verdict = float("nan"), "RAISED (not finite)"
        finally:
            dma.BICGSTAB_RAISE = old
    return res, verdict


def step_errors(A, x):
    """Where the check loses accuracy: the session's A x (the same call and
    settings check_bicgstab uses: maxm uncapped, cutoff 0) against the dense
    product of v3's own MPO, and the session's sum A x - b against the dense
    difference of its own two inputs, both relative to ||b||."""
    from dmrgpy import vev
    old, oldcut = s3c.maxm, s3c.cutoff
    s3c.maxm, s3c.cutoff = vev._UNCAPPED, 0.
    try:
        ax = dma.applyoperator_dmrg(s3c, A, x)
        r = dma.summps_dmrg(s3c, ax, b3*(-1.))
    finally:
        s3c.maxm, s3c.cutoff = old, oldcut
        s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
    Ad = mpo_dense_arrays(s3.mpo_site_arrays(s3.build_operator(A.to_terms())))
    vx = vec_arrays(s3.mps_site_arrays(x.cpp_handle))
    vb = vec_arrays(s3.mps_site_arrays(b3.cpp_handle))
    vax = vec_arrays(s3.mps_site_arrays(ax.cpp_handle))
    vr = vec_arrays(s3.mps_site_arrays(r.cpp_handle))
    nb = np.linalg.norm(vb)
    return (np.linalg.norm(vax - Ad @ vx)/nb,
            np.linalg.norm(vr - (vax - vb))/nb)


for shift in (1e-4, 1e-6, 1e-7):
    e = E0 - shift
    A3 = h3 - e*multioperator.identity()
    for delta in (1e-4, 1e-5, 1e-6):
        s3.set_sweep_params(s3c.maxm, s3c.nsweeps, s3c.cutoff, s3c.noise)
        s3.set_mpomaxm(max(s3c.maxm, s3c.mpomaxm))
        x3 = dmps.MPS(MBO=s3c, cpp_handle=s3.apply_inverse(
            A3.to_terms(), b3.cpp_handle, delta*nb3, 400))
        rd3, ratio = dense_residual3(A3, x3, b3)
        r3, v3 = run_check(A3, x3, delta)
        apply_err, sum_err = step_errors(A3, x3)
        print("STEP shift=%.0e delta=%.0e apply_err=%.2e sum_err=%.2e"
              % (shift, delta, apply_err, sum_err), flush=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                s3c.applyinverse(A3, b3, delta=delta, maxn=400)
                pub = ("warned" if any("BiCGSTAB" in str(m.message) for m in w)
                       else "passed")
            except RuntimeError:
                pub = "RAISED"
        honest = "dense<=10delta" if rd3 <= 10.*delta else "dense>10delta"
        print("ROW shift=%.0e delta=%.0e ||x||/||b||=%.2e res_dense_v3=%.2e "
              "check_v3=%.2e [%s] public=%s %s"
              % (shift, delta, ratio, rd3, r3, v3, pub, honest), flush=True)
print("done")
```

Observed:

```
E0 (v3, Rayleigh quotient) = -3.374932598687892
STEP shift=1e-04 delta=1e-04 apply_err=4.37e-08 sum_err=3.60e-12
ROW shift=1e-04 delta=1e-04 ||x||/||b||=8.94e+03 res_dense_v3=7.25e-05 check_v3=7.25e-05 [passed] public=passed dense<=10delta
STEP shift=1e-04 delta=1e-05 apply_err=8.99e-08 sum_err=4.49e-11
ROW shift=1e-04 delta=1e-05 ||x||/||b||=8.94e+03 res_dense_v3=5.39e-06 check_v3=5.40e-06 [passed] public=passed dense<=10delta
STEP shift=1e-04 delta=1e-06 apply_err=4.74e-08 sum_err=1.47e-09
ROW shift=1e-04 delta=1e-06 ||x||/||b||=8.94e+03 res_dense_v3=5.28e-07 check_v3=5.30e-07 [passed] public=passed dense<=10delta
STEP shift=1e-06 delta=1e-04 apply_err=2.58e-03 sum_err=9.75e-13
ROW shift=1e-06 delta=1e-04 ||x||/||b||=8.94e+05 res_dense_v3=1.19e-04 check_v3=2.58e-03 [warned] public=warned dense<=10delta
STEP shift=1e-06 delta=1e-05 apply_err=2.56e-03 sum_err=4.12e-13
ROW shift=1e-06 delta=1e-05 ||x||/||b||=8.94e+05 res_dense_v3=1.10e-04 check_v3=2.58e-03 [warned] public=warned dense>10delta
STEP shift=1e-06 delta=1e-06 apply_err=1.42e-03 sum_err=8.59e-13
ROW shift=1e-06 delta=1e-06 ||x||/||b||=8.94e+05 res_dense_v3=1.11e-04 check_v3=1.44e-03 [warned] public=warned dense>10delta
STEP shift=1e-07 delta=1e-04 apply_err=1.51e-01 sum_err=6.97e-15
ROW shift=1e-07 delta=1e-04 ||x||/||b||=8.95e+06 res_dense_v3=4.19e-03 check_v3=1.51e-01 [warned] public=warned dense>10delta
STEP shift=1e-07 delta=1e-05 apply_err=7.78e-02 sum_err=9.59e-15
ROW shift=1e-07 delta=1e-05 ||x||/||b||=8.95e+06 res_dense_v3=4.19e-03 check_v3=7.78e-02 [warned] public=warned dense>10delta
STEP shift=1e-07 delta=1e-06 apply_err=1.45e-01 sum_err=2.35e-15
ROW shift=1e-07 delta=1e-06 ||x||/||b||=8.95e+06 res_dense_v3=4.19e-03 check_v3=1.45e-01 [warned] public=warned dense>10delta
done
[exit 0, 8 s]
```

**Reviewer (CONFIRMED)**: Confirmed by the counter-repro's rule, which restores the sub-claim of a
false warning on the public API. The reviewer graded it LOW while that
sub-claim was struck; the lead raised it to MEDIUM once the sub-claim came
back, the grade of the parent finding by the parent's own rule (a real solve
gets a different verdict from the default check than from the exact one),
on the default backend. Struck: the size law
`eps (||x||/||b||)^2`, since the measured floor is 20 to 60 times above it,
closer to `eps (||A|| ||x||/||b||)^2` with `||A||` about 5. The rows where
the dense residual is above 10 delta show v3's own solve stopping above
delta, which bears on the open divergence lead, not on this finding.

**Suggested fix**: Give the v3 session an exact residual: contract the MPO into the MPS site by
site with no truncation (bond dimension `D_K D_x`, as v2's `exactApplyMPO`
does), take the direct sum with -b, make one orthogonalizing sweep at cutoff
0 and read the norm off the orthogonality centre. ITensor v3's `applyMPO`
accepts only `DensityMatrix` and `Fit`, so it cannot be a Method flag. Do it
together with the `"python"` fix, so "summed exactly" holds on every
backend. Whether v3's `bicstab` recurrence should take the exact apply too
is a separate pick.

### 18. The per-bond MPO cutoff does nothing when every bond is crossed by an O(1) term, so a weak link of 1e-7 in a ring, or a 1e-7 XY part on an Ising chain, reads 0.000 of itself on v2, v3 and `"python"`

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `pyitensor, v3-cpp` &middot; candidate `pyitensor-3` and `v3-cpp-3`, one defect

**Where**: `src/dmrgpy/mpscpp3/mo_terms.h:371-400 (bond_cutoff_factor: m = min over bonds of the strongest crossing term)`; `src/dmrgpy/mpscpp2/mo_terms.h:88-117 (same rule)`; `src/dmrgpy/pyitensor/mpobuilder.py:295-323 and 366-372 (_bond_cutoff_factor, lowering skipped when f = 1)`; `src/dmrgpy/pyitensor/chain.py:46 (_BUILD_CUTOFF = 1e-14, relative)`; `src/dmrgpy/ysr.py:900-931 (_metal_layout, the comb layout)`

The per-bond cutoff of `5a53e40` lowers the cutoff of the MPO build by
`m^2`, with m the minimum over bonds of the strongest coefficient crossing
each bond. So it does nothing when every bond of the operator is crossed by
an O(1) term. A weak multi-site term whose channel on such a bond is its own
is then cut exactly as before the fix: at an absolute 1e-13 on squared
singular values on v2 and v3, and at a relative 1e-14 on `"python"`. On 6
S=1/2 sites at s = 1e-7, three cases read 0.000000 of the weak part on v2,
v3 and `"python"`, against 1.000000 on ED: a weak link J' = s in a periodic
ring, the same open-chain link next to one strong `Sz0*Sz5`, and an Ising
chain plus `s*(XX+YY)`. The open-chain weak link the fix was built for reads
1.000000, and every case reads 1 at s = 1e-5. A geometric chain jumped by a
strong `S_0.S_11` loses its hoppings of 1e-7 and 1e-8 the same way, on
`"python"` and v3. The counter-repro reaches it through a public builder:
the `comb` layout of the YSR chain on a metal, a spin after the Wilson
chain, loses the hoppings of 1e-8 and below on both backends, which needs
Wilson chains of about 29 sites or more at Lambda = 4.

**Expected**: every term above the build's relative cutoff kept, as for the open chain.

Repro (v2, v3, python: ring, strong end-to-end term, Ising plus XY):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/v3-cpp/v3cpp_1_bond_cutoff_bypass.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# The per-bond MPO cutoff (mo_terms.h bond_cutoff_factor, 5a53e40) lowers
# svdMPO's absolute Cutoff by m^2, m = min over bonds of the STRONGEST
# coefficient crossing that bond, over the largest coefficient. A strong
# term that only passes through a weak bond (a periodic closing bond, an
# end-to-end coupling) sets m = 1 there, so the weak bond's own channels
# meet the old absolute 1e-13 again. Same for a weak channel on a bond that
# also holds a strong one (Ising ZZ plus eps*XY).
# Probe: ratio = [vev(X) - vev(X without the small part)] / s / vev(small
# part), on one fixed state. Exact answer: 1.0 for every s (linear in s).
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext

N = 6


def dot(sc, i, j):
    return sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]


def xy(sc, i, j):
    return sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j]


def chain(version):
    np.random.seed(1)
    v = "python" if version == "ED" else version
    sc = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=v)
    sc.maxm = 30
    sc.nsweeps = 10
    return sc


def heis(sc, bonds):
    """sum of J*S_i.S_j over bonds = [(i, j, J)]"""
    h = 0
    for (i, j, J) in bonds:
        h = h + J*dot(sc, i, j)
    return h


OBC = [(i, i+1) for i in range(N-1)]
WEAK = (2, 3)


def cases(sc, s):
    """name -> (X, big, small): X = big + s*small"""
    out = {}
    strong_obc = [(i, j, 1.0) for (i, j) in OBC if (i, j) != WEAK]
    # control: the open-chain weak link the per-bond fix was built for
    big = heis(sc, strong_obc)
    out["obc_weak_link(control)"] = (big + s*dot(sc, *WEAK), big, dot(sc, *WEAK))
    # periodic ring with one weak link: the closing bond (0,5) passes
    # through the weak bond (2,3)
    big = heis(sc, strong_obc + [(0, N-1, 1.0)])
    out["pbc_weak_link"] = (big + s*dot(sc, *WEAK), big, dot(sc, *WEAK))
    # open chain with the weak link and one strong end-to-end Sz0*Sz5
    big = heis(sc, strong_obc) + 1.0*sc.Sz[0]*sc.Sz[N-1]
    out["obc_weak_link+Sz0Sz5"] = (big + s*dot(sc, *WEAK), big, dot(sc, *WEAK))
    # every bond holds ZZ at 1 and XY at s (weak channel on a strong bond)
    big = 0
    small = 0
    for (i, j) in OBC:
        big = big + sc.Sz[i]*sc.Sz[j]
        small = small + xy(sc, i, j)
    out["ising+s*XY"] = (big + s*small, big, small)
    return out


backends = ["ED", "python"]
for v in (3, 2):
    if cppext.available(v):
        backends.append(v)
    else:
        print("backend v%d = not compiled, skipped" % v)

print("# ratio exact = 1.0; a dropped small part reads 0.0 (or 1/3, 2/3)")
for s in (1e-5, 1e-7):
    for b in backends:
        sc = chain(b)
        mode = "ED" if b == "ED" else "DMRG"
        # fixed state: ground state of the uniform periodic ring
        ring = heis(sc, [(i, j, 1.0) for (i, j) in OBC] + [(0, N-1, 1.0)])
        sc.set_hamiltonian(ring)
        if mode == "DMRG":
            wf = sc.get_gs()
        for name, (X, big, small) in cases(sc, s).items():
            if mode == "ED":
                vx = np.real(sc.vev(X, mode="ED"))
                vb = np.real(sc.vev(big, mode="ED"))
                vs = np.real(sc.vev(small, mode="ED"))
            else:
                vx = np.real(sc.vev(X, wf=wf))
                vb = np.real(sc.vev(big, wf=wf))
                vs = np.real(sc.vev(small, wf=wf))
            r = (vx - vb)/s/vs
            print("ratio s=%.0e backend=%-6s case=%-24s = %.6f" % (s, b, name, r))

# ground-state energy of the ring with one weak link J'=s, against ED
s = 1e-7
print("# gs: (E_dmrg - E_ed)/s should be ~0 (+-1e-3); a dropped weak link "
      "reads -<B23>_0 (printed as ref_if_dropped)")
for b in backends:
    if b == "ED":
        continue
    sc = chain(b)
    strong = [(i, j, 1.0) for (i, j) in OBC if (i, j) != WEAK] + [(0, N-1, 1.0)]
    big = heis(sc, strong)
    sc.set_hamiltonian(big)
    e0ed = float(np.real(sc.gs_energy(mode="ED")))
    b23 = float(np.real(sc.vev(dot(sc, *WEAK), mode="ED")))
    X = big + s*dot(sc, *WEAK)
    sc.set_hamiltonian(X)
    e = float(np.real(sc.gs_energy()))
    eed = float(np.real(sc.gs_energy(mode="ED")))
    print("gs backend=%-6s (E_dmrg-E_ed)/s = %.6f  (E_ed-E0_ed)/s = %.6f  "
          "ref_if_dropped = %.6f" % (b, (e - eed)/s, (eed - e0ed)/s, -b23))
```

Observed:

```
# ratio exact = 1.0; a dropped small part reads 0.0 (or 1/3, 2/3)
ratio s=1e-05 backend=ED     case=obc_weak_link(control)   = 1.000000
ratio s=1e-05 backend=ED     case=pbc_weak_link            = 1.000000
ratio s=1e-05 backend=ED     case=obc_weak_link+Sz0Sz5     = 1.000000
ratio s=1e-05 backend=ED     case=ising+s*XY               = 1.000000
ratio s=1e-05 backend=python case=obc_weak_link(control)   = 1.000000
ratio s=1e-05 backend=python case=pbc_weak_link            = 1.000000
ratio s=1e-05 backend=python case=obc_weak_link+Sz0Sz5     = 1.000000
ratio s=1e-05 backend=python case=ising+s*XY               = 1.000000
ratio s=1e-05 backend=3      case=obc_weak_link(control)   = 1.000000
ratio s=1e-05 backend=3      case=pbc_weak_link            = 1.000000
ratio s=1e-05 backend=3      case=obc_weak_link+Sz0Sz5     = 1.000000
ratio s=1e-05 backend=3      case=ising+s*XY               = 1.000000
ratio s=1e-05 backend=2      case=obc_weak_link(control)   = 1.000000
ratio s=1e-05 backend=2      case=pbc_weak_link            = 1.000000
ratio s=1e-05 backend=2      case=obc_weak_link+Sz0Sz5     = 1.000000
ratio s=1e-05 backend=2      case=ising+s*XY               = 1.000000
ratio s=1e-07 backend=ED     case=obc_weak_link(control)   = 1.000000
ratio s=1e-07 backend=ED     case=pbc_weak_link            = 1.000000
ratio s=1e-07 backend=ED     case=obc_weak_link+Sz0Sz5     = 1.000000
ratio s=1e-07 backend=ED     case=ising+s*XY               = 1.000000
ratio s=1e-07 backend=python case=obc_weak_link(control)   = 1.000000
ratio s=1e-07 backend=python case=pbc_weak_link            = 0.000000
ratio s=1e-07 backend=python case=obc_weak_link+Sz0Sz5     = -0.000000
ratio s=1e-07 backend=python case=ising+s*XY               = -0.000000
ratio s=1e-07 backend=3      case=obc_weak_link(control)   = 1.000000
ratio s=1e-07 backend=3      case=pbc_weak_link            = -0.000000
ratio s=1e-07 backend=3      case=obc_weak_link+Sz0Sz5     = -0.000000
ratio s=1e-07 backend=3      case=ising+s*XY               = -0.000000
ratio s=1e-07 backend=2      case=obc_weak_link(control)   = 1.000000
ratio s=1e-07 backend=2      case=pbc_weak_link            = -0.000000
ratio s=1e-07 backend=2      case=obc_weak_link+Sz0Sz5     = -0.000000
ratio s=1e-07 backend=2      case=ising+s*XY               = -0.000000
# gs: (E_dmrg - E_ed)/s should be ~0 (+-1e-3); a dropped weak link reads -<B23>_0 (printed as ref_if_dropped)
gs backend=python (E_dmrg-E_ed)/s = 0.140460  (E_ed-E0_ed)/s = -0.140460  ref_if_dropped = 0.140460
gs backend=3      (E_dmrg-E_ed)/s = 0.140460  (E_ed-E0_ed)/s = -0.140460  ref_if_dropped = 0.140460
gs backend=2      (E_dmrg-E_ed)/s = 0.140460  (E_ed-E0_ed)/s = -0.140460  ref_if_dropped = 0.140460
[exit 0, 4 s]
```

Repro (python, v3: a geometric chain jumped by a strong term):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/pyitensor/pyit_3_mpo_cutoff_jump.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: the per-bond MPO cutoff of 5a53e40 (pyitensor
# mpobuilder._bond_cutoff_factor, mo_terms.h bond_cutoff_args) lowers the
# build cutoff by m^2 with m the STRONGEST coefficient crossing the weakest
# bond. A weak term on a bond that a strong term also crosses (a coupling
# that jumps over a geometrically decaying chain, as a Kondo or spin bond
# jumps a Wilson chain in the multi-lead layouts) keeps m = 1, so the weak
# channels are cut at the old relative 1e-14 ("python") / absolute 1e-13
# (v3) and the deep hoppings vanish from the MPO. Probe: a chain 1..10 with
# t_k = 10^-(k-1) on bond (k, k+1), O(1) bonds (0,1) and (10,11), and a
# jump J S_0.S_11. ratio_k = [vev(O) - vev(O - t_k B_k)] / (t_k vev(B_k))
# must be 1; with J=0 the per-bond factor sees the weak bonds.
from dmrgpy import spinchain, cppext
from dmrgpy.pyitensor.mpscontainer import _link_at

L = 12


def build(version):
    np.random.seed(4)
    sc = spinchain.Spin_Chain(["S=1/2"] * L)
    if version == "python":
        sc.setup_python()
    else:
        sc.setup_cpp(version=3)
    hu = 0
    for i in range(L - 1):
        hu = hu + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(hu)       # the state: a uniform Heisenberg ground state
    sc.maxm = 64
    sc.nsweeps = 8
    sc.gs_energy()
    return sc


def B(sc, i, j):
    return sc.Sx[i] * sc.Sx[j] + sc.Sy[i] * sc.Sy[j] + sc.Sz[i] * sc.Sz[j]


def t(k):
    return 10.0 ** (-(k - 1))


def operator(sc, J, skip=None):
    op = B(sc, 0, 1) + B(sc, L - 2, L - 1)
    for k in range(1, L - 2):          # bonds (1,2) .. (9,10)
        if k != skip:
            op = op + t(k) * B(sc, k, k + 1)
    if J != 0:
        op = op + J * B(sc, 0, L - 1)
    return op



def print_dims(sc):
    s = sc._session
    for J in (0.0, 1.0):
        M = s._mpo(operator(sc, J).to_terms())
        dims = [_link_at(M, i, i + 1).dim for i in range(1, L)]
        print("python MPO link dims, J = %.1f: %s" % (J, dims), flush=True)


deep = [5, 6, 7, 8, 9]                 # t = 1e-4 .. 1e-8
print("deep bonds k =", deep, " t_k =", ["%.0e" % t(k) for k in deep], flush=True)
runs = [("python", "DMRG"), ("python", "ED")]
if cppext.available(3):
    runs.append((3, "DMRG"))
chains = {}
for version, mode in runs:
    if version not in chains:
        chains[version] = build(version)
        if version == "python":
            print_dims(chains[version])  # before any v3 call
    sc = chains[version]
    for J in (0.0, 1.0):
        full = sc.vev(operator(sc, J), mode=mode).real
        ratios = []
        for k in deep:
            red = sc.vev(operator(sc, J, skip=k), mode=mode).real
            bk = sc.vev(B(sc, k, k + 1), mode=mode).real
            ratios.append((full - red) / (t(k) * bk))
        print("backend = %-6s mode = %-4s J = %.1f  ratio_k = %s"
              % (version, mode, J, " ".join("%.4f" % r for r in ratios)), flush=True)
```

Observed:

```
deep bonds k = [5, 6, 7, 8, 9]  t_k = ['1e-04', '1e-05', '1e-06', '1e-07', '1e-08']
python MPO link dims, J = 0.0: [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4]
python MPO link dims, J = 1.0: [4, 8, 8, 8, 8, 8, 8, 8, 5, 5, 4]
backend = python mode = DMRG J = 0.0  ratio_k = 1.0000 1.0000 1.0000 1.0000 1.0000
backend = python mode = DMRG J = 1.0  ratio_k = 1.0000 1.0000 1.0000 -0.0000 -0.0000
backend = python mode = ED   J = 0.0  ratio_k = 1.0000 1.0000 1.0000 1.0000 1.0000
backend = python mode = ED   J = 1.0  ratio_k = 1.0000 1.0000 1.0000 1.0000 1.0000
backend = 3      mode = DMRG J = 0.0  ratio_k = 1.0000 1.0000 1.0000 1.0000 1.0000
backend = 3      mode = DMRG J = 1.0  ratio_k = 1.0000 1.0000 1.0000 -0.0000 -0.0000
[exit 0, 4 s]
```

Reviewer's counter-repro, the YSR comb and sides layouts:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g09/g09_native_jw_span.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis (by reading, 86ec994): on Spinful_Fermionic_Chain_Native every
# fermionic factor carries its own Jordan-Wigner string F_0..F_{i-1} in the
# term list (jordanwigner_spinful.one_fermion; Sx and Sy are Cdag*C there),
# and AutoMPO.from_terms / HTerm::add keep the paired F factors, so
# mpobuilder._bond_cutoff_factor and mo_terms.h bond_cutoff_factor see every
# hopping and every spin-spin term spanning from site 0. A Wilson chain whose
# deep end is not the right end of the MPS then has every deep bond "crossed"
# by an O(1) term, m stays about 1/32 and its weak hoppings are cut as before
# 5a53e40:
#   order="sides", contact in the left half: the chain is reversed on the
#     left end, crossed only through the strings (the JW-span defect);
#   order="comb", a spin after the chain: Js genuinely jumps the chain
#     (pyitensor-3's jump, through a public builder);
#   order="linear": deep end = right end of the MPS, the control.
# Probe: ratio_m = [vev(H) - vev(H - t_m hop_m)] / (t_m vev(hop_m)) on one
# fixed state (the ground state with every hopping at 1); exact answer 1.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import ysr, cppext
from dmrgpy.pyitensor import task, mpobuilder

NS = 2
T = [0.5, 1e-4, 1e-8, 1e-9, 1e-10]   # Wilson hoppings t_m, bond m = (m, m+1)
N = len(T) + 1                       # 6 Wilson sites, 8 native sites in all
EPS = [0.0]*N


def chain(order, version):
    np.random.seed(11)
    fc = ysr.spin_chain_on_metal(NS, N, site=[0], order=order,
                                 wilson=(EPS, T), itensor_version=version)
    fc.maxm = 32
    fc.nsweeps = 6
    return fc


def link_dims(W, w):
    # the link between MPS positions p and p+1 is the right link of site p
    return [int(np.shape(W[min(w[m], w[m + 1])])[-1]) for m in range(len(w) - 1)]


def ratios(fc, H, hops, t, **kw):
    vH = np.real(fc.vev(H, **kw))
    out = []
    for m in range(len(t)):
        vr = np.real(fc.vev(H - t[m]*hops[m], **kw))
        vb = np.real(fc.vev(hops[m], **kw))
        out.append((vH - vr)/(t[m]*vb))
    return " ".join("%.4f" % r for r in out)


print("t_m =", ["%.0e" % x for x in T], flush=True)
backends = ["python"] + ([3] if cppext.available(3) else [])
for version in backends:
    for order in ("linear", "sides", "comb"):
        fc = chain(order, version)
        w = list(fc.bath_sites[0])
        t = np.asarray(fc.lead_chains[0][1])
        H = fc.hamiltonian
        terms = H.to_terms()
        s = fc._session
        if version == "python":
            f = mpobuilder._bond_cutoff_factor(s._ampo(terms))
            dims = link_dims(task.mpo_arrays(s._mpo(terms)), w)
            print("python order=%-6s wilson sites=%s bond_cutoff_factor=%.3e "
                  "link dims per Wilson bond=%s" % (order, w, f, dims), flush=True)
        else:
            dims = link_dims(s.mpo_site_arrays(s.build_operator(terms)), w)
            print("v3     order=%-6s wilson sites=%s link dims per Wilson bond=%s"
                  % (order, w, dims), flush=True)
        hops = [ysr._hop(fc, w[m], w[m + 1]) for m in range(N - 1)]
        Hst = H
        for m in range(N - 1):   # the state: every Wilson hopping at 1
            Hst = Hst + (1.0 - t[m])*hops[m]
        fc.set_hamiltonian(Hst)
        wf = fc.get_gs()
        print("ratio backend=%-6s mode=DMRG order=%-6s ratio_m=%s"
              % (version, order, ratios(fc, H, hops, t, wf=wf)), flush=True)
        if version == "python":
            print("ratio backend=ED     mode=ED   order=%-6s ratio_m=%s"
                  % (order, ratios(fc, H, hops, t, mode="ED")), flush=True)
```

Observed:

```
t_m = ['5e-01', '1e-04', '1e-08', '1e-09', '1e-10']
python order=linear wilson sites=[2, 3, 4, 5, 6, 7] bond_cutoff_factor=1.563e-22 link dims per Wilson bond=[6, 6, 6, 6, 5]
ratio backend=python mode=DMRG order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
ratio backend=ED     mode=ED   order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
python order=sides  wilson sites=[5, 4, 3, 2, 1, 0] bond_cutoff_factor=3.906e-03 link dims per Wilson bond=[6, 5, 1, 1, 1]
ratio backend=python mode=DMRG order=sides  ratio_m=1.0000 1.0000 -0.0000 0.0000 -0.0000
ratio backend=ED     mode=ED   order=sides  ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
python order=comb   wilson sites=[1, 2, 3, 4, 5, 6] bond_cutoff_factor=1.000e+00 link dims per Wilson bond=[9, 9, 5, 5, 5]
ratio backend=python mode=DMRG order=comb   ratio_m=1.0000 1.0000 -0.0000 -0.0000 -0.0000
ratio backend=ED     mode=ED   order=comb   ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
v3     order=linear wilson sites=[2, 3, 4, 5, 6, 7] link dims per Wilson bond=[7, 7, 7, 7, 6]
ratio backend=3      mode=DMRG order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
v3     order=sides  wilson sites=[5, 4, 3, 2, 1, 0] link dims per Wilson bond=[8, 9, 6, 5, 4]
ratio backend=3      mode=DMRG order=sides  ratio_m=1.0000 1.0000 0.2500 -0.0000 -0.0000
v3     order=comb   wilson sites=[1, 2, 3, 4, 5, 6] link dims per Wilson bond=[11, 10, 6, 6, 6]
ratio backend=3      mode=DMRG order=comb   ratio_m=1.0000 1.0000 0.0000 -0.0000 -0.0000
[exit 0, 65 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Both reviewers confirmed it. Struck: "CLAUDE.md says a hierarchy inside H is
covered", since CLAUDE.md states the rule exactly; and the headline "a weak
bond that a strong term merely passes through", which needs the qualifier
that m is a global minimum, so a single other bond crossed by weak terms
alone protects the whole operator. The YSR reach, first struck as unmeasured,
is shown by the counter-repro for the `comb` layout. The `sides` layout of
the same counter-repro is a separate mechanism, the reviewer-found finding
on Jordan-Wigner spans below. LOW, since it needs a weak term below about
3e-7 of the largest one with every bond crossed by a strong term.

**Suggested fix**: Taking m as the minimum over all crossing terms is the blanket lowering that
`.claude/rules/cpp-backend.md` forbids for long-range MPOs: on the 1/r^2
chain it takes the v3 cutoff to 1e-21 at 100 sites, and that MPO would no
longer be the same to the bit. The reviewer's rule keeps both constraints:
lower the cutoff by the square of the weakest merged multi-site coefficient
whose ratio to the largest lies between 1e-14 and a risk threshold of about
1e-5, whatever else crosses its bonds. It leaves the six reference chains
bit-identical and 1/r^2 compressible up to about 300 sites, and covers all
the cases above. The trade-off on long tails is the author's to decide.

### 19. On native spinful sites the per-bond cutoff reads every fermionic term as spanning from site 0, so a Wilson chain reversed on the left end loses its hoppings of 1e-8 and below, silently

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `pyitensor, v3-cpp` &middot; candidate found by the reviewer of `pyitensor-3`

**Where**: `src/dmrgpy/pyitensor/mpobuilder.py:295-323 (_bond_cutoff_factor reads a term's span from its first listed site)`; `src/dmrgpy/mpscpp3/mo_terms.h:371-400 (bond_cutoff_factor, the same)`; `src/dmrgpy/multioperatortk/jordanwigner_spinful.py (one_fermion dresses every factor with F from site 0)`; `src/dmrgpy/ysr.py:900-931 and 973-974 (_metal_layout, the sides order)`

On `Spinful_Fermionic_Chain_Native` every fermionic factor reaches the MPO
builder with its own Jordan-Wigner string `F_0 ... F_{i-1}`, and the paired
`F` factors below the leftmost fermion are not cancelled, although
`F_k F_k = 1`. Both `mpobuilder._bond_cutoff_factor` and `mo_terms.h`'s
`bond_cutoff_factor` read a term's span from its first to its last listed
site, so every fermionic term, the spin terms included, counts as crossing
every bond to its left. For a Wilson chain reversed on the left end of the
MPS, the YSR `sides` order with no contact in the right half, the factor is
then set by the spin bonds: 3.9e-3 on `"python"` against 1.6e-22 for the
`linear` order. The Wilson hoppings of 1e-8, 1e-9 and 1e-10 then vanish from
the MPO: ratio 0, 0, 0 on `"python"` and 0.25, 0, 0 on v3, against 1 on ED
and for `linear`, and the deep Wilson link dimensions on `"python"` fall to
1 against 6. Reaching it needs hoppings below about 1e-8, so Wilson chains of
about 27 sites or more at Lambda = 4.

**Expected**: the span read from the leftmost factor that is not a string, so the factor follows the deepest hopping as for `linear`.

Repro (the same script as the counter-repro of the per-bond cutoff finding above; its sides rows):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g09/g09_native_jw_span.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis (by reading, 86ec994): on Spinful_Fermionic_Chain_Native every
# fermionic factor carries its own Jordan-Wigner string F_0..F_{i-1} in the
# term list (jordanwigner_spinful.one_fermion; Sx and Sy are Cdag*C there),
# and AutoMPO.from_terms / HTerm::add keep the paired F factors, so
# mpobuilder._bond_cutoff_factor and mo_terms.h bond_cutoff_factor see every
# hopping and every spin-spin term spanning from site 0. A Wilson chain whose
# deep end is not the right end of the MPS then has every deep bond "crossed"
# by an O(1) term, m stays about 1/32 and its weak hoppings are cut as before
# 5a53e40:
#   order="sides", contact in the left half: the chain is reversed on the
#     left end, crossed only through the strings (the JW-span defect);
#   order="comb", a spin after the chain: Js genuinely jumps the chain
#     (pyitensor-3's jump, through a public builder);
#   order="linear": deep end = right end of the MPS, the control.
# Probe: ratio_m = [vev(H) - vev(H - t_m hop_m)] / (t_m vev(hop_m)) on one
# fixed state (the ground state with every hopping at 1); exact answer 1.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import ysr, cppext
from dmrgpy.pyitensor import task, mpobuilder

NS = 2
T = [0.5, 1e-4, 1e-8, 1e-9, 1e-10]   # Wilson hoppings t_m, bond m = (m, m+1)
N = len(T) + 1                       # 6 Wilson sites, 8 native sites in all
EPS = [0.0]*N


def chain(order, version):
    np.random.seed(11)
    fc = ysr.spin_chain_on_metal(NS, N, site=[0], order=order,
                                 wilson=(EPS, T), itensor_version=version)
    fc.maxm = 32
    fc.nsweeps = 6
    return fc


def link_dims(W, w):
    # the link between MPS positions p and p+1 is the right link of site p
    return [int(np.shape(W[min(w[m], w[m + 1])])[-1]) for m in range(len(w) - 1)]


def ratios(fc, H, hops, t, **kw):
    vH = np.real(fc.vev(H, **kw))
    out = []
    for m in range(len(t)):
        vr = np.real(fc.vev(H - t[m]*hops[m], **kw))
        vb = np.real(fc.vev(hops[m], **kw))
        out.append((vH - vr)/(t[m]*vb))
    return " ".join("%.4f" % r for r in out)


print("t_m =", ["%.0e" % x for x in T], flush=True)
backends = ["python"] + ([3] if cppext.available(3) else [])
for version in backends:
    for order in ("linear", "sides", "comb"):
        fc = chain(order, version)
        w = list(fc.bath_sites[0])
        t = np.asarray(fc.lead_chains[0][1])
        H = fc.hamiltonian
        terms = H.to_terms()
        s = fc._session
        if version == "python":
            f = mpobuilder._bond_cutoff_factor(s._ampo(terms))
            dims = link_dims(task.mpo_arrays(s._mpo(terms)), w)
            print("python order=%-6s wilson sites=%s bond_cutoff_factor=%.3e "
                  "link dims per Wilson bond=%s" % (order, w, f, dims), flush=True)
        else:
            dims = link_dims(s.mpo_site_arrays(s.build_operator(terms)), w)
            print("v3     order=%-6s wilson sites=%s link dims per Wilson bond=%s"
                  % (order, w, dims), flush=True)
        hops = [ysr._hop(fc, w[m], w[m + 1]) for m in range(N - 1)]
        Hst = H
        for m in range(N - 1):   # the state: every Wilson hopping at 1
            Hst = Hst + (1.0 - t[m])*hops[m]
        fc.set_hamiltonian(Hst)
        wf = fc.get_gs()
        print("ratio backend=%-6s mode=DMRG order=%-6s ratio_m=%s"
              % (version, order, ratios(fc, H, hops, t, wf=wf)), flush=True)
        if version == "python":
            print("ratio backend=ED     mode=ED   order=%-6s ratio_m=%s"
                  % (order, ratios(fc, H, hops, t, mode="ED")), flush=True)
```

Observed:

```
t_m = ['5e-01', '1e-04', '1e-08', '1e-09', '1e-10']
python order=linear wilson sites=[2, 3, 4, 5, 6, 7] bond_cutoff_factor=1.563e-22 link dims per Wilson bond=[6, 6, 6, 6, 5]
ratio backend=python mode=DMRG order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
ratio backend=ED     mode=ED   order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
python order=sides  wilson sites=[5, 4, 3, 2, 1, 0] bond_cutoff_factor=3.906e-03 link dims per Wilson bond=[6, 5, 1, 1, 1]
ratio backend=python mode=DMRG order=sides  ratio_m=1.0000 1.0000 -0.0000 0.0000 -0.0000
ratio backend=ED     mode=ED   order=sides  ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
python order=comb   wilson sites=[1, 2, 3, 4, 5, 6] bond_cutoff_factor=1.000e+00 link dims per Wilson bond=[9, 9, 5, 5, 5]
ratio backend=python mode=DMRG order=comb   ratio_m=1.0000 1.0000 -0.0000 -0.0000 -0.0000
ratio backend=ED     mode=ED   order=comb   ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
v3     order=linear wilson sites=[2, 3, 4, 5, 6, 7] link dims per Wilson bond=[7, 7, 7, 7, 6]
ratio backend=3      mode=DMRG order=linear ratio_m=1.0000 1.0000 1.0000 1.0000 1.0000
v3     order=sides  wilson sites=[5, 4, 3, 2, 1, 0] link dims per Wilson bond=[8, 9, 6, 5, 4]
ratio backend=3      mode=DMRG order=sides  ratio_m=1.0000 1.0000 0.2500 -0.0000 -0.0000
v3     order=comb   wilson sites=[1, 2, 3, 4, 5, 6] link dims per Wilson bond=[11, 10, 6, 6, 6]
ratio backend=3      mode=DMRG order=comb   ratio_m=1.0000 1.0000 0.0000 -0.0000 -0.0000
[exit 0, 65 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: "hoppings below about 5e-9 lose their channels", since `"python"`
loses 1e-8 as well and v3 keeps a quarter of it, so the depth is about 27
sites and not 29; and "that layout is the default for several spins", which
holds by reading only, since the repro forced `order="sides"` with one
contact. MEDIUM, since it is silent and the YSR builders reach it.

**Suggested fix**: Cancel the paired strings rather than lowering the cutoff further: in
`product2jw`, or in a pass after it, drop `F_k F_k = 1` for every k below the
leftmost fermionic factor of a parity-even term; every factor contributes
the same `F_k` there, so they always pair. The span then runs from the
leftmost factor that is not a string, and the terms handed to both builders
get shorter. Pin `order="sides"` with one left contact and hoppings to
1e-10 on `"python"` and v3, next to the `linear` control.

### 20. On `"python"`, `Chain.sum_mps_weighted` takes the density-matrix sum, so a weighted sum of nearly parallel states 5e-6 to 5e-9 of their norms comes back off by 4.9e-4 to 0.97

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `pyitensor` &middot; candidate `pyitensor-1`

**Where**: `src/dmrgpy/pyitensor/chain.py:910-922 (sum_mps_weighted: sum_many with no method)`; `src/dmrgpy/pyitensor/mpsalgebra.py:396 (SUM_METHOD = "densitymatrix")`; `src/dmrgpy/pyitensor/mpsalgebra.py:401-460 (_sum_densitymatrix and its documented floor)`; `src/dmrgpy/nhlevels.py:163-166 (_combine, the consumer; reach not shown)`

`sum_mps_weighted` (`4ed4891`) serves the Ritz vectors of `nhlevels.py`,
whose coefficients cancel between nearly parallel states, as its docstring
says. It calls `sum_many` with no method, so since `c0c4d2d` it takes the
density-matrix sum, whose own docstring states a floor of about 1e-8 of the
input norms. Fourteen hours earlier `17017a9` had moved every other
cancelling caller to `method="exact"`. With two nearly parallel inputs whose
weighted sum is 5e-6, 5e-7, 5e-8 and 5e-9 of their norms (coefficients up to
1e8), the result is off by 4.9e-4, 4.7e-3, 0.19 and 0.97 against a dense
anchor, while `method="exact"` stays at the construction error of the inputs
(6e-10 to 6e-7). The reach into `get_nh_levels(vectors=True)` was measured
and not found: on the 10-site PT-symmetric XX chain the accepted pairs have
`sum|c|` of at most 3.8, and the density-matrix and exact routes give the
same vectors to 1e-6.

**Expected**: `sum_mps_weighted` resolving a sum that cancels, as the exact route does.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/pyitensor/pyit_1_dm_sum_cancellation.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Hypothesis: on itensor_version="python" the public MPS difference a - b,
# H*psi - e*psi and Chain.sum_mps_weighted take the density-matrix sum
# (mpsalgebra.SUM_METHOD = "densitymatrix", default since c0c4d2d), whose
# floor is ~1e-8 of the INPUT norms. A difference of nearly equal states
# (and a weighted sum with cancelling coefficients, which nhlevels and the
# ESR KRYLOV route form) then carries a relative error ~1e-8/kappa, kappa
# = ||result||/||inputs||, while the exact route is right to roundoff.
# Reference: dense vectors of the MPS (8 sites, 256 amplitudes), no code
# shared with either route.
import time
from dmrgpy import spinchain
from dmrgpy import mps as dmps
from dmrgpy.pyitensor import mpsalgebra as ma

np.random.seed(1)
n = 8
sc = spinchain.Spin_Chain(["S=1/2"] * n)
sc.setup_python()
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
sc.set_hamiltonian(h)
sc.maxm = 32      # 8 sites: the exact bond dimension is 16
sc.nsweeps = 12
sc.cutoff = 1e-12
e_gs = sc.gs_energy()
a = sc.get_gs()
a = a * (1. / np.sqrt(abs(a.dot(a))))
s = sc._session
s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
sites = s.sites
print("SUM_METHOD default =", ma.SUM_METHOD, "APPLY_METHOD default =", ma.APPLY_METHOD, flush=True)
print("e_gs =", e_gs, flush=True)


def vec(h_):
    T = h_.A(1)
    for i in range(2, h_.length() + 1):
        T = T * h_.A(i)
    arr = T.transpose_to([sites.si(i) for i in range(1, h_.length() + 1)])
    return np.asarray(arr).reshape(-1)


def mpo_dense(mpo):
    T = mpo.A(1)
    for i in range(2, mpo.length() + 1):
        T = T * mpo.A(i)
    order = ([sites.si(i).prime(1) for i in range(1, n + 1)]
             + [sites.si(i) for i in range(1, n + 1)])
    return np.asarray(T.transpose_to(order)).reshape(2**n, 2**n)


def exact_sum(hs, cs):
    terms = [x * complex(c) for x, c in zip(hs, cs)]
    return ma.sum_many(terms, cutoff=0.0, maxdim=10**6, method="exact")


def with_routes(method, fn):
    old = (ma.SUM_METHOD, ma.APPLY_METHOD)
    ma.SUM_METHOD, ma.APPLY_METHOD = method, method
    try:
        return fn()
    finally:
        ma.SUM_METHOD, ma.APPLY_METHOD = old


# c = Sz_0|a>, a fixed MPS; b = a + eps c is built by the exact sum, so
# a - b = -eps c exactly and sum_mps_weighted([b, a], [1/eps, -1/eps]) = c.
c = sc.Sz[0] * a
va, vc = vec(a.cpp_handle), vec(c.cpp_handle)
nvc = np.linalg.norm(vc)
print("||c|| =", nvc, flush=True)
for eps in (1e-2, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
    b_h = exact_sum([a.cpp_handle, c.cpp_handle], [1.0, eps])
    b = dmps.MPS(MBO=sc, cpp_handle=b_h)
    vb = vec(b_h)
    berr = np.linalg.norm(vb - va - eps * vc) / (eps * nvc)
    target = -eps * vc
    out = {}
    for route in ("densitymatrix", "exact"):
        t0 = time.time()
        d = with_routes(route, lambda: a - b)
        vd = vec(d.cpp_handle)
        rel_d = np.linalg.norm(vd - target) / (eps * nvc)
        norm_ratio = np.linalg.norm(vd) / (eps * nvc)
        y_h = with_routes(route, lambda: s.sum_mps_weighted(
            [b.cpp_handle, a.cpp_handle], [1.0 / eps, -1.0 / eps]))
        rel_y = np.linalg.norm(vec(y_h) - vc) / nvc
        out[route] = (rel_d, norm_ratio, rel_y, time.time() - t0)
    print("eps = %.0e  b_build_err = %.2e" % (eps, berr), flush=True)
    for route in ("densitymatrix", "exact"):
        rel_d, norm_ratio, rel_y, dt = out[route]
        print("  route = %-13s  a-b rel_err = %.3e  ||a-b||/true = %.6f  "
              "sum_mps_weighted rel_err = %.3e" % (route, rel_d, norm_ratio, rel_y), flush=True)

# H*psi - e*psi for the (near-)eigenstate a: its squared norm is the
# variance; dense reference from the MPO of h contracted to a matrix.
Hd = mpo_dense(s._mpo(h.to_terms()))
e = np.vdot(va, Hd @ va).real
var_dense = np.linalg.norm(Hd @ va - e * va)**2
print("variance dense reference =", var_dense, flush=True)
for route in ("densitymatrix", "exact"):
    r = with_routes(route, lambda: h * a - e * a)
    print("  route = %-13s  ||H a - e a||^2 public = %.3e  (dense %.3e)"
          % (route, abs(r.dot(r)), var_dense), flush=True)
```

Observed:

```
SUM_METHOD default = densitymatrix APPLY_METHOD default = densitymatrix
e_gs = -3.37493259868789
||c|| = 0.4999999999999989
eps = 1e-02  b_build_err = 5.48e-13
  route = densitymatrix  a-b rel_err = 4.130e-11  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 9.299e-11
  route = exact          a-b rel_err = 5.976e-13  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 5.589e-13
eps = 1e-04  b_build_err = 5.86e-11
  route = densitymatrix  a-b rel_err = 2.511e-06  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 1.375e-06
  route = exact          a-b rel_err = 6.460e-11  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 6.266e-11
eps = 1e-05  b_build_err = 5.26e-10
  route = densitymatrix  a-b rel_err = 6.463e-04  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 4.886e-04
  route = exact          a-b rel_err = 5.368e-10  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 5.984e-10
eps = 1e-06  b_build_err = 4.48e-09
  route = densitymatrix  a-b rel_err = 1.678e-02  ||a-b||/true = 0.999859  sum_mps_weighted rel_err = 4.693e-03
  route = exact          a-b rel_err = 7.735e-09  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 4.393e-09
eps = 1e-07  b_build_err = 6.46e-08
  route = densitymatrix  a-b rel_err = 8.646e-02  ||a-b||/true = 0.996255  sum_mps_weighted rel_err = 1.850e-01
  route = exact          a-b rel_err = 6.940e-08  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 6.614e-08
eps = 1e-08  b_build_err = 5.25e-07
  route = densitymatrix  a-b rel_err = 9.216e-01  ||a-b||/true = 0.388030  sum_mps_weighted rel_err = 9.716e-01
  route = exact          a-b rel_err = 5.876e-07  ||a-b||/true = 1.000000  sum_mps_weighted rel_err = 5.791e-07
variance dense reference = 2.795079500163924e-24
  route = densitymatrix  ||H a - e a||^2 public = 6.574e-29  (dense 2.795e-24)
  route = exact          ||H a - e a||^2 public = 2.790e-24  (dense 2.795e-24)
[exit 0, 4 s]
```

Reviewer's counter-repro, the reach into get_nh_levels:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g07/g07_p1_nhlevels_reach.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Group 07 counter-repro P1 (the "reached" sub-claims of pyitensor-1 and
# v3-cpp-2). Is the sum defect reached by get_nh_levels(vectors=True) on
# the documented PT-symmetric XX chain, i 0.3 (-1)^i Sz_i, 10 sites? The
# script runs levels()'s own steps (nhlevels._right_states,
# _rayleigh_ritz, _accepted, _combine), and for every accepted Ritz pair
# prints the coefficient size, the deviation `dev` of the vector _combine
# returns from the exact combination sum_k c_k |s_k> (dense, from the MPS
# site arrays of the snapshots), and the dense residual of the exact
# combination (`res_combo`) and of the returned vector (`res_returned`),
# both relative as levels() reports them. "python" on both sum routes, v3
# on its own. v3's random start is not seeded by numpy, so the number of
# levels can differ from run to run.
# Run with MKL_NUM_THREADS=1 OMP_NUM_THREADS=1; about 1 to 2 minutes.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext
from dmrgpy import nhlevels as nl
from dmrgpy.pyitensor import mpsalgebra as ma
from dmrgpy.pyitensor import task


def dense_state(s, handle):
    arrs = [np.asarray(x) for x in s.mps_site_arrays(handle)]
    T = arrs[0][0]                      # (d, chi)
    for A in arrs[1:]:                  # (chi_l, d, chi_r)
        T = np.einsum("xa,asb->xsb", T, A).reshape(-1, A.shape[2])
    return T[:, 0]


def dense_mpo(arrs):
    M = np.asarray(arrs[0])[0]          # (d_out, d_in, w)
    for A in arrs[1:]:                  # (w_l, d_out, d_in, w_r)
        A = np.asarray(A)
        M = np.einsum("xyw,wabv->xaybv", M, A)
        M = M.reshape(M.shape[0]*M.shape[1], M.shape[2]*M.shape[3], M.shape[4])
    return M[:, :, 0]


def mpo_arrays(sc, M):
    if sc.itensor_version == "python":
        return task.mpo_arrays(M)
    return sc._session.mpo_site_arrays(M)


def routes(v):
    return ("densitymatrix", "exact") if v == "python" else ("native",)


def with_route(v, route, fn):
    if v != "python":
        return fn()
    old = ma.SUM_METHOD
    ma.SUM_METHOD = route
    try:
        return fn()
    finally:
        ma.SUM_METHOD = old


versions = ["python"]
if cppext.available(3):
    versions.insert(0, 3)
else:
    print("backend v3 = not compiled, skipped")

for v in versions:
    np.random.seed(11)
    n = 10
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm = 40       # 10 sites: the exact bond dimension is 32
    sc.nsweeps = 10
    sc.cutoff = 1e-12
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        h = h + 0.3j*((-1)**i)*sc.Sz[i]
    sc.set_hamiltonian(h)
    H = sc.hamiltonian
    snaps = nl._right_states(sc, H, None, nl.NSTEPS, nl.NSNAP)
    s = sc._session
    V = np.array([dense_state(s, x.cpp_handle) for x in snaps]).T
    Hd = dense_mpo(mpo_arrays(sc, s.build_operator(H.to_terms())))
    m = sc.get_ED_obj().MO2matrix(h)
    m = m.todense() if hasattr(m, "todense") else m
    e_mpo = np.min(np.linalg.eigvals(Hd).real)
    e_ed = np.min(np.linalg.eigvals(np.asarray(m)).real)
    print("P1: v=%s min Re E dense-MPO = %.10f  ED = %.10f  |diff| = %.2e"
          % (v, e_mpo, e_ed, abs(e_mpo - e_ed)), flush=True)
    w = np.linalg.eigvalsh(V.conj().T @ V)
    kept = w > nl.GRAM_FLOOR*w.max()
    print("P1: v=%s Gram eigenvalues kept = %d of %d, smallest kept/largest = %.2e"
          % (v, int(np.sum(kept)), len(w), np.min(w[kept])/w.max()), flush=True)
    acc = nl._accepted(nl._rayleigh_ritz(snaps, H), H, nl.RITZ_TOL)
    rscale, e_id = nl._scale(H)
    print("P1: v=%s accepted levels = %d" % (v, len(acc)), flush=True)
    for k, (t, r, c) in enumerate(acc):
        scl = rscale + abs(t - e_id)
        yt = V @ c
        nyt = np.linalg.norm(yt)
        res_combo = np.linalg.norm(Hd @ yt - t*yt)/nyt/scl
        print("P1: v=%s level %d E = %.8f%+.8fi  reported = %.3e  res_combo = %.3e  "
              "sum|c| = %.3e  max|c| = %.3e"
              % (v, k, t.real, t.imag, r/scl, res_combo,
                 np.sum(np.abs(c)), np.max(np.abs(c))), flush=True)
        for route in routes(v):
            y = with_route(v, route, lambda: nl._combine(sc, snaps, c))
            vy = dense_state(s, y.cpp_handle)
            dev = np.linalg.norm(vy - yt)/nyt
            res_ret = np.linalg.norm(Hd @ vy - t*vy)/np.linalg.norm(vy)/scl
            print("P1: v=%s level %d route=%-13s dev = %.3e  res_returned = %.3e"
                  % (v, k, route, dev, res_ret), flush=True)
```

Observed:

```
P1: v=3 min Re E dense-MPO = -2.5667821672  ED = -2.5667821672  |diff| = 4.44e-15
P1: v=3 Gram eigenvalues kept = 6 of 24, smallest kept/largest = 3.24e-09
P1: v=3 accepted levels = 4
P1: v=3 level 0 E = -2.56678219+0.26409558i  reported = 4.420e-07  res_combo = 5.882e-07  sum|c| = 1.669e+00  max|c| = 1.182e-01
P1: v=3 level 0 route=native        dev = 9.496e-07  res_returned = 1.035e-06
P1: v=3 level 1 E = -2.56678218-0.00000000i  reported = 5.982e-07  res_combo = 6.106e-07  sum|c| = 2.477e+00  max|c| = 1.941e-01
P1: v=3 level 1 route=native        dev = 1.166e-06  res_returned = 1.507e-06
P1: v=3 level 2 E = -2.56678213-0.26409567i  reported = 8.066e-07  res_combo = 9.421e-07  sum|c| = 2.711e+00  max|c| = 1.918e-01
P1: v=3 level 2 route=native        dev = 1.405e-06  res_returned = 1.354e-06
P1: v=3 level 3 E = -2.56678119+0.00000045i  reported = 3.850e-05  res_combo = 3.850e-05  sum|c| = 1.170e+02  max|c| = 1.308e+01
P1: v=3 level 3 route=native        dev = 1.350e-05  res_returned = 3.585e-05
P1: v=python min Re E dense-MPO = -2.5667821672  ED = -2.5667821672  |diff| = 1.51e-14
P1: v=python Gram eigenvalues kept = 3 of 24, smallest kept/largest = 8.52e-02
P1: v=python accepted levels = 3
P1: v=python level 0 E = -2.56678219-0.26409559i  reported = 3.070e-07  res_combo = 4.680e-07  sum|c| = 1.556e+00  max|c| = 8.925e-02
P1: v=python level 0 route=densitymatrix dev = 6.000e-07  res_returned = 7.221e-07
P1: v=python level 0 route=exact         dev = 6.000e-07  res_returned = 7.221e-07
P1: v=python level 1 E = -2.56678217-0.00000000i  reported = 1.559e-07  res_combo = 3.774e-07  sum|c| = 1.293e+00  max|c| = 7.385e-02
P1: v=python level 1 route=densitymatrix dev = 7.033e-07  res_returned = 1.039e-06
P1: v=python level 1 route=exact         dev = 7.033e-07  res_returned = 1.039e-06
P1: v=python level 2 E = -2.56678212+0.26409572i  reported = 7.954e-07  res_combo = 1.151e-06  sum|c| = 3.765e+00  max|c| = 2.159e-01
P1: v=python level 2 route=densitymatrix dev = 9.737e-07  res_returned = 1.131e-06
P1: v=python level 2 route=exact         dev = 9.737e-07  res_returned = 1.131e-06
[exit 0, 58 s]
```

Reviewer's counter-repro, the two-term floor on v3 and python:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g07/g07_p2_two_term_floor.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Group 07 counter-repro P2 (both candidates, fix opinion of v3-cpp-2 and
# the a - b sub-claim of pyitensor-1). Does a TWO-term sum on v3, ITensor's
# sum(a,b) = plusEq -> addAssumeOrth -> MPS::orthogonalize (no pairwise
# tree, Nt = 2), carry the same eps*||inputs||^2 floor as python's
# density-matrix route? Inputs a (8-site Heisenberg ground state, exact at
# maxm=32) and b = a + eps Sz_0 a, built by the backend's own public sum.
# The anchor is dense(a) - dense(b), the inputs as they are, read off the
# MPS site arrays, so it does not matter how b was built.
# Run with MKL_NUM_THREADS=1 OMP_NUM_THREADS=1; well under a minute.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext
from dmrgpy.pyitensor import mpsalgebra as ma


def dense_state(s, handle):
    arrs = [np.asarray(x) for x in s.mps_site_arrays(handle)]
    T = arrs[0][0]                      # (d, chi)
    for A in arrs[1:]:                  # (chi_l, d, chi_r)
        T = np.einsum("xa,asb->xsb", T, A).reshape(-1, A.shape[2])
    return T[:, 0]


def routes(v):
    return ("densitymatrix", "exact") if v == "python" else ("native",)


def with_route(v, route, fn):
    if v != "python":
        return fn()
    old = ma.SUM_METHOD
    ma.SUM_METHOD = route
    try:
        return fn()
    finally:
        ma.SUM_METHOD = old


versions = ["python"]
if cppext.available(3):
    versions.insert(0, 3)
else:
    print("backend v3 = not compiled, skipped")

for v in versions:
    np.random.seed(1)
    n = 8
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm = 32       # 8 sites: exact bond 16, a sum of two at most 32
    sc.nsweeps = 12
    sc.cutoff = 1e-12
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    print("P2: v=%s e_gs = %.12f" % (v, sc.gs_energy()), flush=True)
    a = sc.get_gs()
    a = a*(1./np.sqrt(abs(a.dot(a))))
    s = sc._session
    s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    c = sc.Sz[0]*a
    va = dense_state(s, a.cpp_handle)
    for eps in (1e-4, 1e-5, 1e-6, 1e-8):
        b = a + eps*c              # this backend's own public sum
        vb = dense_state(s, b.cpp_handle)
        tgt = va - vb
        nt = np.linalg.norm(tgt)
        print("P2: v=%s eps=%.0e ||a-b||_dense/eps = %.6f" % (v, eps, nt/eps), flush=True)
        for route in routes(v):
            d = with_route(v, route, lambda: a - b)
            rel_d = np.linalg.norm(dense_state(s, d.cpp_handle) - tgt)/nt
            yh = with_route(v, route, lambda: s.sum_mps_weighted(
                [b.cpp_handle, a.cpp_handle], [1.0/eps + 0j, -1.0/eps + 0j]))
            rel_y = np.linalg.norm(dense_state(s, yh) - (vb - va)/eps)/(nt/eps)
            print("P2: v=%s eps=%.0e route=%-13s a-b rel_err = %.3e  "
                  "sum_mps_weighted rel_err = %.3e" % (v, eps, route, rel_d, rel_y),
                  flush=True)
```

Observed:

```
P2: v=3 e_gs = -3.374932598688
P2: v=3 eps=1e-04 ||a-b||_dense/eps = 0.500000
P2: v=3 eps=1e-04 route=native        a-b rel_err = 3.023e-08  sum_mps_weighted rel_err = 2.643e-08
P2: v=3 eps=1e-05 ||a-b||_dense/eps = 0.500000
P2: v=3 eps=1e-05 route=native        a-b rel_err = 3.572e-06  sum_mps_weighted rel_err = 2.288e-06
P2: v=3 eps=1e-06 ||a-b||_dense/eps = 0.500000
P2: v=3 eps=1e-06 route=native        a-b rel_err = 5.733e-04  sum_mps_weighted rel_err = 6.229e-04
P2: v=3 eps=1e-08 ||a-b||_dense/eps = 0.500000
P2: v=3 eps=1e-08 route=native        a-b rel_err = 9.988e-01  sum_mps_weighted rel_err = 4.892e-01
P2: v=python e_gs = -3.374932598688
P2: v=python eps=1e-04 ||a-b||_dense/eps = 0.500000
P2: v=python eps=1e-04 route=densitymatrix a-b rel_err = 2.072e-06  sum_mps_weighted rel_err = 9.190e-07
P2: v=python eps=1e-04 route=exact         a-b rel_err = 2.147e-11  sum_mps_weighted rel_err = 2.044e-11
P2: v=python eps=1e-05 ||a-b||_dense/eps = 0.500000
P2: v=python eps=1e-05 route=densitymatrix a-b rel_err = 8.621e-04  sum_mps_weighted rel_err = 9.801e-04
P2: v=python eps=1e-05 route=exact         a-b rel_err = 1.573e-10  sum_mps_weighted rel_err = 2.126e-10
P2: v=python eps=1e-06 ||a-b||_dense/eps = 0.500000
P2: v=python eps=1e-06 route=densitymatrix a-b rel_err = 1.545e-02  sum_mps_weighted rel_err = 3.080e-03
P2: v=python eps=1e-06 route=exact         a-b rel_err = 3.021e-09  sum_mps_weighted rel_err = 2.352e-09
P2: v=python eps=1e-08 ||a-b||_dense/eps = 0.500000
P2: v=python eps=1e-08 route=densitymatrix a-b rel_err = 7.524e-01  sum_mps_weighted rel_err = 9.900e-01
P2: v=python eps=1e-08 route=exact         a-b rel_err = 2.267e-07  sum_mps_weighted rel_err = 2.164e-07
[exit 0, 5 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: the public MPS arithmetic `a - b` and `H*psi - e*psi` on
`"python"`, since the floor is the documented and chosen behaviour of the
default and v3 has the same floor (the second counter-repro: 5.7e-4 at
5e-7 relative and 1.0 at 5e-9 on v3's `a - b`); the secondary sign on the
energy variance, since both numbers are roundoff of an exact eigenstate;
and the reach into `get_nh_levels` and the ESR Krylov eigenstates, not
shown. LOW, since no public caller was shown to land in the regime.

**Suggested fix**: Pass `method="exact"` from `nhlevels._combine`, for example through a
`method=` keyword on `sum_mps_weighted`. Do not set it unconditionally
inside `sum_mps_weighted`: the ESR Krylov caller passes the whole basis (273
vectors on 12 spins), and the exact route's dense block placement would then
need gigabytes per site.

### 21. On v3, `Chain::sum_mps_weighted` is a tree of truncated pairwise sums, not the one direct sum its comment promises: [A, B, A] with [M, 1, -M] comes back off B by 1.1e-3 at M=1e4 and 1.4e-2 at M=1e5

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `v3-cpp` &middot; candidate `v3-cpp-2`

**Where**: `src/dmrgpy/mpscpp3/chain_session.h:1439-1456 (Chain::sum_mps_weighted and its comment)`; `vendored mpscpp3/ITensor/itensor/mps/mps.cc:1471-1502 (sum(std::vector<MPS>): pairwise tree)`; `vendored mps.cc:739-790 (MPS::orthogonalize: density matrix)`; `src/dmrgpy/nhlevels.py:163-166 (_combine, the consumer)`

The comment of `Chain::sum_mps_weighted` says it is one direct sum, "so a
chain of pairwise sums would truncate the partial sums first". It calls
ITensor v3's `sum(std::vector<MPS>, args)`, which adds in pairs and
recurses, and each pair ends in `MPS::orthogonalize`, a density-matrix
truncation relative to that partial sum. With cancelling coefficients the
error grows with the coefficients: for [A, B, A] with [M, 1, -M] on 10 sites
at cutoff 1e-12, whose exact answer is B, the result is off B by 8.0e-5,
1.1e-3 and 1.4e-2 at M = 1e2, 1e4 and 1e5, and `1-<R|R>` is 1.9e-4 at 1e5.
The reach into `get_nh_levels(vectors=True)` exists but is small: on the
10-site PT-symmetric XX chain one accepted level has `sum|c| = 117`, and its
returned vector is 1.35e-5 off its own Ritz combination.

**Expected**: the direct sum of the comment, truncated once.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/v3-cpp/v3cpp_3_sum_mps_weighted.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Chain::sum_mps_weighted (mpscpp3, 4ed4891) is "one ITensor
# sum(std::vector<MPS>) call", meant to avoid truncating partial sums of
# coefficients that cancel (the Ritz vectors of nhlevels.py). Hypothesis:
# ITensor's sum(std::vector<MPS>, args) is itself a tree of pairwise
# sum(a, b, args), each truncated at Cutoff relative to the PARTIAL sum, so
# with cancelling coefficients the result carries an absolute error of
# about |c| sqrt(cutoff) and loses the small Schmidt weights of the
# survivor. "python" (sum_many: one K-way direct sum, then one truncation
# relative to the result) is the control.
# A. [A, B, A] with coefficients [M, 1, -M]: the exact result is B.
# B. get_nh_levels(vectors=True): the reported residual of each Ritz pair
#    (from Gram matrices) against the true residual of the returned MPS.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext, mps
from dmrgpy.nhdmrg import _residual_scale

N = 10


def dot(sc, i, j):
    return sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]


versions = ["python"]
if cppext.available(3): versions.insert(0, 3)
else: print("backend v3 = not compiled, skipped")

# ---------------------------------------------------------------- part A
for v in versions:
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*N, itensor_version=v)
    sc.maxm = 80 # holds A + B exactly (bond dimension <= 32 each)
    sc.nsweeps = 12
    sc.cutoff = 1e-12
    # A: ground state of an XX chain in a field, B: Heisenberg ground state
    hA = 0
    for i in range(N-1):
        hA = hA + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(N):
        hA = hA + (0.3*(1.0 if i % 3 == 0 else -0.5))*sc.Sz[i]
    sc.set_hamiltonian(hA)
    A = sc.get_gs().copy()
    hB = 0
    for i in range(N-1):
        hB = hB + dot(sc, i, i+1)
    sc.set_hamiltonian(hB)
    B = sc.get_gs().copy()
    nA = np.sqrt(A.dot(A).real); nB = np.sqrt(B.dot(B).real)
    A = A*(1./nA); B = B*(1./nB)
    s = sc._session
    s.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    print("A: v=%s |<A|B>| = %.3e" % (v, abs(A.dot(B))))
    for M in (1.0, 1e2, 1e4, 1e5):
        h = s.sum_mps_weighted([A.cpp_handle, B.cpp_handle, A.cpp_handle],
                               [complex(M), 1.0+0j, complex(-M)])
        R = mps.MPS(sc, cpp_handle=h)
        rr = R.dot(R).real; rb = R.dot(B); bb = B.dot(B).real
        err = np.sqrt(max(rr - 2*rb.real + bb, 0.0))
        print("A: v=%s M=%.0e ||R - B|| = %.3e  <R|R> = %.10f" % (v, M, err, rr))

# ---------------------------------------------------------------- part B
n = 8
for v in versions:
    np.random.seed(11)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm = 40
    sc.nsweeps = 10
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        h = h + 0.3j*((-1)**i)*sc.Sz[i]
    sc.set_hamiltonian(h)
    out = sc.get_nh_levels(vectors=True)
    es, rs, vecs = out
    rscale, e_id = _residual_scale(h.to_terms())
    try:
        m = sc.get_ED_obj().MO2matrix(h)
        m = m.todense() if hasattr(m, "todense") else m
        remin = "%.8f" % np.min(np.linalg.eigvals(np.asarray(m)).real)
    except Exception as exc: # the reference is optional
        remin = "n/a (%s)" % type(exc).__name__
    print("B: v=%s levels found = %d, exact min Re E = %s" % (v, len(es), remin))
    for e, r, y in zip(es, rs, vecs):
        Hy = h*y
        yy = y.dot(y).real
        r2 = Hy.dot(Hy).real - 2*(np.conj(e)*y.dot(Hy)).real + abs(e)**2*yy
        true = np.sqrt(max(r2, 0.0)/yy)/(rscale + abs(e - e_id))
        print("B: v=%s E = %.8f%+.8fi  reported_rel_residual = %.3e  "
              "true_rel_residual = %.3e  ratio = %.2e"
              % (v, e.real, e.imag, r, true, true/max(r, 1e-300)))
```

Observed:

```
A: v=3 |<A|B>| = 8.878e-01
A: v=3 M=1e+00 ||R - B|| = 0.000e+00  <R|R> = 1.0000000000
A: v=3 M=1e+02 ||R - B|| = 8.016e-05  <R|R> = 0.9999999883
A: v=3 M=1e+04 ||R - B|| = 1.103e-03  <R|R> = 0.9999988861
A: v=3 M=1e+05 ||R - B|| = 1.358e-02  <R|R> = 0.9998139120
A: v=python |<A|B>| = 8.878e-01
A: v=python M=1e+00 ||R - B|| = 9.105e-07  <R|R> = 1.0000000000
A: v=python M=1e+02 ||R - B|| = 9.107e-07  <R|R> = 1.0000000000
A: v=python M=1e+04 ||R - B|| = 9.104e-07  <R|R> = 1.0000000000
A: v=python M=1e+05 ||R - B|| = 9.102e-07  <R|R> = 1.0000000000
B: v=3 levels found = 3, exact min Re E = -1.99537545
B: v=3 E = -1.99537545-0.00000000i  reported_rel_residual = 0.000e+00  true_rel_residual = 0.000e+00  ratio = 0.00e+00
B: v=3 E = -1.99537545+0.24463506i  reported_rel_residual = 9.900e-09  true_rel_residual = 0.000e+00  ratio = 0.00e+00
B: v=3 E = -1.99537545-0.24463506i  reported_rel_residual = 0.000e+00  true_rel_residual = 0.000e+00  ratio = 0.00e+00
B: v=python levels found = 4, exact min Re E = -1.99537545
B: v=python E = -1.99537545+0.24463506i  reported_rel_residual = 0.000e+00  true_rel_residual = 0.000e+00  ratio = 0.00e+00
B: v=python E = -1.99537545-0.24463506i  reported_rel_residual = 0.000e+00  true_rel_residual = 0.000e+00  ratio = 0.00e+00
B: v=python E = -1.99537545-0.00000000i  reported_rel_residual = 0.000e+00  true_rel_residual = 1.407e-08  ratio = 1.41e+292
B: v=python E = -1.99537538-0.00000003i  reported_rel_residual = 4.185e-07  true_rel_residual = 3.649e-07  ratio = 8.72e-01
[exit 0, 12 s]
```

Reviewer's counter-repro, the reach into get_nh_levels (shared with the previous finding):

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g07/g07_p1_nhlevels_reach.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Group 07 counter-repro P1 (the "reached" sub-claims of pyitensor-1 and
# v3-cpp-2). Is the sum defect reached by get_nh_levels(vectors=True) on
# the documented PT-symmetric XX chain, i 0.3 (-1)^i Sz_i, 10 sites? The
# script runs levels()'s own steps (nhlevels._right_states,
# _rayleigh_ritz, _accepted, _combine), and for every accepted Ritz pair
# prints the coefficient size, the deviation `dev` of the vector _combine
# returns from the exact combination sum_k c_k |s_k> (dense, from the MPS
# site arrays of the snapshots), and the dense residual of the exact
# combination (`res_combo`) and of the returned vector (`res_returned`),
# both relative as levels() reports them. "python" on both sum routes, v3
# on its own. v3's random start is not seeded by numpy, so the number of
# levels can differ from run to run.
# Run with MKL_NUM_THREADS=1 OMP_NUM_THREADS=1; about 1 to 2 minutes.
import warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext
from dmrgpy import nhlevels as nl
from dmrgpy.pyitensor import mpsalgebra as ma
from dmrgpy.pyitensor import task


def dense_state(s, handle):
    arrs = [np.asarray(x) for x in s.mps_site_arrays(handle)]
    T = arrs[0][0]                      # (d, chi)
    for A in arrs[1:]:                  # (chi_l, d, chi_r)
        T = np.einsum("xa,asb->xsb", T, A).reshape(-1, A.shape[2])
    return T[:, 0]


def dense_mpo(arrs):
    M = np.asarray(arrs[0])[0]          # (d_out, d_in, w)
    for A in arrs[1:]:                  # (w_l, d_out, d_in, w_r)
        A = np.asarray(A)
        M = np.einsum("xyw,wabv->xaybv", M, A)
        M = M.reshape(M.shape[0]*M.shape[1], M.shape[2]*M.shape[3], M.shape[4])
    return M[:, :, 0]


def mpo_arrays(sc, M):
    if sc.itensor_version == "python":
        return task.mpo_arrays(M)
    return sc._session.mpo_site_arrays(M)


def routes(v):
    return ("densitymatrix", "exact") if v == "python" else ("native",)


def with_route(v, route, fn):
    if v != "python":
        return fn()
    old = ma.SUM_METHOD
    ma.SUM_METHOD = route
    try:
        return fn()
    finally:
        ma.SUM_METHOD = old


versions = ["python"]
if cppext.available(3):
    versions.insert(0, 3)
else:
    print("backend v3 = not compiled, skipped")

for v in versions:
    np.random.seed(11)
    n = 10
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm = 40       # 10 sites: the exact bond dimension is 32
    sc.nsweeps = 10
    sc.cutoff = 1e-12
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
    for i in range(n):
        h = h + 0.3j*((-1)**i)*sc.Sz[i]
    sc.set_hamiltonian(h)
    H = sc.hamiltonian
    snaps = nl._right_states(sc, H, None, nl.NSTEPS, nl.NSNAP)
    s = sc._session
    V = np.array([dense_state(s, x.cpp_handle) for x in snaps]).T
    Hd = dense_mpo(mpo_arrays(sc, s.build_operator(H.to_terms())))
    m = sc.get_ED_obj().MO2matrix(h)
    m = m.todense() if hasattr(m, "todense") else m
    e_mpo = np.min(np.linalg.eigvals(Hd).real)
    e_ed = np.min(np.linalg.eigvals(np.asarray(m)).real)
    print("P1: v=%s min Re E dense-MPO = %.10f  ED = %.10f  |diff| = %.2e"
          % (v, e_mpo, e_ed, abs(e_mpo - e_ed)), flush=True)
    w = np.linalg.eigvalsh(V.conj().T @ V)
    kept = w > nl.GRAM_FLOOR*w.max()
    print("P1: v=%s Gram eigenvalues kept = %d of %d, smallest kept/largest = %.2e"
          % (v, int(np.sum(kept)), len(w), np.min(w[kept])/w.max()), flush=True)
    acc = nl._accepted(nl._rayleigh_ritz(snaps, H), H, nl.RITZ_TOL)
    rscale, e_id = nl._scale(H)
    print("P1: v=%s accepted levels = %d" % (v, len(acc)), flush=True)
    for k, (t, r, c) in enumerate(acc):
        scl = rscale + abs(t - e_id)
        yt = V @ c
        nyt = np.linalg.norm(yt)
        res_combo = np.linalg.norm(Hd @ yt - t*yt)/nyt/scl
        print("P1: v=%s level %d E = %.8f%+.8fi  reported = %.3e  res_combo = %.3e  "
              "sum|c| = %.3e  max|c| = %.3e"
              % (v, k, t.real, t.imag, r/scl, res_combo,
                 np.sum(np.abs(c)), np.max(np.abs(c))), flush=True)
        for route in routes(v):
            y = with_route(v, route, lambda: nl._combine(sc, snaps, c))
            vy = dense_state(s, y.cpp_handle)
            dev = np.linalg.norm(vy - yt)/nyt
            res_ret = np.linalg.norm(Hd @ vy - t*vy)/np.linalg.norm(vy)/scl
            print("P1: v=%s level %d route=%-13s dev = %.3e  res_returned = %.3e"
                  % (v, k, route, dev, res_ret), flush=True)
```

Observed:

```
P1: v=3 min Re E dense-MPO = -2.5667821672  ED = -2.5667821672  |diff| = 4.44e-15
P1: v=3 Gram eigenvalues kept = 6 of 24, smallest kept/largest = 3.24e-09
P1: v=3 accepted levels = 4
P1: v=3 level 0 E = -2.56678219+0.26409558i  reported = 4.420e-07  res_combo = 5.882e-07  sum|c| = 1.669e+00  max|c| = 1.182e-01
P1: v=3 level 0 route=native        dev = 9.496e-07  res_returned = 1.035e-06
P1: v=3 level 1 E = -2.56678218-0.00000000i  reported = 5.982e-07  res_combo = 6.106e-07  sum|c| = 2.477e+00  max|c| = 1.941e-01
P1: v=3 level 1 route=native        dev = 1.166e-06  res_returned = 1.507e-06
P1: v=3 level 2 E = -2.56678213-0.26409567i  reported = 8.066e-07  res_combo = 9.421e-07  sum|c| = 2.711e+00  max|c| = 1.918e-01
P1: v=3 level 2 route=native        dev = 1.405e-06  res_returned = 1.354e-06
P1: v=3 level 3 E = -2.56678119+0.00000045i  reported = 3.850e-05  res_combo = 3.850e-05  sum|c| = 1.170e+02  max|c| = 1.308e+01
P1: v=3 level 3 route=native        dev = 1.350e-05  res_returned = 3.585e-05
P1: v=python min Re E dense-MPO = -2.5667821672  ED = -2.5667821672  |diff| = 1.51e-14
P1: v=python Gram eigenvalues kept = 3 of 24, smallest kept/largest = 8.52e-02
P1: v=python accepted levels = 3
P1: v=python level 0 E = -2.56678219-0.26409559i  reported = 3.070e-07  res_combo = 4.680e-07  sum|c| = 1.556e+00  max|c| = 8.925e-02
P1: v=python level 0 route=densitymatrix dev = 6.000e-07  res_returned = 7.221e-07
P1: v=python level 0 route=exact         dev = 6.000e-07  res_returned = 7.221e-07
P1: v=python level 1 E = -2.56678217-0.00000000i  reported = 1.559e-07  res_combo = 3.774e-07  sum|c| = 1.293e+00  max|c| = 7.385e-02
P1: v=python level 1 route=densitymatrix dev = 7.033e-07  res_returned = 1.039e-06
P1: v=python level 1 route=exact         dev = 7.033e-07  res_returned = 1.039e-06
P1: v=python level 2 E = -2.56678212+0.26409572i  reported = 7.954e-07  res_combo = 1.151e-06  sum|c| = 3.765e+00  max|c| = 2.159e-01
P1: v=python level 2 route=densitymatrix dev = 9.737e-07  res_returned = 1.131e-06
P1: v=python level 2 route=exact         dev = 9.737e-07  res_returned = 1.131e-06
[exit 0, 58 s]
```

**Reviewer (CONFIRMED, NARROWED)**: Struck: the `"python"` control "below 1e-6", since `"python"`'s
`sum_mps_weighted` takes the density-matrix route of the previous finding
and is only exact on this probe by its construction; and "the Ritz vectors
of `get_nh_levels` would be far worse than the residuals it reports", since
the hunter's own Part B shows true residuals at or below the reported ones.
The size is narrowed to the measured one and depends on how the tree pairs
the operands.

**Suggested fix**: Pairwise sums with `Cutoff=-1` do not reach roundoff, since every ITensor
pairwise sum ends in a density-matrix orthogonalize; the second
counter-repro of the previous finding shows a single v3 pair sum already
carries the floor (6.2e-4 at 5e-7 relative). Place the K-way direct sum by
hand (the `plussers` loop of `addAssumeOrth` generalized to K, without the
closing orthogonalize), then `position(1)` losslessly and
`position(N, {MaxDim, Cutoff})`, which is a real SVD. Correct the comment in
any case.

### 22. The NRG tail is capped at an absolute 1e-30, so for a band of D of about 4e4 or more the default tail runs into the broken end of the Wilson chain and moves the curve by 7.3e-3 (D=1e5) and 9.2e-2 (D=1e8) of its peak

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `nrg-ysr-metals` &middot; candidate `nrg-ysr-metals-2`

**Where**: `src/dmrgpy/nrgmps.py:603-604 (N_MAX, OMEGA_MIN absolute)`; `src/dmrgpy/nrgmps.py:629-634 (max_sites)`; `src/dmrgpy/nrgmps.py:637-654 (tail_chain default)`; `src/dmrgpy/nrgmps.py:690-695 (_check_prefix, absolute 1e-12)`; `src/dmrgpy/nrgmps.py:784-793 (head_state solve before _check_prefix)`; `src/dmrgpy/nrgmps.py:992-1002 (head="exact": broken tail, no check)`; `src/dmrgpy/wilsonchain.py:159-196 (scale-covariant Lanczos)`

`nrgmps.max_sites` stops the tail where `D Lambda^(-n/2)` reaches
`OMEGA_MIN = 1e-30`, an absolute energy. The Wilson chain itself breaks at a
fixed site index whatever D is, since `wilsonchain.wilson_chain` is a
Lanczos on levels proportional to D: at Lambda = 4 the first site with
`|eps|/t > 1e-3` is 114 at D = 1e5 and at 1e8, the same as on the D = 1
chain of that length, and the hoppings over sites 0 to 112 equal the D = 1
ones times D to 7e-16. So the default tail is 99 sites at D = 1, 109 at
1e3, 116 at 1e5 and 120 at 1e8, and from D of about 4e4 on it ends in broken
sites, with `|eps|/t` up to 2.3e8 at D = 1e5. On a 3-spin chain with the
metal under the middle spin (`head="exact"`, keep 256), the default tail
moves the curve by 7.3e-3 of its peak at D = 1e5 and by 9.2e-2 at D = 1e8,
against a tail cut at 113 sites. With `head="dmrg"`, the default, the
hand-over check `_check_prefix` compares the head and tail chains to an
absolute 1e-12 while they agree to about 7e-17 D, so a correct hand-over
raises `ValueError` from D of about 1e4 on, and the raise comes only after
the DMRG solve of the head.

**Expected**: a tail and a prefix check that are the same at every D, relative to the scale of the band.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/nrg-ysr-metals/nrg_2_tail_chain_absolute_floor.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# nrgmps.max_sites caps the tail chain where D Lambda^(-n/2) reaches the
# ABSOLUTE OMEGA_MIN = 1e-30, while wilsonchain.wilson_chain is a Lanczos
# on levels proportional to D, so it breaks at a fixed site (113 at
# Lambda = 4, |eps|/t up to 8e8 per the max_sites docstring) whatever D is.
# For D = 1 the cap (99 sites) keeps out of the broken region; for a band
# in other units (D = 1e5, 1e8) the default tail includes broken sites.
# No DMRG solve here: only the chains the NRG route would be given.
from dmrgpy import nrgmps, wilsonchain, ysr

L = 4.0
for D in (1.0, 1e3, 1e5, 1e8):
    n = nrgmps.max_sites(L, D)
    eps, t = wilsonchain.wilson_chain(n, Lambda=L, band=(-D, D))
    ratio = np.abs(eps[1:])/np.abs(t) # onsite energy of site k against t_{k-1}
    bad = np.where(ratio > 1e-3)[0]
    print("D = %g: max_sites = %d" % (D, n))
    print("D = %g: t_last/D = %.3e" % (D, abs(t[-1])/D))
    print("D = %g: max |eps_k|/t_(k-1) over the chain = %.3e" % (D, np.max(ratio)))
    print("D = %g: first site with |eps|/t > 1e-3 = %s" % (D,
          (int(bad[0]) + 1) if len(bad) else "none"))
    # the chain is scale covariant: the same n sites on [-1, 1], times D
    e1, t1 = wilsonchain.wilson_chain(n, Lambda=L, band=(-1.0, 1.0))
    dev = np.max(np.abs(t/D - t1)/np.abs(t1))
    print("D = %g: max relative |t/D - t(D=1, same n)| = %.3e" % (D, dev))
    print("D = %g: the D = 1 chain of the same %d sites has max |eps|/t = %.3e"
          % (D, n, np.max(np.abs(e1[1:])/np.abs(t1))))
    # the hand-over prefix check of head_tail_poles, absolute 1e-12
    fc = ysr.spin_chain_on_metal(1, 12, JK=0.3*D, D=D, Lambda=L,
                                 itensor_version="python")
    eh, th = nrgmps.wilson_parameters(fc)
    et, tt = nrgmps.tail_chain(fc)
    diff = max(np.max(np.abs(et[:12] - eh)), np.max(np.abs(tt[:11] - th)))
    try:
        nrgmps._check_prefix(et, tt, eh, th)
        ok = "passes"
    except ValueError as err:
        ok = "raises (%s)" % str(err)[:40]
    print("D = %g: head/tail prefix max abs diff = %.3e, _check_prefix %s"
          % (D, diff, ok))
```

Observed:

```
D = 1: max_sites = 99
D = 1: t_last/D = 3.944e-30
D = 1: max |eps_k|/t_(k-1) over the chain = 6.952e-16
D = 1: first site with |eps|/t > 1e-3 = none
D = 1: max relative |t/D - t(D=1, same n)| = 0.000e+00
D = 1: the D = 1 chain of the same 99 sites has max |eps|/t = 6.952e-16
D = 1: head/tail prefix max abs diff = 2.746e-17, _check_prefix passes
D = 1000: max_sites = 109
D = 1000: t_last/D = 3.852e-33
D = 1000: max |eps_k|/t_(k-1) over the chain = 3.357e-16
D = 1000: first site with |eps|/t > 1e-3 = none
D = 1000: max relative |t/D - t(D=1, same n)| = 8.882e-16
D = 1000: the D = 1 chain of the same 109 sites has max |eps|/t = 5.107e-16
D = 1000: head/tail prefix max abs diff = 6.661e-14, _check_prefix passes
D = 100000: max_sites = 116
D = 100000: t_last/D = 3.009e-35
D = 100000: max |eps_k|/t_(k-1) over the chain = 2.251e+08
D = 100000: first site with |eps|/t > 1e-3 = 114
D = 100000: max relative |t/D - t(D=1, same n)| = 1.421e-15
D = 100000: the D = 1 chain of the same 116 sites has max |eps|/t = 3.917e+02
D = 100000: head/tail prefix max abs diff = 7.276e-12, _check_prefix raises (the tail chain does not start with the W)
D = 1e+08: max_sites = 120
D = 1e+08: t_last/D = 5.151e-08
D = 1e+08: max |eps_k|/t_(k-1) over the chain = 2.249e+24
D = 1e+08: first site with |eps|/t > 1e-3 = 114
D = 1e+08: max relative |t/D - t(D=1, same n)| = 1.944e+14
D = 1e+08: the D = 1 chain of the same 120 sites has max |eps|/t = 1.060e+22
D = 1e+08: head/tail prefix max abs diff = 7.451e-09, _check_prefix raises (the tail chain does not start with the W)
[exit 0, 3 s]
```

Reviewer's counter-repro, the break relative to D and its effect on the curve:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=<repo>/src python3 <scratch>/reviews/g04/nrg2_counter_broken_tail_reached.py
```
```python
import sys; sys.path.insert(0, "<repo>/src"); import numpy as np
# Counter-repro for nrg-ysr-metals-2 (reviewer, group 04).
# A: is the breakdown of wilson_chain at a fixed site index for every D
#    (relative), with the chain clean before it?
# B: does a T = 0 head="exact" NRG on a marginal flow (3 spins, metal under
#    the middle spin, an effectively ferromagnetic exchange) reach the broken
#    end of the default tail at D = 1e5 and 1e8, and does that change the
#    curve against the same run on a tail cut at 113 sites (all clean)?
# No DMRG solve, no random start: head="exact" is deterministic.
import time
from dmrgpy import nrgmps, wilsonchain, ysr

L = 4.0


def first_broken(eps, t, tol=1e-3):
    r = np.abs(eps[1:])/np.abs(t)
    bad = np.where(r > tol)[0]
    return (int(bad[0]) + 1) if len(bad) else None


for D in (1e5, 1e8):
    n = nrgmps.max_sites(L, D)
    e, t = wilsonchain.wilson_chain(n, Lambda=L, band=(-D, D))
    e1, t1 = wilsonchain.wilson_chain(n, Lambda=L, band=(-1.0, 1.0))
    k = 112 # hoppings t[0..111], Wilson sites 0..112
    print("A: D = %g, n = %d: first broken site (|eps|/t > 1e-3) = %s, of "
          "the D = 1 chain of the same n = %s"
          % (D, n, first_broken(e, t), first_broken(e1, t1)))
    print("A: D = %g: max relative |t/D - t(D=1)| over sites 0..112 = %.3e"
          % (D, np.max(np.abs(t[:k]/D - t1[:k])/np.abs(t1[:k]))))
    print("A: D = %g: max |eps_k|/t_(k-1) over sites 1..112 = %.3e (D), "
          "%.3e (D = 1)" % (D, np.max(np.abs(e[1:k+1])/np.abs(t[:k])),
                            np.max(np.abs(e1[1:k+1])/np.abs(t1[:k]))))


def run(D, n_tail, keep=256): # keep 256, not the default 800, for time
    fc = ysr.spin_chain_on_metal(3, 4, site=1, Js=0.5*D, JK=0.5*D,
                                 Lambda=L, D=D, itensor_version="python")
    x = np.logspace(-10, -1, 37)
    es = D*np.concatenate([-x[::-1], x])
    t0 = time.perf_counter()
    _, A, res = ysr.get_stm_spectrum(fc, es, submode="NRG", head="exact",
                                     keep=keep, n_tail=n_tail,
                                     return_poles=True)
    return es, np.asarray(A), res[0], time.perf_counter() - t0


for D in (1e5, 1e8):
    try:
        es, A, r, dt = run(D, None)
    except Exception as err:
        print("B: D = %g default tail: raises %s: %s"
              % (D, type(err).__name__, str(err)[:80]))
        continue
    # the start shell holds Wilson site 0, so the last site is nshells - 1
    print("B: D = %g default tail (%d sites): stop = %s, nshells = %d, "
          "last Wilson site = %d, omega_N/D = %.3e, %.1f s"
          % (D, nrgmps.max_sites(L, D), r.stop, r.nshells, r.nshells - 1,
             r.omega_n[-1]/D, dt))
    if r.nshells - 1 < 114:
        print("B: D = %g: the run stopped before Wilson site 114, the broken "
              "sites were not reached" % D)
        continue
    es2, A2, r2, dt2 = run(D, 113)
    print("B: D = %g tail of 113 sites: stop = %s, nshells = %d, "
          "omega_N/D = %.3e, %.1f s" % (D, r2.stop, r2.nshells,
                                        r2.omega_n[-1]/D, dt2))
    phys = np.abs(es) >= 1e-8*D
    scale = np.max(np.abs(A2))
    print("B: D = %g: max|A_default - A_113|/max|A_113| over |es| >= 1e-8 D "
          "= %.3e" % (D, np.max(np.abs(A[phys] - A2[phys]))/scale))
    print("B: D = %g: max|A_default - A_113|/max|A_113| over all es = %.3e"
          % (D, np.max(np.abs(A - A2))/scale))
```

Observed:

```
A: D = 100000, n = 116: first broken site (|eps|/t > 1e-3) = 114, of the D = 1 chain of the same n = 114
A: D = 100000: max relative |t/D - t(D=1)| over sites 0..112 = 7.105e-16
A: D = 100000: max |eps_k|/t_(k-1) over sites 1..112 = 5.195e-13 (D), 3.109e-15 (D = 1)
A: D = 1e+08, n = 120: first broken site (|eps|/t > 1e-3) = 114, of the D = 1 chain of the same n = 115
A: D = 1e+08: max relative |t/D - t(D=1)| over sites 0..112 = 7.105e-16
A: D = 1e+08: max |eps_k|/t_(k-1) over sites 1..112 = 1.121e-14 (D), 5.329e-16 (D = 1)
B: D = 100000 default tail (116 sites): stop = chain, nshells = 116, last Wilson site = 115, omega_N/D = 3.009e-35, 5.8 s
B: D = 100000 tail of 113 sites: stop = chain, nshells = 113, omega_N/D = 2.407e-34, 4.5 s
B: D = 100000: max|A_default - A_113|/max|A_113| over |es| >= 1e-8 D = 7.292e-03
B: D = 100000: max|A_default - A_113|/max|A_113| over all es = 7.292e-03
B: D = 1e+08 default tail (120 sites): stop = chain, nshells = 120, last Wilson site = 119, omega_N/D = 5.151e-08, 4.4 s
B: D = 1e+08 tail of 113 sites: stop = chain, nshells = 113, omega_N/D = 2.407e-34, 5.9 s
B: D = 1e+08: max|A_default - A_113|/max|A_113| over |es| >= 1e-8 D = 9.181e-02
B: D = 1e+08: max|A_default - A_113|/max|A_113| over all es = 1.275e-01
[exit 0, 25 s]
```

**Reviewer (CONFIRMED, NARROWED)**: The hunter's line "max relative |t/D - t(D=1)|" at D = 1e8 prints 1.9e14,
which met its own refute rule; it was computed over the broken region too.
The counter-repro restricts it to sites 0 to 112, where it is 7e-16, which
confirms the mechanism, and it measures the effect on the curve, which the
hunter had not. LOW, since it needs a band in units where D is above about
1e4.

**Suggested fix**: Drop D from `max_sites`, `n = 2 ln(1/OMEGA_MIN)/ln Lambda`, and call
`OMEGA_MIN` relative to D in the docstring and in `.claude/rules/kondo.md`;
`test_tail_chain_stays_above_the_broken_end`, which asserts 99 at D = 1,
stays valid. Compare `_check_prefix` relative to `max|t_h|` of the head
chain, and run it before `head_state` so a refused hand-over does not pay
for a DMRG solve first. The docstrings of `tail_chain` and `head_tail_poles`
name `N_MAX` as the default, which contradicts the code.

## Left open, and new leads

- Executed, not reviewed: `_check_basis` returns at once for a name outside
  `canonical._DIAGONAL`, which lists `N` but not its diagonal aliases `n`,
  `projOcc`, `projEmp` and `FermiPhase`. So `basis_ops=("n",)` or
  `("projOcc",)` on a number-conserving spinless chain is not refused, while
  `("N",)` is. On 4 sites at T=1 (ED `<Ntot>` = 1.854) the samples froze at
  `<Ntot>` = 3, 1, 1 on `"python"` and 4, 0, 2 on v3, at stderr below 1e-16,
  with the `("N",)` guard firing on both. The fix of the METTS guard-gap
  finding covers it if it tests charges whatever the names are.
- The cost of the spin average of `ysr.get_stm_spectrum(mode="ED")` on a
  pure multiplet (refuted as a finding, since it is recorded as the price of
  the 2026-11-20 change): measured 1.89x on `submode="CVM"` and 2.10x on
  `submode="ED"` for a doublet, with the curves equal to 2e-11 (ratios
  measured while other sessions held the CPU). A shortcut
  would need `[H, S^a_tot] = 0` checked on the ED matrices, since the
  canonical-form proof cannot show it, and it touches a "Do NOT go back to
  the single member" rule, so it is the author's pick.
- v3's `bicstab` recurrence takes ITensor's pairwise `sum()` at every step,
  the same sum the v3 `sum_mps_weighted` finding measures, while its
  `"python"` twin moved to exact sums in `17017a9`. `check_bicgstab` now
  catches a bad answer, so this was not raised as a candidate.
- The BiCGSTAB check may in principle also pass a bad solve, since its
  density-matrix floor drops real directions; not measured. On
  `julia_live` the check sums with `add(alg="densitymatrix")` and would warn
  falsely in the same way at `delta` below about 1e-9; not measured (JIT).
- `evolve_imaginary_time` and `evolve_imaginary_time_and_measure` on a
  1-site `"python"` chain return the input state unevolved:
  `imaginarytime._check_dmrg` guards v3 below 3 sites and `julia_live`
  below 2, not `"python"` below 2, and `_galerkin_sweep` loops over no
  bond. Read, not run.
- The "keep the better pair" logic of `ef9a885` compares certificate
  residuals only, never real parts, so an NH-DMRG pair that converged onto a
  non-lowest eigenpair can replace a right imaginary-time pair without a
  warning. Read, not run.
- `thermal.auto_anneal_method` reads the unit off the simplified H with its
  identity offset and one-site fields, so a Hamiltonian written as
  `sum(S_i.S_j - 1/4)` or in a strong field takes the v3 implicit
  purification at T of order J, where it was measured less accurate than
  Taylor (4.1e-6 against 3.1e-8 at T=1 on 4 sites). `metts_plan` compares a
  relative implicit grid with an absolute TDVP step (`dbeta_half_step=0.05`),
  so its TDVP default is not scale covariant. Read, not run.
- `julia_live` `TDVP_GSE` caps each Krylov vector at `maxlinkdim(psi)+1`
  with ITensorMPS's projection cutoff of 1.5e-8, where `"python"` and v3
  apply at `maxm` with 1e-12. `metts.jl` calls `tdvp_step` without
  `herm=true`, so its exponentiation runs Arnoldi; a possible speedup. Read,
  not run.
- The ESR Krylov routes check only that the final states of the lines
  converge, not that the out-rates of the upper states are complete, so on a
  chain with no conserved quantity and a small `nkry` the widths may come
  out low. Saturated cascades of exactly equal frequencies drop the
  two-photon coherence that `method="floquet"` keeps, unmeasured. At T=0
  `rates.stationary` takes one ground state when two exactly degenerate ones
  are not coupled. Read, not run.
- The NRG route: `twin_chain` copies only `maxm`, `nsweeps`, `cutoff` and
  `noise` onto the other twists; `head="exact"` silently ignores `m`, `k0`,
  `nkry`, `symmetric` and `join`; a negative-U level with a charge doublet
  ground level is handed over in one charge only. Read, not run.
- `task.real_mpo` keeps a real basis of rank `sum(S > 1e-12 S[0])`, which
  overlaps the known depth floor near Wilson site 37 at Lambda = 4.
  `rootndmrg`'s breakdown test at 1e-10 sits below the density-matrix floor
  on `"python"`, which costs spurious poles, not wrong numbers. Read, not
  run.
- Documentation drift: `docs/development.md` 4.6 still says the
  `julia_live` defaults "stay where they were" (BiCGSTAB, TDVP for METTS,
  Taylor), which `0323513` and `b107ace` changed; the docstrings of
  `Many_Body_Chain.metts_vev` and `metts_dynamical_correlator` still say
  they are implemented for `"python"` and 3 only; the 2026-11-15 ledger
  entry on the imaginary-time route still says the smallest real part is
  taken.
