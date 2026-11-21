"""ESR-STM spectra: get_esr_spectrum, as a mixin for the chains whose sites
carry a spin (Spin_Chain and the two spinful fermion chains).

The signal is the DC change of the tunnel current under a continuous RF
drive, as a function of the RF frequency. The drive is Omega V cos(w t),
V = m.S_tip by default (exchange-type driving, m the tip moment), and the
readout is the magnetoresistance of the tip site, D = m.S_tip. The steady
state is that of the rate equation between eigenstates with the Korringa
relaxation of eq. 42 of arXiv:2606.27830 (esrtk/rates.py), and the current
has a DC part, the change of the populations read through D, and a homodyne
part, the coherence of the driven transition read through the RF modulation
of the conductance (esrtk/lineshape.py). For one spin 1/2 it is the Bloch
formula with the asymmetry of eq. 12 of the review of Delgado and Lorente
(arXiv:2208.14105). The derivation and the choice of the method are in
docs/esr_stm_plan.md.
"""
import warnings

import numpy as np

from .esrtk import floquet, lineshape, spectrum

H_EV_GHZ = 4.135667696e-6 # Planck's constant in eV/GHz
KB_EV_K = 8.617333262e-5 # Boltzmann's constant in eV/K
# The relative Boltzmann weight below which a state is not populated, and
# so roughly the accuracy of the DC part. 1e-6 by the author's pick of
# 2026-11-21 (it was 1e-8): about half the states (29 against 56 on 10
# Heisenberg spins at J = 240 GHz, 1.2 K; 56 against 73 on 12 spins at
# J = 30 GHz, 50 mK) for a DC part 1e-5 of its peak off a converged curve
# instead of 1e-7
P_CUT = 1e-6
# submode="KRYLOV" keeps a Ritz state when its residual is below RITZ_TOL of
# the half width 1/(2 pi T2) of a free-spin line (in energy, h/(2 pi T2)):
# its energy is then within that of an eigenvalue
RITZ_TOL = 0.1
# The operator norm of the drive, for the estimate of the saturation of the
# strongest line (lineshape.band_margin), is taken from NORM_ITER steps of
# power iteration on its ED matrix
NORM_ITER = 50


def _ed_operator_norm(chain, V):
    """An estimate of ||V||, the largest |eigenvalue| of the Hermitian V, by
    power iteration on its ED matrix from a fixed random start: from below,
    and exact after two steps for m.S of a spin 1/2 or a spin 1"""
    import scipy.sparse as sp
    M = chain.get_ED_obj().get_operator(V)
    if not sp.issparse(M):
        M = np.asarray(M)
    x = np.random.default_rng(0).standard_normal(M.shape[0])
    x = x/np.linalg.norm(x)
    nv = 0.
    for _ in range(NORM_ITER):
        y = M @ x
        nv = float(np.linalg.norm(y))
        if nv == 0.: break
        x = y/nv
    return nv


def _mps_operator_norm(chain, V, spec):
    """||V|| on the DMRG route: from the local matrix of a one-site V
    (spectrum.local_operator, v3 and "python"), otherwise the largest
    |V_nm| between the kept states, which is only a lower estimate"""
    loc = spectrum.local_operator(chain, V)
    if loc is not None:
        return float(np.linalg.norm(loc[1], 2))
    Vm = np.abs(np.asarray(spec.V))
    return float(np.max(Vm)) if Vm.size else 0.


def _ed_band(chain, V, D, f, h, kT, T1, T2, rabi, p_cut, f_ref, tip_weight,
             ethermal, site_weights):
    """(spec, lines) on mode="ED" with the kept band following the drive.
    The first cut keeps the states up to h(f_max + band) above the
    populated ones (ethermal above the lowest), band the power-broadened
    margin of lineshape.band_margin. Under a saturating drive the upper
    state of a line carries population too, and a second photon reaches
    from it a state above that cut (on a spin 1 at 50 mK that cut leaves
    out the rung |0> to |+1>, and the line is 4.3e-2 low at the default
    drive, 50% at rabi=0.05), so emax is set again to h(f_max + band)
    above the highest state that can carry population (Lines.carry_top),
    and the spectrum is rebuilt until that state stops rising, one pass
    per photon. It ends:
    each pass needs a higher state of the spectrum, and it stops once every
    state the builder solved for is kept."""
    band = lineshape.band_margin(T1, T2, rabi, _ed_operator_norm(chain, V))
    fmax = float(np.max(f))
    emax = h*(fmax + band) + ethermal
    while True:
        spec = spectrum.ed_spectrum(chain, V, D, emax, site_weights)
        lines = lineshape.Lines(spec, f, h, kT, T1, T2, rabi, p_cut, f_ref,
                                tip_weight=tip_weight)
        if spec.reached <= emax: # every state is kept
            break
        need = lines.carry_top + h*(fmax + band)
        # the tolerance is that of a degeneracy, above the rounding of the
        # energies between two eigsh calls on the sparse route
        if not need > emax + lines.deg:
            break
        emax = need
    return spec, lines


def _warn_band(chain, V, spec, lines, emax, f, h, T1, T2, rabi):
    """On mode="DMRG" the kept band is not grown, since more states cost
    more solves: warn when it holds less than lineshape.BAND_WARN of the
    band the drive needs above the highest state that can carry
    population, the band _ed_band keeps on ED"""
    band = lineshape.band_margin(T1, T2, rabi,
                                 _mps_operator_norm(chain, V, spec))
    fmax = float(np.max(f))
    top = lines.carry_top/h # GHz above the lowest state
    kept = emax/h - top - fmax # GHz kept above the window, from top
    if kept >= lineshape.BAND_WARN*band:
        return
    warnings.warn(
        "get_esr_spectrum(mode=\"DMRG\"): the states up to %.4g GHz above the "
        "lowest can take part under this drive (%.4g GHz above the highest "
        "state the drive and the temperature populate, %.4g GHz, an "
        "estimate with the power broadening of the lines), but the states "
        "kept end at %.4g GHz, so the lines into the states above them, and "
        "their tails in the window, are left out. mode=\"ED\" keeps them; a "
        "weaker drive (rabi) needs fewer" % (top + fmax + band, fmax + band,
                                              top, emax/h),
        RuntimeWarning, stacklevel=4)


def _tip_operator(chain, site, tip):
    m = np.asarray(tip, dtype=float)
    if m.shape != (3,) or not np.all(np.isfinite(m)):
        raise ValueError("tip must be a 3-vector, the direction of the tip "
                         "moment, got %r" % (tip,))
    nm = np.linalg.norm(m)
    if nm == 0.:
        raise ValueError("tip must not be the zero vector")
    m = m/nm
    S = (chain.Sx[site], chain.Sy[site], chain.Sz[site])
    return sum(m[a]*S[a] for a in range(3) if m[a] != 0.)


def get_esr_spectrum(chain, f, tip=None, site=0, drive=None, detection=None,
                     T=1.0, T1=100., T2=50., f_ref=None, rabi=1e-3,
                     vrf_over_vdc=0.,
                     phase=0., site_weights=None, mode="ED", submode=None,
                     nstates=16, nkry=None, method="rates", nharm=2,
                     p_cut=P_CUT, h=H_EV_GHZ, kB=KB_EV_K, return_parts=False,
                     return_poles=False):
    """See ESRSpectrumMixin.get_esr_spectrum"""
    if mode not in ("ED", "DMRG"):
        raise ValueError("mode must be \"ED\" or \"DMRG\"")
    if submode not in (None, "KRYLOV"):
        raise ValueError("submode must be None or \"KRYLOV\"")
    if nkry is not None and (isinstance(nkry, bool)
                             or not isinstance(nkry, (int, np.integer))
                             or nkry < 1):
        raise ValueError("nkry must be a positive integer, got %r" % (nkry,))
    if method not in ("rates", "floquet"):
        raise ValueError("method must be \"rates\" or \"floquet\"")
    if method == "floquet" and mode != "ED":
        raise ValueError("method=\"floquet\" is the reference of small "
                         "systems and runs on mode=\"ED\" only")
    if tip is None:
        raise ValueError("get_esr_spectrum needs tip=, the direction of the "
                         "tip moment (it sets the readout D = m.S_tip and, by "
                         "default, the drive)")
    if not 0 <= site < len(chain.Sx):
        raise ValueError("site must be in [0, %d), got %r"
                         % (len(chain.Sx), site))
    if not T >= 0.:
        raise ValueError("T must be >= 0, got %r" % (T,))
    if not rabi >= 0.:
        raise ValueError("rabi must be >= 0, got %r" % (rabi,))
    if not (0. < p_cut < 1.):
        raise ValueError("p_cut must be in (0, 1), got %r" % (p_cut,))
    f = np.asarray(f, dtype=float)
    if f.ndim != 1 or len(f) == 0 or not np.all(np.isfinite(f)):
        raise ValueError("f must be a one-dimensional array of frequencies "
                         "in GHz")
    if not np.min(f) > 0.:
        raise ValueError("the frequencies must be positive (GHz)")
    # required (the author's pick of 2026-11-20): with the centre of the
    # window as the reference, zooming onto a line moved it by 2 to 8%
    if f_ref is None:
        raise ValueError("get_esr_spectrum needs f_ref=, the frequency in GHz "
                         "at which T1 and T2 are given (the Zeeman frequency "
                         "of the tip spin where they were measured)")
    if not f_ref > 0.:
        raise ValueError("f_ref must be a positive frequency in GHz, got %r"
                         % (f_ref,))
    m_op = _tip_operator(chain, site, tip)
    V = m_op if drive is None else drive
    D = m_op if detection is None else detection
    # the formula reads V_nm = conj(V_mn) and a real diagonal of D, and the
    # DMRG route fills both matrices by conjugation
    for name, op in (("drive", drive), ("detection", detection)):
        if op is not None and not chain.is_hermitian(op):
            raise ValueError("%s= must be a Hermitian operator" % name)
    kT = kB*T
    # the states that can take part: populated above p_cut (kT ln(1/p_cut)
    # above the lowest), and a transition inside the window and its margin
    # (lineshape.WINDOW_MARGIN full widths of a free spin with this T2).
    # The states above that, which the upper state of a line reaches only
    # by a thermal upward rate, are left out, and p_cut is then the
    # accuracy of the DC part: it changes by 1.3e-4 at p_cut=1e-4 and by
    # 1.8e-8 at p_cut=1e-8 against a converged curve (10 Heisenberg spins,
    # J = 240 GHz, T = 1.2 K), with 6.6 times fewer states than a second
    # kT ln(1/p_cut) kept at the default p_cut. The 1% window dependence
    # that second term was added for was the homodyne tail of far lines
    # (lineshape.homodyne_taper), not the relaxation: zooming onto a line
    # moves the DC part by 1.7e-5 on 6 spins with one term. This is the
    # band the DMRG and the two Krylov routes keep. On mode="ED" the margin
    # is the power-broadened one and the band follows the drive (_ed_band),
    # and the DMRG eigenstate route warns when its band is short of that
    # (_warn_band).
    if not (T1 > 0. and T2 > 0.):
        raise ValueError("T1 and T2 must be positive, got %r and %r" % (T1, T2))
    # T1 and T2 are those of the tip spin: the Korringa scale is divided by
    # the weight of the tip site (rates.korringa_scale); this also checks
    # site_weights
    tip_weight = float(spectrum._weights(chain, site_weights)[site])
    margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
    emax = h*(np.max(f) + margin)
    if kT > 0.:
        emax += kT*np.log(1./p_cut)
    if mode == "DMRG":
        mode = chain.get_mode(mode="DMRG") # the automatic fallbacks to ED
    if nkry is None: # blocks of npop vectors on ED, single steps on DMRG
        nkry = 100 if mode == "ED" else 600
    if mode == "ED" and submode == "KRYLOV":
        epop = kT*np.log(1./p_cut) if kT > 0. else 0.
        ritz_tol = RITZ_TOL*h/(2.*np.pi*T2)
        window = (h*(np.min(f) - margin), h*(np.max(f) + margin))
        spec = spectrum.ed_krylov_spectrum(chain, V, D, emax, epop, window,
                                           kT, nkry, ritz_tol, site_weights)
    elif mode == "ED":
        spec, lines = _ed_band(chain, V, D, f, h, kT, T1, T2, rabi, p_cut,
                               f_ref, tip_weight,
                               kT*np.log(1./p_cut) if kT > 0. else 0.,
                               site_weights)
    elif submode == "KRYLOV":
        rm = lineshape.READ_MARGIN/(2.*np.pi*T2)
        epop = kT*np.log(1./p_cut) if kT > 0. else 0.
        ritz_tol = RITZ_TOL*h/(2.*np.pi*T2)
        window = (h*(np.min(f) - margin), h*(np.max(f) + margin))
        spec = spectrum.mps_krylov_spectrum(chain, V, D, emax, epop, window,
                                            kT, nkry, ritz_tol, nstates,
                                            site_weights,
                                            read=(h*(np.min(f) - rm),
                                                  h*(np.max(f) + rm)))
    else:
        spec = spectrum.mps_spectrum(chain, V, D, emax, nstates, site_weights)
    if not (mode == "ED" and submode is None): # _ed_band made its own
        lines = lineshape.Lines(spec, f, h, kT, T1, T2, rabi, p_cut, f_ref,
                                tip_weight=tip_weight)
    if mode == "DMRG" and submode is None:
        _warn_band(chain, V, spec, lines, emax, f, h, T1, T2, rabi)
    if method == "rates":
        dc, hom = lines.signal(f, vrf_over_vdc, phase)
    else:
        dc, hom = floquet.signal(lines, f, vrf_over_vdc, phase, nharm)
    lines.warn_positions(h)
    out = [f, dc + hom]
    if return_parts:
        out.append(dict(dc=dc, homodyne=hom))
    if return_poles:
        out.append(lines.poles(h))
    return tuple(out)


class ESRSpectrumMixin():
    """get_esr_spectrum for a chain with Sx/Sy/Sz[site] operators"""
    def get_esr_spectrum(self, f, tip=None, site=0, drive=None,
                         detection=None, T=1.0, T1=100., T2=50., f_ref=None,
                         rabi=1e-3, vrf_over_vdc=0., phase=0.,
                         site_weights=None,
                         mode="ED", submode=None, nstates=16, nkry=None,
                         method="rates", nharm=2,
                         p_cut=P_CUT, h=H_EV_GHZ, kB=KB_EV_K,
                         return_parts=False, return_poles=False):
        """ESR-STM spectrum: the DC change of the current under a CW RF
        drive, against the RF frequency.

        The Hamiltonian of the chain must hold everything static, the
        external field and the field of the tip on the tip site included,
        so every tip position is its own Hamiltonian. Returns (f, dI), dI in
        units of G_j a V_DC, that is the change of <D> for the DC part.

        Arguments:
          f: RF frequencies in GHz.
          tip: direction of the tip moment m (a 3-vector, normalized here).
              The readout is D = m.S_site and the default drive V = m.S_site.
          site: the site under the tip.
          drive, detection: MultiOperators replacing V and D (for instance
              the transverse spin for a drive through a field gradient, or
              an exchange term for a remote spin driven by a nearby atom).
          T: temperature in K (kB= changes the unit).
          T1, T2: relaxation and coherence times in ns of the tip spin, a
              free spin 1/2 on the tip site with Zeeman frequency f_ref,
              whatever the weight of that site in site_weights. T1 fixes the
              scale of the Korringa rates of every pair of states, through
              every site with the weights site_weights, and 1/T2 - 1/(2 T1)
              is the pure dephasing of every coherence.
          site_weights: the weight of each site in the relaxation, one
              non-negative entry per site (equal by default). The scale is
              divided by the weight of the tip site, so only the ratios
              matter; with a weight of 0 on the tip site (relaxation through
              the other sites only), T1 and T2 refer to a site of weight 1.
          f_ref: the frequency in GHz at which T1 and T2 are given, the
              Zeeman frequency of the tip spin where they were measured.
              Required: the Korringa rates grow as f coth(hf/2kT), so a
              reference that followed the window would make every rate
              depend on it.
          rabi: Rabi frequency in GHz of a transition with |V_nm| = 1/2;
              the Rabi frequency of a pair is rabi*2|V_nm|. Small values give
              the linear response, larger ones saturation and power
              broadening.
          vrf_over_vdc, phase: the homodyne part, V_RF/V_DC and the phase
              between the RF voltage and the drive; 0 gives the DC part only.
              Each line enters it only within 30 free-spin widths of its
              resonance (lineshape.homodyne_taper): its dispersive tail
              falls as 1/detuning, and the far background is left out.
          mode: "ED" (full diagonalization up to 2000 states, the lowest
              states of the sparse H above; the states kept reach 30
              power-broadened widths of the strongest line above the window,
              counted from the highest state the temperature or the drive
              populates, so a state a second photon reaches is kept) or
              "DMRG" (the nstates lowest DMRG eigenstates, kept up to 30
              free-spin widths above the window and the thermal states;
              raises if they do not reach every state the window and the
              temperature need, warns when the drive needs more states than
              that, and warns when the energy variance of a strong line's
              states bounds its position only to more than its half width).
          submode: None, or "KRYLOV": the final states of the lines come
              from a Krylov space started from V|m> for every populated m,
              instead of from the eigenstates up to the top of the window.
              On mode="ED" the populated eigenstates and one block Krylov
              space of nkry blocks (default 100). On mode="DMRG"
              (itensor_version 3 and "python") the nstates lowest DMRG
              states must reach the populated states and every state below
              the window, and the Krylov space of nkry steps (default 600)
              runs in the tangent space of each populated state (TASK).
              Faster when few states are populated and many lie in the
              window, slower when most are populated; raises when nkry is
              too small to converge the final states of the window lines.
          method: "rates" (the rate equation, default) or "floquet" (the
              driven Lindblad equation in the same eigenbasis without the
              secular and rotating-wave approximations, nharm harmonics,
              mode="ED" and up to 64 kept states, as a reference).
          p_cut: relative population below which a state is empty, from
              the temperature or, on mode="ED", from the drive, and so
              roughly the relative accuracy of the DC part (1e-6 by
              default: 1e-5 of the DC peak on 10 and 12 spins).
          h, kB: Planck's constant in energy/GHz and Boltzmann's in
              energy/K, eV by default (pass h=1 and kB=20.8366 for H in GHz).
          return_parts: also return dict(dc=..., homodyne=...).
          return_poles: also return the driven transitions: frequency
              (GHz), upper and lower state, weight |V_nm|^2, form_factor
              D_nn - D_mm, population_difference p_m - p_n, width (FWHM in
              GHz at weak drive) and, on mode="DMRG", position_bound (GHz).

        Two states degenerate within the kept spectrum and coupled by V, D or
        the relaxation raise, since their populations depend on the basis
        inside the level (a zero field or a level crossing)."""
        return get_esr_spectrum(self, f, tip=tip, site=site, drive=drive,
                                detection=detection, T=T, T1=T1, T2=T2,
                                f_ref=f_ref, rabi=rabi, vrf_over_vdc=vrf_over_vdc,
                                phase=phase, site_weights=site_weights,
                                mode=mode, submode=submode, nstates=nstates,
                                nkry=nkry, method=method,
                                nharm=nharm, p_cut=p_cut, h=h, kB=kB,
                                return_parts=return_parts,
                                return_poles=return_poles)
