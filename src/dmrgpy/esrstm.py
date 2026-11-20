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
import numpy as np

from .esrtk import floquet, lineshape, spectrum

H_EV_GHZ = 4.135667696e-6 # Planck's constant in eV/GHz
KB_EV_K = 8.617333262e-5 # Boltzmann's constant in eV/K
# The relative Boltzmann weight below which a state is not populated
P_CUT = 1e-8


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
                     phase=0., site_weights=None, mode="ED", nstates=16,
                     method="rates", nharm=2, p_cut=P_CUT, h=H_EV_GHZ,
                     kB=KB_EV_K, return_parts=False, return_poles=False):
    """See ESRSpectrumMixin.get_esr_spectrum"""
    if mode not in ("ED", "DMRG"):
        raise ValueError("mode must be \"ED\" or \"DMRG\"")
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
    # moves the DC part by 1.7e-5 on 6 spins with one term.
    if not (T1 > 0. and T2 > 0.):
        raise ValueError("T1 and T2 must be positive, got %r and %r" % (T1, T2))
    margin = lineshape.WINDOW_MARGIN/(np.pi*T2)
    emax = h*(np.max(f) + margin)
    if kT > 0.:
        emax += kT*np.log(1./p_cut)
    if mode == "DMRG":
        mode = chain.get_mode(mode="DMRG") # the automatic fallbacks to ED
    if mode == "ED":
        spec = spectrum.ed_spectrum(chain, V, D, emax, site_weights)
    else:
        spec = spectrum.mps_spectrum(chain, V, D, emax, nstates, site_weights)
    lines = lineshape.Lines(spec, f, h, kT, T1, T2, rabi, p_cut, f_ref)
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
                         mode="ED", nstates=16, method="rates", nharm=2,
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
          T1, T2: relaxation and coherence times in ns of a free spin 1/2 on
              the tip site with Zeeman frequency f_ref. T1 fixes the scale
              of the Korringa rates of every pair of states, through every
              site with the weights site_weights (equal by default), and
              1/T2 - 1/(2 T1) is the pure dephasing of every coherence.
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
              states of the sparse H above) or "DMRG" (the nstates lowest
              DMRG eigenstates; raises if they do not reach every state the
              window and the temperature need, and warns when the energy
              variance of a strong line's states bounds its position only to
              more than its half width).
          method: "rates" (the rate equation, default) or "floquet" (the
              driven Lindblad equation in the same eigenbasis without the
              secular and rotating-wave approximations, nharm harmonics,
              mode="ED" and up to 64 kept states, as a reference).
          p_cut: relative Boltzmann weight below which a state is empty,
              and so roughly the relative accuracy of the DC part.
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
                                mode=mode, nstates=nstates, method=method,
                                nharm=nharm, p_cut=p_cut, h=h, kB=kB,
                                return_parts=return_parts,
                                return_poles=return_poles)
