import numpy as np

# get_kondo_spectrum as a mixin, so that every chain whose sites carry a
# spin (Spin_Chain, and the two spinful fermion chains of fermionchain.py)
# gets the same method. Nothing below reads the statistics of the sites:
# the spectrum objects (edkondo.KondoSpectrum, taskkondo.TaskKondoSpectrum)
# and the correlator routes only read the chain's Sx/Sy/Sz[site], its
# Hamiltonian and its session, so on a Hubbard chain the impurity spin is
# the local electron spin c^dagger sigma c/2, charge fluctuations included.

# The dynamical-correlator submodes known to measure from the chain's own
# ground state, the ones get_kondo_spectrum(n_gs>1) can average by setting
# each member of the manifold as that state (measured with set_gs on the
# two members of the 3-site Heisenberg doublet: each of these separates
# them as mode="ED" does). maxent and CVMimag need modules this package
# does not ship, and SECTOR measures from a per-sector solve of its own.
_N_GS_SUBMODES = ("KPM", "CVM", "CVM_explicit", "ROOTN", "TD", "TDZ", "EX",
                  "TASK")


def _task_backends():
    """The itensor_version values that have TASK, as the text of an error"""
    from ..taskdmrg import BACKENDS
    return ", ".join(repr(b) for b in BACKENDS)

class KondoSpectrumMixin():
    """get_kondo_spectrum for a chain with Sx/Sy/Sz[site] operators"""
    def get_kondo_spectrum(self, eV, site=0, Jrho_s=0.0, U=0.0, T=1.0,
                            T0=1.0, omega0=20e-3, Gamma0=5e-6, order=3,
                            kB=8.617333262e-5, mode="ED", **kwargs):
        """Third-order STM/Kondo perturbation-theory tunneling spectrum
        dI/dV(eV) for a single impurity site under the tip, following
        Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.

        The impurity spin is the chain's own Sx/Sy/Sz[site]: a spin site on
        Spin_Chain, the electron spin c^dagger sigma c/2 of an orbital on
        Spinful_Fermionic_Chain and Spinful_Fermionic_Chain_Native. On the
        fermion chains mode="ED" diagonalizes the whole Fock space unless a
        conserved sector is set; Nf is exact at T=0, since S_site conserves
        the electron number, and canonical at T>0, while an Sz sector raises
        (S+- leaves it). submode="TASK" is exact at saturation on both
        classes; on Spinful_Fermionic_Chain the orbital is two MPS sites, so
        S_site|GS> is a two-site change of the ground state, and at a
        truncated maxm it measured less accurate than the native class.

        Parameters (see kondospectrumtk.conductance for the underlying
        equations):
          eV: array of bias energies (same units as the Hamiltonian, e.g.
              eV if built with eV-scale couplings)
          site: chain site index coupled to the tip
          Jrho_s: dimensionless Kondo exchange coupling (J*rho_sample)
          U: dimensionless potential-scattering ratio (eq. "Matrix1")
          T: temperature in Kelvin (kB below is eV/K by default; pass a
             matching kB if your Hamiltonian is in different energy units).
             T=0 is a valid, exact limit, not an error.
          T0: overall tunneling-strength scale (only sets the absolute
              scale of the returned dI/dV, in units of 2*pi*e^2*T0^2/hbar)
          omega0, Gamma0: band cutoff and lifetime broadening for the
              third-order Kondo function F(eps,T)
          order: 2 for the second-order (Fermi golden rule) term alone, 3
              to add the third-order Kondo and (if U!=0) potential-
              interference terms. Every term is the full, bidirectional
              net-current derivative: the two tunneling directions add for
              the second-order and third-order Kondo terms (both even in
              eV) and subtract for the potential-interference term (odd in
              eV, the origin of the bias asymmetry) -- see
              kondospectrumtk/conductance.py's module docstring.
              "LL" sums the leading logs (Jrho_s)^n ln^n to all orders
              (kondospectrumtk/leadinglog.py); its first order in Jrho_s is
              exactly order=3 at the same T. mode="ED" only; at T>0
              the cutoff kernel is thermally averaged with the weight of
              the third order (kondospectrumtk/thermalkernel.py); see
              docs/new_methods/kondo_leading_log/README.md. mem_gb
              (order="LL" only, default None) cuts the biases in chunks
              so that one flow stays below about mem_gb GB; with None
              every bias is one flow on the CPU, which at T>0 holds a
              table per occupied state per bias (leadinglog.memory_model
              gives the size, and the call warns when it is above the
              memory the process can still take).
          mode: "ED" (default) always uses full ED diagonalization of this
              chain's Hamiltonian (every eigenstate is needed as a
              possible virtual intermediate state), via
              kondospectrumtk.edkondo.KondoSpectrum -- independent of this
              chain's own itensor_version/DMRG-vs-ED mode setting, and
              valid at any T>=0. It reads the named parameters only, and
              any other keyword raises TypeError: it is either one of
              the mode="DMRG" parameters below or a misspelling, which
              used to be dropped silently and leave that parameter at
              its default (Jrho= for Jrho_s= returned the order=2
              curve).
              mode="ED" with submode="KRYLOV" (order 2 and 3, any T)
              does not diagonalize the whole Hamiltonian: the initial
              states are its lowest eigenstates from a sparse solver,
              down to a relative Boltzmann weight p_cut (default 1e-8),
              and the final and intermediate states of each one are the
              Ritz states of one block Krylov space started from
              S_x|i>, S_y|i>, S_z|i> of the tip site, of at most 3*nkry
              states (nkry, default 100). The sums are those of
              mode="ED", and they are exact once each Krylov space is
              exhausted. su2=None (the default) takes one member per
              spin multiplet with the weight 2S+1 when the Hamiltonian
              commutes with the total spin (checked) and every state
              otherwise; su2=True raises instead when it does not
              commute, and su2=False takes every state. order="LL" needs the full eigenbasis and
              raises there. See
              docs/new_methods/kondo_krylov_thermal/README.md.
              "DMRG" instead uses this chain's own itensor_version
              (3 or "python", i.e. either the compiled ITensor v3
              extension or the pure-Python pyitensor backend -- both
              expose the identical Chain method surface this feature's
              DMRG-side modules call through, so neither is hardcoded
              anywhere in kondospectrumtk/dmrgtwotime.py or
              secondorder_dc.py) real time evolution throughout, never
              diagonalizing beyond the ground state -- see
              kondospectrumtk/twotime.py's module docstring for the
              construction, and
              examples/kondo_third_order_timing_ED_v3_pyitensor for a
              three-way ED/v3/pyitensor timing comparison. The correlator
              and two-time routes support T=0 only (the T>0 excited-state
              Boltzmann sum this feature was originally scoped around was
              never built for DMRG -- T=0 turned out to admit a cleaner,
              diagonalization-free construction instead, which is what
              shipped); submode="TASK" takes T>0 as the temperature of
              the electrodes with the sample in its ground state, see
              below. The potential-interference term (U!=0, order=3)
              is also supported for DMRG, via
              kondospectrumtk.potentialdc.third_order_potential_dIdV_dc --
              like the second-order term, its T=0 limit collapses to a
              dynamical-correlator convolution (against the F0 kernel
              instead of a Theta0-weighted cumulative sum), so it needs
              no excited-state enumeration either; it carries the same
              general-S-extrapolation caveat as
              conductance.third_order_potential_dIdV (see that
              function's docstring). With mode="DMRG" the extra
              kwargs are forwarded: `submode` (default "TASK" on
              itensor_version="python", 3 and "julia_live", see below,
              and "KPM" on v2, where only order=2 is available),
              `delta`,
              `es` to kondospectrumtk.secondorder_dc.second_order_dIdV_dc
              (second-order term) and
              kondospectrumtk.potentialdc.third_order_potential_dIdV_dc
              (potential-interference term, if U!=0). `es` has no safe
              default and must always be supplied, and it is ONE grid
              shared by both terms, so with U!=0 and order=3 it has to
              meet the potential term's stricter requirement: it must
              reach past the top of the site's S_k spectrum (every
              transition from the ground state that carries weight in
              sum_k |<m|S_k|GS>|^2, plus several delta, and w=0 itself),
              where the second-order term alone only needs the
              transitions below max|eV| -- see both functions'
              docstrings. Any further kwargs go on to
              chain.get_dynamical_correlator via both of the above, and
              are whatever that submode's correlator takes: for
              submode="KPM" that is `kernel`, `hodc_order` and
              `hodc_eta`. The number of KPM moments is not a call
              argument (`n=` raises TypeError): it follows from `delta`
              and the chain's own attributes (kpm_n_scale, kpm_scale,
              kpmmaxm, ...), see kpmdmrg.get_dynamical_correlator.
              `n_gs` (default 1) is the size of the degenerate
              ground-state manifold to average over. mode="ED" defines
              T=0 as the T->0+ limit, the equal-weight average over the
              degenerate ground manifold (edkondo.KondoSpectrum); the
              DMRG route's T=0 is the single state the solver converged
              to unless n_gs is given, which at an accidental degeneracy
              (a level crossing) is whichever member, or superposition
              of members, the random start reached, so the spectrum can
              land anywhere between those of the members themselves.
              With n_gs>1 every term (second order, third-order Kondo,
              potential) is averaged with equal weight over the states
              of get_excited_states(n=n_gs, purify=True), each set as
              the chain's ground state in turn; every term is linear in
              the ground-state density matrix, so this is exactly
              mode="ED"'s average whichever orthonormal basis of the
              manifold DMRG returns. The members are taken to be
              degenerate with the ground state (every submode measures a
              member's transitions from the member's own energy), so
              n_gs is the caller's
              statement of the degeneracy, the same contract as
              `dex` in the ED dynamical correlator: a warning is issued
              when a member lies more than `delta` away from the
              ground-state energy. It costs n_gs times the n_gs=1 run
              plus one excited-state solve, and needs an MPS backend
              with its own session (itensor_version 2, 3 or "python";
              v2 behaves like v3 in every measured run). Each member is
              measured exactly as returned, unswept, and gs_energy() is
              that member's own energy while it is measured (what the
              two-time Kondo term reads, and the ED references'
              convention); the chain's state is put back afterwards,
              gs_energy() included. Only the submodes that
              read the chain's state are accepted with n_gs>1 (KPM, CVM,
              CVM_explicit, ROOTN, TD, TDZ, EX and TASK); SECTOR, which
              measures from its own per-sector solve, and the others
              raise NotImplementedError.
              `dt2`, `n_t2_half`,
              `dtau`, `n_tau_half`
              to kondospectrumtk.dmrgtwotime.two_time_kondo_term_dmrg for
              the third-order Kondo term (order=3) -- these four also have
              no safe default (see that function's docstring: a grid
              fine/wide enough for the default omega0/Gamma0 is
              computationally infeasible to pick automatically, while a
              small fast default silently returns a badly wrong result
              instead of erroring) and must be supplied explicitly.
              Validated against ITensor v3 once a
              compiled backend became available (see
              test_kondo_spectrum_dmrgtwotime.py and
              kondospectrumtk/dmrgtwotime.py's module docstring for what
              that surfaced and fixed): the third-order Kondo term's
              G(t2,tau) matches the ED reference to ~1e-9-1e-10, and the
              swept second-order term (KPM, delta=2e-5) agrees with the
              exact sum to 0.2% at every bias point (it was "a few tens
              of percent at thresholds" until the cumulative integral in
              secondorder_dc.py was made a trapezoid rule on 2026-09-12;
              the error was never KPM's).

          submode="TASK" (itensor_version="python", 3 and
              "julia_live", where it is the default) replaces all of the above by the exact sums of mode="ED" over the tangent-space
              Krylov spectrum of the ground state
              (kondospectrumtk/taskkondo.py, arXiv:2510.07279): no `es`,
              `delta` or time grid, which raise if passed, only `nkry`
              (default 100) and `n_gs`. Exact when maxm saturates the
              Hilbert space (cutoff=0) and nkry exhausts the Krylov space.
              At T=0 every pole is a sharp step, and a warning says how
              much of the weight below max|eV| is in poles the Krylov
              space has not converged (TaskKondoSpectrum.convergence).
              T>0 broadens the steps and the Kondo log with the
              temperature of the electrodes, as mode="ED" does, and keeps
              the sample in its ground state (or in the n_gs manifold),
              which is the thermal state while kT is far below the lowest
              excitation: the call raises ValueError when the Boltzmann
              factor of the lowest pole above the ground state exceeds
              taskkondo.THERMAL_OCCUPATION (1e-3), and mode="ED" is then
              the route. The check sees the states the spin of the site
              reaches from the ground state, not the others. It also
              raises when that lowest pole is not converged in nkry (a
              Ritz energy is an upper bound, so the excitation could be
              lower), when more states lie within 1e-3 kT of the ground
              state than n_gs declares, and when the spin reaches a state
              below the one DMRG converged to.

        Returns (eV, dIdV)."""
        if order not in (2, 3, "LL"):
            raise ValueError("order must be 2, 3 or \"LL\"")
        eV = np.asarray(eV, dtype=float)
        # the mode is validated before anything qualifies it: a mistyped
        # name must not be told that LL is unsupported on it
        if mode not in ("ED", "DMRG"):
            raise ValueError("mode must be \"ED\" or \"DMRG\"")
        if order == "LL" and mode != "ED":
            # the flow sums over every intermediate state and needs
            # products of spin operators that the TASK Krylov space built
            # from S_k|GS> does not hold, so a DMRG spectrum would be a
            # silently truncated one
            raise NotImplementedError(
                "order=\"LL\" is implemented on mode=\"ED\" only; see "
                "docs/new_methods/kondo_leading_log/README.md")
        if mode == "ED":
            # The ED route reads nothing beyond the named parameters, and
            # a bare **kwargs used to swallow the rest silently: since the
            # defaults are physical values, a misspelled one (Jrho= for
            # Jrho_s=, u= for U=, t= for T=) returned a plausible spectrum
            # at the default instead, e.g. exactly the order=2 curve for
            # the Jrho typo (development.md 4.10's "**kwargs with no
            # consumer"). No allow-list for the mode="DMRG" keywords
            # either: delta=/es= would then be inert here, the same hole
            # one level down. The one submode of mode="ED" is "KRYLOV",
            # and its keywords are read only under it; any other submode
            # is a mode="DMRG" keyword and stays in the list below.
            submode = None
            krylov = {}
            if kwargs.get("submode") == "KRYLOV":
                submode = kwargs.pop("submode")
                for k in ("nkry", "p_cut", "su2"):
                    if k in kwargs:
                        krylov[k] = kwargs.pop(k)
            # the memory budget of the leading-log flow, without which a
            # T>0 call on the CPU runs every bias in one flow
            mem_gb = kwargs.pop("mem_gb", None) if (
                order == "LL" and submode is None) else None
            if kwargs:
                raise TypeError(
                    "get_kondo_spectrum(mode=\"ED\"%s) got unexpected keyword "
                    "argument(s): " % (", submode=\"KRYLOV\"" if submode else "")
                    +", ".join(sorted(kwargs))+". They are "
                    "either mode=\"DMRG\" parameters (submode other than "
                    "\"KRYLOV\", delta, es, "
                    "dt2, n_t2_half, dtau, n_tau_half, n_gs and the "
                    "dynamical-correlator keywords), which the exact "
                    "full-spectrum route does not use -- its T=0 already "
                    "averages a degenerate ground manifold -- or "
                    "parameters of submode=\"KRYLOV\" (nkry, p_cut, su2), "
                    "which need that submode, or mem_gb, which needs "
                    "order=\"LL\" without a submode, or not "
                    "parameters of this method at all (check the spelling "
                    "of Jrho_s, U, T, T0, omega0, Gamma0, order, kB)")
            if submode == "KRYLOV" and order == "LL":
                # the flow sums products of spin operators between any two
                # states, so it needs the full eigenbasis, not the Krylov
                # spaces of the thermal states
                raise NotImplementedError(
                    "order=\"LL\" needs the full eigenbasis and is not "
                    "available with submode=\"KRYLOV\"; use mode=\"ED\" "
                    "without a submode")
            if submode == "KRYLOV":
                return self._get_kondo_spectrum_krylov(eV, site, Jrho_s, U, T,
                                                       T0, omega0, Gamma0, order,
                                                       kB, **krylov)
            return self._get_kondo_spectrum_ed(eV, site, Jrho_s, U, T, T0,
                                                omega0, Gamma0, order, kB,
                                                mem_gb=mem_gb)
        elif mode == "DMRG":
            if not T >= 0.:
                raise ValueError("T must be >= 0, got %r" % (T,))
            return self._get_kondo_spectrum_dmrg(eV, site, Jrho_s, U, T0,
                                                  omega0, Gamma0, order,
                                                  T=T, kB=kB, **kwargs)
    def _get_kondo_spectrum_krylov(self, eV, site, Jrho_s, U, T, T0, omega0,
                                    Gamma0, order, kB, nkry=100, p_cut=None,
                                    su2=None):
        from . import krylovkondo
        if not T >= 0.:
            raise ValueError("T must be >= 0, got %r" % (T,))
        p_cut = krylovkondo.P_CUT if p_cut is None else p_cut
        kks = krylovkondo.KrylovKondoSpectrum(self, site, T, kB=kB, nkry=nkry,
                                              p_cut=p_cut, su2=su2)
        return eV, krylovkondo.krylov_kondo_dIdV(kks, eV, Jrho_s, U=U, T0=T0,
                                                 omega0=omega0, Gamma0=Gamma0,
                                                 order=order)
    def _get_kondo_spectrum_ed(self, eV, site, Jrho_s, U, T, T0, omega0,
                                Gamma0, order, kB, mem_gb=None):
        from .edkondo import KondoSpectrum
        from . import conductance
        from .stepfunctions import cached_fbuilder
        ks = KondoSpectrum(self, site, T, kB=kB)
        if order == "LL":
            from .leadinglog import leading_log_dIdV
            return eV, leading_log_dIdV(ks, eV, Jrho_s, U=U, T0=T0,
                                        omega0=omega0, Gamma0=Gamma0,
                                        mem_gb=mem_gb)
        dIdV = conductance.second_order_dIdV(ks, eV, T0=T0, U=U)
        if order == 3:
            # shared between both calls below, and kept across calls at
            # the same T, omega0, Gamma0 and kB: building it tabulates F
            # (0.2 to 0.6 s, see FBuilder)
            Fb = (cached_fbuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB)
                  if T > 0. else None)
            dIdV = dIdV + conductance.third_order_kondo_dIdV(
                    ks, eV, Jrho_s, T0=T0, omega0=omega0, Gamma0=Gamma0, Fb=Fb)
            if U != 0.0:
                dIdV = dIdV + conductance.third_order_potential_dIdV(
                        ks, eV, Jrho_s, U, T0=T0, omega0=omega0, Fb=Fb,
                        Gamma0=Gamma0)
        return eV, dIdV
    def _get_kondo_spectrum_dmrg(self, eV, site, Jrho_s, U, T0, omega0,
                                  Gamma0, order, submode=None, delta=None,
                                  es=None, dt2=None, n_t2_half=None,
                                  dtau=None, n_tau_half=None, n_gs=1,
                                  T=0., kB=8.617333262e-5, **dc_kwargs):
        # first, before the submode is resolved or anything is built
        if (isinstance(n_gs, bool) or not isinstance(n_gs, (int, np.integer))
                or n_gs < 1):
            raise ValueError("n_gs must be a positive integer, the size of "
                             "the degenerate ground-state manifold to "
                             "average over, got %r" % (n_gs,))
        implicit = submode is None
        if submode is None:
            # TASK wherever it exists (taskdmrg.BACKENDS: "python", 3 and
            # "julia_live"): exact at saturation, no grids to choose, and
            # orders of magnitude faster than the two-time TDVP Kondo term
            # on the same chain. v2 has no TASK, so it gets the correlator
            # route, which covers order=2 only: the two-time Kondo term
            # needs the TDVP of itensor_version="python", 3 or "julia_live".
            from ..taskdmrg import BACKENDS
            submode = "TASK" if self.itensor_version in BACKENDS else "KPM"
        if submode == "TASK":
            # No frequency or time grid at all: the exact sums of the ED
            # route (conductance.py) on the tangent-space Krylov spectrum
            # (kondospectrumtk/taskkondo.py). The grid keywords of the other
            # DMRG routes would be inert here, so they raise.
            nkry = dc_kwargs.pop("nkry", 100)
            # checked here and not only in TaskKondoSpectrum, which a chain
            # too short for the session never builds
            if (isinstance(nkry, bool)
                    or not isinstance(nkry, (int, np.integer)) or nkry < 1):
                raise ValueError("nkry must be a positive integer, got %r"
                                 % (nkry,))
            unused = sorted(k for k, v in (("delta", delta), ("es", es),
                    ("dt2", dt2), ("n_t2_half", n_t2_half), ("dtau", dtau),
                    ("n_tau_half", n_tau_half)) if v is not None)
            unused += sorted(dc_kwargs)
            if unused:
                raise TypeError(
                    "get_kondo_spectrum(mode=\"DMRG\", submode=\"TASK\") "
                    "takes nkry and n_gs only, it sums the poles exactly "
                    "with no frequency or time grid; got "+", ".join(unused)
                    +(". TASK is the default on itensor_version in %s; "
                      "pass submode=\"KPM\" (or another correlator submode) "
                      "for the correlator and two-time route these grids "
                      "belong to" % (_task_backends(),) if implicit else ""))
            Fb = None
            if T > 0. and order == 3:
                # the thermal Kondo log, tabulated once for every member
                # and kept across calls at the same T, omega0, Gamma0, kB
                from .stepfunctions import cached_fbuilder
                Fb = cached_fbuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB)
            def terms():
                from .taskkondo import TaskKondoSpectrum
                from . import conductance
                if self.get_mode(mode="DMRG") == "DMRG":
                    # T>0 is the temperature of the electrodes, with the
                    # sample in its ground state; TaskKondoSpectrum raises
                    # when that is not the thermal state; a member of an
                    # n_gs manifold may sit up to Gamma0 off its partner,
                    # the split the n_gs average accepts below
                    ks = TaskKondoSpectrum(self, site, nkry=nkry, kB=kB,
                                           T=T, n_gs=n_gs,
                                           manifold_window=Gamma0 if n_gs > 1 else 0.)
                    # the poles the Krylov space has not converged are
                    # steps at energies that move with nkry
                    if len(eV): ks.warn_unconverged(np.max(np.abs(eV)))
                else: # a chain too short for the session (mode.py): the
                    # exact spectrum, as every other DMRG route falls back
                    from .edkondo import KondoSpectrum
                    ks = KondoSpectrum(self, site, T=T, kB=kB)
                dIdV = conductance.second_order_dIdV(ks, eV, T0=T0, U=U)
                if order == 3:
                    dIdV = dIdV + conductance.third_order_kondo_dIdV(
                            ks, eV, Jrho_s, T0=T0, omega0=omega0,
                            Gamma0=Gamma0, Fb=Fb)
                    if U != 0.0:
                        dIdV = dIdV + conductance.third_order_potential_dIdV(
                                ks, eV, Jrho_s, U, T0=T0, omega0=omega0,
                                Gamma0=Gamma0, Fb=Fb)
                return dIdV
            split_tol = Gamma0 # no broadening: the F kernel's own width
        else:
            if T != 0.:
                # the correlator and two-time routes are built on the T=0
                # limit of the step functions and have no thermal form
                raise ValueError(
                    "get_kondo_spectrum(mode=\"DMRG\", submode=%r) only "
                    "supports T=0; T>0 is implemented for submode=\"TASK\" "
                    "(itensor_version=\"python\", 3 and \"julia_live\") "
                    "and for mode=\"ED\"" % (submode,))
            if delta is None: delta = 2e-6
            split_tol = delta
            short = self.get_mode(mode="DMRG") != "DMRG"
            if (order == 3 and not short
                    and self.itensor_version not in ("python", 3,
                                                     "julia_live")):
                raise NotImplementedError(
                    "get_kondo_spectrum(mode=\"DMRG\", order=3) with "
                    "submode=%r: the two-time Kondo term runs on the TDVP "
                    "of itensor_version=\"python\", 3 and \"julia_live\" only "
                    "(got "
                    "itensor_version=%r); use order=2, mode=\"ED\", or one "
                    "of those backends" % (submode, self.itensor_version))
            def terms():
                """Every requested term on the state that is this chain's
                ground state right now"""
                from .secondorder_dc import second_order_dIdV_dc
                dIdV = second_order_dIdV_dc(self, site, eV, T0=T0, U=U,
                                             mode="DMRG", submode=submode,
                                             delta=delta, es=es, **dc_kwargs)
                if order == 3 and short:
                    # a chain too short for the session (mode.py): the
                    # correlators above fell back to ED, and so does the
                    # Kondo term, which has no session to evolve
                    from .edkondo import KondoSpectrum
                    from . import conductance
                    dIdV = dIdV + conductance.third_order_kondo_dIdV(
                            KondoSpectrum(self, site, T=0.), eV, Jrho_s,
                            T0=T0, omega0=omega0, Gamma0=Gamma0)
                elif order == 3:
                    from .dmrgtwotime import two_time_kondo_term_dmrg
                    term = two_time_kondo_term_dmrg(
                            self, site, eV, omega0=omega0, Gamma0=Gamma0,
                            dt2=dt2, n_t2_half=n_t2_half, dtau=dtau,
                            n_tau_half=n_tau_half)
                    dIdV = dIdV + 4*np.pi*T0**2*Jrho_s*term
                if order == 3 and U != 0.0:
                    from .potentialdc import third_order_potential_dIdV_dc
                    dIdV = dIdV + third_order_potential_dIdV_dc(
                            self, site, eV, Jrho_s, U, T0=T0, omega0=omega0,
                            Gamma0=Gamma0, mode="DMRG", submode=submode,
                            delta=delta, es=es, **dc_kwargs)
                return dIdV
        if n_gs == 1: return eV, terms() # the single converged state
        # n_gs>1: equal-weight average over an orthonormal basis of the
        # manifold, mode="ED"'s T->0+ limit (edkondo.KondoSpectrum). Every
        # term is linear in the ground-state density matrix, so the
        # average is Tr[P0 X]/n_gs whichever basis DMRG returns, while a
        # single state gives <psi|X|psi>, anywhere between the members'
        # own values: on an S=1 impurity at its |0>/|-1> crossing the
        # zero-bias dI/dV/2pi is 1-<Sz_0> of the converged state, anywhere
        # in [1.0, 2.0] against 1.5 (2026-09-24 hole hunt, finding 12).
        # Detecting the degeneracy from energies instead would need a
        # tolerance that no DMRG run can honour, since nothing separates
        # "exactly degenerate" from "split below what the sweep resolved",
        # hence the caller-supplied n_gs.
        if submode not in _N_GS_SUBMODES:
            # an allow-list, not a deny-list (development.md 4.10): a
            # submode is averaged only if it is known to measure from the
            # chain's own state. SECTOR measures from its own per-sector
            # solve on a clone and returned the n_gs=1 value for every
            # member; EX did the same from its cached basis until it was
            # reprojected (2026-09-24b hole hunt, finding 15).
            raise NotImplementedError(
                "get_kondo_spectrum(n_gs=%d) averages over the manifold by "
                "setting each member as the chain's ground state, and "
                "submode=%r does not read that state (it measures from a "
                "state of its own), so every member would give the same "
                "number. Use one of %s." % (n_gs, submode,
                                             ", ".join(_N_GS_SUBMODES)))
        session = getattr(self, "_session", None)
        if (self.get_mode(mode="DMRG") != "DMRG" or session is None
                or not hasattr(session, "set_wavefunction")):
            raise NotImplementedError(
                "n_gs>1 sets each member of the ground-state manifold as "
                "the ground state of this chain's MPS session in turn, and "
                "this chain has none (itensor_version=%r, or a fallback to "
                "ED). mode=\"ED\" averages the degenerate manifold at T=0 "
                "on its own." % (self.itensor_version,))
        import warnings
        from .. import groundstate
        gs0 = self.get_gs() # the solved state, put back afterwards
        e0 = self.gs_energy()
        energies, members = self.get_excited_states(n=n_gs, purify=True)
        if len(members) < n_gs:
            raise RuntimeError("n_gs=%d: get_excited_states returned only %d "
                               "independent states" % (n_gs, len(members)))
        split = np.max(np.abs(np.real(np.asarray(energies)) - e0))
        if split > split_tol:
            warnings.warn(
                "get_kondo_spectrum(n_gs=%d): a member of the averaged "
                "manifold lies %.3g away from the ground-state energy, more "
                "than %.3g (delta, or Gamma0 under submode=\"TASK\"). Every member is weighted equally, as in "
                "the T->0+ average over a degenerate manifold, so this is "
                "not a T=0 average unless those states are meant to be "
                "degenerate." % (n_gs, split, split_tol), RuntimeWarning,
                stacklevel=3)
        # set_gs() marks each member as injected, so the first correlator
        # of terms() hands it to the session unswept, where KPM and TD read
        # it, and gives it its own energy <wf|H|wf> as self.e0, which is
        # what the two-time Kondo term and the ED references measure each
        # initial state from. The band edges KPM rescales with were filled
        # by get_excited_states() above from the solved state, so every
        # member is measured on the same window. This loop used to install
        # each member with set_gs() plus session.set_wavefunction(), which
        # dropped the session's energy, and the correlator's own ground-
        # state re-verification then ran a real DMRG sweep from the
        # member, which relaxed the upper one onto the lower one in 5 of 6
        # runs on v3 at a split below delta (2026-09-24b hole hunt,
        # finding 13); the dropped energy was not, as this comment said,
        # something nothing read again.
        #
        # The chain's state is restored as a unit afterwards: wf0, e0,
        # computed_gs, the solver key and the injection mark, plus the
        # session's own state. Restoring wf0 alone left gs_energy() on the
        # last member's energy (finding 14).
        snapshot = (self.wf0, self.e0, self.computed_gs,
                    getattr(self, "_gs_solver_key", None),
                    getattr(self, "_gs_injected", None),
                    getattr(self, "_gs_supplied", False))
        total = 0.
        try:
            for wf in members:
                self.set_gs(wf)
                total = total + terms()
        finally:
            (self.wf0, self.e0, self.computed_gs, self._gs_solver_key,
             self._gs_injected, self._gs_supplied) = snapshot
            session.set_wavefunction(groundstate.detached_copy(gs0).cpp_handle)
        return eV, total/len(members)
