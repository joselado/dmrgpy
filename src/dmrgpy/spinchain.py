from .manybodychain import Many_Body_Chain
import numpy as np
from .algebra import algebra
from . import effectivehamiltonian
from . import pychainwrapper
from . import multioperator

class Coupling():
  def __init__(self,i,j,g):
    """Store a two-site coupling constant g between sites i and j"""
    self.i = i
    self.j = j
    self.g = g

Spin_Chain = Many_Body_Chain

# dictionary for the sites, with a more readable nomenclature
label2site = dict() # dictionary
label2site["1/2"] = 2
label2site["S=1/2"] = 2
label2site[2] = 2
label2site["1"] = 3
label2site["S=1"] = 3
label2site[3] = 3
label2site["3/2"] = 4
label2site["S=3/2"] = 4
label2site[4] = 4
label2site["2"] = 5
label2site["S=2"] = 5
label2site[5] = 5
label2site["5/2"] = 6
label2site["S=5/2"] = 6
label2site["S=3"] = 7
label2site[6] = 6


def get_logdimension(self):
    """Return the logarithm of the dimension"""
    return np.sum(np.log(np.array(self.sites))) # return dimension



def get_site(label):
    if label in label2site: return label2site[label]
    else: return None


class Spin_Chain(Many_Body_Chain):
    """Class for spin Hamiltonians"""
    def __init__(self,sites,**kwargs):
        """Build a spin chain from a list of site labels (e.g. "1/2",
        "1", "3/2", ..., see label2site)"""
        sites = [label2site[s] for s in sites]
        Many_Body_Chain.__init__(self,sites,**kwargs)
        # default exchange constants
        self.use_ampo_hamiltonian = True # use ampo
        self.pychain_object = None # pychain object
        self.Sx = [self.get_operator("Sx",i) for i in range(self.ns)]
        self.Sy = [self.get_operator("Sy",i) for i in range(self.ns)]
        self.Sz = [self.get_operator("Sz",i) for i in range(self.ns)]
        self.Si = [self.Sx,self.Sy,self.Sz]
    def SS(self,i,j):
        """Return the Heisenberg dot product S_i . S_j"""
        return self.Sx[i]*self.Sx[j] + self.Sy[i]*self.Sy[j] + self.Sz[i]*self.Sz[j]
    def set_fields(self,fun):
        """Set the local magnetic field term of the Hamiltonian"""
        h = 0
        for i in range(self.ns):
            b = fun(i)
            for j in range(3):  h = h + b[j]*self.Si[j][i]
        self.fields = h
        self.hamiltonian = self.exchange + self.fields # update Hamiltonian
    def test(self,ntries=3,**kwargs):
        """Check the anticommunation relations"""
        Sx = self.Sx
        Sy = self.Sy
        Sz = self.Sz
        for ii in range(ntries):
            i = np.random.randint(self.ns)
            j = np.random.randint(self.ns)
            op = Sx[i]*Sy[j] - Sy[j]*Sx[i]
            if i==j: op = op - 1j*Sz[i]
            if not self.is_zero_operator(op,**kwargs): raise
    def get_logdimension(self):
        """Return the logarithm of the Hilbert space dimension"""
        return get_logdimension(self)
    def set_exchange(self,fun):
        """Set the exchange coupling between sites"""
        def terms():
            for i in range(self.ns): # loop
              for j in range(self.ns):  # loop
                g = fun(i,j).real # call the function
                if np.sum(np.abs(fun(i,j)-fun(j,i)))>1e-5: raise # something wrong
                one = np.identity(3) # identity matrix
                g = g*one # multiply by the identity
                for ii in range(3):
                  for jj in range(3):
                      yield g[ii,jj]*self.Si[ii][i]*self.Si[jj][j]
        h = multioperator.msum(terms())
        self.exchange = h # exchange matrix
        self.hamiltonian = self.exchange + self.fields # update Hamiltonian
    def get_ED_obj(self):
        """Return the ED object (pychain wrapper), building it if not
        already cached"""
        if self.has_ED_obj:
            return self.ED_obj
        else:
            self.ED_obj = pychainwrapper.get_pychain(self)
            self.has_ED_obj = True # store
            return self.ED_obj
    def get_pychain(self):
        """Return the underlying pychain object"""
        return pychainwrapper.get_pychain(self)
    def get_full_hamiltonian(self):
        """Return the full Hamiltonian"""
        from . import pychainwrapper
        return pychainwrapper.get_full_hamiltonian(self)
    def get_magnetization(self,**kwargs):
        """Return the magnetization on each site, and save it to
        MAGNETIZATION.OUT"""
        mx = [self.vev(self.Sx[i],**kwargs) for i in range(self.ns)]
        my = [self.vev(self.Sy[i],**kwargs) for i in range(self.ns)]
        mz = [self.vev(self.Sz[i],**kwargs) for i in range(self.ns)]
        np.savetxt("MAGNETIZATION.OUT",np.array([mx,my,mz]).T)
        return np.array([mx,my,mz]).real
    def get_full_SS_correlator(self,**kwargs):
        """Return the full spin correlator"""
        from .dynamicstk import spincorrelators
        return spincorrelators.get_full_SS_correlator(self,**kwargs)
    def get_effective_hamiltonian(self,**kwargs):
        """Return the effective Hamiltonian"""
        return effectivehamiltonian.get_effective_hamiltonian(self,
                    name="XX",**kwargs)
    def get_hamiltonian(self):
        """Return Hamiltonian as a multioperator"""
        if self.hamiltonian is not None: return self.hamiltonian
        else: # conventional way
            sxs = [self.get_operator("Sx",i) for i in range(self.ns)]
            sys = [self.get_operator("Sy",i) for i in range(self.ns)]
            szs = [self.get_operator("Sz",i) for i in range(self.ns)]
            ss = [sxs,sys,szs]
        out = multioperator.msum(c.g[i,j]*ss[i][c.i]*ss[j][c.j]
                for c in self.exchange # exchange coupling
                for i in range(3) for j in range(3))
        out.clean()
        if len(self.fields)>0:
            fieldterms = multioperator.msum(b[j]*ss[j][i]
                    for i,b in enumerate(self.fields) for j in range(3))
            out = out + fieldterms
        # still have to add the fields!!
        return out # return multioperator
    def get_kondo_spectrum(self, eV, site=0, Jrho_s=0.0, U=0.0, T=1.0,
                            T0=1.0, omega0=20e-3, Gamma0=5e-6, order=3,
                            kB=8.617333262e-5, mode="ED", **kwargs):
        """Third-order STM/Kondo perturbation-theory tunneling spectrum
        dI/dV(eV) for a single impurity site under the tip, following
        Ternes, New J. Phys. 17 063016 (2015), arXiv:1505.04430.

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
              interference terms. The third-order terms are d(I^{t->s})/dV
              only -- see third_order_kondo_dIdV's docstring for why (the
              paper never gives a general t<->s formula for them); the
              second-order term is the full, bidirectional net-current
              derivative.
          mode: "ED" (default) always uses full ED diagonalization of this
              chain's Hamiltonian (every eigenstate is needed as a
              possible virtual intermediate state), via
              kondospectrumtk.edkondo.KondoSpectrum -- independent of this
              chain's own itensor_version/DMRG-vs-ED mode setting, and
              valid at any T>=0.
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
              three-way ED/v3/pyitensor timing comparison. Only T=0 is
              supported (the T>0 excited-state
              Boltzmann sum this feature was originally scoped around was
              never built for DMRG -- T=0 turned out to admit a cleaner,
              diagonalization-free construction instead, which is what
              shipped). The potential-interference term (U!=0, order=3)
              is also supported for DMRG, via
              kondospectrumtk.potentialdc.third_order_potential_dIdV_dc --
              like the second-order term, its T=0 limit collapses to a
              dynamical-correlator convolution (against the F0 kernel
              instead of a Theta0-weighted cumulative sum), so it needs
              no excited-state enumeration either; it carries the same
              general-S-extrapolation caveat as
              conductance.third_order_potential_dIdV (see that
              function's docstring). Extra
              kwargs are forwarded: `submode` (default "KPM"), `delta`,
              `es` to kondospectrumtk.secondorder_dc.second_order_dIdV_dc
              (second-order term) and
              kondospectrumtk.potentialdc.third_order_potential_dIdV_dc
              (potential-interference term, if U!=0) (`es` has no safe
              default -- see either function's docstring -- and must be
              supplied when order requests it). Any further kwargs (e.g.
              `n`, the number of KPM moments -- its own default of 1000
              is accurate but expensive; a compiled-backend smoke test
              may want far fewer) are forwarded on to
              chain.get_dynamical_correlator via both of the above, for
              whichever mode/submode combination is in use. `dt2`, `n_t2_half`,
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
              swept second-order term (KPM) agrees to within a few tens
              of percent at thresholds, consistent with the expected
              delta-broadening/moment-truncation error.

        Returns (eV, dIdV)."""
        if order not in (2, 3): raise ValueError("order must be 2 or 3")
        eV = np.asarray(eV, dtype=float)
        if mode == "ED":
            return self._get_kondo_spectrum_ed(eV, site, Jrho_s, U, T, T0,
                                                omega0, Gamma0, order, kB)
        elif mode == "DMRG":
            if T != 0.:
                raise ValueError("mode=\"DMRG\" only supports T=0")
            return self._get_kondo_spectrum_dmrg(eV, site, Jrho_s, U, T0,
                                                  omega0, Gamma0, order,
                                                  **kwargs)
        else: raise ValueError("mode must be \"ED\" or \"DMRG\"")
    def _get_kondo_spectrum_ed(self, eV, site, Jrho_s, U, T, T0, omega0,
                                Gamma0, order, kB):
        from .kondospectrumtk.edkondo import KondoSpectrum
        from .kondospectrumtk import conductance
        from .kondospectrumtk.stepfunctions import FBuilder
        ks = KondoSpectrum(self, site, T, kB=kB)
        dIdV = conductance.second_order_dIdV(ks, eV, T0=T0, U=U)
        if order == 3:
            # shared between both calls below: building it tabulates an
            # expensive adaptive-quadrature integral (see FBuilder)
            Fb = FBuilder(T, omega0=omega0, Gamma0=Gamma0, kB=kB) if T>0. else None
            dIdV = dIdV + conductance.third_order_kondo_dIdV(
                    ks, eV, Jrho_s, T0=T0, omega0=omega0, Gamma0=Gamma0, Fb=Fb)
            if U != 0.0:
                dIdV = dIdV + conductance.third_order_potential_dIdV(
                        ks, eV, Jrho_s, U, T0=T0, omega0=omega0, Fb=Fb,
                        Gamma0=Gamma0)
        return eV, dIdV
    def _get_kondo_spectrum_dmrg(self, eV, site, Jrho_s, U, T0, omega0,
                                  Gamma0, order, submode="KPM", delta=2e-6,
                                  es=None, dt2=None, n_t2_half=None,
                                  dtau=None, n_tau_half=None, **dc_kwargs):
        from .kondospectrumtk.secondorder_dc import second_order_dIdV_dc
        dIdV = second_order_dIdV_dc(self, site, eV, T0=T0, U=U, mode="DMRG",
                                     submode=submode, delta=delta, es=es,
                                     **dc_kwargs)
        if order == 3:
            from .kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg
            term = two_time_kondo_term_dmrg(
                    self, site, eV, omega0=omega0, Gamma0=Gamma0, dt2=dt2,
                    n_t2_half=n_t2_half, dtau=dtau, n_tau_half=n_tau_half)
            dIdV = dIdV + 4*np.pi*T0**2*Jrho_s*term
            if U != 0.0:
                from .kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
                dIdV = dIdV + third_order_potential_dIdV_dc(
                        self, site, eV, Jrho_s, U, T0=T0, omega0=omega0,
                        Gamma0=Gamma0, mode="DMRG", submode=submode,
                        delta=delta, es=es, **dc_kwargs)
        return eV, dIdV

Spin_Hamiltonian = Spin_Chain # backwards compatibility
