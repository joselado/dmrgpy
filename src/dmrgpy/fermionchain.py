from .manybodychain import Many_Body_Chain
import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .fermionchaintk import staticcorrelator
from .fermionchaintk import hamiltonian
from . import funtk
from . import gap
from . import multioperator

class Fermionic_Chain(Many_Body_Chain):
    """Class for fermionic Hamiltonians"""
    def __init__(self,n,**kwargs):
        self.F = [self.get_operator("F",i) for i in range(n)] # Fermi string
        self.C = [self.get_operator("C",i) for i in range(n)]
        self.Cdag = [self.get_operator("Cdag",i) for i in range(n)]
        self.A = [self.get_operator("A",i) for i in range(n)]
        self.Adag = [self.get_operator("Adag",i) for i in range(n)]
        self.N = [self.get_operator("N",i) for i in range(n)]
#        self.N = [self.Cdag[i]*self.C[i] for i in range(n)]
        self.Id = self.get_operator("Id",1)
        Many_Body_Chain.__init__(self,[0 for i in range(n)],**kwargs)
        self.fermionic = True
        self.use_ampo_hamiltonian = True # use ampo
    def get_charge_gap(self,**kwargs):
        """Return the charge gap"""
        return gap.sector_gap(self,sum(self.N),**kwargs)
    def set_hoppings(self,fun):
        """Add the spin independent hoppings"""
        self.set_hoppings_MB(fun)
    def get_logdimension(self):
        return len(self.C)*np.log(2) # log dimension
    def get_density_spinless(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_spinless(self,**kwargs)
    def get_density(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_spinless(self,**kwargs)
    def set_hubbard_spinless(self,fun):
        """ Hubbard term """
        hamiltonian.set_hubbard_spinless(self,fun)
    def set_hubbard(self,fun):
        """ Hubbard term """
        hamiltonian.set_hubbard_spinless(self,fun)
    def get_density_fluctuation_spinless(self,**kwargs):
        """Return the electronic density fluctuations"""
        return staticcorrelator.get_density_fluctuation_spinless(self,**kwargs)
    def get_density_fluctuation(self,**kwargs):
        """Return the electronic density fluctuations"""
        return staticcorrelator.get_density_fluctuation_spinless(self,**kwargs)
    def get_pairing(self):
        """
        Return the superfluid density
        """
        return staticcorrelator.get_pairing_spinless(self,**kwargs)
    def vev_spinless(self,MO,**kwargs):
        """ Return a vacuum expectation value"""
        return self.vev(MO,**kwargs)
    def excited_vev_spinless(self,MO,mode="DMRG",**kwargs):
        """ Return a vaccum expectation value"""
        if mode=="DMRG": return self.excited_vev_MB(MO,**kwargs)
        elif mode=="ED": return self.get_ED_obj().excited_vev(MO,**kwargs) 
    def excited_vev(self,MO,**kwargs): 
        return self.excited_vev_spinless(MO,**kwargs)
    def hamiltonian_free(self,pairs=[[]]):
        """
        Return the free part of the fermionic Hamiltonian
        """
        m = np.zeros((self.ns,self.ns),dtype=np.complex128) # matrix
        for key in self.hoppings:
              t = self.hoppings[key]
              m[t.i,t.j] = t.g
        return m
#    def get_excited(self,mode="DMRG",**kwargs):
#          """
#          Wrapper for static correlator
#          """
#          if mode=="DMRG": # using DMRG
#            return Many_Body_Chain.get_excited(self,**kwargs)
#          elif mode=="ED":
#            MBF = self.get_ED_obj() # get the object
#            return algebra.lowest_eigenvalues(MBF.h,**kwargs)
    def get_correlator_free(self,pairs=[[]]):
          """Get the correlator for free fermions"""
          m = self.hamiltonian_free() # get the single body matrix
          (es,vs) = lg.eigh(m) # diagonalize
          vs = vs.transpose()
          out = []
          for p in pairs:
              o = 0.0 # initialize
              for (e,v) in zip(es,vs):
                  if e<=0.0: 
                      if self.spinful: # spinful Hamiltonian
                          for i in range(2):
                            o += v[2*p[0]+i]*np.conjugate(v[2*p[1]+i]) # add
                      else: 
                            o += v[p[0]]*np.conjugate(v[p[1]]) # add
              out.append(o)
          return np.array(out) # return
    def gs_energy_free(self):
        """Get the energy for free fermions"""
        m = self.hamiltonian_free() # get the single body matrix
        es = lg.eigvalsh(m) # get the energies
        return np.sum(es[es<0.0]) # return energies
    def get_gr(self,**kwargs):
        return get_gr(self,**kwargs)
    def get_gr_free(self,**kwargs):
        return get_gr_free(self,**kwargs)
    def gs_energy(self,mode="DMRG",**kwargs):
        """Compute ground state energy, overrriding the method"""
        if mode=="DMRG": 
            return Many_Body_Chain.gs_energy(self,**kwargs)
        elif mode=="ED":
            MBF = self.get_ED_obj()
            return algebra.lowest_eigenvalues(MBF.h,n=1)[0]
        else: raise # unrecognised
    def get_ED_obj(self):
        """
        Return the many body fermion object
        """
        if self.has_ED_obj: # if the ED object has been computed
            return self.ED_obj # return the stored object
        else:
            MBf = mbfermion.MBFermion(self.ns) # create object
            MBf.add_multioperator(self.hamiltonian) # add the Hamiltonian
            self.ED_obj = MBf # store the object
            self.has_ED_obj = True # set to True
            return self.ED_obj # return the object
    def execute(self,f):
        """
        This is a temporal fix to use the C operators in Julia ITensor
        """
        if self.itensor_version=="julia": # use the fermionic representation
            from . import multioperator
            multioperator.use_jordan_wigner = False
        return Many_Body_Chain.execute(self,f)






def get_gr_free(self,es=np.linspace(-10.,10.,800),delta=0.1,i=0,j=0):
    m = self.hamiltonian_free() # get the single body matrix
#    print(m)
    y = np.zeros(es.shape[0],dtype=np.complex128) # output
    iden = np.identity(m.shape[0])
    for ii in range(len(es)):
        yi = np.matrix(m-(es[ii]+1j*delta)*iden).I[i,j]
        y[ii] = yi # store
    return es,y













def get_gr(self,delta=0.002,es=np.linspace(-10.0,10.0,800),i=0,j=0):
    """Compute the advanced Green's function"""
    from . import kpmdmrg
    (x1,y1) = kpmdmrg.get_dynamical_correlator(self,es=es,i=i,j=j,
            name="cdc",delta=delta)
    (x2,y2) = kpmdmrg.get_dynamical_correlator(self,es=es,i=i,j=j,
            name="ccd",delta=delta)
    x1 = x1 + self.e0 # shift by the fermi energy
    x2 = x2 + self.e0 # shift by the fermi energy
    # define interpolating function
    from scipy.interpolate import interp1d
    f1r = interp1d(x1, y1.real,fill_value=0.0,bounds_error=False)
    f1i = interp1d(x1, y1.imag,fill_value=0.0,bounds_error=False)
    f2r = interp1d(x2, y2.real,fill_value=0.0,bounds_error=False)
    f2i = interp1d(x2, y2.imag,fill_value=0.0,bounds_error=False)
    # compute the result
    yr = f1r(es) #+ f2r(-es) # real part
    yi = f1i(es) #- f2i(-es) # imaginary part
    # now add the imaginary part
    from scipy.signal import hilbert
    y = yr + 1j*hilbert(yr) + 1j*yi - hilbert(yi)
#    y = 1j*y
#    y = 1j*yr
    return (es,y)



class Majorana_Chain(Fermionic_Chain):
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        nf = (n+1)//2 # number of conventional fermions
        nf = max([2,nf]) # fix
        super().__init__(nf) # initialize the Hamiltonian
        # define the Majorana operators
        G = [0 for i in range(nf*2)] # empty list
        for jf in range(nf): # loop over fermions
            G[2*jf] = (self.C[jf] + self.Cdag[jf])/np.sqrt(2)
            G[2*jf+1] = 1j*(self.C[jf] - self.Cdag[jf])/np.sqrt(2)
        self.G = [G[i] for i in range(n)] # sotre those operators
        del self.C  # clean
        del self.Cdag  # clean
        del self.N # clean





class Spinful_Fermionic_Chain(Fermionic_Chain):
    """
    Class to deal with fermionic Hamiltonians with
    spin degree of freedom
    """
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        super().__init__(2*n) # initialize the Hamiltonian
        self.Sx = [0.5*self.Cdag[2*i]*self.C[2*i+1] +
                0.5*self.Cdag[2*i+1]*self.C[2*i] for i in range(n)]
        self.Sy = [-0.5*1j*self.Cdag[2*i]*self.C[2*i+1] +
                1j*0.5*self.Cdag[2*i+1]*self.C[2*i] for i in range(n)]
        self.Sz = [0.5*self.N[2*i] +
                (-1)*0.5*self.N[2*i+1] for i in range(n)]
        self.Delta = [0.5*self.C[2*i]*self.C[2*i+1] for i in range(n)]
        self.Cup = [self.C[2*i] for i in range(n)]
        self.Cdagup = [self.Cdag[2*i] for i in range(n)]
        self.Cdn = [self.C[2*i+1] for i in range(n)]
        self.Cdagdn = [self.Cdag[2*i+1] for i in range(n)]
        self.Nup = [self.N[2*i] for i in range(n)]
        self.Ndn = [self.N[2*i+1] for i in range(n)]
        self.Ntot = [self.Nup[i]+self.Ndn[i] for i in range(n)]
        self.use_ampo_hamiltonian = True # use ampo
    def get_density_spinful(self,**kwargs):
        """
        Return the density in each site, summing over spin channels
        """
        return staticcorrelator.get_density_spinful(self,**kwargs)
    def get_density(self,**kwargs):
        """
        Return the density in each site, summing over spin channels
        """
        return staticcorrelator.get_density_spinful(self,**kwargs)
    def get_magnetization(self,**kwargs):
        """Return magnetization"""
        return staticcorrelator.get_magnetization_spinful(self,**kwargs)
    def get_onsite_pairing(self,**kwargs):
        """
        Return the expectation value of the onsite pairing
        """
        return staticcorrelator.get_onsite_pairing_spinful(self,**kwargs)
    def set_hubbard_spinful(self,fun):
        """
        Add Hubbard interation in a spinful manner
        The Hubbard term will be defined as
        n_i n_j, with n_i = n_{i,up} + n_{i,,down}
        """
        hamiltonian.set_hubbard_spinful(self,fun)
    def set_hubbard(self,fun):
        """
        Add Hubbard interation in a spinful manner
        The Hubbard term will be defined as
        n_i n_j, with n_i = n_{i,up} + n_{i,,down}
        """
        hamiltonian.set_hubbard_spinful(self,fun)
    def set_swave_pairing(self,fun):
        """
        Add onsite swave pairing to a spinful Hamiltonian
        The pairing term is of the form
        Delta_i c_{i,up} c_{i,down} + h.c.
        """
        hamiltonian.set_swave_pairing_spinful(self,fun)
    def get_density_fluctuation_spinful(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_fluctuation_spinful(self,**kwargs)
    def get_density_fluctuation(self,**kwargs):
        """Return the electronic density"""
        return staticcorrelator.get_density_fluctuation_spinful(self,**kwargs)
    def set_hoppings_spinful(self,fun):
        """
        Function to Add hopping in a spinful manner
        """
        def fun2(i,j):
            if i%2==j%2: return fun(i//2,j//2)
            return 0.0
        self.set_hoppings(fun2)



class Spinful_Fermionic_Chain_Native(Many_Body_Chain):
    """
    Spinful fermionic (Hubbard-like) chain built directly on native
    spinful sites: each orbital is a *single* tensor-network site with a
    local Hilbert space of dimension 4 (Empty, Up, Down, UpDn -- ITensor
    v3's own "Electron"/"Hubbard" site, mpscpp3/get_sites.h's site-type
    code 1: HubbardSite, an alias for ElectronSite -- see
    mpscpp3/ITensor/itensor/mps/sites/electron.h), instead of the two
    separate spinless-fermion sites per orbital that
    Spinful_Fermionic_Chain uses (site-type code 0, interleaved
    up/down). Jordan-Wigner strings are still threaded explicitly at the
    Python level (multioperatortk/jordanwigner_spinful.py), exactly as
    done for the spinless case in multioperatortk/jordanwigner.py -- this
    backend does not rely on ITensor AutoMPO's own automatic fermionic
    sign insertion (autompo.cc's isFermionic()/fermionicTerm()), even
    though it exists and would also work here; keeping the Jordan-Wigner
    transform in Python keeps every backend (ED, the C++ extensions, the
    pure-Python engine, Julia) fed from the exact same term list.

    Only itensor_version=3 wires up a native spinful site on the DMRG
    side today: mpscpp2/get_sites.h never registered a Hubbard/Electron
    site type, and neither pyitensor (itensor_version="python") nor
    mpsjulialive (itensor_version="julia_live") implement one. A DMRG
    call with any other itensor_version simply is not available for this
    class -- unlike an uncompiled C++ extension, there is no automatic
    ED fallback for a site type mode.py doesn't know how to build, so
    such calls raise from mpscpp2/pyitensor/mpsjulialive's own site-type
    dispatch. Cross-checking against mode="ED" always works (see
    get_ED_obj() below): ED never builds tensor-network sites at all, it
    diagonalizes a plain occupation-number Fock space directly.

    Performance (measured, not just argued from theory): halving the
    number of tensor-network sites (n instead of 2n) does keep every
    nearest-neighbor hopping term nearest-neighbor and every
    Jordan-Wigner string half as long, so the a priori expectation was
    that this class would be faster. Directly benchmarked against
    Spinful_Fermionic_Chain instead (same Hubbard-chain Hamiltonian,
    same itensor_version=3, same maxm/nsweeps/noise, n=6..14): it is
    NOT. At matched bond dimension, Spinful_Fermionic_Chain's DMRG run
    is both faster (roughly 2-3x at maxm>=40, growing with n) *and*
    converges to a slightly lower/more accurate energy -- this class
    never wins on either axis in the sizes tested. The reason is that
    ITensor v3's two-site dmrg() sweep cost is governed by the LOCAL
    physical dimension of the two-site block being diagonalized/SVD'd,
    which scales with the *square* of each site's local dimension: a
    two-site block here spans two dimension-4 sites (16 combined local
    states) versus two dimension-2 sites (4 combined local states) on
    the interleaved chain. That per-step cost increase more than
    outweighs sweeping over half as many sites -- the entanglement
    entropy across a cut (what actually sets the bond dimension needed
    for a given accuracy) does not improve just from repackaging two
    flavors already at the same physical location into one site instead
    of two, so there is no compensating win on the memory/bond-dimension
    side either.

    The same disadvantage was checked, and confirmed, in every other
    regime that could plausibly have favored native sites instead:
    strong on-site coupling (U/t=20, where the two flavors' entanglement
    is purely local -- doubled still wins, at every bond dimension
    tried); long-range/power-law-decaying hopping (where native sites'
    shorter Jordan-Wigner strings should shrink the MPO the most --
    doubled still wins at matched bond dimension, from maxm=20 up to
    160); real-time evolution via two-site TDVP (same combined-local-
    dimension penalty as two-site DMRG -- doubled wins by a wider margin
    once an actual entangling quench is run, not just a near-static
    check); one-site TDVP with global subspace expansion
    ("TDVP_GSE", the one method whose cost is linear rather than
    quadratic in local dimension, so the best a priori candidate -- near
    parity on a weakly-entangling test, but doubled again slightly ahead
    once a real quench is run); and the KPM dynamical correlator itself
    (doubled faster at every kpmmaxm tried, by a smaller margin than the
    ground-state case -- roughly 1.2-1.6x rather than 2-3x, since KPM's
    per-moment cost is an MPO-MPS application/truncation rather than a
    two-site diagonalization+SVD, but still not a win).

    One case DOES flip in this class's favor: the 4-point correlator
    tensor <Cdag_i C_j Cdag_k C_l> (mps.MPS.get_four_correlation_tensor(),
    entropytk/correlationentropy.py). This class exposes self.C/self.Cdag
    as flat, single-flavor-per-entry lists (mode 2*i=up, 2*i+1=down,
    matching Spinful_Fermionic_Chain's own indexing exactly, so the two
    classes' tensors are directly comparable index for index) purely so
    the existing, backend-agnostic ctmode="explicit" implementation works
    unchanged. That implementation is a Python loop of independent
    static overlaps <wf|Op|wf> (one per (i,j,k,l), each a single-shot
    MPO-MPS-MPS contraction, not an iterative two-site variational
    search), so it does not pay the combined-local-dimension penalty
    that dominates two-site DMRG/TDVP. Measured (n=3,4,5,6,12 orbitals,
    same Hubbard-chain Hamiltonian): this class's ctmode="explicit" beats
    not only Spinful_Fermionic_Chain's own ctmode="explicit", but also
    its specialized ctmode="full" (mpscpp3/chain_session.h's
    C++-accelerated AutoMPO implementation) at n=3..6, by a margin that
    grows with n there (n=6: ~13-22s vs ~16-35s depending on run) --
    but that growth does NOT continue indefinitely: at n=12 (24 flat
    modes) the two are back to essentially tied, ~700s each, with
    ctmode="full" very slightly ahead (697s vs 707s) -- see
    examples/four_correlation_tensor_spinful_native. So the win is real
    but bounded to smaller/moderate sizes, not an asymptotic advantage;
    treat "native wins here" as size-dependent, not as a rule that
    extrapolates to larger n. ctmode="full" is not available for this
    class at all: it hardcodes the literal "Cdag"/"C" operator names,
    which ITensor's ElectronSite does not define (only
    Cup/Cdn/Cdagup/Cdagdn are), so there was nothing to lose by
    comparing against it here.

    Kept as an alternative backend (correctness cross-checked exactly
    against ED and against Spinful_Fermionic_Chain, for the ground-state
    energy, the KPM dynamical correlator, and the 4-point correlator
    tensor) -- a genuine, if narrow and size-bounded, performance edge
    for the one static-overlap calculation checked so far, not a general-purpose
    speedup for anything iterative.
    """
    def __init__(self,n,**kwargs):
        """Create the sites"""
        self.Cup = [self.get_operator("Cup",i) for i in range(n)]
        self.Cdagup = [self.get_operator("Cdagup",i) for i in range(n)]
        self.Cdn = [self.get_operator("Cdn",i) for i in range(n)]
        self.Cdagdn = [self.get_operator("Cdagdn",i) for i in range(n)]
        self.Nup = [self.get_operator("Nup",i) for i in range(n)]
        self.Ndn = [self.get_operator("Ndn",i) for i in range(n)]
        self.Ntot = [self.Nup[i]+self.Ndn[i] for i in range(n)]
        self.Sx = [0.5*self.Cdagup[i]*self.Cdn[i] +
                0.5*self.Cdagdn[i]*self.Cup[i] for i in range(n)]
        self.Sy = [-0.5*1j*self.Cdagup[i]*self.Cdn[i] +
                0.5*1j*self.Cdagdn[i]*self.Cup[i] for i in range(n)]
        self.Sz = [0.5*self.Nup[i] + (-1)*0.5*self.Ndn[i] for i in range(n)]
        self.Delta = [0.5*self.Cup[i]*self.Cdn[i] for i in range(n)]
        # Flat, single-flavor-per-entry operator lists (2*n entries: mode
        # 2*i = up, 2*i+1 = down at physical site i), matching exactly
        # the same flat-mode convention used throughout
        # jordanwigner_spinful.py and the interleaved
        # Spinful_Fermionic_Chain's own self.C/self.Cdag (there, C[2*i]
        # and C[2*i+1] literally *are* two separate spinless-fermion
        # sites; here they are the same Cup[i]/Cdn[i] MultiOperators
        # already built above, just re-exposed under the generic
        # name/indexing entropytk/correlationentropy.py's
        # get_four_correlation_tensor()/get_correlation_matrix() expect
        # (via wf.MBO.C/wf.MBO.Cdag) -- this is what makes those
        # backend-agnostic functions work unchanged for this class too,
        # and keeps the resulting tensor/matrix indices directly
        # comparable, index for index, against Spinful_Fermionic_Chain's.
        self.C = []
        self.Cdag = []
        for i in range(n):
            self.C.append(self.Cup[i]); self.C.append(self.Cdn[i])
            self.Cdag.append(self.Cdagup[i]); self.Cdag.append(self.Cdagdn[i])
        self.Id = self.get_operator("Id",1)
        Many_Body_Chain.__init__(self,[1 for i in range(n)],**kwargs)
        self.fermionic = True
        self.use_ampo_hamiltonian = True # use ampo
    def get_charge_gap(self,**kwargs):
        """Return the charge gap"""
        return gap.sector_gap(self,sum(self.Ntot),**kwargs)
    def set_hoppings_spinful(self,fun):
        """
        Add a spin-conserving hopping fun(i,j), for both flavors
        """
        h = self.generate_bilinear(fun,self.Cdagup,self.Cup)
        h = h + self.generate_bilinear(fun,self.Cdagdn,self.Cdn)
        self.hopping = h # store
        self.update_hamiltonian()
    def set_hoppings(self,fun):
        """Add the spin-conserving hoppings"""
        self.set_hoppings_spinful(fun)
    def set_hubbard_spinful(self,fun):
        """
        Add Hubbard interaction in a spinful manner
        The Hubbard term will be defined as
        n_i n_j, with n_i = n_{i,up} + n_{i,down}
        """
        self.hubbard = self.generate_bilinear(fun,self.Ntot,self.Ntot)
        self.update_hamiltonian()
    def set_hubbard(self,fun):
        """
        Add Hubbard interaction in a spinful manner
        """
        self.set_hubbard_spinful(fun)
    def set_swave_pairing(self,fun):
        """
        Add onsite swave pairing to a spinful Hamiltonian
        The pairing term is of the form
        Delta_i c_{i,up} c_{i,down} + h.c.
        """
        h = multioperator.msum(fun(i)*self.Delta[i]
                for i in range(len(self.Delta)))
        self.pairing = h + h.get_dagger()
        self.update_hamiltonian()
    def get_density(self,**kwargs):
        """
        Return the density in each site, summing over spin channels
        """
        return np.array([self.vev(Ni,**kwargs).real for Ni in self.Ntot])
    def get_magnetization(self,**kwargs):
        """Return magnetization"""
        mx = np.array([self.vev(self.Sx[i],**kwargs).real
                for i in range(len(self.Sx))])
        my = np.array([self.vev(self.Sy[i],**kwargs).real
                for i in range(len(self.Sy))])
        mz = np.array([self.vev(self.Sz[i],**kwargs).real
                for i in range(len(self.Sz))])
        return np.array([mx,my,mz]).T # return magnetization
    def get_onsite_pairing(self,**kwargs):
        """
        Return the expectation value of the onsite pairing
        """
        return np.array([self.vev(Di,**kwargs) for Di in self.Delta])
    def get_density_fluctuation(self,**kwargs):
        """Return the electronic density fluctuations"""
        d = self.get_density(**kwargs) # total density
        d2 = np.array([self.vev(Ni*Ni,**kwargs).real for Ni in self.Ntot])
        return d2 - d**2 # return density fluctuations
    def get_ED_obj(self):
        """
        Return the many body fermion object, built on the 2*n flat
        fermionic modes (up,down per orbital, see
        jordanwigner_spinful.py's mode-ordering convention) that this
        chain's n native spinful sites represent -- unlike
        Fermionic_Chain.get_ED_obj(), self.ns is the number of *sites*
        (n), not the number of fermionic modes (2*n), so it cannot be
        used directly here.
        """
        if self.has_ED_obj: return self.ED_obj
        from .pyfermion import mbfermion
        MBf = mbfermion.MBFermion(2*len(self.Cup))
        MBf.add_multioperator(self.hamiltonian)
        self.ED_obj = MBf
        self.has_ED_obj = True
        return self.ED_obj



class Spinon_Chain(Spinful_Fermionic_Chain):
    """Class for spinon chains"""
    def get_gs(self,**kwargs):
        """Redefine the ground state method"""
        if self.computed_gs: return self.wf0 # return the wavefunction
        P = 1. # parton projector
        for i in range(len(self.Sx)): 
            P = P*(-2*self.Nup[i]*self.Ndn[i] + self.Ndn[i] + self.Nup[i])
            #P = P*(1.- self.Nup[i]*self.Ndn[i]) #  Gutzwiller projection
        from .mpsalgebra import mpsarnoldi
        super().gs_energy(**kwargs) # get the GS
        wf = self.wf0
#        wf = None # no initial guess
#        print("Projection",wf.dot(P*wf).real)
        wf = mpsarnoldi(self,self.hamiltonian,mode="GS",P=P,wf=wf)
        print("Projection",wf.dot(P*wf).real)
        self.computed_gs = True # computed GS
        self.wf0 = wf # store ground state
        return wf
    def gs_energy(self,**kwargs):
        wf = self.get_gs(**kwargs) # ground state
        return wf.dot(self.hamiltonian*wf).real # return energy










Fermionic_Hamiltonian = Fermionic_Chain
Spinful_Fermionic_Hamiltonian = Spinful_Fermionic_Chain


def isfermion(self):
    """Function to determine if an object is a valid fermionic object"""
    from .pyfermion.mbfermion import MBFermion
    if type(self)==Fermionic_Chain: return True
    if type(self)==Spinful_Fermionic_Chain: return True
    if type(self)==Spinful_Fermionic_Chain_Native: return True
    if type(self)==Spinful_F_Fermionic_Chain: return True
    if type(self)==MBFermion: return True
    else: return False
    



class Spinful_F_Fermionic_Chain(Fermionic_Chain):
    """
    Class to deal with fermionic Hamiltonians with
    spin degree of freedom
    """
    def __init__(self,n):
        """ Rewrite the init method to take twice as many sites"""
        super().__init__(3*n) # initialize the Hamiltonian
        self.Sx = [0.5*self.Cdag[3*i]*self.C[3*i+1] +
                0.5*self.Cdag[3*i+1]*self.C[3*i] for i in range(n)]
        self.Sy = [-0.5*1j*self.Cdag[3*i]*self.C[3*i+1] +
                1j*0.5*self.Cdag[3*i+1]*self.C[3*i] for i in range(n)]
        self.Sz = [0.5*self.N[3*i] +
                (-1)*0.5*self.N[3*i+1] for i in range(n)]
        self.Cup = [self.C[3*i] for i in range(n)]
        self.Cdagup = [self.Cdag[3*i] for i in range(n)]
        self.Cdn = [self.C[3*i+1] for i in range(n)]
        self.Cdagdn = [self.Cdag[3*i+1] for i in range(n)]
        self.F = [self.C[3*i+2] for i in range(n)]
        self.Fdag = [self.Cdag[3*i+2] for i in range(n)]
        self.Nup = [self.N[3*i] for i in range(n)]
        self.Ndn = [self.N[3*i+1] for i in range(n)]
        self.NF = [self.N[3*i+2] for i in range(n)]










