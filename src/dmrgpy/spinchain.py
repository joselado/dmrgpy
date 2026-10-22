from .manybodychain import Many_Body_Chain
import numpy as np
from .algebra import algebra
from . import effectivehamiltonian
from . import pychainwrapper
from . import multioperator
from .kondospectrumtk.kondomethod import KondoSpectrumMixin

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


class Spin_Chain(KondoSpectrumMixin,Many_Body_Chain):
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
    # set_fields(fun) used to live here: it built sum_i b(i).S_i and then
    # assigned self.hamiltonian = self.exchange + self.fields directly,
    # bypassing set_hamiltonian(). Once set_exchange() was removed,
    # self.exchange was permanently the integer 0, so this silently
    # REPLACED whatever Hamiltonian the caller had built with the field
    # term alone rather than adding to it -- and, going around
    # set_hamiltonian(), never told the backend session about it either.
    # Write the field into the Hamiltonian directly instead
    # (h = h + b[j]*sc.Si[j][i], then set_hamiltonian(h)), which is what
    # its only caller (meanfield.py) now does.
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
            if not self.is_zero_operator(op,**kwargs): raise RuntimeError("Spin commutation relation [Sx,Sy] failed for sites %d %d" % (i,j))
    def get_logdimension(self):
        """Return the logarithm of the Hilbert space dimension"""
        return get_logdimension(self)
    # set_exchange(fun) used to live here: it built
    # sum_{i,j} fun(i,j) S_i.S_j -- over both orderings of every pair, so
    # a nearest-neighbour fun=1 gave 2*sum_<ij> S_i.S_j, not the
    # conventional Heisenberg sum -- stored it as self.exchange and set
    # self.hamiltonian = self.exchange + self.fields. It has been removed
    # in favour of writing the Hamiltonian out with SS(i,j) and
    # set_hamiltonian(), which is what every other model in this codebase
    # already does and what its callers now do (see
    # tests/test_benchmarks.py's all-pairs exchange test, which pins the
    # same reference energy the old builder produced).
    def get_sector_charge_operators(self):
        """A spin chain conserves total Sz, in ITensor's integer 2*Sz
        units (Sz=1 is one spin-1/2's worth), which is the unit
        set_conserved_sector takes on every backend"""
        return {"Sz":2*sum(self.Sz)}
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
        return effectivehamiltonian.get_effective_hamiltonian(self,**kwargs)
    # get_hamiltonian() used to live here, as an override with a
    # "conventional way" fallback that rebuilt the Hamiltonian from
    # self.exchange/self.fields when none had been set. Those two
    # attributes were only ever populated by set_exchange()/set_fields(),
    # both removed (see the comments above), so Many_Body_Chain.__init__
    # leaves them as the integer 0 and that branch could only ever die
    # with "TypeError: 'int' object is not iterable" -- the same corpse
    # 43d1a35 removed from meanfield.py, on any chain with no
    # set_hamiltonian() call, and inherited by the public
    # gs_energy_fluctuation(), which calls get_hamiltonian()
    # unconditionally. The live half (return self.hamiltonian) is exactly
    # what Many_Body_Chain.get_hamiltonian already does, so the override
    # is gone rather than reduced to a duplicate.

Spin_Hamiltonian = Spin_Chain # backwards compatibility
