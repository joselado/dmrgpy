import numpy as np
import scipy.linalg as lg
from .pyfermion import mbfermion
from .algebra import algebra
from .manybodychain import Many_Body_Chain
from .pyboson import boson


# The only boson site mpscpp2 (ITensor v2) knows is ITensor's own
# BosonFourSite, spelled as the site-type code 104 (mpscpp2/get_sites.h).
# Anything else falls through to that header's Error(), and ITensor's
# Error() calls abort(), so a non-default maxnb on that backend used to
# take the whole interpreter down with SIGABRT before any Python code
# could see it. The check therefore has to happen here, before the
# session is built -- nothing thrown from C++ is catchable.
V2_BOSON_DIM = 4 # the single local dimension itensor_version=2 supports


def _check_backend_supports_maxnb(version,maxnb):
    """Raise before a session is built for a (backend, local dimension)
    combination the backend cannot represent. maxnb entries that are None
    belong to non-bosonic sites and are skipped."""
    dims = [m for m in maxnb if m is not None]
    if version!=2: return # every other backend builds an arbitrary dimension
    bad = sorted(set(m for m in dims if m!=V2_BOSON_DIM))
    if len(bad)==0: return
    raise ValueError(
        "itensor_version=2 only implements the fixed %d-level boson site "
        "(ITensor's BosonFourSite), but this chain asks for local boson "
        "dimension(s) %s. Use itensor_version=3, itensor_version=\"python\" "
        "or mode=\"ED\", which all build a boson site of the requested "
        "dimension."%(V2_BOSON_DIM,str(bad)))


class Bosonic_Chain(Many_Body_Chain):
    """Bosonic Hamiltonian"""
    def __init__(self,n,maxnb=None,**kwargs):
        if maxnb is None: maxnb = [4 for i in range(n)] # maximum # of bosons
        elif len(maxnb)!=n:
            raise ValueError("maxnb has length %d, but n=%d sites were requested"%(len(maxnb),n))
        self.maxnb = maxnb # local Hilbert space dimension of each site
        # Boson type codes are 100+dim (see mpscpp3/get_sites.h), so the
        # DMRG session actually gets a site of the requested dimension
        # instead of always the fixed 4-level BosonFourSite. This works on
        # itensor_version=3, on "python" (pyitensor/sites/boson.py builds
        # an arbitrary-dimension site) and under mode="ED"; mpscpp2 and
        # the Julia backend only understand the plain 104 code, and v2 is
        # refused outright in initialize() below rather than aborting the
        # process inside ITensor. See docs/user_guide.md.
        Many_Body_Chain.__init__(self,[100+m for m in maxnb],**kwargs)
        self.use_ampo_hamiltonian = True # use ampo
        self.N = [self.get_operator("N",i) for i in range(self.ns)]
        self.A = [self.get_operator("A",i) for i in range(self.ns)]
        self.Adag = [self.get_operator("Adag",i) for i in range(self.ns)]
        # occupation-number projectors N0..N{maxnb[i]-1}, per site
        self.D = [[self.get_operator("N"+str(k),i) for k in range(self.maxnb[i])]
                   for i in range(self.ns)]
        # D0..D3 convenience aliases, kept for backwards compatibility,
        # only defined when every site actually has >=4 levels
        if all(m>=4 for m in self.maxnb):
            self.D0 = [self.D[i][0] for i in range(self.ns)]
            self.D1 = [self.D[i][1] for i in range(self.ns)]
            self.D2 = [self.D[i][2] for i in range(self.ns)]
            self.D3 = [self.D[i][3] for i in range(self.ns)]
    def initialize(self,**kwargs):
        """Build the backend session, refusing a (backend, local
        dimension) pair the backend cannot represent first. This hook
        covers both routes onto such a pair: itensor_version=2 passed to
        the constructor (Many_Body_Chain.__init__ calls initialize) and
        setup_cpp(version=2) on an existing chain (which does too)."""
        if self.mode=="ED": return # no session is built at all, nothing to check
        _check_backend_supports_maxnb(self.itensor_version,self.maxnb)
        return Many_Body_Chain.initialize(self,**kwargs)
    def get_density(self,**kwargs):
        """Return the average boson occupation in each site"""
        out = [self.vev(self.N[i],**kwargs) for i in range(self.ns)]
        return np.array(out).real
    def get_density_fluctuation(self,**kwargs):
        """Return the occupation-number fluctuations <N^2>-<N>^2 in each site"""
        d = self.get_density(**kwargs) # get the density
        d2 = np.array([self.vev(self.N[i]*self.N[i],**kwargs) for i in range(self.ns)])
        return d2.real-d**2
    def get_sector_charge_operators(self):
        """A bosonic chain conserves the total boson number"""
        return {"Nb":sum(self.N)}
    def get_ED_obj(self):
        """Return the associated ED object"""
        # Same has_ED_obj/ED_obj caching protocol every other chain class
        # uses (see Fermionic_Chain.get_ED_obj): restart() clears the
        # flag, and set_hamiltonian()/set_conserved_sector() both call
        # restart(), so the cached object can never outlive the
        # Hamiltonian it was built for. Without it every single
        # vev(...,mode="ED") rebuilt the sparse operator dictionary and
        # redid the full ground-state solve.
        if self.has_ED_obj: # if the ED object has been computed
            return self.ED_obj # return the stored object
        else:
            dim = np.exp(np.sum(np.log(self.maxnb))) # kept in floating
                # point: the product overflows int64 on a long chain
            if dim>10000:
                raise ValueError("this bosonic chain has Hilbert space "
                        "dimension %g, too large for the ED backend"%dim)
            out = boson.bosonchain(self.maxnb)
            out.hamiltonian = self.hamiltonian
            self._apply_sector_to_ed(out) # conserved sector, if any
            self.ED_obj = out # store object
            self.has_ED_obj = True # set to True
            return self.ED_obj


from .spinchain import get_site as get_site_spin
from .spinchain import label2site as spin_label2site

# boson site labels accepted by SpinBoson_Chain. "B" is the historical
# spelling and means the default 4-level site; "B<k>" names the local
# dimension explicitly (so "B4" == "B", and "B6" is a 6-level site).
BOSON_LABELS = "\"B\" (4 levels) or \"B<k>\" for a k-level site, e.g. B4, B6"


def get_site(label,maxnb=None):
    """Site-type code for one label of a SpinBoson_Chain. maxnb, when
    given, is the local dimension asked for at this position and applies
    only to bosonic labels (a spin label carries its own dimension)."""
    out = get_site_spin(label) # get the spin site
    if out is None: # this is not a spin, try a boson
        dim = None
        if label=="B": dim = 4 # historical spelling, the default site
        elif isinstance(label,str) and label.startswith("B") and label[1:].isdigit():
            dim = int(label[1:]) # "B4", "B6", ... names the dimension
        if dim is None:
            raise ValueError("unknown site label "+repr(label)+"; this "
                    "chain accepts the spin labels "
                    +str(sorted(k for k in spin_label2site if isinstance(k,str)))
                    +" and the boson labels "+BOSON_LABELS)
        if maxnb is not None: dim = maxnb # an explicit maxnb wins
        if dim<2:
            raise ValueError("a boson site needs at least 2 levels, got "
                    +repr(dim))
        return 100+dim # boson type code, see mpscpp3/get_sites.h
    if maxnb is not None:
        raise ValueError("maxnb was given for site label "+repr(label)
                +", which is a spin site and carries its own dimension; "
                "pass None at that position")
    return out # otherwise


def is_boson(n):
    """Check if a certain site is a boson"""
    if 1<n<10: return False # nope, this is a spin
    elif 101<n<200: return True # yes
    else:
        raise ValueError("unrecognized site type code "+repr(n))

class SpinBoson_Chain(Many_Body_Chain):
    """Bosonic Hamiltonian"""
    def __init__(self,sitesin,n=None,maxnb=None,**kwargs):
        if n is not None and n!=len(sitesin):
            # n is redundant -- sitesin already gives the length -- and it
            # used to be accepted and dropped on the floor
            raise ValueError("n=%d contradicts the %d site labels given"
                    %(n,len(sitesin)))
        if maxnb is None: maxnb = [None for s in sitesin]
        elif len(maxnb)!=len(sitesin):
            raise ValueError("maxnb has length %d, but %d site labels were "
                    "given; pass one entry per site (None at the spin "
                    "sites)"%(len(maxnb),len(sitesin)))
        sites = [get_site(sitesin[i],maxnb[i]) for i in range(len(sitesin))]
        # local boson dimension per site, None at the spin sites -- this
        # is what initialize()'s backend check and get_ED_obj() read, and
        # it is set before Many_Body_Chain.__init__ because that calls
        # initialize()
        self.maxnb = [sites[i]-100 if is_boson(sites[i]) else None
                      for i in range(len(sites))]
        Many_Body_Chain.__init__(self,sites,**kwargs) # initialize
        self.use_ampo_hamiltonian = True # use ampo
        self.N = [self.get_operator("N",i) for i in range(self.ns)]
        self.A = [self.get_operator("A",i) for i in range(self.ns)]
        self.Adag = [self.get_operator("Adag",i) for i in range(self.ns)]
        self.Sx = [self.get_operator("Sx",i) for i in range(self.ns)]
        self.Sy = [self.get_operator("Sy",i) for i in range(self.ns)]
        self.Sz = [self.get_operator("Sz",i) for i in range(self.ns)]
        # occupation-number projectors per site, empty at the spin sites
        self.D = [[self.get_operator("N"+str(k),i)
                   for k in range(self.maxnb[i])] if self.maxnb[i] is not None
                  else [] for i in range(self.ns)]
        # now depurate the operators
        for i in range(self.ns): # loop over sites
            if is_boson(sites[i]): # for bosonic sites
                self.Sx[i] = 0 # set to zero
                self.Sy[i] = 0 # set to zero
                self.Sz[i] = 0 # set to zero
            else: # for spin sites
                self.N[i] = 0 # set to zero
                self.A[i] = 0 # set to zero
                self.Adag[i] = 0 # set to zero
        # D0..D3 convenience aliases, kept for backwards compatibility.
        # 0 at a spin site, as the other per-site operators are, and only
        # defined when every boson site actually has >=4 levels.
        if all(m is None or m>=4 for m in self.maxnb):
            for k in range(4):
                setattr(self,"D"+str(k),
                        [self.D[i][k] if self.maxnb[i] is not None else 0
                         for i in range(self.ns)])
    def initialize(self,**kwargs):
        """See Bosonic_Chain.initialize -- same backend check, same reason."""
        if self.mode=="ED": return # no session is built at all
        _check_backend_supports_maxnb(self.itensor_version,self.maxnb)
        return Many_Body_Chain.initialize(self,**kwargs)
    def get_ED_obj(self):
        """Return the associated ED object"""
        # same caching protocol as Bosonic_Chain.get_ED_obj above
        if self.has_ED_obj: # if the ED object has been computed
            return self.ED_obj # return the stored object
        else:
            out = boson.SpinBosonChain(self.sites) # not implemented yet
            out.hamiltonian = self.hamiltonian
            self._apply_sector_to_ed(out) # conserved sector, if any
            self.ED_obj = out # store object
            self.has_ED_obj = True # set to True
            return self.ED_obj
