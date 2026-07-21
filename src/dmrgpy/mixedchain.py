from .manybodychain import Many_Body_Chain
from .spinchain import label2site

_fermion_labels = ("F", "f", "fermion", "Fermion")


def is_fermion_label(label):
    """Check if a mixed-chain site label denotes a spinful-fermion
    location (as opposed to a spin location)"""
    return label in _fermion_labels


def get_site_label(label):
    """Classify a mixed-chain site label as "spin" or "fermion" """
    if is_fermion_label(label): return "fermion"
    elif label in label2site: return "spin"
    else: raise ValueError("Unrecognized mixed-chain site label "+str(label))


class Mixed_Spin_Fermion_Chain(Many_Body_Chain):
    """
    Chain mixing genuine spin sites and spinful-fermion locations in the
    same chain (e.g. for Kondo-lattice-like models coupling a literal
    local moment to a conduction-electron site, instead of the
    large-U two-fermion-site trick used by
    fermionchain.Spinful_Fermionic_Chain/spinfermionchain.py).

    `sitesin` is a list with one entry per logical location, each either
    a spin label recognized by spinchain.label2site ("1/2","1","3/2",...)
    or a fermion marker ("F"/"f"/"fermion"). A fermion location is built,
    exactly like fermionchain.Spinful_Fermionic_Chain, as a pair of
    physical spinless-fermion (type-code 0) sites for the up/down
    channels; a spin label is a single physical site of the
    corresponding spin type. This reuses the existing C/Cdag/A/Adag/F
    operator vocabulary and multioperatortk/jordanwigner.py unchanged.

    Operator lists are indexed by *logical* location (length
    len(sitesin)), not by physical site (self.ns, which counts each
    fermion location twice) -- the same convention
    Spinful_Fermionic_Chain already uses for its own Sx/Sy/Sz/Cup/Cdn/...

    Sx/Sy/Sz are defined at every location: at a spin location they are
    the native spin operators, at a fermion location they are the usual
    spin-density bilinear (0.5*Cdagup*Cdn + h.c., etc, same formula as
    Spinful_Fermionic_Chain). The fermionic operators (Cup/Cdagup/Cdn/
    Cdagdn/Nup/Ndn/Ntot/Delta) are only meaningful at fermion locations
    and are set to the literal int 0 at spin locations (same masking
    pattern bosonchain.py::SpinBoson_Chain uses for its boson/spin mix).

    For itensor_version=3, a Jordan-Wigner string that has to cross a
    spin location (e.g. a hopping term between two fermion locations
    with a spin location in between) is handled correctly because
    ITensor v3's SiteSet::op() (mpscpp3/ITensor/itensor/mps/siteset.h)
    resolves an unrecognized "F" request -- which is what a spin site
    receives, since it defines no "F" operator -- to the identity,
    rather than erroring; see tests/test_mixed_chain.py for a
    regression test exercising this directly. itensor_version=2 has no
    such fallback and hard-aborts the whole process on the same input,
    so it is rejected up front in __init__; this class only supports
    itensor_version in (3,"python").

    There is no ED backend for this chain type yet (get_ED_obj() raises
    NotImplementedError): building one needs to combine genuine
    fermionic Fock-space statistics with tensor-product spin operators,
    which no existing EDchain subclass does today (see
    pyboson.boson.SpinBosonChain, the closest analogue for spin+boson,
    which is itself an unfinished stub). Because of this,
    itensor_version=3's automatic ns<3 -> ED fallback (mode.py, worked
    around for every other chain type) will raise that same
    NotImplementedError for a mixed chain with fewer than 3 physical
    sites -- use a larger chain or itensor_version="python" instead.

    The generic Many_Body_Chain._MB Hamiltonian-building helpers
    (set_hoppings_MB/set_hubbard_MB/set_pairings_MB) are not usable
    here: they reference self.C/self.Cdag/self.N, which this class does
    not define (a single-species "C" is not well-defined at a spinful
    fermion location) -- build the Hamiltonian directly from
    Cup/Cdagup/Cdn/Cdagdn/Nup/Ndn/Ntot/Delta and Sx/Sy/Sz instead, as
    tests/test_mixed_chain.py and examples/mixed_spin_fermion_chain do.
    """
    def __init__(self, sitesin, **kwargs):
        itensor_version = kwargs.get("itensor_version", None)
        if itensor_version is None:
            from .cppext import DEFAULT_ITENSOR_VERSION
            itensor_version = DEFAULT_ITENSOR_VERSION
        if itensor_version not in (3, "python"):
            raise ValueError("Mixed_Spin_Fermion_Chain only supports "
                    "itensor_version in (3,'python'), got "
                    +str(itensor_version)+" -- itensor_version=2 has no "
                    "fallback for a Jordan-Wigner string crossing a spin "
                    "site and hard-aborts the process (see this class's "
                    "docstring)")
        kinds = [get_site_label(s) for s in sitesin]
        phys_types = [] # flat physical type-code list
        phys_index = [] # per logical location: int (spin) or (up,dn) tuple (fermion)
        for label, kind in zip(sitesin, kinds):
            if kind=="fermion":
                up = len(phys_types)
                phys_types.append(0)
                phys_types.append(0)
                phys_index.append((up, up+1))
            else: # spin
                i = len(phys_types)
                phys_types.append(label2site[label])
                phys_index.append(i)
        Many_Body_Chain.__init__(self, phys_types, **kwargs)
        self.kind = kinds # "spin" or "fermion" per logical location
        self.phys_index = phys_index
        self.use_ampo_hamiltonian = True # use ampo
        n = len(sitesin)
        # per-physical-site fermionic operators (only meaningful at the
        # physical sites belonging to a fermion location)
        Cp = [self.get_operator("C",i) for i in range(self.ns)]
        Cdagp = [self.get_operator("Cdag",i) for i in range(self.ns)]
        Np = [self.get_operator("N",i) for i in range(self.ns)]
        # Masked-out entries are real (numerically zero) MultiOperator
        # instances, not the literal int 0: a plain 0 would degrade to a
        # bare Python int the moment two masked entries at spin locations
        # get multiplied together (0*0==0, an int with none of
        # MultiOperator's API, e.g. no .get_dagger()), and would also
        # break vev() on a masked entry, since MultiOperator.obj2MO's
        # bare-number branch does not use the target vev name (see
        # multioperator.py). "0*self.Id" builds a fresh, independent,
        # properly-typed zero operator per list entry.
        self.Sx = [0 for i in range(n)]
        self.Sy = [0 for i in range(n)]
        self.Sz = [0 for i in range(n)]
        self.Cup = [0*self.Id for i in range(n)]
        self.Cdagup = [0*self.Id for i in range(n)]
        self.Cdn = [0*self.Id for i in range(n)]
        self.Cdagdn = [0*self.Id for i in range(n)]
        self.Nup = [0*self.Id for i in range(n)]
        self.Ndn = [0*self.Id for i in range(n)]
        self.Ntot = [0*self.Id for i in range(n)]
        self.Delta = [0*self.Id for i in range(n)]
        for loc in range(n):
            if kinds[loc]=="spin":
                i = phys_index[loc]
                self.Sx[loc] = self.get_operator("Sx",i)
                self.Sy[loc] = self.get_operator("Sy",i)
                self.Sz[loc] = self.get_operator("Sz",i)
            else: # fermion
                up,dn = phys_index[loc]
                self.Cup[loc] = Cp[up]
                self.Cdagup[loc] = Cdagp[up]
                self.Cdn[loc] = Cp[dn]
                self.Cdagdn[loc] = Cdagp[dn]
                self.Nup[loc] = Np[up]
                self.Ndn[loc] = Np[dn]
                self.Ntot[loc] = Np[up] + Np[dn]
                self.Delta[loc] = 0.5*Cp[up]*Cp[dn]
                self.Sx[loc] = 0.5*Cdagp[up]*Cp[dn] + 0.5*Cdagp[dn]*Cp[up]
                self.Sy[loc] = -0.5*1j*Cdagp[up]*Cp[dn] + 1j*0.5*Cdagp[dn]*Cp[up]
                self.Sz[loc] = 0.5*Np[up] - 0.5*Np[dn]
        self.Si = [self.Sx, self.Sy, self.Sz]
    def SS(self,i,j):
        """Spin-spin coupling between two locations (either kind)"""
        return self.Sx[i]*self.Sx[j] + self.Sy[i]*self.Sy[j] + self.Sz[i]*self.Sz[j]
    def setup_cpp(self,version=3):
        """Switch the C++ DMRG backend -- only version=3 is supported
        here (see this class's docstring)."""
        if version!=3:
            raise ValueError("Mixed_Spin_Fermion_Chain only supports "
                    "itensor_version=3, got "+str(version))
        Many_Body_Chain.setup_cpp(self,version)
    def get_ED_obj(self):
        """
        There is no ED backend for a mixed spin/spinful-fermion Hilbert
        space yet (see this class's docstring): unlike the pure-spin or
        pure-fermion chains, that needs a chain-of-tensor-product spin
        operators combined with genuine Fock-space fermion statistics,
        which no existing EDchain subclass implements. Raising here
        (instead of leaving get_ED_obj undefined) turns mode.py's
        automatic DMRG->ED fallback for itensor_version==3 chains with
        self.ns<3 into a clear error instead of an unrelated
        AttributeError.
        """
        raise NotImplementedError("Mixed_Spin_Fermion_Chain has no ED "
                "backend yet; use itensor_version=3 or \"python\" DMRG "
                "on a chain with self.ns>=3 physical sites instead of "
                "mode=\"ED\"")
