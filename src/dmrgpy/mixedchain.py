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

    mode="ED" is exact diagonalization on the tensor product of the
    physical sites (edtk/mixededchain.py): a fermion mode is a two-state
    site, a spin S site has 2S+1 states, and C/Cdag carry the
    Jordan-Wigner string over the fermion modes to their left only, the
    same convention as the DMRG backends (F is the identity on a spin
    site). It is what itensor_version=3's automatic ns<3 fallback and
    "python"'s ns<2 fallback (mode.py) answer with. Checked in
    tests/test_mixed_chain_ed.py against Spin_Chain's and
    Spinful_Fermionic_Chain's own ED backends (all-spin and all-fermion
    sitesin, identical spectra), against the dense oracle
    tests/_dense_jw_reference.py, and against DMRG on both backends.

    Conserved sectors (set_conserved_sector): the total Sz (2*Sz units,
    spins and electrons together) and, when there is a fermion location,
    Nf. The DMRG site sets cannot carry them on a genuinely mixed chain:
    a spinless fermion mode carries Nf and no Sz, a spin site Sz and no
    Nf, so Nf or Sz alone is rejected by both backends and Nf+Sz together
    would grade only the spin sites' Sz, a different quantum number. A
    sector on a chain with both kinds is therefore ED only: it is refused
    unless the chain is pinned to mode="ED" (as the interleaved
    Spinful_Fermionic_Chain does for Sz). An all-spin sitesin conserves
    Sz on DMRG and an all-fermion one Nf, exactly as the pure chains.

    The Hamiltonian is built directly from
    Cup/Cdagup/Cdn/Cdagdn/Nup/Ndn/Ntot/Delta and Sx/Sy/Sz, as
    tests/test_mixed_chain.py, examples/fermion_models/mixed_spin_fermion_chain
    and examples/fermion_models/mixed_spin_fermion_chain_VS_ED do.
    """
    def __init__(self, sitesin, **kwargs):
        itensor_version = kwargs.get("itensor_version", None)
        if itensor_version is None:
            from .cppext import default_backend
            itensor_version = default_backend()
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
        """The ED object (edtk/mixededchain.py) on this chain's physical
        sites, built once and cached until the next restart()"""
        if self.has_ED_obj: return self.ED_obj
        from .edtk.mixededchain import MixedSpinFermionED
        ed = MixedSpinFermionED(self.sites)
        ed.hamiltonian = self.hamiltonian
        self._apply_sector_to_ed(ed) # conserved sector, if any
        self.ED_obj = ed
        self.has_ED_obj = True
        return self.ED_obj
    def get_sector_charge_operators(self):
        """Total Sz in 2*Sz units (2 Sz on a spin site, Nup-Ndn on a fermion
        location) and, when the chain has a fermion location, the electron
        number Nf. Both are diagonal in the ED product basis."""
        sz = 0*self.Id
        for loc in range(len(self.kind)):
            if self.kind[loc]=="spin": sz = sz + 2*self.Sz[loc]
            else: sz = sz + self.Nup[loc] + (-1)*self.Ndn[loc]
        out = {"Sz":sz}
        if "fermion" in self.kind: out["Nf"] = sum(self.Ntot)
        return out
    def _dmrg_quantum_numbers(self):
        """The quantum numbers the DMRG site set carries with the meaning
        of get_sector_charge_operators(): Sz on an all-spin chain, Nf on
        an all-fermion chain, none on a chain with both kinds (see the
        class docstring)"""
        if "fermion" not in self.kind: return {"Sz"}
        if "spin" not in self.kind: return {"Nf"}
        return set()
    def _ed_only_quantum_numbers(self):
        return set(self.get_sector_charge_operators()) - self._dmrg_quantum_numbers()
    def _apply_conserved_sector(self):
        """A sector naming a quantum number the DMRG sites cannot carry
        never reaches the session: on a chain pinned to mode="ED" the ED
        object targets it alone, otherwise it is refused here. Without
        this, Nf+Sz on a chain with both kinds would be accepted by the
        session with Sz graded on the spin sites only."""
        sector = self.conserved_sector or {}
        if not sector or not hasattr(self,"kind"): # nothing to decide yet
            return Many_Body_Chain._apply_conserved_sector(self)
        ed_only = set(sector) & self._ed_only_quantum_numbers()
        if len(ed_only)==0:
            return Many_Body_Chain._apply_conserved_sector(self)
        if self.mode!="ED":
            raise NotImplementedError(
                "conserved sector %s: the DMRG sites of this chain cannot "
                "carry %s (a fermion mode carries Nf and no Sz, a spin site "
                "Sz and no Nf), only its ED backend can. Pin the chain to "
                "ED with mode=\"ED\" before setting the sector. For the "
                "same reason the sector-resolved dynamical correlator "
                "(submode=\"SECTOR\", which sets this sector on an "
                "internal copy) is not available on a chain with both "
                "spin and fermion sites; submode=\"KPM\", \"TD\" or "
                "\"CVM\" work on it, with mode=\"DMRG\" or mode=\"ED\"."
                %(dict(sector),sorted(ed_only)))
        self._sector_on_session = False
