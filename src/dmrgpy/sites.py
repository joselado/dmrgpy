# initialize the sites for the C++ calculation
import os

def initialize(self,**kwargs):
    # self.path is kept only as a label for a few legacy helpers
    # (Many_Body_Chain.clone()/to_folder()/to_origin()) inherited from the
    # old file-based DMRG backend -- the actual calculation is entirely
    # in-process now (see cppext.py/chain_session.h), so nothing writes to
    # or reads from this folder anymore, and it is no longer created on
    # disk.
    self.path = os.getcwd()+"/.mpsfolder/" # folder of the calculations
    self.inipath = os.getcwd() # original folder
    # build the in-process extension session (mpscpp2/mpscpp3's
    # chain_session.h Chain, or pyitensor.chain.Chain, depending on
    # itensor_version). If the extension isn't compiled, self._session
    # stays None and mode.py's get_mode() falls back to ED for this chain
    # -- there is no file-based DMRG backend left to fall back to. (The
    # "python" backend has no compiled-extension precondition, so this
    # never happens for it -- see cppext.py.)
    if self.itensor_version in (2,3,"python"):
        from . import cppext
        backend = cppext.get_backend(self.itensor_version)
        if backend is not None:
            self._session = backend.Chain(self.sites)
            # A conserved sector is a property of the chain, but it lives in
            # the session's site set (chain_session.h rebuilds sites_ with
            # QN-carrying indices), so a session built or rebuilt here --
            # setup_cpp()/setup_python() switching backend, clone() making a
            # fresh one -- has to be told about it again.
            # _apply_conserved_sector() is what decides what that means for
            # this backend: it refuses outright for one with no quantum
            # numbers at all (which would silently answer with the global
            # ground state), unless the chain says mode="ED", in which case
            # ED targets the sector by itself and the session is simply left
            # out of it.
            if getattr(self,"conserved_sector",None):
                self._apply_conserved_sector()


# The settings of a chain: the solver parameters a constructor keyword may
# set, every one of them read by a solver. An allowlist, where the rule used
# to be "any public attribute the chain has that is not a method, less
# STATE": that rule admitted every attribute STATE forgot, the Hamiltonian
# accumulators below among them, so Fermionic_Chain(4, hubbard=2.0) put the
# constant 2*Id into the Hamiltonian of the next set_hoppings() (-0.2360679775
# against the free-fermion -2.2360679775), and five attributes nothing reads
# were accepted and stored (2026-09-25b hole hunt, finding 4). A setting
# left out of this list fails loudly, with a TypeError naming it; a state
# name left out of STATE is refused all the same, with the generic message.
# itensor_version is not here: it is the constructor's own named argument.
SETTINGS = frozenset((
    "maxm","nsweeps","noise","cutoff","mpomaxm","verbose","mode",
    "bond_ramp","bond_ramp_start","bond_ramp_fraction",
    "bond_ramp_noise_decay",
    "kpmmaxm","kpmcutoff","kpm_scale","kpm_accelerate","kpm_n_scale",
    "kpm_extrapolate","kpm_extrapolate_factor","kpm_extrapolate_mode",
    "kpm_energy_truncate","kpm_truncate_dK","kpm_truncate_nsweeps",
    "kpm_truncate_threshold",
    "cvm_tol","cvm_nit","cvm_patience","cvm_blowup","cvm_solver",
    "cvm_nsweeps","cvm_maxm",
    "tevol_method","tevol_custom_exp","tdvp_gse_sweeps",
    "tdvp_gse_krylov_order","tdvp_gse_cutoff",
    "excited_gram_schmidt",
))


# Attributes a chain holds that a constructor keyword may not set, each with
# the reason the message gives: its state, each with its own entry point
# (the site list is the positional argument), the Hamiltonian accumulators
# update_hamiltonian() sums (a keyword would be summed in as value*Id), and
# the two flags the model class fixes after Many_Body_Chain.__init__ has run,
# which would overwrite the keyword. They are refused because they are not
# in SETTINGS; this table only names the way to set each.
STATE = {
    "sites":"the site list is the constructor's positional argument",
    "ns":"the site list is the constructor's positional argument",
    "Id":"the site list is the constructor's positional argument",
    "hamiltonian":"use set_hamiltonian()",
    "conserved_sector":"use set_conserved_sector()",
    "wf0":"use set_gs() or set_initial_wf()",
    "e0":"use set_gs() or set_initial_wf()",
    "computed_gs":"use set_gs() or set_initial_wf()",
    "gs_from_file":"use set_gs() or set_initial_wf()",
    "skip_dmrg_gs":"use set_gs() or set_initial_wf()",
    "sites_from_file":"use set_gs() or set_initial_wf()",
    "excited_from_file":"use set_gs() or set_initial_wf()",
    "has_ED_obj":"the ED object is built on demand",
    "ED_obj":"the ED object is built on demand",
    "fermionic":"the chain class sets it",
    "use_ampo_hamiltonian":"the chain class sets it",
    "hopping":"use set_hoppings() or set_hamiltonian()",
    "hubbard":"use set_hubbard() or set_hamiltonian()",
    "pairing":"use set_pairings_MB() or set_hamiltonian()",
    "exchange":"use set_hamiltonian()",
}


def check_settings(self,settings,model=()):
    """Refuse every constructor keyword that is not a setting of this chain,
    naming all of them, sorted; return the ones that are, as a dict.

    Many_Body_Chain.__init__(**kwargs) used to hand its keywords to
    initialize(), which ignores them, so Spin_Chain(sites, maxm=50) built a
    chain at the default maxm=30 and a misspelled keyword was accepted
    without a word, on every model chain, since they all forward to it. A
    setting is a name in SETTINGS above; anything else is refused, the
    state in STATE and `model`, the attributes the model class built before
    calling Many_Body_Chain.__init__ (Fermionic_Chain's N,
    Parafermionic_Chain's Sig, ...), with a message saying so, since a
    keyword would overwrite them with a number. The check reads nothing off
    `self` but its class name, which names the constructor in the message,
    so a wrapper such as Thermal_Spin_Chain runs it on itself before it
    builds its chain. Called before any session is built, so a bad keyword
    costs nothing; mode= is checked for its value here too, as sc.mode is
    when it is read."""
    def why(k):
        if k in STATE: return STATE[k]
        return "the %s class builds it" % type(self).__name__
    state = sorted(k for k in settings if k in STATE or k in model)
    if state:
        raise TypeError(
            "%s() cannot set %s at construction: %s is the chain's state, "
            "not a setting (%s)" % (type(self).__name__, ", ".join(state),
            "that" if len(state)==1 else "each",
            "; ".join("%s: %s" % (k,why(k)) for k in state)))
    unknown = sorted(k for k in settings if k not in SETTINGS)
    if unknown:
        raise TypeError(
            "%s() got unexpected keyword argument(s) %s. A constructor "
            "keyword names a setting of the chain (maxm, nsweeps, noise, "
            "cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes "
            "effect as if assigned right after construction; any other "
            "name would be stored where nothing reads it, or overwrite the "
            "chain's own state or behaviour"
            % (type(self).__name__, ", ".join(unknown)))
    if settings.get("mode") is not None:
        from .mode import _check_mode
        _check_mode(settings["mode"],"chain mode (mode=)")
    return dict(settings)
