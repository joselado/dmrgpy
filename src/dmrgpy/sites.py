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


# Attributes a chain holds that a constructor keyword may not set: its
# state, each with its own entry point (the site list is the positional
# argument), and the two flags the model class fixes after
# Many_Body_Chain.__init__ has run, which would overwrite the keyword.
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
}


def check_settings(self,settings,model=()):
    """Refuse every constructor keyword that is not a setting of this chain,
    naming all of them, sorted; return the ones that are, as a dict.

    Many_Body_Chain.__init__(**kwargs) used to hand its keywords to
    initialize(), which ignores them, so Spin_Chain(sites, maxm=50) built a
    chain at the default maxm=30 and a misspelled keyword was accepted
    without a word, on every model chain, since they all forward to it. A
    setting is what the check in Infinite_Many_Body_Chain.kpm_finite
    accepts for window_chain_kwargs: a public attribute the chain already
    has that is not a method (a missing name would be stored where nothing
    reads it, and a private one or a method would overwrite the chain's own
    state or behaviour), less the state in STATE above and less `model`,
    the attributes the model class built before calling
    Many_Body_Chain.__init__ (Fermionic_Chain's N, Parafermionic_Chain's
    Sig, ...), which a keyword would overwrite with a number. Called
    before any session is built, so a bad keyword costs nothing; mode= is
    checked for its value here too, as sc.mode is when it is read."""
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
    unknown = sorted(k for k in settings
                     if k.startswith("_") or not hasattr(self,k)
                     or callable(getattr(self,k)))
    if unknown:
        raise TypeError(
            "%s() got unexpected keyword argument(s) %s. A constructor "
            "keyword names a setting of the chain (maxm, nsweeps, noise, "
            "cutoff, kpmmaxm, kpm_scale, tevol_method, mode, ...) and takes "
            "effect as if assigned right after construction; a name the "
            "chain does not have would be stored where nothing reads it"
            % (type(self).__name__, ", ".join(unknown)))
    if settings.get("mode") is not None:
        from .mode import _check_mode
        _check_mode(settings["mode"],"chain mode (mode=)")
    return dict(settings)
