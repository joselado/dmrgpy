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




