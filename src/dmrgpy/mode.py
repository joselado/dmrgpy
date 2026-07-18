import os

def run(self,automatic=False):
    """
    Run the DMRG calculation (Julia backend only -- the C++ backend is
    handled in-process via self._session, see chain_session.h, and never
    reaches this function)
    """
    if get_mode(self)=="ED":
        return # do nothing
    if self.itensor_version in ["julia","Julia","jl"]:
        from . import juliarun
        juliarun.run(self)
        return
    else:
        print("Unrecognized mode",get_mode(self))
        raise



def get_mode(self,mode="DMRG"):
    """Return the mode of the calculation"""
    # if the in-process C++ extension isn't available for a C++ chain,
    # fall back to ED (Julia chains are unaffected by this check). The
    # "python" backend always passes this check (cppext.available("python")
    # is unconditionally True -- no compiled extension to be missing).
    if self.itensor_version in (2,3,"python"):
        from . import cppext
        if not cppext.available(self.itensor_version):
            print("C++ extension not compiled, using default ED routines")
            return "ED" # use exact diagonalization
    # ITensor v3's dmrg() always does two-site updates (see chain_session.h's
    # dmrg_args()) and its sweep loop aborts the whole process (SIGABRT,
    # "LocalOp is default constructed" deep inside ITensor's own
    # DMRGWorker/davidson sweep loop, mps/localop.h) rather than raising a
    # catchable Python exception, for chains with too few sites to have a
    # well-defined sweep range. Confirmed empirically: n=1 and n=2 both
    # abort (regardless of model or Hamiltonian -- tried spin, fermionic,
    # diagonal-only and hopping terms), n=3 is fine. That sweep loop is
    # vendored ITensor code (mpscpp3/ITensor/), out of scope to patch here,
    # so fall back to ED for these unsupported sizes, same as the
    # extension-not-compiled fallback above. v2 and the pure-Python backend
    # both handle n=1/n=2 correctly and are unaffected.
    if self.itensor_version==3 and self.ns<3:
        print("ITensor v3's two-site DMRG can't handle a chain this short "
              "(n=%d < 3 sites), using default ED routines"%self.ns)
        return "ED" # use exact diagonalization
    # if there is an enforced mode, then use that one
    if self.mode is not None: return self.mode # use the enforced mode
    else:
        if mode in ["ED","DMRG"]: return mode # use default
        else:
            print("Unrecognized mode",mode)
            raise


