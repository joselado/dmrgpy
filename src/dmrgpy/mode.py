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
    # fall back to ED (Julia chains are unaffected by this check)
    if self.itensor_version==2:
        from . import cppext
        if not cppext.available():
            print("C++ extension not compiled, using default ED routines")
            return "ED" # use exact diagonalization
    # if there is an enforced mode, then use that one
    if self.mode is not None: return self.mode # use the enforced mode
    else:
        if mode in ["ED","DMRG"]: return mode # use default
        else:
            print("Unrecognized mode",mode)
            raise


