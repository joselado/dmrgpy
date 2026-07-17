# initialize the sites for the C++ calculation
import subprocess
import os

def initialize(self,**kwargs):
    self.path = os.getcwd()+"/.mpsfolder/" # folder of the calculations
    self.clean() # clean calculation
    self.inipath = os.getcwd() # original folder
    subprocess.run(["mkdir","-p",self.path]) # create folder
    # build the in-process extension session (mpscpp2/mpscpp3's
    # chain_session.h Chain, depending on itensor_version). If the
    # extension isn't compiled, self._session stays None and mode.py's
    # get_mode() falls back to ED for this chain -- there is no file-based
    # DMRG backend left to fall back to.
    if self.itensor_version in (2,3):
        from . import cppext
        backend = cppext.get_backend(self.itensor_version)
        if backend is not None:
            self._session = backend.Chain(self.sites)




