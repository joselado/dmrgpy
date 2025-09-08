from .manybodychain import dmrgpath
import os
import subprocess

def run(self):
    if self.itensor_version in [2,"2","v2","C++","cpp","c","C"]:
        mpscpp = dmrgpath+"/mpscpp2/mpscpp.x"
    else: raise
    def fexec():
        with open("status.txt", "w") as f:
            subprocess.run(mpscpp.split(), stdout=f)
    self.execute(fexec) # execute function

