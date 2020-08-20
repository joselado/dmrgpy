from .manybodychain import dmrgpath
import os

def run(self):
    if self.itensor_version in [2,"2","v2","C++","cpp","c","C"]:
        mpscpp = dmrgpath+"/mpscpp2/mpscpp.x"
    else: raise
    self.execute(lambda : os.system(mpscpp+" > status.txt"))

