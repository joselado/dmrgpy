#!/usr/bin/python3
import os
import sys


mpath = os.path.dirname(os.path.realpath(__file__))


# install python dependences
os.system("pip install juliacall") # bridge to the live Julia session
                                    # (mpsjulialive/), replacing the old
                                    # PyJulia ("pip install julia") bridge
os.system("pip install pmdarima")


# Trigger juliacall/juliapkg to provision Julia (via juliaup, reusing any
# existing Julia install juliapkg finds; no manual download/tarball dance
# needed here any more) and resolve+precompile the packages declared in
# src/dmrgpy/juliapkg.json (ITensors, ITensorMPS, ITensorNHDMRG) into its
# own managed project. This can take several minutes the first time
# (registry update + precompilation); later imports reuse the same
# managed project and are fast.
print("Setting up the Julia backend (juliacall) -- this may take a few "
      "minutes the first time")
sys.path.insert(0,mpath+"/src")
try:
    import dmrgpy.mpsjulialive.juliasession
    print("Julia backend ready")
except Exception as e:
    print("Could not set up the Julia backend:",e)
    print("setup_julia() will not be usable until this is resolved")


##################
# and add dmrgpy to the python path

def get_lp():
    for p in sys.path:
        if "/site-packages" in p: return p
    print("Path to add the Python library not found")
    exit()


p = get_lp()+"/dmrgpy" # get the path
os.system("rm -f "+p) # remove the symbolic link
os.system("ln -s "+mpath+"/src/dmrgpy  "+p) # remove the symbolic link





