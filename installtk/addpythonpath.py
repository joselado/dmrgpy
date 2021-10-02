import os
import sys

def addpath():
    mpath = os.path.dirname(os.path.realpath(__file__))
    def get_lp():
        for p in sys.path:
            if "/site-packages" in p: return p
        print("Path to add the Python library not found")
        exit()
    p = get_lp()+"/dmrgpy" # get the path
    os.system("rm -f "+p) # remove the symbolic link
    os.system("ln -s "+mpath+"/../src/dmrgpy  "+p) # remove the symbolic link

