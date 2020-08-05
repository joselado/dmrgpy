#!/usr/bin/python3

# compile C++
import os
os.system("python install/install2.py") # execute installation script


# and add dmrpy to the python path
##############################################
import os
import sys
import subprocess


mpath = os.path.dirname(os.path.realpath(__file__))


def get_lp():
    for p in sys.path:
        if "/site-packages" in p: return p
    print("Path to add the Python library not found")
    exit()


p = get_lp()+"/dmrgpy" # get the path
os.system("rm -f "+p) # remove the symbolic link
os.system("ln -s "+mpath+"/src/dmrgpy  "+p) # remove the symbolic link


# install two required python libraries
#os.system("pip install julia")
os.system("pip install pmdarima")



