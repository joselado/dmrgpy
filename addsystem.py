#!/usr/bin/python
from __future__ import print_function
import os


pwd = os.getcwd() # get the current location

# different files for Linux and Mac
import platform
if platform.system()=="Linux":
  bashrc = os.environ["HOME"]+"/.bashrc" # path to .bashrc
  print("Detected Linux system")
else:
  bashrc = os.environ["HOME"]+"/.bash_profile" # path to .bashrc
  print("Detected Mac system")
print("Adding DMRGROOR to your ",bashrc)

route = "\n###############################\n"
route += "# Added by dmrgpy\n"
route += "###############################\n"
route += "  export DMRGROOT=\""+pwd+"/src\"\n"
route += "###############################\n"


f = open(bashrc).read() # read bashrc
open(bashrc,"w").write(f+route) # write in the bashrc

print("Added \n"+route+"\n route to ",bashrc)




