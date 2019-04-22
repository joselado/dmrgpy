import os

# compile Itensor library
path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
os.chdir(path+"/src/dmrgpy/mpscpp/ITensor-master")
os.system("cp options.save options.mk")
import cppversion
if cppversion.correct_version():
    print("C++ compiler is fine")
else: # C++ compiler is not the right one
    import autoinstallubuntu
    autoinstallubuntu.install() # install dependencies
    # replace the compiler in options.mk
    out = open("options.mk").read().replace("CCCOM=g++","CCCOM=g++-6")
    open("options.mk","w").write(out) # write new file
# now do the compilation
os.system("make clean ")
os.system("make")
os.chdir(path) # original directory

# compile DMRG program
os.chdir(path+"/src/dmrgpy/mpscpp")
os.system("make clean")
os.system("make")
os.system("mv mpscpp mpscpp.x")
os.chdir(path) # go to the main path

import addsystem
addsystem.addbashrc() # add to the .bashrc
