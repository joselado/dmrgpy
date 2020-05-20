import os



# compile Itensor library
path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
os.chdir(path+"/src/dmrgpy/mpscpp2/ITensor")
import platform
if platform.system()=="Linux": # Linux system
  os.system("cp options.save options.mk")
  print("Detected Linux system")
  import cppversion
  if cppversion.correct_version():
      print("C++ compiler is OK")
  else: # C++ compiler is not the right one
      import autoinstallubuntu
      autoinstallubuntu.install() # install dependencies
      # replace the compiler in options.mk
      out = open("options.mk").read().replace("CCCOM=g++","CCCOM=g++-6")
      open("options.mk","w").write(out) # write new file
else:
  print("Detected Mac system")
  print("You may need to install the GNU C++ compiler")
  os.system("cp options.mac options.mk")
# now do the compilation
os.system("make clean ")
os.system("make")
os.chdir(path) # original directory

# compile DMRG program
os.chdir(path+"/src/dmrgpy/mpscpp2")
os.system("make clean")
os.system("make")
os.system("mv mpscpp mpscpp.x")
os.chdir(path) # go to the main path

import addsystem
addsystem.addbashrc() # add to the .bashrc
