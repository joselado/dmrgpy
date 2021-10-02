import os


def compile(gpp="g++"):
    # compile Itensor library
    path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
    os.chdir(path+"/src/dmrgpy/mpscpp2/ITensor")
    import platform
    if platform.system()=="Linux" or platform.system()=="Darwin": # System
      os.system("cp options.save options.mk")
      print("Detected Unix system")
      from . import cppversion
      writemk(gpp=gpp) # write options.mk
      if cppversion.correct_version(gpp=gpp):
          print("C++ compiler "+gpp+" is ok and will be used")
      else: # C++ compiler is not the right one
          print("You may need to install the GNU C++ compiler")
          print("In Ubuntu systems, run the following commands")
          print("#######################################################")
          print("sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y")
          print("sudo apt-get update -y")
          print("sudo sudo apt-get install g++-6 -y")
          print("sudo apt-get install liblapack-dev -y")
          print("#######################################################")
          print("and afterwards rerun the installation script with")
          print("python install.py --cpp=g++-6")
#          import autoinstallubuntu
#          autoinstallubuntu.install() # install dependencies
#          # replace the compiler in options.mk
#          out = open("options.mk").read().replace("CCCOM=g++","CCCOM=g++-6")
#          open("options.mk","w").write(out) # write new file
    else:
      print("Not a Unix system. The code only works in Unix systems. Stopping")
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


def writemk(gpp="g++"):
    """Write options.mk"""
    path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
    path = path+"/src/dmrgpy/mpscpp2/ITensor"
    out = open(path+"/options.save").read().replace("CCCOM=g++","CCCOM="+gpp)
    open(path+"/options.mk","w").write(out) # write file


#import addsystem
#addsystem.addbashrc() # add to the .bashrc
