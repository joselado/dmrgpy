import os
import re


def compile(gpp="g++",check_gpp=True,**kwargs):
    # compile Itensor library
    path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
    os.chdir(path+"/src/dmrgpy/mpscpp2/ITensor")
    import platform
    if platform.system()=="Linux" or platform.system()=="Darwin": # System
      os.system("cp options.save options.mk")
      print("Detected Unix system")
      from . import cppversion
      writemk(gpp=gpp,**kwargs) # write options.mk
      if check_gpp: # check the compiler
        print("Checking if the compilar is ok")
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
          exit()
#          import autoinstallubuntu
#          autoinstallubuntu.install() # install dependencies
#          # replace the compiler in options.mk
#          out = open("options.mk").read().replace("CCCOM=g++","CCCOM=g++-6")
#          open("options.mk","w").write(out) # write new file
    else:
      print("Not a Unix system. The code only works in Unix systems. Stopping")
      exit()
    # now do the compilation
    os.system("make clean ")
    os.system("make")
    os.chdir(path) # original directory

    # compile DMRG program
    os.chdir(path+"/src/dmrgpy/mpscpp2")
    os.system("make clean")
    os.system("make")
    os.system("mv mpscpp mpscpp.x")
    # build the in-process pybind11 extension alongside the old
    # subprocess-based executable (see the file-I/O-to-in-memory migration
    # plan); kept optional/best-effort for now since it's not yet used by
    # the Python side (that starts once mpscpp2/*.h are refactored)
    if cppversion.has_pybind11():
        os.system("make pybind")
    else:
        print("pybind11 not found, skipping the in-process extension build")
        print("Install it with: pip install pybind11")
    os.chdir(path) # go to the main path


def writemk(gpp="g++",openblas=False,openblas_libdir=None,openblas_includedir=None):
    """Write options.mk"""
    path = os.path.dirname(os.path.realpath(__file__))+"/../" # main path
    path = path+"/src/dmrgpy/mpscpp2/ITensor"
    out = open(path+"/options.save").read().replace("CCCOM=g++","CCCOM="+gpp)
    if openblas:
        # PLATFORM matters beyond linking: itensor/tensor/lapack_wrap.h
        # branches on -DPLATFORM_$(PLATFORM) for the LAPACK_COMPLEX layout
        # and several function signatures, so just swapping the link flags
        # (as before) risks a silent ABI mismatch against libopenblas.
        libflags = "-lpthread"
        if openblas_libdir:
            libflags += " -L"+openblas_libdir
        libflags += " -lopenblas"
        includeflags = "-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE"
        if openblas_includedir:
            includeflags = "-I"+openblas_includedir+" "+includeflags
        out = out.replace("PLATFORM=lapack","PLATFORM=openblas")
        out = re.sub(r"(?m)^BLAS_LAPACK_LIBFLAGS=.*$",
                "BLAS_LAPACK_LIBFLAGS="+libflags+
                "\nBLAS_LAPACK_INCLUDEFLAGS="+includeflags,
                out,count=1)
    open(path+"/options.mk","w").write(out) # write file


#import addsystem
#addsystem.addbashrc() # add to the .bashrc
