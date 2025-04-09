#!/usr/bin/python3

# create the optional arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--gpp",default="g++",
        help='C++ Compiler to be used')
parser.add_argument("--check_gpp",default="True",
        help='C++ Compiler to be used')
parser.add_argument("--openblas",default="False",
        help='Use openblas, instead of usual lapack')
args = parser.parse_args() # get the arguments



from installtk import install2
check_gpp = "True"==args.check_gpp # if gpp is checked
openblas = "True"==args.openblas # if openblas should be used
install2.compile(gpp=args.gpp,check_gpp=check_gpp) # compile the C++ library
from installtk import addpythonpath
addpythonpath.addpath() # add the library to the Python path
from installtk import addsystem
addsystem.addrc() # add the library to the bash file
