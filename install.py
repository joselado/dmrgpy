#!/usr/bin/python3

# create the optional arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--gpp",default="g++",
        help='C++ Compiler to be used')
args = parser.parse_args() # get the arguments



from installtk import install2
install2.compile(gpp=args.gpp) # compile the C++ library
from installtk import addpythonpath
addpythonpath.addpath() # add the library to the Python path
from installtk import addsystem
addsystem.addrc() # add the library to the bash file
