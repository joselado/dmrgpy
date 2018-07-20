#!/usr/bin/python


import os


pwd = os.getcwd()
for d in os.walk("."): # loop over subdirectories
  os.chdir(d[0])
  os.system("rm -f *.o")
  os.system("rm -f *.so")
  os.system("rm -rf __pycache__")
  os.chdir(pwd)

