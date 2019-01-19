#!/usr/bin/python


import os


pwd = os.getcwd()
for d in os.walk("."): # loop over subdirectories
  os.chdir(d[0])
  os.system("rm -rf .mpsfolder")
  os.system("rm -rf .pychainfolder")
  os.system("rm -rf .dmrgfolder")
  os.chdir(pwd)

