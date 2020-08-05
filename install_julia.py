#!/usr/bin/python3
import os
import sys
import subprocess


mpath = os.path.dirname(os.path.realpath(__file__)) 


# main path
try:
    out,err = subprocess.Popen(['julia', '--version'],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT).communicate()
except: out = ""
# check if JUlia has the correct version
hasjulia = "julia version 1.4" in str(out)
#print(hasjulia) ; exit()

if not hasjulia: # if the correct Julia version is not present
    print("Julia not present in path, downloading")
    os.system("mkdir "+mpath+"/src/julia") # create a subfodler for julia
    os.chdir(mpath+"/src/julia") # go to the subfolder
# download julia
    os.system("wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz")
    os.system("tar -xvf julia-1.4.2-linux-x86_64.tar.gz") # untar the file
    os.system("rm julia-1.4.2-linux-x86_64.tar.gz") # rm the file
    julia = mpath+"/src/julia/julia-1.4.2/bin/julia" # path for julia
    os.chdir(mpath) # go back

else: 
    julia = "julia"
    print("Correct Julia version found in path")

# install Itensor julia

os.system(julia+" "+mpath+"/install/juliainstall.jl") # julia installation
# install python dependences
os.system("pip install julia")
os.system("pip install pmdarima")






##################
# and add dmrgpy to the python path

def get_lp():
    for p in sys.path:
        if "/site-packages" in p: return p
    print("Path to add the Python library not found")
    exit()


p = get_lp()+"/dmrgpy" # get the path
os.system("rm -f "+p) # remove the symbolic link
os.system("ln -s "+mpath+"/src/dmrgpy  "+p) # remove the symbolic link





