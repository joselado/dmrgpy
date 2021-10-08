# routines to run the code with Julia
import os
import subprocess

dmrgpath = os.path.dirname(os.path.realpath(__file__))




# check the system image
sysimage = os.environ["HOME"]+"/.julia/sysimages/sys_itensors.so"
if not os.path.isfile(sysimage):
    print("No ITensors system image found, this may take some time")
    sysimage = None

try: # create the executable
    from julia.api import Julia
    jlsession = Julia(compiled_modules=False,
            sysimage=sysimage) # start the Julia session
    jlsession.eval("using Suppressor") # suppress output
except:
    print("Julia cannot be executed, you probably need to install the Julia version")
    class JLdummy(): 
        def eval(self,c):
            print("Julia cannot be executed, you probably need to install the Julia version")
            exit()
    jlsession = JLdummy()


def run(self):
    """Execute the Julia program"""
    import contextlib
    c = "@suppress_out include(\""+dmrgpath+"/mpsjulia/mpsjulia.jl\");"
    self.execute(lambda: jlsession.eval(c)) # evaluate Julia



def install():
    """Install Julia and ITensor"""
    julia = "julia" # julia command
    os.system(julia+" --eval "+"\"import Pkg; Pkg.add(\\\"ITensors\\\")\"")
    os.system(julia+" --eval "+"\"import Pkg; Pkg.add(\\\"PyCall\\\")\"")
    os.system(julia+" --eval "+"\"import Pkg; Pkg.add(\\\"Suppressor\\\")\"")
    precompile()


def precompile():
    """Precompile Julia"""
    julia = "julia" # julia command
    os.system(julia+" --eval  \"using ITensors; ITensors.compile()\"") 

