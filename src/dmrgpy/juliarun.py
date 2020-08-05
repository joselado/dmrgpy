# routines to run the code with Julia
import os
import subprocess

dmrgpath = os.path.dirname(os.path.realpath(__file__))

try:
    out,err = subprocess.Popen(['julia', '--version'],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT).communicate()
except: out = ""
# check if Julia has the correct version
hasjulia = "julia version 1.4" in str(out)

if not hasjulia: # no julia present, add manually to the path
  jpath = dmrgpath + "/../julia/julia-1.4.2/bin/"
  os.environ["PATH"] = jpath+":"+os.environ["PATH"]



precompile = False

# check the system image
sysimage = os.environ["HOME"]+"/.julia/sysimages/sys_itensors.so"
if not os.path.isfile(sysimage):
    print("No ITensors system image found, this may take some time")
    sysimage = None


from julia.api import Julia
jlsession = Julia(compiled_modules=False,
        sysimage=sysimage) # start the Julia session
jlsession.eval("using Suppressor") # suppress output


def run(self):
    if precompile: # precompile the ITensor code
        jlsession.eval("using ITensors; ITensors.compile()")
#        mpsjl = dmrgpath+"/mpsjulia/mpsjulia.jl"
#        txt = "using Fezzik\n" # Fezzik package
#        txt += "Fezzik.auto_trace()\n"
#        txt += "include(\""+mpsjl+"\")\n"
#        txt += "Fezzik.brute_build_julia()\n"
#        def fun():
#            f = open("run_and_compile.jl","w")
#            f.write(txt)
#            f.close()
#        self.execute(fun) # write the compilation script
#   #     exit()
#        print("Precompiling Julia, this may take a while")
#        self.execute(lambda: os.system("julia run_and_compile.jl > status.jl"))
#        print("FINISHED precompiling Julia")
    else:
        import contextlib
        c = "@suppress_out include(\""+dmrgpath+"/mpsjulia/mpsjulia.jl\");"
        self.execute(lambda: jlsession.eval(c)) # evaluate Julia


