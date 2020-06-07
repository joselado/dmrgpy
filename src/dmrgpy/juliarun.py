# routines to run the code with Julia
import os

dmrgpath = os.path.dirname(os.path.realpath(__file__))

precompile = False

from julia.api import Julia
jlsession = Julia(compiled_modules=False) # start the Julia session
jlsession.eval("using Suppressor") # suppress output


def run(self):
    if precompile:
        mpsjl = dmrgpath+"/mpsjulia/mpsjulia.jl"
        txt = "using Fezzik\n" # Fezzik package
        txt += "Fezzik.auto_trace()\n"
        txt += "include(\""+mpsjl+"\")\n"
        txt += "Fezzik.brute_build_julia()\n"
        def fun():
            f = open("run_and_compile.jl","w")
            f.write(txt)
            f.close()
        self.execute(fun) # write the compilation script
   #     exit()
        print("Precompiling Julia, this may take a while")
        self.execute(lambda: os.system("julia run_and_compile.jl > status.jl"))
        print("FINISHED precompiling Julia")
    else:
        import contextlib
        c = "@suppress_out include(\""+dmrgpath+"/mpsjulia/mpsjulia.jl\");"
        self.execute(lambda: jlsession.eval(c)) # evaluate Julia





