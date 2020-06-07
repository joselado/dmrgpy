using Pkg # using the pre-compiled package
#pkg"add https://github.com/TsurHerman/Fezzik"
#Pkg.add("JLD")
Pkg.add("ITensors")
Pkg.add("PyCall")
Pkg.add("Suppressor")
#using Fezzik
#Fezzik.auto_trace()
#@time include("setup.jl") # this executes all the basic orders
#Fezzik.brute_build_julia() # and precompile everything
