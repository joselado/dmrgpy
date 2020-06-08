using Pkg # using the pre-compiled package
#pkg"add https://github.com/TsurHerman/Fezzik"
#Pkg.add("JLD")
#pkg"add https://github.com/JuliaComputing/MKL.jl" # MKL for linear algebra
Pkg.add("ITensors") # install ITensor
Pkg.add("PyCall") # communication between Python and Julia
Pkg.add("Suppressor") # to hide Julia outputs
#using Fezzik
#Fezzik.auto_trace()
#@time include("setup.jl") # this executes all the basic orders
#Fezzik.brute_build_julia() # and precompile everything
