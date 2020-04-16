#include("common.jl")
push!(LOAD_PATH, @__DIR__) # include this path
using common

let
	if get_bool("GS") # ground state energy
		get_gs()
	end
	if get_bool("vev") # vacuum expectation value
		get_vev()
	end
	if get_bool("many_vev") # vacuum expectation value
		get_many_vev()
	end
	if get_bool("dynamical_correlator") # vacuum expectation value
		dynamical_correlator_kpm()
	end
end
