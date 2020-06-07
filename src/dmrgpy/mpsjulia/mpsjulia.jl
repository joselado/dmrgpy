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
	if get_bool("applyoperator") # apply an operator
		applyoperator()
	end
	if get_bool("general_kpm") # apply an operator
		general_kpm()
	end
	if get_bool("overlap") # apply an operator
		overlap()
	end
	if get_bool("summps") # sum two mps
		summps()
	end
	if get_bool("exponential_eMwf") # apply an operator
		exponential()
	end
end
