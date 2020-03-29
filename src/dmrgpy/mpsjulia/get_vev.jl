using ITensors
include("common.jl")

function get_vev()
    sites = get_sites()
    wf = get_gs(sites) # get the ground state wavefunction
    Am = read_operator("vev_multioperator.in") # read operator
    A = MPO(Am,sites)
    out = dot(wf,dot(A,wf))
    print(out)
end

get_vev()
