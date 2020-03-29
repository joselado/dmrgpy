__precompile__()
module common
using ITensors
export get_gs,get_bool,get_vev
include("read_operator.jl")
include("read_wf.jl")
include("get_input.jl")
include("get_sites.jl")
include("get_gs.jl")
include("write_in_file.jl")
include("get_vev.jl")
end
