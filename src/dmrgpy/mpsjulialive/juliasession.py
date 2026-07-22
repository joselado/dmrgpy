## Bridge to a live, in-process Julia session (ITensors.jl).
#
# Uses juliacall/PythonCall.jl rather than PyJulia: no libpython-matching
# precondition, no compiled_modules=False startup workaround, and the
# required Julia packages (ITensors, ITensorMPS, ITensorNHDMRG) are
# declared in ../juliapkg.json and auto-installed/precompiled into a
# juliacall-managed project the first time this module is imported --
# reusing whatever Julia install juliapkg finds (or provisioning its own
# via juliaup if none is found).

import os
from juliacall import Main
from juliacall import convert as jlconvert


def to_julia_strvec(ls):
    """Convert a Python list of str into a genuine Julia Vector{String}.
    juliacall wraps a plain Python list as a lazy PyList{Any} by default,
    which fails to dispatch against Julia functions declared with a typed
    Vector{String} argument (sites.jl's get_sites, read_operator.jl's
    read_operator/toMPO) -- unlike PyJulia, which converted a homogeneous
    list of str to Vector{String} implicitly."""
    return jlconvert(Main.Vector[Main.String], ls)


# current path
path = os.path.dirname(os.path.realpath(__file__))


# now execute the Julia code

# list of all the source files

def initialize():
    files = ["sites.jl"]
    files += ["get_gs.jl"]
    files += ["read_operator.jl"]
    files += ["mpsalgebra.jl"]
    files += ["kpm.jl"] # KPM moment recursion; calls mpsalgebra.jl's own
                         # applyoperator/summps, so must load after it
    files += ["tdvp.jl"] # real-time TDVP evolution; also calls
                          # mpsalgebra.jl's applyoperator, must load after it
    files += ["excited.jl"] # orthogonality-penalty excited-state dmrg()
    files += ["densitymatrix.jl"] # single-site reduced density matrix
    for name in files: # loop over files
        Main.seval(open(path+"/"+name).read()) # execute this file


initialize() # initialize all the functions
