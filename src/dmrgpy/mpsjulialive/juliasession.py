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
    files = ["sites_custom.jl"] # site types ITensors.jl lacks, before sites.jl
    files += ["sites.jl"]
    files += ["get_gs.jl"]
    files += ["read_operator.jl"]
    files += ["mpsalgebra.jl"]
    files += ["kpm.jl"] # KPM moment recursion; calls mpsalgebra.jl's own
                         # apply_op/summps, so must load after it
    files += ["tdvp.jl"] # real-time TDVP evolution; also calls
                          # mpsalgebra.jl's apply_op, must load after it
    files += ["metts.jl"] # METTS finite-temperature sampling; reuses
                           # tdvp.jl's tdvp_step for imaginary-time
                           # evolution, must load after it
    files += ["tebd.jl"] # real-time TEBD (Trotter-gate) evolution; reuses
                          # tdvp.jl's apply_clean for observable
                          # application, must load after it
    files += ["excited.jl"] # orthogonality-penalty excited-state dmrg()
    files += ["nhdmrg.jl"] # biorthogonal non-Hermitian DMRG pair; uses
                            # get_gs.jl's make_sweeps/run_quiet, so must
                            # load after it
    files += ["generalized.jl"] # generalized-eigenvalue (Hermitian and
                                 # non-Hermitian) DMRG; uses get_gs.jl's
                                 # make_sweeps/run_quiet and nhdmrg.jl's
                                 # nh_biorthogonal_pair, so must load
                                 # after both
    files += ["densitymatrix.jl"] # single-site reduced density matrix
    files += ["entropy.jl"] # SVD-based bond entanglement entropy
    files += ["fermion_ops.jl"] # bare A/Adag on Fermion sites, for the
                                # Jordan-Wigner-dressed parity-odd MPOs
    files += ["mpo_ops.jl"] # MPO product, sum, dagger, trace, <a|A|b>
    files += ["bicgstab.jl"] # MPS BiCGSTAB for applyinverse; uses
                              # mpsalgebra.jl's apply_op/summps/conj_mps
    files += ["nhkpm.jl"] # non-Hermitian KPM moments; uses mpsalgebra.jl's
                           # apply_op/summps
    files += ["arrays.jl"] # MPS/MPO site tensors as arrays, for TASK
    files += ["correlations.jl"] # native <Cdag_i C_j> correlation matrix
    files += ["gs_extra.jl"] # exact <a|M|b>, uncapped MPO application
                              # and the four-point tensor
    files += ["sites_qn.jl"] # QN site sets of conserved-sector mode
    for name in files: # loop over files
        Main.seval(open(path+"/"+name).read()) # execute this file


initialize() # initialize all the functions
