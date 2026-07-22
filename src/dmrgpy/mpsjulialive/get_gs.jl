
using ITensorMPS
using ITensors
using ITensorNHDMRG


# Shared by every mpsjulialive/*.jl DMRG-family entry point (get_gs_dmrg/
# get_gs_nhdmrg here, excited_state_dmrg in excited.jl) -- factored out
# after the same Sweeps-construction + stdout-silencing pair had been
# copy-pasted a third time in excited.jl. Defined here (loaded before
# excited.jl, see juliasession.py's initialize() file list) but usable
# from any later-loaded file: Julia resolves a global function call at
# call time, not at the calling function's definition time, so load
# order only needs to put this before the point of first *use*, which is
# well after every mpsjulialive/*.jl file has finished loading.
function make_sweeps(nsweeps,maxm,cutoff)
  sweeps = Sweeps(nsweeps)
  maxdim!(sweeps, maxm)
  cutoff!(sweeps, cutoff)
  return sweeps
end

function run_quiet(f)
  open("/dev/null", "w") do devnull
    redirect_stdout(devnull) do
      return f()
    end
  end
end


function get_gs_dmrg(H,psi0; nsweeps = 10,cutoff=1e-8,maxm=80,
	ishermitian=true)
  sweeps = make_sweeps(nsweeps,maxm,cutoff)
  run_quiet() do
    dmrg(H,psi0,sweeps,ishermitian=ishermitian)
  end
end




function get_gs_nhdmrg(H,psi0; nsweeps = 10,cutoff=1e-8,maxm=80,
        alg="onesided",biorthoalg = "biorthoblock",
	eigsolve_krylovdim=30,eigsolve_maxite=3)
  sweeps = make_sweeps(nsweeps,maxm,cutoff)
  run_quiet() do
    e, wfl0, wfr0 = nhdmrg(
        H,
        psi0,
        psi0,
        sweeps;
        alg=alg,
        biorthoalg=biorthoalg,
        outputlevel=1,
        eigsolve_krylovdim=eigsolve_krylovdim,
        eigsolve_maxiter=eigsolve_maxite,
    )
    e,wfl0,wfr0
  end
end





#### test the function
#function main() 
#    N = 20
#    sites = siteinds("S=1", N)
#    os = OpSum()
#    for j in 1:(N - 1)
#        os += "Sz", j, "Sz", j + 1
#        os += 0.5, "S+", j, "S-", j + 1
#        os += 0.5, "S-", j, "S+", j + 1
#    end
#    H = MPO(os, sites)
#
#    # Create an initial random matrix product state
#    psi0 = random_mps(sites)
#    nsweeps = 5
#    maxdim = 20
#    cutoff = 1.0e-10
#    # Run the DMRG algorithm, returning energy
#    # (dominant eigenvalue) and optimized MPS
#    energy, psi = get_gs_dmrg(H, psi0; nsweeps=nsweeps, maxm=maxdim, cutoff=cutoff)
#    println("Final energy = $energy")
#end
#
#
#main()
