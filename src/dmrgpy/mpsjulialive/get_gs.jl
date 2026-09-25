
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
function make_sweeps(nsweeps,maxm,cutoff; noise=0.0, taper=true)
  sweeps = Sweeps(nsweeps)
  maxdim!(sweeps, maxm)
  cutoff!(sweeps, cutoff)
  # DMRG's density-matrix noise term (Many_Body_Chain.noise, default 1e-7,
  # manybodychain.py) -- what lets a sweep escape a local minimum. Every
  # other backend forwards it (mpscpp3/pyitensor both take it through
  # their own set_sweep_params), so a julia_live caller who leaves
  # self.noise at its default must get the same algorithm here, not a
  # silently noise-free one. Applied only when nonzero so the callers that
  # don't pass it (the non-Hermitian ones, see nhdmrg.jl's note) keep
  # their exact previous behaviour.
  #
  # taper=true puts it on the first half of the schedule only, the
  # sweeps 1..div(nsweeps,2) -- the schedule mpscpp2/mpscpp3's
  # make_sweeps and pyitensor's _make_sweeps build, whose own comment
  # records why: a flat noise returns the output of a noisy sweep as the
  # state and the energy. So a one-sweep schedule gets no noise here;
  # generalized.jl, which runs one sweep per outer iteration, tapers over
  # its iterations instead and passes taper=false.
  #
  # Nothing forwarded the noise at all until 2026-09-25 (hole hunt 25b,
  # finding 7): get_gs_dmrg called this without it, so with the
  # bond-dimension-1 start the Julia solves used to take, a Hamiltonian
  # coupling two sites across a decoupled one (every Thermal_Spin_Chain's
  # physical H on the doubled chain) stopped at a product state.
  nnoisy = taper ? div(nsweeps,2) : nsweeps
  if noise != 0.0 && nnoisy > 0
    noise!(sweeps, fill(float(noise), nnoisy)..., 0.0)
  end
  return sweeps
end

function run_quiet(f)
  open("/dev/null", "w") do devnull
    redirect_stdout(devnull) do
      return f()
    end
  end
end


# psi0 is dmrg_start_state's (mpsalgebra.jl) for a fresh solve and the
# caller's state for a warm start; noise is Many_Body_Chain.noise, see
# make_sweeps above.
function get_gs_dmrg(H,psi0; nsweeps = 10,cutoff=1e-8,maxm=80,
	ishermitian=true,noise=0.0)
  sweeps = make_sweeps(nsweeps,maxm,cutoff; noise=noise)
  run_quiet() do
    dmrg(H,psi0,sweeps,ishermitian=ishermitian)
  end
end




# Delegates to nhdmrg.jl's nhdmrg_raw rather than calling ITensorNHDMRG
# directly: the two used to be verbatim copies of the same invocation in
# two files, so any upstream keyword/default change would have been
# applied to one and silently missed by the other, leaving
# chain.gs_energy() and chain.nhdmrg() running differently-configured
# solvers on the same chain. nhdmrg.jl loads *after* this file, which is
# fine -- Julia resolves a global function call at call time, and every
# .jl file has finished loading well before the first call (see
# make_sweeps' own note above).
#
# No tie-break here, deliberately: this entry point returns only the right
# eigenvector (groundstate.py's plain non-Hermitian gs_energy() path
# discards wfl0), and the tie-break exists solely to fix the LEFT vector.
# Running it would cost an extra solve to correct something this caller
# throws away.
function get_gs_nhdmrg(H,psi0; nsweeps = 10,cutoff=1e-8,maxm=80,
        alg="onesided",biorthoalg = "biorthoblock",
	eigsolve_krylovdim=30,eigsolve_maxite=3)
  sweeps = make_sweeps(nsweeps,maxm,cutoff)
  return nhdmrg_raw(H,psi0,psi0,sweeps; alg=alg,biorthoalg=biorthoalg,
      eigsolve_krylovdim=eigsolve_krylovdim,eigsolve_maxite=eigsolve_maxite)
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
