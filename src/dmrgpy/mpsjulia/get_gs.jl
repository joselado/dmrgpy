using ITensors

function get_gs(sites)
  ampo = read_operator("hamiltonian.in")
  H = MPO(ampo,sites)
  psi0 = randomMPS(sites)
  sweeps = Sweeps(get_int("nsweeps"))
  maxdim!(sweeps, get_int("maxm"))
  cutoff!(sweeps, get_float("cutoff"))
  energy, psi = dmrg(H,psi0, sweeps)
  return psi
end

