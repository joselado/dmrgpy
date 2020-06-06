
function get_gs()
  sites = get_sites()
  return get_gs(sites)
end

function get_gs(sites)
  ampo = read_operator("hamiltonian.in")
  H = MPO(ampo,sites)
  if get_bool("gs_from_file")
	  psi0 = load_mps("psi_GS.mps")
	  if get_bool("skip_dmrg_gs") 
		  return psi0
	  end
  else
          psi0 = randomMPS(sites)
  end
  sweeps = Sweeps(get_int("nsweeps"))
  maxdim!(sweeps, get_int("maxm"))
  cutoff!(sweeps, get_float("cutoff"))
  energy, psi = dmrg(H,psi0, sweeps)
  write_in_file("GS_ENERGY.OUT", energy)
  save_mps("psi_GS.mps",psi)
  return psi
end

