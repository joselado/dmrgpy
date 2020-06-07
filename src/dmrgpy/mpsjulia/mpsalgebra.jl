function applyoperator()
	sites = get_sites()
	nameop = get_string("applyoperator_multioperator")
	A = read_mpo(nameop) # read the operator
	psi0 = load_mps(get_string("applyoperator_wf0"))
	iden_op = identity_mpo() # identity operator
	maxdim = get_int("maxm")
	cutoff = get_float("cutoff")
	psi1 = contract(A,psi0,maxdim=maxdim,cutoff=cutoff)
	psi1 = contract(iden_op,psi1,maxdim=maxdim) # this fixes a bug
	save_mps(get_string("applyoperator_wf1"),psi1)
end


function summps()
	sites = get_sites()
	psi1 = load_mps("summps_wf1.mps")
	psi2 = load_mps("summps_wf2.mps")
	maxdim = get_int("maxm")
	psi3 = add(psi1,psi2,maxdim=maxdim)
	save_mps("summps_wf3.mps",psi3)
end



function exponential()
	dtr = get_float("tevol_dt_real")
	dti = get_float("tevol_dt_imag")
	maxdim = get_int("maxm")
	cutoff = get_float("cutoff")
	H = read_mpo("hamiltonian.in")  # input wavefunction
	psi1 = load_mps("input_wavefunction.mps") # Operator to exponentiate
	nt0 = get_int("tevol_n") # number of steps
	tau = dtr + im*dti # evolution time
	taui = tau/nt0 # small time step
	iden_op = identity_mpo() # identity operator
        for i=1:nt0
          	psi2 = taui*contract(H,1*psi1;maxdim=maxdim,cutoff=cutoff) # sum
          	psi2 = contract(iden_op,psi2;maxdim=maxdim,cutoff=cutoff) # fix
          	psi1 = add(1*psi1,1*psi2;maxdim=maxdim,cutoff=cutoff) # sum
		truncate!(psi1;maxdim=maxdim,cutoff=cutoff)
	end
	save_mps("output_wavefunction.mps",psi1) # write result
end




function overlap()
	wf1 = load_mps("overlap_wf1.mps")
	wf2 = load_mps("overlap_wf2.mps")
	c = inner(wf1,wf2)
	write_in_file("OVERLAP.OUT",convert(Complex,c),"w")
end



